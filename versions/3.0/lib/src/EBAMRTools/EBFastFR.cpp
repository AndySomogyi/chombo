#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EBArith.H"
#include "EBFastFR.H"
#include "LayoutIterator.H"
#include "BaseIVFactory.H"
#include "EBCoarToCoarRedist.H"
#include "EBCoarToFineRedist.H"
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "EBIndexSpace.H"
#include "parstream.H"
#include "EBAlias.H"
#include "NamespaceHeader.H"
#include "CH_Timer.H"
/*******************/
int
EBFastFR::
index(int dir, Side::LoHiSide side)
{
  CH_assert(dir >= 0);
  CH_assert(dir < SpaceDim);
  CH_assert((side == Side::Lo) || (side == Side::Hi));

  int ioffset;
  if(side == Side::Lo)
    ioffset = 0;
  else
    ioffset = 1;

  return ioffset*SpaceDim+dir;
}
/*******************/
EBFastFR::
EBFastFR(const EBLevelGrid& a_eblgFine,
         const EBLevelGrid& a_eblgCoar,
         const int&         a_refRat,
         const int&         a_nvar)
{
  m_isDefined = false;
  define(a_eblgFine, a_eblgCoar, a_refRat, a_nvar);
}
/*******************/
void
EBFastFR::
define(const EBLevelGrid& a_eblgFine,
       const EBLevelGrid& a_eblgCoar,
       const int&         a_refRat,
       const int&         a_nvar)
{
  CH_TIME("EBFastFR::define");


  m_refRat   = a_refRat;
  m_nComp    = a_nvar;
  m_eblgFine = a_eblgFine;
  m_eblgCoar = a_eblgCoar;
#if (CH_SPACEDIM == 2)
  m_nrefdmo = m_refRat;
#elif (CH_SPACEDIM == 3)
  m_nrefdmo = m_refRat*m_refRat;
#else
  bogus spacedim;
#endif

  m_levelFluxReg.define(a_eblgFine.getDBL(),
                        a_eblgCoar.getDBL(),
                        a_eblgFine.getDomain(),
                        a_refRat, a_nvar);

  //find the fine C-F cells on the coarse level and put them into an invectset
  LayoutData<IntVectSet> ivsCoarsenedFine;
  EBArith::defineCFIVS(ivsCoarsenedFine, m_eblgFine.getDBL(), m_eblgFine.getDomain());
  IntVectSet cfivstotloc;
  for(DataIterator dit = m_eblgFine.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      cfivstotloc |= ivsCoarsenedFine[dit()];
    }
  cfivstotloc.coarsen(a_refRat);

  Vector<IntVectSet> allTags;

  const int destProc = uniqueProc(SerialTask::compute);

  gather(allTags,cfivstotloc,destProc);

  IntVectSet cfivstot;
  if (procID() == uniqueProc(SerialTask::compute))
    {
      for (int i = 0; i < allTags.size(); ++i)
        {
          cfivstot |= allTags[i];
        }
    }
  broadcast(cfivstot, destProc);

  //these are the objects needed for caching the old solution
  m_vofitCoar.define(m_eblgCoar.getDBL());
  m_saveoldc. define(m_eblgCoar.getDBL());
  for(DataIterator dit = a_eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      IntVectSet ivsIrreg = m_eblgCoar.getEBISL()[dit()].getIrregIVS(a_eblgCoar.getDBL().get(dit()));
      ivsIrreg &= cfivstot;
      m_vofitCoar[dit()].define(ivsIrreg, m_eblgCoar.getEBISL()[dit()].getEBGraph());
      m_saveoldc [dit()].define(ivsIrreg, m_eblgCoar.getEBISL()[dit()].getEBGraph(), m_nComp);
    }

  m_hasEBCF = false;
  for(int idir=0 ; idir < SpaceDim; ++idir)
    {
      for(SideIterator side; side.ok(); ++side)
        {
          // step one, build fineBoxes, flux register boxes
          // indexed by the fine level but in the coar index
          // space
          int iindex = index(idir, side());
          DisjointBoxLayout& regBoxesFine = m_bufGridsFine[iindex];
          regBoxesFine.deepCopy(m_eblgFine.getDBL());

          for(LayoutIterator lit = regBoxesFine.layoutIterator();
              lit.ok(); ++lit)
            {
              Box& abutbox = regBoxesFine.ref(lit());
              abutbox.coarsen(a_refRat);
              if(side() == Side::Lo)
                abutbox = adjCellLo(abutbox, idir);
              else
                abutbox = adjCellHi(abutbox, idir);

              // its ok for adjacent cells to hang off the domain, that
              // gets resolved by the time refluxng happens.
              //because the invectset is intersected with the domain
            }
          regBoxesFine.close();

          //now figure out the register boxes on the coar level.
          DisjointBoxLayout& regBoxesCoar = m_bufGridsCoar[iindex];
          regBoxesCoar = DisjointBoxLayout();
          for(LayoutIterator litcf = regBoxesFine.layoutIterator();
              litcf.ok(); ++litcf)
            {
              const Box& cfRegBox = regBoxesFine.get(litcf());
              for(LayoutIterator litco = m_eblgCoar.getDBL().layoutIterator();
                  litco.ok(); ++litco)
                {
                  const Box& gridCoar = m_eblgCoar.getDBL().get(litco());
                  if(cfRegBox.intersectsNotEmpty(gridCoar))
                    {
                      unsigned int proc = m_eblgCoar.getDBL().procID(litco());
                      Box regBox = cfRegBox;
                      regBox &= gridCoar;
                      regBoxesCoar.addBox(regBox,proc);
                    }
                }
            }
          regBoxesCoar.close();

          //figure out the mapping---which coar boxes
          //intersect each of the coar  registers.
          //remember that the mapping goes from the coar input
          //layout to the coar registers.
          const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
          LayoutData<Vector<DataIndex> >& mapsLD = m_coarIndexMap[iindex];
          mapsLD.define(m_eblgCoar.getDBL());
          for(DataIterator dit = m_eblgCoar.getDBL().dataIterator();
              dit.ok(); ++dit)
            {
              const Box& coarGrid = m_eblgCoar.getDBL().get(dit());
              Vector<DataIndex>& localMap = mapsLD[dit()];
              localMap.resize(0);
              for(LayoutIterator litcr = regBoxesCoar.layoutIterator();
                  litcr.ok(); ++litcr)
                {
                  const Box& coarReg = regBoxesCoar.get(litcr());
                  if(coarGrid.intersectsNotEmpty(coarReg))
                    {
                      CH_assert(coarGrid.contains(coarReg));
                      DataIndex regIndex = DataIndex(litcr());
                      localMap.push_back(regIndex);
                    }
                }
            }

          //make the ebisl for the buffers.
          //both live on the coar index space
          EBISLayout& ebislIndexCoar = m_ebislBufCoar[iindex];
          EBISLayout& ebislIndexFine = m_ebislBufFine[iindex];

          //need one ghost so that it can check the vofs on either side
          //of the faces
          //remember that they are both in the coarse index space
          int nghost = 1;
          ebisPtr->fillEBISLayout(ebislIndexCoar, regBoxesCoar,
                                  m_eblgCoar.getDomain(), nghost);
          ebisPtr->fillEBISLayout(ebislIndexFine, regBoxesFine,
                                  m_eblgCoar.getDomain(),   nghost);
          //the fine one has to be able to refine by refRat
          ebislIndexFine.setMaxRefinementRatio(m_refRat, ebisPtr);

          // now, for every buffer,
          // we need an IntVectSet that indicates
          // which parts are ACTUALLY coar-fine interfaces, and which
          // parts are covered by another fine box.
          LayoutData<IntVectSet> & ivsldFine = m_cfivsFine[iindex];
          LayoutData<IntVectSet> & ivsldCoar = m_cfivsCoar[iindex];

          LayoutData<bool>       & emptyFine = m_noEBCFine[iindex];
          LayoutData<bool>       & emptyCoar = m_noEBCCoar[iindex];
          //could also define with m_gridsFine
          emptyFine.define(regBoxesFine);
          emptyCoar.define(regBoxesCoar);
          ivsldFine.define(regBoxesFine);
          ivsldCoar.define(regBoxesCoar);
          for(DataIterator ditrbf = regBoxesFine.dataIterator();
              ditrbf.ok(); ++ditrbf)
            {
              IntVectSet& ivsBuf = ivsldFine[ditrbf()];
              ivsBuf = IntVectSet(regBoxesFine.get(ditrbf()));
              //remember this is all in the coarse index space.
              ivsBuf &= m_eblgCoar.getDomain().domainBox();
              for(LayoutIterator litf = m_eblgFine.getDBL().layoutIterator();
                  litf.ok(); ++litf)
                {
                  ivsBuf -= coarsen(m_eblgFine.getDBL().get(litf()), m_refRat);
                }

              IntVectSet ivsIrreg = m_ebislBufFine[iindex][ditrbf()].getIrregIVS(regBoxesFine.get(ditrbf()));
              ivsBuf &=    ivsIrreg;
              emptyFine[ditrbf()] = ivsBuf.isEmpty();
              //m_hasEBCF goes to true if true anywhere
              bool haveEBCFHere = !emptyFine[ditrbf()];
              m_hasEBCF = (m_hasEBCF || haveEBCFHere);
            }
          
          for(DataIterator ditrbc = regBoxesCoar.dataIterator();
              ditrbc.ok(); ++ditrbc)
            {
              IntVectSet& ivsBuf = ivsldCoar[ditrbc()];
              ivsBuf = IntVectSet(regBoxesCoar.get(ditrbc()));
              ivsBuf &= m_eblgCoar.getDomain().domainBox();
              for(LayoutIterator litf = m_eblgFine.getDBL().layoutIterator();
                  litf.ok(); ++litf)
                {
                  ivsBuf -= coarsen(m_eblgFine.getDBL().get(litf()), m_refRat);
                }
              IntVectSet ivsIrreg  =   m_ebislBufCoar[iindex][ditrbc()].getIrregIVS(regBoxesCoar.get(ditrbc()));
              ivsBuf &= ivsIrreg;
              emptyCoar[ditrbc()] = ivsBuf.isEmpty();
              //m_hasEBCF goes to true if true anywhere
              bool haveEBCFHere = !emptyCoar[ditrbc()];
              m_hasEBCF = (m_hasEBCF || haveEBCFHere);
            }

   //in the case of parallel, need to check if ANY of the procs
  //have ebcrossing
#ifdef CH_MPI

  int gatherint = 0;
  if(m_hasEBCF) gatherint = 1;     
  int idoebcf;
  MPI_Allreduce(&gatherint, &idoebcf, 1, MPI_INT,
                MPI_MAX, Chombo_MPI::comm);
  
  m_hasEBCF = (idoebcf==1);

#endif
         

          // now define  the registers
          BaseIVFactory<Real> factFine(ebislIndexFine, ivsldFine);
          BaseIVFactory<Real> factCoar(ebislIndexCoar, ivsldCoar);
          m_regsFine[iindex].define(regBoxesFine, m_nComp,
                                    IntVect::Zero, factFine);
          m_regsCoar[iindex].define(regBoxesCoar, m_nComp,
                                    IntVect::Zero, factCoar);
          //and scratch space for refluxing
          m_scratchc[iindex].define(regBoxesCoar, m_nComp,
                                    IntVect::Zero, factCoar);

          LayoutData<VoFIterator>& vofitFine = m_vofitBufFine[iindex];
          LayoutData<VoFIterator>& vofitCoar = m_vofitBufCoar[iindex];
          vofitFine.define(regBoxesFine);
          vofitCoar.define(regBoxesCoar);
          for(DataIterator ditrbf = regBoxesFine.dataIterator();  ditrbf.ok(); ++ditrbf)
            {
              IntVectSet& ivsBuf = ivsldFine[ditrbf()];
              if(!emptyFine[ditrbf()])
                {
                  vofitFine[ditrbf()].define(ivsBuf, ebislIndexFine[ditrbf()].getEBGraph());
                }
            }
          for(DataIterator ditrbc = regBoxesCoar.dataIterator();  ditrbc.ok(); ++ditrbc)
            {
              IntVectSet& ivsBuf = ivsldCoar[ditrbc()];
              if(!emptyCoar[ditrbc()])
                {
                  vofitCoar[ditrbc()].define(ivsBuf, ebislIndexCoar[ditrbc()].getEBGraph());
                }
            }
        }
    }
  //set all registers to zero to start.
  m_isDefined = true;
  setToZero();
}
/*******************/
EBFastFR::
EBFastFR()
{
  m_isDefined = false;
  m_nComp = -1;
  m_refRat = -1;
}
/*******************/
EBFastFR::~EBFastFR()
{
}
/*******************/
bool
EBFastFR::
isDefined() const
{
  return m_isDefined;
}
/*******************/
void
EBFastFR::
setToZero()
{
  CH_TIME("EBFastFR::setToZero");
  m_levelFluxReg.setToZero();

  if(m_hasEBCF)
    {
      irregSetToZero();
    }
}
void
EBFastFR::
irregSetToZero()
{
  CH_TIME("EBFastFR::irregSetToZero");
  for(int idir=0 ; idir<SpaceDim; ++idir)
    {
      for(SideIterator side; side.ok(); ++side)
        {
          int iindex = index(idir, side());
          LevelData<BaseIVFAB<Real> >& fineReg = m_regsFine[index(idir, side())];
          for(DataIterator dit = fineReg.dataIterator(); dit.ok(); ++dit)
            {
              if(!m_noEBCFine[iindex][dit()])
                {
                  fineReg[dit()].setVal(0.0);
                }
            }

          LevelData<BaseIVFAB<Real> >& coarReg = m_regsCoar[index(idir, side())];
          LevelData<BaseIVFAB<Real> >& scratch = m_scratchc[index(idir, side())];
          for(DataIterator dit = coarReg.dataIterator(); dit.ok(); ++dit)
            {
              if(!m_noEBCCoar[iindex][dit()])
                {
                  coarReg[dit()].setVal(0.0);
                  scratch[dit()].setVal(0.0);
                }
            }
        }
    }
}
/*******************/
void
EBFastFR::
incrementCoarseBoth(const EBFaceFAB& a_coarFlux,
                    const Real&      a_scale,
                    const DataIndex& a_coarDatInd,
                    const Interval&  a_variables,
                    const int&       a_dir,
                    const Side::LoHiSide& a_sd)
{
  CH_TIME("EBFastFR::incrementCoarseBoth");
  FArrayBox& coarFluxFAB = (FArrayBox&)(a_coarFlux.getFArrayBox());

  //increment  as if there were no EB
  m_levelFluxReg.incrementCoarse(coarFluxFAB, a_scale, a_coarDatInd,
                                 a_variables, a_variables, a_dir, a_sd);

  //fill buffers at irregular cells
  if(m_hasEBCF)
    {
      incrementCoarIrreg(a_coarFlux, a_scale, a_coarDatInd, a_variables, a_dir, a_sd);
    } //end if we have any EBCF in the whole class
}
/*******************/
void
EBFastFR::
incrementCoarIrreg(const EBFaceFAB&       a_coarFlux,
                   const Real&            a_scale,
                   const DataIndex&       a_coarDatInd,
                   const Interval&        a_variables,
                   const int&             a_dir,
                   const Side::LoHiSide& a_sd)
{
  CH_TIME("EBFastFR::incrementCoarseIrreg");


  int iindex = index(a_dir, a_sd);
  LevelData<BaseIVFAB<Real> >& registerLD = m_regsCoar[iindex];
  const EBISLayout& ebislBuf = m_ebislBufCoar[iindex];
  const Vector<DataIndex>& vecIndex =  m_coarIndexMap[iindex][a_coarDatInd];
  for(int iind = 0; iind < vecIndex.size(); iind++)
    {
      const DataIndex& regDatInd = vecIndex[iind];
      const LayoutData<bool>& boolLD = m_noEBCCoar[iindex];
      if(!(boolLD[regDatInd]))
        {
          BaseIVFAB<Real>& coarReg = registerLD[regDatInd];
          const EBISBox& ebisBox = ebislBuf[regDatInd];
          //increment coar buffer with the area-weighted sum
          //of the fluxes on the faces.
          //buf += sum(areaFrac*flux)
          VoFIterator vofit = m_vofitBufCoar[iindex][regDatInd];
          for(vofit.reset();    vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              for(int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++)
                {
                  //flip because the buffer is on the
                  //other side of the CF interface
                  Vector<FaceIndex> faces = ebisBox.getFaces(vof, a_dir, flip(a_sd));

                  Real areaSum = 0.0;
                  for(int iface = 0; iface < faces.size(); iface++)
                    {
                      const FaceIndex& face = faces[iface];
                      Real area = ebisBox.areaFrac(face);
                      Real flux = a_coarFlux(face, ivar);
                      areaSum += area*flux;
                    }
                  //consistent with non-eb version
                  areaSum *= -a_scale;
                  coarReg(vof, ivar) += areaSum;
                }//end iteration over vars
            }//end iteration over vofs in buffer
        } //end if(have ebcf)
    }//end loop over buffers
}
/*******************/
void
EBFastFR::
incrementFineBoth(const EBFaceFAB&      a_fineFlux,
                  const Real&           a_scale,
                  const DataIndex&      a_fineDatInd,
                  const Interval&       a_variables,
                  const int&            a_dir,
                  const Side::LoHiSide& a_sd)
{
  CH_TIME("EBFastFR::incrementFineBoth");
  //increment  as if there were no EB
  FArrayBox& fineFluxFAB = (FArrayBox&)(a_fineFlux.getFArrayBox());
  m_levelFluxReg.incrementFine(fineFluxFAB, a_scale, a_fineDatInd,
                               a_variables, a_variables, a_dir, a_sd);

  //fill the buffers at irregular cells
  if(m_hasEBCF)
    {
      Real newScale = a_scale/m_nrefdmo;
      incrementFineIrreg(a_fineFlux, newScale, a_fineDatInd, a_variables, a_dir, a_sd);
    } //end if (have EBCF in whole class)
}
/*******************/
void
EBFastFR::
incrementFineIrreg(const EBFaceFAB&      a_fineFlux,
                   const Real&           a_newScale,
                   const DataIndex&      a_fineDatInd,
                   const Interval&       a_variables,
                   const int&            a_dir,
                   const Side::LoHiSide& a_sd)
{
  CH_TIME("EBFastFR::incrementFineIrreg");
  int iindex = index(a_dir, a_sd);
  if(!m_noEBCFine[iindex][a_fineDatInd])
    {
      BaseIVFAB<Real>&  regFab = m_regsFine[iindex][a_fineDatInd];
      const EBISLayout& ebislBuf = m_ebislBufFine[iindex];

      const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[a_fineDatInd];
      VoFIterator& vofit = m_vofitBufFine[iindex][a_fineDatInd];
      for(vofit.reset(); vofit.ok(); ++vofit)
        {
          //remember the registers live at the coar level
          const VolIndex& coarVoF = vofit();
          //box to weed out which cells are on the  cf interface
          const IntVect& coarIndex = coarVoF.gridIndex();
          Box cfbox(coarIndex, coarIndex);
          cfbox.refine(m_refRat);
          if(a_sd == Side::Lo) cfbox.growLo(a_dir,-(m_refRat-1));
          else                   cfbox.growHi(a_dir,-(m_refRat-1));
          Vector<VolIndex> fineVoFs = ebislBuf.refine(coarVoF,m_refRat, a_fineDatInd);
          Vector<Real> areaSum(m_nComp, 0.0);

          //loop over fine vofs at the coar VoFs cf interface.
          for(int ifinev = 0; ifinev < fineVoFs.size(); ifinev++)
            {
              const VolIndex& fineVoF = fineVoFs[ifinev];
              const IntVect& fineIV = fineVoF.gridIndex();
              //the data for the fine cells is on the other side of the
              //CF interface so need to check for irregularity there
              //if(vof is on coar-fine interface)
              if(cfbox.contains(fineIV))
                {
                  Vector<FaceIndex> faces = ebisBoxFine.getFaces(fineVoF, a_dir, flip(a_sd));
                  for(int iface = 0; iface < faces.size(); iface++)
                    {
                      const FaceIndex& fineFace = faces[iface];
                      for(int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++)
                        {
                          Real area = ebisBoxFine.areaFrac(fineFace);
                          Real flux = a_fineFlux(fineFace, ivar);
                          //consistent with non-eb version.
                          areaSum[ivar] += a_newScale*area*flux;
                        }
                    } //end loop over faces at cf interface
                }//end check if the vof is at the coarfine interface
            }//end loop over fine vofs that compose coar vof
          //add area weighted sum to register
          for(int ivar = a_variables.begin();
              ivar <= a_variables.end(); ivar++)
            {
              regFab(coarVoF, ivar) += areaSum[ivar];
            }
        } //end loop over coar vofs in register
    } //end if we have EBCF in this box
}// mmm deep loops
/*******************/
void
EBFastFR::
reflux(LevelData<EBCellFAB>& a_uCoar,
       const Interval&       a_variables,
       const Real&           a_scale)
{
  CH_TIME("EBFastFR::reflux");
  LevelData<FArrayBox> uCoarLDF;
  //save initial values of ucoar because the non-eb flux  reg will
  //change it in its ignorance
  if(m_hasEBCF)
    {
      cacheOldSolution(a_uCoar, a_variables);
    }


  //reflux as if there were no EB
  aliasEB(uCoarLDF, a_uCoar);
  m_levelFluxReg.reflux(uCoarLDF, a_variables, a_variables, a_scale);

  //correct at irregular cells
  if(m_hasEBCF)
    {
      restoreOldSolution(a_uCoar, a_variables);
      irregReflux(a_uCoar, a_variables, a_scale);
    }
}
/*******************/
void
EBFastFR::
irregReflux(LevelData<EBCellFAB>& a_uCoar,
            const Interval&       a_variables,
            const Real&           a_scale)
{
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      for(SideIterator side; side.ok(); ++side)
        {
          int iindex = index(idir, side());
          const LevelData<BaseIVFAB<Real> >& fineReg     = m_regsFine[iindex];
          const LevelData<BaseIVFAB<Real> >& coarReg     = m_regsCoar[iindex];
          LevelData<BaseIVFAB<Real> >      & scratch     = m_scratchc[iindex];
          LayoutData<VoFIterator>          & vofitldCoar = m_vofitBufCoar[iindex];
          // make a local copy in the coar layout space of
          // the fine register increments

          fineReg.copyTo(a_variables, scratch, a_variables);

          for(DataIterator dit = m_eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit)
            {
              EBCellFAB& ufab = a_uCoar[dit()];
              Vector<DataIndex>&  datIndVec = m_coarIndexMap[iindex][dit()];
              for(int ivec=0; ivec < datIndVec.size(); ++ivec)
                {
                  const DataIndex& datInd = datIndVec[ivec];
                  if(!(m_noEBCCoar[iindex][datInd]))
                    {
                      const BaseIVFAB<Real>& coarRegFAB = coarReg[datInd];
                      const BaseIVFAB<Real>& fineRegFAB = scratch[datInd];
                      VoFIterator& vofit = vofitldCoar[datInd];
                      for(vofit.reset(); vofit.ok(); ++vofit)
                        {
                          const VolIndex& vof = vofit();
                          Real newScale = sign(side())*a_scale;
                          for(int icomp = a_variables.begin();
                              icomp <= a_variables.end(); icomp++)
                            {
                              //the non-eb version has these added and has
                              //the coar one multiplied by the negative of  the scale.
                              //so that is what we do.
                              //the -= is to cope with the fact that sign(side()) gives
                              //the wrong sign
                              ufab(vof, icomp) -= newScale*(fineRegFAB(vof, icomp) +
                                                            coarRegFAB(vof, icomp));
                            } //end loop over components
                        } //end loop over irreg vofs in buffer
                    }//end if(we have irreg cells in this buffer)
                }//end loop over buffers in this coarse box
            } //end data iterator loop
        } //end loop over sides
    } //end loop over directions
} //now that was a lot of curly braces
/*******************/
void
EBFastFR::
restoreOldSolution(LevelData<EBCellFAB>&       a_uCoar,
                   const Interval&             a_variables)
{

  CH_TIME("EBFastFR::saveOldSolution");
  for(DataIterator dit = m_eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      VoFIterator& vofit = m_vofitCoar[dit()];
      for(vofit.reset(); vofit.ok(); ++vofit)
        {
          for(int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++)
            {
              a_uCoar[dit()](vofit(), icomp) =  m_saveoldc[dit()](vofit(), icomp) ;
            }
        }
    }
}

/*******************/
void
EBFastFR::
cacheOldSolution(const LevelData<EBCellFAB>& a_uCoar,
                 const Interval&             a_variables)
{

  CH_TIME("EBFastFR::saveOldSolution");
  for(DataIterator dit = m_eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      VoFIterator& vofit = m_vofitCoar[dit()];
      for(vofit.reset(); vofit.ok(); ++vofit)
        {
          for(int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++)
            {
              m_saveoldc[dit()](vofit(), icomp) = a_uCoar[dit()](vofit(), icomp);
            }
        }
    }
}
/*******************/
void
EBFastFR::
incrementDensityArray(LevelData<EBCellFAB>& a_coarDense,
                      const Interval&       a_variables,
                      const Real&           a_scale)
{
  //only do stuff on irregular cells because on non-irregular cells 1-kappa=0
  //does not need to be fast because this is a debugging tool
  if(m_hasEBCF)
    {
      for(int idir=0 ; idir<SpaceDim; ++idir)
        {
          for(SideIterator side; side.ok(); ++side)
            {
              int iindex = index(idir, side());
              const LevelData<BaseIVFAB<Real> >& fineReg = m_regsFine[iindex];
              const LevelData<BaseIVFAB<Real> >& coarReg = m_regsCoar[iindex];
              LevelData<BaseIVFAB<Real> >& scratch =       m_scratchc[iindex];
              const LayoutData<IntVectSet>& ivsldCoar  =   m_cfivsCoar[iindex];
              const EBISLayout& ebislCoar = m_ebislBufCoar[iindex];
              // make a local copy in the coar layout space of
              // the fine register increments

              fineReg.copyTo(a_variables, scratch, a_variables);

              for(DataIterator dit = m_eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit)
                {
                  EBCellFAB& denfab = a_coarDense[dit()];
                  Vector<DataIndex>&  datIndVec = m_coarIndexMap[iindex][dit()];
                  for(int ivec=0; ivec < datIndVec.size(); ++ivec)
                    {
                      const DataIndex& datInd = datIndVec[ivec];
                      const BaseIVFAB<Real>& coarRegFAB = coarReg[datInd];
                      const BaseIVFAB<Real>& fineRegFAB = scratch[datInd];
                      const IntVectSet& ivs  = ivsldCoar[datInd];
                      const EBISBox& ebisBox = ebislCoar[datInd];
                      for(VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
                        {
                          const VolIndex& vof = vofit();
                          Real volfrac = ebisBox.volFrac(vof);
                          //multiply this by volfrac
                          //to increment a mass instead of a density
                          Real newScale = (1.0-volfrac)*sign(side())*a_scale;
                          for(int icomp = a_variables.begin();
                              icomp <= a_variables.end(); icomp++)
                            {
                              //the non-eb version has these added and has
                              //the coar one multiplied by the negative of  the scale.
                              //so that is what we do
                              //the -= is to cope with the fact that sign(side()) gives
                              //the wrong sign
                              denfab(vof, icomp) -= newScale*(fineRegFAB(vof, icomp) +
                                                              coarRegFAB(vof, icomp));
                            }
                        }
                    }
                }
            }
        }
    }
}
/*************************/
#include "NamespaceFooter.H"
