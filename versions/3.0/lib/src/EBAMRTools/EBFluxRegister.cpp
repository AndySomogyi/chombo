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
#include "EBFluxRegister.H"
#include "LayoutIterator.H"
#include "BaseIVFactory.H"
#include "EBCoarToCoarRedist.H"
#include "EBCoarToFineRedist.H"
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "EBIndexSpace.H"
#include "parstream.H"
#include "NamespaceHeader.H"
#include "CH_Timer.H"
/*******************/
void
EBFluxRegister::setDefaultValues()
{
  m_isDefined = false;
  m_nComp = -1;
  m_refRat = -1;
}
/*******************/
int EBFluxRegister::index(int dir, Side::LoHiSide side)
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
EBFluxRegister::EBFluxRegister(const DisjointBoxLayout& a_dblFine,
                               const DisjointBoxLayout& a_dblCoar,
                               const EBISLayout&        a_ebislFine,
                               const EBISLayout&        a_ebislCoar,
                               const Box&               a_domainCoar,
                               const int&               a_refRat,
                               const int&               a_nvar,
                               const EBIndexSpace*      ebisPtr)
{
  setDefaultValues();
  define(a_dblFine, a_dblCoar, a_ebislFine, a_ebislCoar,
         a_domainCoar, a_refRat, a_nvar, ebisPtr);
}
/*******************/
void
EBFluxRegister::define(const DisjointBoxLayout& a_gridsFine,
                       const DisjointBoxLayout& a_gridsCoar,
                       const EBISLayout&        a_ebislFine,
                       const EBISLayout&        a_ebislCoar,
                       const ProblemDomain&     a_domainCoar,
                       const int&               a_refRat,
                       const int&               a_nvar,
                       const EBIndexSpace*      ebisPtr)
{
  CH_TIME("EBFluxRegister::define");
  CH_assert(a_refRat > 0);
  CH_assert(a_nvar > 0);
  CH_assert(!a_domainCoar.isEmpty());

  m_nComp = a_nvar;
  m_refRat  = a_refRat;

  m_isDefined = true;
  m_domainCoar = a_domainCoar;
  //  Box domainFine = refine(m_domainCoar, m_refRat);
  m_gridsFine = a_gridsFine;
  m_gridsCoar = a_gridsCoar;
  m_ebislFine = a_ebislFine;
  m_ebislCoar = a_ebislCoar;

  for(int idir=0 ; idir < SpaceDim; ++idir)
    {
      for(SideIterator side; side.ok(); ++side)
        {
          // step one, build fineBoxes, flux register boxes
          // indexed by the fine level but in the coar index
          // space
          int iindex = index(idir, side());
          DisjointBoxLayout& regBoxesFine = m_bufGridsFine[iindex];
          regBoxesFine.deepCopy(m_gridsFine);

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
              for(LayoutIterator litco = m_gridsCoar.layoutIterator();
                  litco.ok(); ++litco)
                {
                  const Box& gridCoar = m_gridsCoar.get(litco());
                  if(cfRegBox.intersectsNotEmpty(gridCoar))
                    {
                      unsigned int proc = m_gridsCoar.procID(litco());
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
          LayoutData<Vector<DataIndex> >& mapsLD = m_coarIndexMap[iindex];
          mapsLD.define(m_gridsCoar);
          for(DataIterator dit = m_gridsCoar.dataIterator();
              dit.ok(); ++dit)
            {
              const Box& coarGrid = m_gridsCoar.get(dit());
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

          // now, for every buffer,
          // we need an IntVectSet that indicates
          // which parts are ACTUALLY coar-fine interfaces, and which
          // parts are covered by another fine box.
          LayoutData<IntVectSet>& ivsldFine = m_cfivsFine[iindex];
          LayoutData<IntVectSet>& ivsldCoar = m_cfivsCoar[iindex];
          //could also define with m_gridsFine
          ivsldFine.define(regBoxesFine);
          ivsldCoar.define(regBoxesCoar);
          for(DataIterator ditrbf = regBoxesFine.dataIterator();
              ditrbf.ok(); ++ditrbf)
            {
              IntVectSet& ivsBuf = ivsldFine[ditrbf()];
              ivsBuf = IntVectSet(regBoxesFine.get(ditrbf()));
              //remember this is all in the coarse index space.
              ivsBuf &= m_domainCoar;
              for(LayoutIterator litf = m_gridsFine.layoutIterator();
                  litf.ok(); ++litf)
                {
                  ivsBuf -= coarsen(m_gridsFine.get(litf()), m_refRat);
                }
            }
          for(DataIterator ditrbc = regBoxesCoar.dataIterator();
              ditrbc.ok(); ++ditrbc)
            {
              IntVectSet& ivsBuf = ivsldCoar[ditrbc()];
              ivsBuf = IntVectSet(regBoxesCoar.get(ditrbc()));
              ivsBuf &= m_domainCoar;
              for(LayoutIterator litf = m_gridsFine.layoutIterator();
                  litf.ok(); ++litf)
                {
                  ivsBuf -= coarsen(m_gridsFine.get(litf()), m_refRat);
                }
            }

          //make the ebisl for the buffers.
          //both live on the coar index space
          EBISLayout ebislIndexCoar = m_ebislBufCoar[iindex];
          EBISLayout ebislIndexFine = m_ebislBufFine[iindex];

          //need one ghost so that it can check the vofs on either side
          //of the faces
          //remember that they are both in the coarse index space
          int nghost = 1;
          ebisPtr->fillEBISLayout(ebislIndexCoar, regBoxesCoar,
                                  m_domainCoar, nghost);
          ebisPtr->fillEBISLayout(ebislIndexFine, regBoxesFine,
                                  m_domainCoar,   nghost);
          //the fine one has to be able to refine by refRat
          ebislIndexFine.setMaxRefinementRatio(m_refRat, ebisPtr);
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
        }
    }
  //set all registers to zero to start.
  setToZero();
}
/*******************/
EBFluxRegister::EBFluxRegister()
{
  setDefaultValues();
}
/*******************/
EBFluxRegister::~EBFluxRegister()
{
}
/*******************/
bool
EBFluxRegister::isDefined() const
{
  return m_isDefined;
}
/*******************/
void
EBFluxRegister::setToZero()
{
  CH_TIME("EBFluxRegister::setToZero");
  for(int idir=0 ; idir<SpaceDim; ++idir)
    {
      for(SideIterator side; side.ok(); ++side)
        {
          LevelData<BaseIVFAB<Real> >& fineReg = m_regsFine[index(idir, side())];
          LevelData<BaseIVFAB<Real> >& coarReg = m_regsCoar[index(idir, side())];

          for(DataIterator dit = fineReg.dataIterator(); dit.ok(); ++dit)
            fineReg[dit()].setVal(0.0);

          for(DataIterator dit = coarReg.dataIterator(); dit.ok(); ++dit)
            coarReg[dit()].setVal(0.0);

        }
    }
}
/*******************/
void
EBFluxRegister::incrementCoarseBoth(const EBFaceFAB& a_coarFlux,
                                    const Real&      a_scale,
                                    const DataIndex& a_coarDatInd,
                                    const Interval&  a_variables,
                                    const int&       a_dir)
{
  CH_TIME("EBFluxRegister::incrementCoarseBoth");
  incrementCoarseRegular(  a_coarFlux, a_scale, a_coarDatInd, a_variables, a_dir);
  incrementCoarseIrregular(a_coarFlux, a_scale, a_coarDatInd, a_variables, a_dir);
}
/*******************/

/*******************/
void
EBFluxRegister::incrementCoarseBothRZ(const EBFaceFAB& a_coarFlux,
                                      const Real&      a_scale,
                                      const DataIndex& a_coarDatInd,
                                      const Interval&  a_variables,
                                      const int&       a_dir,
                                      const Real&      a_dx)
{
  incrementCoarseRegulRZ(a_coarFlux, a_scale, a_coarDatInd, a_variables, a_dir, a_dx);
  IntVectSet ivs(a_coarFlux.getCellRegion());
  const EBISBox& ebisBoxCoar = a_coarFlux.getEBISBox();
  BaseIFFAB<Real> irregFlux(ivs, ebisBoxCoar.getEBGraph(), a_dir, a_coarFlux.nComp());

  for(FaceIterator faceit(ivs, ebisBoxCoar.getEBGraph(), a_dir,
                          FaceStop::SurroundingWithBoundary);
      faceit.ok(); ++faceit)
    {
      for(int icomp = 0; icomp < a_coarFlux.nComp(); icomp++)
        {
          irregFlux(faceit(), icomp) = a_coarFlux(faceit(), icomp);
        }
    }
  incrementCoarseIrregulRZ(irregFlux, a_scale, a_coarDatInd, a_variables, a_dir, a_dx);
}
/*******************/
void
EBFluxRegister::incrementFineBoth(const EBFaceFAB&      a_fineFlux,
                                  const Real&           a_scale,
                                  const DataIndex&      a_fineDatInd,
                                  const Interval&       a_variables,
                                  const int&            a_dir,
                                  const Side::LoHiSide& a_sd)
{
  CH_TIME("EBFluxRegister::incrementFineBoth");
  incrementFineRegular(a_fineFlux,   a_scale, a_fineDatInd, a_variables, a_dir, a_sd);
  incrementFineIrregular(a_fineFlux, a_scale, a_fineDatInd, a_variables, a_dir, a_sd);
}
/*******************/
/*******************/
void
EBFluxRegister::incrementFineBothRZ(const EBFaceFAB&      a_fineFlux,
                                    const Real&           a_scale,
                                    const DataIndex&      a_fineDatInd,
                                    const Interval&       a_variables,
                                    const int&            a_dir,
                                    const Side::LoHiSide& a_sd,
                                    const Real&           a_dx)
{
  incrementFineRegulRZ(a_fineFlux, a_scale, a_fineDatInd, a_variables, a_dir, a_sd, a_dx);
  IntVectSet ivs(a_fineFlux.getCellRegion());
  const EBISBox& ebisBoxFine = a_fineFlux.getEBISBox();
  BaseIFFAB<Real> irregFlux(ivs, ebisBoxFine.getEBGraph(), a_dir, a_fineFlux.nComp());

  for(FaceIterator faceit(ivs, ebisBoxFine.getEBGraph(), a_dir,
                          FaceStop::SurroundingWithBoundary);
      faceit.ok(); ++faceit)
    {
      for(int icomp = 0; icomp < a_fineFlux.nComp(); icomp++)
        {
          irregFlux(faceit(), icomp) = a_fineFlux(faceit(), icomp);
        }
    }
  incrementFineIrregulRZ(irregFlux, a_scale, a_fineDatInd, a_variables, a_dir, a_sd, a_dx);
}
/*******************/
void
EBFluxRegister::incrementCoarseRegular(const EBFaceFAB& a_coarFlux,
                                       const Real&      a_scale,
                                       const DataIndex& a_coarDatInd,
                                       const Interval&  a_variables,
                                       const int&       a_dir)
{
  CH_TIME("EBFluxRegister::incrementCoarseRegular");
  CH_assert(isDefined());
  CH_assert(a_dir >= 0);
  CH_assert(a_dir < SpaceDim);
  CH_assert(a_coarFlux.direction() == a_dir);
  CH_assert(a_variables.begin() >= 0);
  CH_assert(a_variables.end() < m_nComp);
  CH_assert(a_variables.end() < a_coarFlux.nComp());

  //iterate over sides of the fine grid.
  for(SideIterator side; side.ok(); ++side)
    {
      int iindex = index(a_dir, side());

      LevelData<BaseIVFAB<Real> >& registerLD = m_regsCoar[iindex];
      const EBISLayout& ebislBuf = m_ebislBufCoar[iindex];
      const LayoutData<IntVectSet>& ldivs = m_cfivsCoar[iindex];

      const Vector<DataIndex>& vecIndex =
        m_coarIndexMap[iindex][a_coarDatInd];

      for(int iind = 0; iind < vecIndex.size(); iind++)
        {
          const DataIndex& regDatInd = vecIndex[iind];
          BaseIVFAB<Real>& coarReg = registerLD[regDatInd];
          const EBISBox& ebisBox = ebislBuf[regDatInd];
          const IntVectSet& cfivs = ldivs[regDatInd];
          //increment coar buffer with the area-weighted sum
          //of the fluxes on the faces.
          //buf += sum(areaFrac*flux)
          for(VoFIterator vofit(cfivs, ebisBox.getEBGraph());
              vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              if(ebisBox.isRegular(vof.gridIndex()))
                {
                  for(int ivar = a_variables.begin();
                      ivar <= a_variables.end(); ivar++)
                    {
                      //flip because the buffer is on the
                      //other side of the CF interface
                      Vector<FaceIndex> faces =
                        ebisBox.getFaces(vof, a_dir, flip(side()));
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
                    }
                } //end if(isRegular())
            } //end iteration over vofs in buffer
        } //end loop over buffers
    } //end iteration over sides
} //end incrementCoar
/*******************/
void
EBFluxRegister::incrementCoarseIrregular(const EBFaceFAB&       a_coarFlux,
                                         const Real&            a_scale,
                                         const DataIndex&       a_coarDatInd,
                                         const Interval&        a_variables,
                                         const int&             a_dir)
{
  CH_TIME("EBFluxRegister::incrementCoarseIrregular 2");
  CH_assert(isDefined());
  CH_assert(a_dir >= 0);
  CH_assert(a_dir < SpaceDim);
  CH_assert(a_coarFlux.direction() == a_dir);
  CH_assert(a_variables.begin() >= 0);
  CH_assert(a_variables.end() < m_nComp);
  CH_assert(a_variables.end() < a_coarFlux.nComp());

  //iterate over sides of the fine grid.
  for(SideIterator side; side.ok(); ++side)
    {
      int iindex = index(a_dir, side());

      LevelData<BaseIVFAB<Real> >& registerLD = m_regsCoar[iindex];
      const EBISLayout& ebislBuf = m_ebislBufCoar[iindex];
      const LayoutData<IntVectSet>& ldivs = m_cfivsCoar[iindex];

      const Vector<DataIndex>& vecIndex =
        m_coarIndexMap[iindex][a_coarDatInd];

      for(int iind = 0; iind < vecIndex.size(); iind++)
        {
          const DataIndex& regDatInd = vecIndex[iind];
          BaseIVFAB<Real>& coarReg = registerLD[regDatInd];
          const EBISBox& ebisBox = ebislBuf[regDatInd];
          const IntVectSet& cfivs = ldivs[regDatInd];
          //increment coar buffer with the area-weighted sum
          //of the fluxes on the faces.
          //buf += sum(areaFrac*flux)
          for(VoFIterator vofit(cfivs, ebisBox.getEBGraph());
              vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              if(ebisBox.isIrregular(vof.gridIndex()))
                {
                  for(int ivar = a_variables.begin();
                      ivar <= a_variables.end(); ivar++)
                    {
                      //flip because the buffer is on the
                      //other side of the CF interface
                      Vector<FaceIndex> faces =
                        ebisBox.getFaces(vof, a_dir, flip(side()));
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
                    }
                } //end if(isIrregularVoF)
            } //end iteration over vofs in buffer
        } //end loop over buffers
    } //end iteration over sides
} //end incrementCoar
/*******************/
/*******************/
void
EBFluxRegister::incrementCoarseIrregular(const BaseIFFAB<Real>& a_coarFlux,
                                         const Real&            a_scale,
                                         const DataIndex&       a_coarDatInd,
                                         const Interval&        a_variables,
                                         const int&             a_dir)
{
  CH_TIME("EBFluxRegister::incrementCoarseIrregular");
  CH_assert(isDefined());
  CH_assert(a_dir >= 0);
  CH_assert(a_dir < SpaceDim);
  CH_assert(a_coarFlux.direction() == a_dir);
  CH_assert(a_variables.begin() >= 0);
  CH_assert(a_variables.end() < m_nComp);
  CH_assert(a_variables.end() < a_coarFlux.nComp());

  //iterate over sides of the fine grid.
  for(SideIterator side; side.ok(); ++side)
    {
      int iindex = index(a_dir, side());

      LevelData<BaseIVFAB<Real> >& registerLD = m_regsCoar[iindex];
      const EBISLayout& ebislBuf = m_ebislBufCoar[iindex];
      const LayoutData<IntVectSet>& ldivs = m_cfivsCoar[iindex];

      const Vector<DataIndex>& vecIndex =
        m_coarIndexMap[iindex][a_coarDatInd];

      for(int iind = 0; iind < vecIndex.size(); iind++)
        {
          const DataIndex& regDatInd = vecIndex[iind];
          BaseIVFAB<Real>& coarReg = registerLD[regDatInd];
          const EBISBox& ebisBox = ebislBuf[regDatInd];
          const IntVectSet& cfivs = ldivs[regDatInd];
          //increment coar buffer with the area-weighted sum
          //of the fluxes on the faces.
          //buf += sum(areaFrac*flux)
          for(VoFIterator vofit(cfivs, ebisBox.getEBGraph());
              vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              if(ebisBox.isIrregular(vof.gridIndex()))
                {
                  for(int ivar = a_variables.begin();
                      ivar <= a_variables.end(); ivar++)
                    {
                      //flip because the buffer is on the
                      //other side of the CF interface
                      Vector<FaceIndex> faces =
                        ebisBox.getFaces(vof, a_dir, flip(side()));
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
                    }
                } //end if(isIrregularVoF)
            } //end iteration over vofs in buffer
        } //end loop over buffers
    } //end iteration over sides
} //end incrementCoar
/*******************/
void
EBFluxRegister::incrementFineRegular(const EBFaceFAB&      a_fineFlux,
                                     const Real&           a_scale,
                                     const DataIndex&      a_fineDatInd,
                                     const Interval&       a_variables,
                                     const int&            a_dir,
                                     const Side::LoHiSide& a_sd)
{
  CH_TIME("EBFluxRegister::incrementFineRegular");
  CH_assert(isDefined());
  CH_assert(a_fineFlux.direction() == a_dir);
  CH_assert(a_variables.begin() >= 0);
  CH_assert(a_variables.end() < a_fineFlux.nComp());
  CH_assert(a_variables.end() < m_nComp);
  CH_assert(a_dir >= 0);
  CH_assert(a_dir < SpaceDim);
  CH_assert((a_sd == Side::Lo)||(a_sd == Side::Hi));

  //denom is the number of fine faces per coar face
  //this is intrinsically dimension-dependent
  //yes, i could do some lame loop here, but that would
  //be more confusing than this
#if (CH_SPACEDIM == 2)
  Real nrefdmo = m_refRat;
#elif (CH_SPACEDIM == 3)
  Real nrefdmo = m_refRat*m_refRat;
#else
  bogus spacedim;
#endif
  Real newScale = a_scale/nrefdmo;

  int iindex = index(a_dir, a_sd);
  BaseIVFAB<Real>&  regFab = m_regsFine[iindex][a_fineDatInd];
  const IntVectSet& ivsReg = m_cfivsFine[iindex][a_fineDatInd];
  const EBISLayout& ebislBuf = m_ebislBufFine[iindex];
  const EBISBox& ebisBoxBuf = ebislBuf[a_fineDatInd];

  const EBISBox& ebisBoxFine = m_ebislFine[a_fineDatInd];
  for(VoFIterator vofit(ivsReg, ebisBoxBuf.getEBGraph()); vofit.ok(); ++vofit)
    {
      //remember the registers live at the coar level
      const VolIndex& coarVoF = vofit();
      //box to weed out which cells are on the  cf interface
      const IntVect& coarIndex = coarVoF.gridIndex();
      Box cfbox(coarIndex, coarIndex);
      cfbox.refine(m_refRat);
      if(a_sd == Side::Lo) cfbox.growLo(a_dir,-(m_refRat-1));
      else                 cfbox.growHi(a_dir,-(m_refRat-1));
      Vector<VolIndex> fineVoFs = ebislBuf.refine(coarVoF,m_refRat, a_fineDatInd);
      Vector<Real> areaSum(m_nComp, 0.0);
      //loop over fine vofs at the coar VoFs cf interface.
      for(int ifinev = 0; ifinev < fineVoFs.size(); ifinev++)
        {
          const VolIndex& fineVoF = fineVoFs[ifinev];
          const IntVect& fineIV = fineVoF.gridIndex();
          int isign = sign(a_sd);
          //the data for the fine cells is on the other side of the
          //CF interface so need to check for irregularity there
          IntVect fineIVInt = fineIV - isign*BASISV(a_dir);
          if((ebisBoxFine.isRegular(fineIVInt))&&(cfbox.contains(fineIV)))
            {
              Vector<FaceIndex> faces =
                ebisBoxFine.getFaces(fineVoF, a_dir, flip(a_sd));
              for(int iface = 0; iface < faces.size(); iface++)
                {
                  const FaceIndex& fineFace = faces[iface];
                  for(int ivar = a_variables.begin();
                      ivar <= a_variables.end(); ivar++)
                    {
                      Real area = ebisBoxFine.areaFrac(fineFace);
                      Real flux = a_fineFlux(fineFace, ivar);
                      //consistent with non-eb version.
                      areaSum[ivar] += newScale*area*flux;
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
}
/*******************/
void
EBFluxRegister::incrementFineIrregular(const BaseIFFAB<Real>& a_fineFlux,
                                       const Real&            a_scale,
                                       const DataIndex&       a_fineDatInd,
                                       const Interval&        a_variables,
                                       const int&             a_dir,
                                       const Side::LoHiSide&  a_sd)
{
  CH_TIME("EBFluxRegister::incrementFineIrregular");
  CH_assert(isDefined());
  CH_assert(a_fineFlux.direction() == a_dir);
  CH_assert(a_variables.begin() >= 0);
  CH_assert(a_variables.end() < a_fineFlux.nComp());
  CH_assert(a_variables.end() < m_nComp);
  CH_assert(a_dir >= 0);
  CH_assert(a_dir < SpaceDim);
  CH_assert((a_sd == Side::Lo)||(a_sd == Side::Hi));

  //nrefdmo is the number of fine faces per coar face
  //this is intrinsically dimension-dependent
  //yes, i could do some lame loop here, but that would
  //be more confusing than this
#if (CH_SPACEDIM == 2)
  Real nrefdmo = m_refRat;
#elif (CH_SPACEDIM == 3)
  Real nrefdmo = m_refRat*m_refRat;
#else
  bogus spacedim;
#endif
  Real newScale = a_scale/nrefdmo;

  int iindex = index(a_dir, a_sd);
  BaseIVFAB<Real>&  regFab = m_regsFine[iindex][a_fineDatInd];
  const IntVectSet& ivsReg = m_cfivsFine[iindex][a_fineDatInd];
  const EBISLayout& ebislBuf = m_ebislBufFine[iindex];
  const EBISBox& ebisBoxBuf = ebislBuf[a_fineDatInd];

  const EBISBox& ebisBoxFine = m_ebislFine[a_fineDatInd];
  for(VoFIterator vofit(ivsReg, ebisBoxBuf.getEBGraph()); vofit.ok(); ++vofit)
    {
      //remember the registers live at the coar level
      const VolIndex& coarVoF = vofit();
      //box to weed out which cells are on the  cf interface
      const IntVect& coarIndex = coarVoF.gridIndex();
      Box cfbox(coarIndex, coarIndex);
      cfbox.refine(m_refRat);
      if(a_sd == Side::Lo) cfbox.growLo(a_dir,-(m_refRat-1));
      else                 cfbox.growHi(a_dir,-(m_refRat-1));
      Vector<VolIndex> fineVoFs = ebislBuf.refine(coarVoF,m_refRat, a_fineDatInd);
      Vector<Real> areaSum(m_nComp, 0.0);
      //loop over fine vofs at the coar VoFs cf interface.
      for(int ifinev = 0; ifinev < fineVoFs.size(); ifinev++)
        {
          const VolIndex& fineVoF = fineVoFs[ifinev];
          const IntVect& fineIV = fineVoF.gridIndex();
          int isign =  sign(a_sd);
          //the data for the fine cells is on the other side of the
          //CF interface so need to check for irregularity there
          IntVect fineIVInt = fineIV - isign*BASISV(a_dir);
          //if(vof is on coar-fine interface)
          if((cfbox.contains(fineIV)) && (ebisBoxFine.isIrregular(fineIVInt)))
            {
              Vector<FaceIndex> faces =
                ebisBoxFine.getFaces(fineVoF, a_dir, flip(a_sd));
              for(int iface = 0; iface < faces.size(); iface++)
                {
                  const FaceIndex& fineFace = faces[iface];
                  for(int ivar = a_variables.begin();
                      ivar <= a_variables.end(); ivar++)
                    {
                      Real area = ebisBoxFine.areaFrac(fineFace);
                      Real flux = a_fineFlux(fineFace, ivar);
                      //consistent with non-eb version.
                      areaSum[ivar] += newScale*area*flux;
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
}
/*******************/
/*******************/
void
EBFluxRegister::incrementFineIrregular(const EBFaceFAB&       a_fineFlux,
                                       const Real&            a_scale,
                                       const DataIndex&       a_fineDatInd,
                                       const Interval&        a_variables,
                                       const int&             a_dir,
                                       const Side::LoHiSide&  a_sd)
{
  CH_TIME("EBFluxRegister::incrementFineIrregular 2");
  CH_assert(isDefined());
  CH_assert(a_fineFlux.direction() == a_dir);
  CH_assert(a_variables.begin() >= 0);
  CH_assert(a_variables.end() < a_fineFlux.nComp());
  CH_assert(a_variables.end() < m_nComp);
  CH_assert(a_dir >= 0);
  CH_assert(a_dir < SpaceDim);
  CH_assert((a_sd == Side::Lo)||(a_sd == Side::Hi));

  //nrefdmo is the number of fine faces per coar face
  //this is intrinsically dimension-dependent
  //yes, i could do some lame loop here, but that would
  //be more confusing than this
#if (CH_SPACEDIM == 2)
  Real nrefdmo = m_refRat;
#elif (CH_SPACEDIM == 3)
  Real nrefdmo = m_refRat*m_refRat;
#else
  bogus spacedim;
#endif
  Real newScale = a_scale/nrefdmo;

  int iindex = index(a_dir, a_sd);
  BaseIVFAB<Real>&  regFab = m_regsFine[iindex][a_fineDatInd];
  const IntVectSet& ivsReg = m_cfivsFine[iindex][a_fineDatInd];
  const EBISLayout& ebislBuf = m_ebislBufFine[iindex];
  const EBISBox& ebisBoxBuf = ebislBuf[a_fineDatInd];

  const EBISBox& ebisBoxFine = m_ebislFine[a_fineDatInd];
  for(VoFIterator vofit(ivsReg, ebisBoxBuf.getEBGraph()); vofit.ok(); ++vofit)
    {
      //remember the registers live at the coar level
      const VolIndex& coarVoF = vofit();
      //box to weed out which cells are on the  cf interface
      const IntVect& coarIndex = coarVoF.gridIndex();
      Box cfbox(coarIndex, coarIndex);
      cfbox.refine(m_refRat);
      if(a_sd == Side::Lo) cfbox.growLo(a_dir,-(m_refRat-1));
      else                 cfbox.growHi(a_dir,-(m_refRat-1));
      Vector<VolIndex> fineVoFs = ebislBuf.refine(coarVoF,m_refRat, a_fineDatInd);
      Vector<Real> areaSum(m_nComp, 0.0);
      //loop over fine vofs at the coar VoFs cf interface.
      for(int ifinev = 0; ifinev < fineVoFs.size(); ifinev++)
        {
          const VolIndex& fineVoF = fineVoFs[ifinev];
          const IntVect& fineIV = fineVoF.gridIndex();
          int isign =  sign(a_sd);
          //the data for the fine cells is on the other side of the
          //CF interface so need to check for irregularity there
          IntVect fineIVInt = fineIV - isign*BASISV(a_dir);
          //if(vof is on coar-fine interface)
          if((cfbox.contains(fineIV)) && (ebisBoxFine.isIrregular(fineIVInt)))
            {
              Vector<FaceIndex> faces =
                ebisBoxFine.getFaces(fineVoF, a_dir, flip(a_sd));
              for(int iface = 0; iface < faces.size(); iface++)
                {
                  const FaceIndex& fineFace = faces[iface];
                  for(int ivar = a_variables.begin();
                      ivar <= a_variables.end(); ivar++)
                    {
                      Real area = ebisBoxFine.areaFrac(fineFace);
                      Real flux = a_fineFlux(fineFace, ivar);
                      //consistent with non-eb version.
                      areaSum[ivar] += newScale*area*flux;
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
}
/*******************/
/*******************/
void
EBFluxRegister::reflux(LevelData<EBCellFAB>& a_uCoar,
                       const Interval&       a_variables,
                       const Real&           a_scale)
{
  //pout() << " reflux 00 " << endl;
  CH_assert(isDefined());
  CH_assert(a_variables.begin() >= 0);
  CH_assert(a_variables.end() < a_uCoar.nComp());
  CH_assert(a_variables.end() < m_nComp);
  //pout() << " reflux 0 " << endl;
  for(int idir=0 ; idir<SpaceDim; ++idir)
    {
      reflux(a_uCoar, a_variables, a_scale, idir);
    }
}
void
EBFluxRegister::reflux(LevelData<EBCellFAB>&       a_uCoar,
                       const Interval&             a_variables,
                       const Real&                 a_scale,
                       const LevelData<EBCellFAB>& a_beta)
{
  //pout() << " reflux 00 " << endl;
  CH_assert(isDefined());
  CH_assert(a_variables.begin() >= 0);
  CH_assert(a_variables.end() < a_uCoar.nComp());
  CH_assert(a_variables.end() < m_nComp);
  //pout() << " reflux 0 " << endl;
  for(int idir=0 ; idir<SpaceDim; ++idir)
    {
      reflux(a_uCoar, a_variables, a_scale, idir,a_beta);
    }
}
/*******************/
/*******************/
void
EBFluxRegister::reflux(LevelData<EBCellFAB>&       a_uCoar,
                       const Interval&             a_variables,
                       const Real&                 a_scale,
                       const int&                  a_idir,
                       const LevelData<EBCellFAB>& a_beta)
{
  CH_TIME("EBFluxRegister::reflux");
  for(SideIterator side; side.ok(); ++side)
    {
      //pout() << " reflux 0. " << sign(side())*a_idir << endl;
      int iindex = index(a_idir, side());
      const LevelData<BaseIVFAB<Real> >& fineReg = m_regsFine[iindex];
      const LevelData<BaseIVFAB<Real> >& coarReg = m_regsCoar[iindex];
      LevelData<BaseIVFAB<Real> >& scratch =       m_scratchc[iindex];
      const LayoutData<IntVectSet>& ivsldCoar  =   m_cfivsCoar[iindex];
      const EBISLayout& ebislCoar = m_ebislBufCoar[iindex];
      // make a local copy in the coar layout space of
      // the fine register increments

      fineReg.copyTo(a_variables, scratch, a_variables);
      //pout() << " reflux 1 " << endl;

      for(DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit)
        {
          EBCellFAB& ufab = a_uCoar[dit()];
          Vector<DataIndex>&  datIndVec = m_coarIndexMap[iindex][dit()];
          for(int ivec=0; ivec < datIndVec.size(); ++ivec)
            {
              const DataIndex& datInd = datIndVec[ivec];
              const BaseIVFAB<Real>& coarRegFAB = coarReg[datInd];
              const BaseIVFAB<Real>& fineRegFAB = scratch[datInd];
              const IntVectSet& ivs  = ivsldCoar[datInd];
              const EBISBox& ebisBox = ebislCoar[datInd];
              //pout() << " reflux 2 " << endl;
              for(VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
                {
                  const VolIndex& vof = vofit();
                  Real betaPt = a_beta[dit()](vof, 1);
                  Real newScale = sign(side())*a_scale*betaPt;
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
                    }
                }
            }
        }
    }
}

/*******************/
void
EBFluxRegister::reflux(LevelData<EBCellFAB>&       a_uCoar,
                       const Interval&             a_variables,
                       const Real&                 a_scale,
                       const int&                  a_idir)
{
  CH_TIME("EBFluxRegister::reflux");
  for(SideIterator side; side.ok(); ++side)
    {
      //pout() << " reflux 0. " << sign(side())*a_idir << endl;
      int iindex = index(a_idir, side());
      const LevelData<BaseIVFAB<Real> >& fineReg = m_regsFine[iindex];
      const LevelData<BaseIVFAB<Real> >& coarReg = m_regsCoar[iindex];
      LevelData<BaseIVFAB<Real> >& scratch =       m_scratchc[iindex];
      const LayoutData<IntVectSet>& ivsldCoar  =   m_cfivsCoar[iindex];
      const EBISLayout& ebislCoar = m_ebislBufCoar[iindex];
      // make a local copy in the coar layout space of
      // the fine register increments

      fineReg.copyTo(a_variables, scratch, a_variables);
      //pout() << " reflux 1 " << endl;

      for(DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit)
        {
          EBCellFAB& ufab = a_uCoar[dit()];
          Vector<DataIndex>&  datIndVec = m_coarIndexMap[iindex][dit()];
          for(int ivec=0; ivec < datIndVec.size(); ++ivec)
            {
              const DataIndex& datInd = datIndVec[ivec];
              const BaseIVFAB<Real>& coarRegFAB = coarReg[datInd];
              const BaseIVFAB<Real>& fineRegFAB = scratch[datInd];
              const IntVectSet& ivs  = ivsldCoar[datInd];
              const EBISBox& ebisBox = ebislCoar[datInd];
              //pout() << " reflux 2 " << endl;
              for(VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
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
                    }
                }
            }
        }
    }
}
/*******************/
void
EBFluxRegister::incrementRedistRegister(EBCoarToCoarRedist& a_register,
                                        const Interval&     a_variables,
                                        const Real&         a_scale)
{
  CH_assert(isDefined());
  CH_assert(a_variables.begin() >= 0);
  CH_assert(a_variables.end() < m_nComp);

  LevelData<BaseIVFAB<Real> >& registerMassLDF = a_register.m_regsCoar;
  LayoutData<IntVectSet>& registerSets = a_register.m_setsCoar;
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

          for(DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit)
            {
              BaseIVFAB<Real>& massfab = registerMassLDF[dit()];
              const IntVectSet& ivsReg = registerSets[dit()];
              Vector<DataIndex>&  datIndVec = m_coarIndexMap[iindex][dit()];
              for(int ivec=0; ivec < datIndVec.size(); ++ivec)
                {
                  const DataIndex& datInd = datIndVec[ivec];
                  const BaseIVFAB<Real>& coarRegFAB = coarReg[datInd];
                  const BaseIVFAB<Real>& fineRegFAB = scratch[datInd];
                  const EBISBox& ebisBox = ebislCoar[datInd];
                  IntVectSet ivs  = ivsldCoar[datInd];
                  ivs &= ivsReg;
                  for(VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
                    {
                      const VolIndex& vof = vofit();
                      Real volfrac = ebisBox.volFrac(vof);
                      //multiplied by volfrac
                      //to increment a mass instead of a density
                      Real newScale = volfrac*(1.0-volfrac)*sign(side())*a_scale;
                      for(int icomp = a_variables.begin();
                          icomp <= a_variables.end(); icomp++)
                        {
                          //the non-eb version has these added and has
                          //the coar one multiplied by the negative of  the scale.
                          //so that is what we do
                          //the -= is to cope with the fact that sign(side()) gives
                          //the wrong sign
                          massfab(vof, icomp) -= newScale*(fineRegFAB(vof, icomp) +
                                                           coarRegFAB(vof, icomp));
                        }
                    }
                }
            }
        }
    }
}
/*******************/
void
EBFluxRegister::incrementRedistRegister(EBCoarToFineRedist& a_register,
                                        const Interval&     a_variables,
                                        const Real&         a_scale)
{
  CH_assert(isDefined());
  CH_assert(a_variables.begin() >= 0);
  CH_assert(a_variables.end() < m_nComp);

  LevelData<BaseIVFAB<Real> >& registerMassLDF = a_register.m_regsCoar;
  LayoutData<IntVectSet>& registerSets = a_register.m_setsCoar;
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

          for(DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit)
            {
              BaseIVFAB<Real>& massfab = registerMassLDF[dit()];
              Vector<DataIndex>&  datIndVec = m_coarIndexMap[iindex][dit()];
              for(int ivec=0; ivec < datIndVec.size(); ++ivec)
                {
                  const DataIndex& datInd = datIndVec[ivec];
                  const BaseIVFAB<Real>& coarRegFAB = coarReg[datInd];
                  const BaseIVFAB<Real>& fineRegFAB = scratch[datInd];
                  const IntVectSet& ivsReg = registerSets[dit()];
                  const EBISBox& ebisBox = ebislCoar[datInd];
                  IntVectSet ivs  = ivsldCoar[datInd];
                  ivs &= ivsReg;
                  for(VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
                    {
                      const VolIndex& vof = vofit();
                      Real volfrac = ebisBox.volFrac(vof);
                      //multiplied by volfrac
                      //to increment a mass instead of a density
                      Real newScale = volfrac*(1.0-volfrac)*sign(side())*a_scale;
                      for(int icomp = a_variables.begin();
                          icomp <= a_variables.end(); icomp++)
                        {
                          //the non-eb version has these added and has
                          //the coar one multiplied by the negative of  the scale.
                          //so that is what we do
                          //the -= is to cope with the fact that sign(side()) gives
                          //the wrong sign
                          massfab(vof, icomp) -= newScale*(fineRegFAB(vof, icomp) +
                                                           coarRegFAB(vof, icomp));
                        }
                    }
                }
            }
        }
    }
}
/*******************/
/*******************/
void
EBFluxRegister::incrementDensityArray(LevelData<EBCellFAB>& a_coarDense,
                                      const Interval&       a_variables,
                                      const Real&           a_scale)
{
  CH_assert(isDefined());
  CH_assert(a_variables.begin() >= 0);
  CH_assert(a_variables.end() < a_coarDense.nComp());
  CH_assert(a_variables.end() < m_nComp);

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

          for(DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit)
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
/*******************/
void
EBFluxRegister::dumpCoar(const int& a_idir,
                         const Side::LoHiSide& a_sd)
{
  pout() << endl << "coarse buffer for idir = " << a_idir
         << " side = " << sign(a_sd) << endl;

  int iindex = index(a_idir, a_sd);
  const DisjointBoxLayout& grids =      m_bufGridsCoar[iindex];
  const LevelData<BaseIVFAB<Real> >& data = m_regsCoar[iindex];
  const EBISLayout& ebisl =             m_ebislBufCoar[iindex];

  for(DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
      pout() << "buffer grid = " << grids.get(dit()) << endl;
      dumpBIVF(data[dit()], ebisl[dit()]);
    }
}
/*******************/
void
EBFluxRegister::dumpFine(const int& a_idir,
                         const Side::LoHiSide& a_sd)
{
  pout() << endl << "fine buffer for idir = " << a_idir
         << " side = " << sign(a_sd) << endl;

  int iindex = index(a_idir, a_sd);
  const DisjointBoxLayout& grids =      m_bufGridsFine[iindex];
  const LevelData<BaseIVFAB<Real> >& data = m_regsFine[iindex];
  const EBISLayout& ebisl =             m_ebislBufFine[iindex];

  for(DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
      pout() << "buffer grid = " << grids.get(dit()) << endl;
      dumpBIVF(data[dit()], ebisl[dit()]);
    }
}
/*******************/
void
EBFluxRegister::dumpBIVF(const BaseIVFAB<Real>& a_reg,
                         const EBISBox& a_ebisBox)
{
  for(VoFIterator vofit(a_reg.getIVS(), a_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      pout() << vofit() << "   ";
      for(int ivar = 0; ivar < a_reg.nComp(); ivar++)
        {
          pout() << a_reg(vofit(), ivar) << "   ";
        }
      pout() << endl;
    }
}
/*******************/
//rz stuff starts here
/*******************/
void
EBFluxRegister::incrementCoarseRegulRZ(const EBFaceFAB& a_coarFlux,
                                       const Real&      a_scale,
                                       const DataIndex& a_coarDatInd,
                                       const Interval&  a_variables,
                                       const int&       a_dir,
                                       const Real&      a_dx)
{
  CH_assert(SpaceDim==2);
  CH_assert(isDefined());
  CH_assert(a_dir >= 0);
  CH_assert(a_dir < SpaceDim);
  CH_assert(a_coarFlux.direction() == a_dir);
  CH_assert(a_variables.begin() >= 0);
  CH_assert(a_variables.end() < m_nComp);
  CH_assert(a_variables.end() < a_coarFlux.nComp());

  Real dxCoar = a_dx;
  //iterate over sides of the fine grid.
  for(SideIterator side; side.ok(); ++side)
    {
      int iindex = index(a_dir, side());

      LevelData<BaseIVFAB<Real> >& registerLD = m_regsCoar[iindex];
      const EBISLayout& ebislBuf = m_ebislBufCoar[iindex];
      const LayoutData<IntVectSet>& ldivs = m_cfivsCoar[iindex];

      const Vector<DataIndex>& vecIndex =
        m_coarIndexMap[iindex][a_coarDatInd];

      for(int iind = 0; iind < vecIndex.size(); iind++)
        {
          const DataIndex& regDatInd = vecIndex[iind];
          BaseIVFAB<Real>& coarReg = registerLD[regDatInd];
          const EBISBox& ebisBox = ebislBuf[regDatInd];
          const IntVectSet& cfivs = ldivs[regDatInd];
          //increment coar buffer with the area-weighted sum
          //of the fluxes on the faces.
          //buf += sum(areaFrac*flux)
          for(VoFIterator vofit(cfivs, ebisBox.getEBGraph());
              vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();

              Real cellVol, kvol;
              EBArith::getKVolRZ(kvol, cellVol, ebisBox, dxCoar, vof);
              Real volFrac = ebisBox.volFrac(vof);
              Real volFactor = dxCoar*dxCoar*volFrac/(kvol*cellVol);

              if(ebisBox.isRegular(vof.gridIndex()))
                {
                  for(int ivar = a_variables.begin();
                      ivar <= a_variables.end(); ivar++)
                    {
                      //flip because the buffer is on the
                      //other side of the CF interface
                      Vector<FaceIndex> faces =
                        ebisBox.getFaces(vof, a_dir, flip(side()));
                      Real areaSum = 0.0;
                      for(int iface = 0; iface < faces.size(); iface++)
                        {
                          const FaceIndex& face = faces[iface];
                          Real area = ebisBox.areaFrac(face);
                          Real flux = a_coarFlux(face, ivar);
                          areaSum += volFactor*area*flux;
                        }
                      //consistent with non-eb version
                      areaSum *= -a_scale;
                      coarReg(vof, ivar) += areaSum;
                    }
                } //end if(isRegular())
            } //end iteration over vofs in buffer
        } //end loop over buffers
    } //end iteration over sides
}
/*************************/
void
EBFluxRegister::incrementCoarseIrregulRZ(const BaseIFFAB<Real>& a_coarFlux,
                                         const Real&            a_scale,
                                         const DataIndex&       a_coarDatInd,
                                         const Interval&        a_variables,
                                         const int&             a_dir,
                                         const Real&            a_dx)
{
  CH_assert(SpaceDim == 2);
  CH_assert(isDefined());
  CH_assert(a_dir >= 0);
  CH_assert(a_dir < SpaceDim);
  CH_assert(a_coarFlux.direction() == a_dir);
  CH_assert(a_variables.begin() >= 0);
  CH_assert(a_variables.end() < m_nComp);
  CH_assert(a_variables.end() < a_coarFlux.nComp());

  Real dxCoar = a_dx;
  //iterate over sides of the fine grid.
  for(SideIterator side; side.ok(); ++side)
    {
      int iindex = index(a_dir, side());

      LevelData<BaseIVFAB<Real> >& registerLD = m_regsCoar[iindex];
      const EBISLayout& ebislBuf = m_ebislBufCoar[iindex];
      const LayoutData<IntVectSet>& ldivs = m_cfivsCoar[iindex];

      const Vector<DataIndex>& vecIndex =
        m_coarIndexMap[iindex][a_coarDatInd];

      for(int iind = 0; iind < vecIndex.size(); iind++)
        {
          const DataIndex& regDatInd = vecIndex[iind];
          BaseIVFAB<Real>& coarReg = registerLD[regDatInd];
          const EBISBox& ebisBox = ebislBuf[regDatInd];
          const IntVectSet& cfivs = ldivs[regDatInd];
          //increment coar buffer with the area-weighted sum
          //of the fluxes on the faces.
          //buf += sum(areaFrac*flux)
          for(VoFIterator vofit(cfivs, ebisBox.getEBGraph());
              vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              Real cellVol, kvol;
              EBArith::getKVolRZ(kvol, cellVol, ebisBox, dxCoar, vof);
              Real volFrac = ebisBox.volFrac(vof);
              Real volFactor = dxCoar*dxCoar*volFrac/(kvol*cellVol);

              if(ebisBox.isIrregular(vof.gridIndex()))
                {
                  for(int ivar = a_variables.begin();
                      ivar <= a_variables.end(); ivar++)
                    {
                      //flip because the buffer is on the
                      //other side of the CF interface
                      Vector<FaceIndex> faces =
                        ebisBox.getFaces(vof, a_dir, flip(side()));
                      Real areaSum = 0.0;
                      for(int iface = 0; iface < faces.size(); iface++)
                        {
                          const FaceIndex& face = faces[iface];
                          Real area = ebisBox.areaFrac(face);
                          Real flux = a_coarFlux(face, ivar);
                          areaSum += volFactor*area*flux;
                        }
                      //consistent with non-eb version
                      areaSum *= -a_scale;
                      coarReg(vof, ivar) += areaSum;
                    }
                } //end if(isIrregularVoF)
            } //end iteration over vofs in buffer
        } //end loop over buffers
    } //end iteration over sides
} //end incrementCoar
/*****/
void
EBFluxRegister::incrementFineIrregulRZ(const BaseIFFAB<Real>& a_fineFlux,
                                       const Real&            a_scale,
                                       const DataIndex&       a_fineDatInd,
                                       const Interval&        a_variables,
                                       const int&             a_dir,
                                       const Side::LoHiSide&  a_sd,
                                       const Real&            a_dx)
{
  CH_assert(SpaceDim==2);
  CH_assert(isDefined());
  CH_assert(a_fineFlux.direction() == a_dir);
  CH_assert(a_variables.begin() >= 0);
  CH_assert(a_variables.end() < a_fineFlux.nComp());
  CH_assert(a_variables.end() < m_nComp);
  CH_assert(a_dir >= 0);
  CH_assert(a_dir < SpaceDim);
  CH_assert((a_sd == Side::Lo)||(a_sd == Side::Hi));

  //nrefdmo is the number of fine faces per coar face
  //this is intrinsically dimension-dependent
  //yes, i could do some lame loop here, but that would
  //be more confusing than this
#if (CH_SPACEDIM == 2)
  Real nrefdmo = m_refRat;
#elif (CH_SPACEDIM == 3)
  Real nrefdmo = m_refRat*m_refRat;
#else
  bogus spacedim;
#endif
  //  Real newScale = a_scale/nrefdmo;
  Real dxCoar = a_dx*m_refRat;
  int iindex = index(a_dir, a_sd);
  BaseIVFAB<Real>&  regFab = m_regsFine[iindex][a_fineDatInd];
  const IntVectSet& ivsReg = m_cfivsFine[iindex][a_fineDatInd];
  const EBISLayout& ebislBuf = m_ebislBufFine[iindex];
  const EBISBox& ebisBoxBuf = ebislBuf[a_fineDatInd];

  const EBISBox& ebisBoxFine = m_ebislFine[a_fineDatInd];
  for(VoFIterator vofit(ivsReg, ebisBoxBuf.getEBGraph()); vofit.ok(); ++vofit)
    {
      //remember the registers live at the coar level
      const VolIndex& coarVoF = vofit();
      Real cellVol, kvol;
      //dx being sent in is already dxcoar
      EBArith::getKVolRZ(kvol, cellVol, ebisBoxBuf, dxCoar, coarVoF);
      Real volFrac = ebisBoxBuf.volFrac(coarVoF);
      Real volFactor = dxCoar*dxCoar*volFrac/(kvol*cellVol*nrefdmo);

      //box to weed out which cells are on the  cf interface
      const IntVect& coarIndex = coarVoF.gridIndex();
      Box cfbox(coarIndex, coarIndex);
      cfbox.refine(m_refRat);
      if(a_sd == Side::Lo) cfbox.growLo(a_dir,-(m_refRat-1));
      else                 cfbox.growHi(a_dir,-(m_refRat-1));
      Vector<VolIndex> fineVoFs = ebislBuf.refine(coarVoF,m_refRat, a_fineDatInd);
      Vector<Real> areaSum(m_nComp, 0.0);
      //loop over fine vofs at the coar VoFs cf interface.
      for(int ifinev = 0; ifinev < fineVoFs.size(); ifinev++)
        {
          const VolIndex& fineVoF = fineVoFs[ifinev];
          const IntVect& fineIV = fineVoF.gridIndex();
          int isign =  sign(a_sd);
          //the data for the fine cells is on the other side of the
          //CF interface so need to check for irregularity there
          IntVect fineIVInt = fineIV - isign*BASISV(a_dir);
          //if(vof is on coar-fine interface)
          if((cfbox.contains(fineIV)) && (ebisBoxFine.isIrregular(fineIVInt)))
            {
              Vector<FaceIndex> faces =
                ebisBoxFine.getFaces(fineVoF, a_dir, flip(a_sd));
              for(int iface = 0; iface < faces.size(); iface++)
                {
                  const FaceIndex& fineFace = faces[iface];
                  for(int ivar = a_variables.begin();
                      ivar <= a_variables.end(); ivar++)
                    {
                      Real area = ebisBoxFine.areaFrac(fineFace);
                      Real flux = a_fineFlux(fineFace, ivar);
                      //consistent with non-eb version.
                      areaSum[ivar] += a_scale*volFactor*area*flux;
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
}
/*****/
void
EBFluxRegister::incrementFineRegulRZ(const EBFaceFAB&      a_fineFlux,
                                     const Real&           a_scale,
                                     const DataIndex&      a_fineDatInd,
                                     const Interval&       a_variables,
                                     const int&            a_dir,
                                     const Side::LoHiSide& a_sd,
                                     const Real&           a_dx)
{
  CH_assert(SpaceDim==2);
  CH_assert(isDefined());
  CH_assert(a_fineFlux.direction() == a_dir);
  CH_assert(a_variables.begin() >= 0);
  CH_assert(a_variables.end() < a_fineFlux.nComp());
  CH_assert(a_variables.end() < m_nComp);
  CH_assert(a_dir >= 0);
  CH_assert(a_dir < SpaceDim);
  CH_assert((a_sd == Side::Lo)||(a_sd == Side::Hi));

  //nrefdmo is the number of fine faces per coar face
  //this is intrinsically dimension-dependent
  //yes, i could do some lame loop here, but that would
  //be more confusing than this
#if (CH_SPACEDIM == 2)
  Real nrefdmo = m_refRat;
#elif (CH_SPACEDIM == 3)
  Real nrefdmo = m_refRat*m_refRat;
#else
  bogus spacedim;
#endif
  //  Real newScale = a_scale/nrefdmo;

  Real dxCoar = a_dx*m_refRat;
  int iindex = index(a_dir, a_sd);
  BaseIVFAB<Real>&  regFab = m_regsFine[iindex][a_fineDatInd];
  const IntVectSet& ivsReg = m_cfivsFine[iindex][a_fineDatInd];
  const EBISLayout& ebislBuf = m_ebislBufFine[iindex];
  const EBISBox& ebisBoxBuf = ebislBuf[a_fineDatInd];

  const EBISBox& ebisBoxFine = m_ebislFine[a_fineDatInd];
  for(VoFIterator vofit(ivsReg, ebisBoxBuf.getEBGraph()); vofit.ok(); ++vofit)
    {
      //remember the registers live at the coar level
      const VolIndex& coarVoF = vofit();
      Real cellVol, kvol;
      //dx being sent in is already dxcoar
      EBArith::getKVolRZ(kvol, cellVol, ebisBoxBuf, dxCoar, coarVoF);
      Real volFrac = ebisBoxBuf.volFrac(coarVoF);
      Real volFactor = dxCoar*dxCoar*volFrac/(kvol*cellVol*nrefdmo);

      //box to weed out which cells are on the  cf interface
      const IntVect& coarIndex = coarVoF.gridIndex();
      Box cfbox(coarIndex, coarIndex);
      cfbox.refine(m_refRat);
      if(a_sd == Side::Lo) cfbox.growLo(a_dir,-(m_refRat-1));
      else                 cfbox.growHi(a_dir,-(m_refRat-1));
      Vector<VolIndex> fineVoFs = ebislBuf.refine(coarVoF,m_refRat, a_fineDatInd);
      Vector<Real> areaSum(m_nComp, 0.0);
      //loop over fine vofs at the coar VoFs cf interface.
      for(int ifinev = 0; ifinev < fineVoFs.size(); ifinev++)
        {
          const VolIndex& fineVoF = fineVoFs[ifinev];
          const IntVect& fineIV = fineVoF.gridIndex();
          int isign = sign(a_sd);
          //the data for the fine cells is on the other side of the
          //CF interface so need to check for irregularity there
          IntVect fineIVInt = fineIV - isign*BASISV(a_dir);
          if((ebisBoxFine.isRegular(fineIVInt))&&(cfbox.contains(fineIV)))
            {
              Vector<FaceIndex> faces =
                ebisBoxFine.getFaces(fineVoF, a_dir, flip(a_sd));
              for(int iface = 0; iface < faces.size(); iface++)
                {
                  const FaceIndex& fineFace = faces[iface];
                  for(int ivar = a_variables.begin();
                      ivar <= a_variables.end(); ivar++)
                    {
                      Real area = ebisBoxFine.areaFrac(fineFace);
                      Real flux = a_fineFlux(fineFace, ivar);
                      //consistent with non-eb version.
                      areaSum[ivar] += a_scale*volFactor*area*flux;
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
}
/*************************/
#include "NamespaceFooter.H"
