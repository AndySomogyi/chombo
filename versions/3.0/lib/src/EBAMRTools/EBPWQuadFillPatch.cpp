#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EBPWQuadFillPatch.H"
#include "EBInterpolateF_F.H"
#include "VoFIterator.H"
#include "DisjointBoxLayout.H"
#include "EBCellFactory.H"
#include "LayoutIterator.H"
#include "EBIndexSpace.H"
#include "NamespaceHeader.H"

/************************************/
/************************************/
void
EBPWQuadFillPatch::setDefaultValues()
{
  m_isDefined = false;
  m_refRat = -1;
  m_nComp = -1;
  m_radius = -1;
}
/************************************/
/************************************/
EBPWQuadFillPatch::EBPWQuadFillPatch()
{
  setDefaultValues();
}
/************************************/
/************************************/
EBPWQuadFillPatch::~EBPWQuadFillPatch()
{
}
/************************************/
/************************************/
EBPWQuadFillPatch::EBPWQuadFillPatch(const DisjointBoxLayout& a_dblFine,
                                     const DisjointBoxLayout& a_dblCoar,
                                     const EBISLayout& a_ebislFine,
                                     const EBISLayout& a_ebislCoar,
                                     const ProblemDomain& a_domainCoar,
                                     const int& a_nref,
                                     const int& a_nvar,
                                     const int& a_radius)
{
  setDefaultValues();

  define(a_dblFine, a_dblCoar, a_ebislFine, a_ebislCoar,
         a_domainCoar, a_nref, a_nvar, a_radius);
}
/************************************/
/************************************/
void
EBPWQuadFillPatch::define(const DisjointBoxLayout& a_dblFine,
                          const DisjointBoxLayout& a_dblCoar,
                          const EBISLayout& a_ebislFine,
                          const EBISLayout& a_ebislCoar,
                          const ProblemDomain& a_domainCoar,
                          const int& a_nref,
                          const int& a_nvar,
                          const int& a_radius,
                          const EBIndexSpace* const ebisPtr)
{
  CH_assert(a_nref > 0);
  CH_assert(a_nvar > 0);
  //first cell is interpolated, rest are extrapolated.
  //if this is 0, there is not even interpolation.
  CH_assert(a_radius > 0);

  m_radius = a_radius;
  m_isDefined = true;
  m_refRat = a_nref;
  m_nComp = a_nvar;

  m_coarDomain = a_domainCoar;
  m_fineEBISL = a_ebislFine;
  if(m_fineEBISL.getMaxCoarseningRatio() < m_refRat)
    {
      m_fineEBISL.setMaxCoarseningRatio(m_refRat);
    }
  m_fineGrids = a_dblFine;
  m_coarGrids = a_dblCoar;
  //needs to be extra because the stencil of
  //the derivs has an extra coarse cell
  const int coarse_slope_radius =
    (m_radius + m_refRat - 1) / m_refRat;
  const int coarse_ghost_radius = coarse_slope_radius + 1;

  m_coarGhostRad = coarse_ghost_radius;
  IntVect ghostiv = m_coarGhostRad*IntVect::Unit;

  CH_assert(ebisPtr->isDefined());
  m_coarsenedFineGrids = DisjointBoxLayout();
  coarsen(m_coarsenedFineGrids, a_dblFine, m_refRat);
  ebisPtr->fillEBISLayout(m_coarsenedFineEBISL,
                          m_coarsenedFineGrids,
                          a_domainCoar.domainBox(), m_coarGhostRad);
  //  m_coarsenedFineEBISL.setMaxCoarseningRatio(m_refRat);
  EBCellFactory ebcellfact(m_coarsenedFineEBISL);
  m_coarOnFDataOld.define(m_coarsenedFineGrids, m_nComp,
                          ghostiv, ebcellfact);
  m_coarOnFDataNew.define(m_coarsenedFineGrids, m_nComp,
                          ghostiv, ebcellfact);
  m_irregRegionsFine.define(m_fineGrids);
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      m_coarCeInterp[idir].define(m_coarsenedFineGrids);
      m_coarLoInterp[idir].define(m_coarsenedFineGrids);
      m_coarHiInterp[idir].define(m_coarsenedFineGrids);
    }
  makeStencils();

  //make a corner copier
  ProblemDomain domainFine = m_coarDomain;
  domainFine.refine(m_refRat);
  m_cornerCopier.define(m_fineGrids, m_fineGrids, domainFine, IntVect::Unit);
}
/************************************/
/************************************/
void
EBPWQuadFillPatch::makeStencils()
{
  ProblemDomain domFine = refine(m_coarDomain, m_refRat);
  for(DataIterator dit = m_coarsenedFineGrids.dataIterator();
      dit.ok(); ++dit)
    {
      //irreg region is grown box - all fine grids.
      IntVectSet& localIrreg   = m_irregRegionsFine[dit()];
      Box bigBox= m_fineGrids.get(dit());
      bigBox.grow(m_radius);
      bigBox &= domFine;
      localIrreg = IntVectSet(bigBox);
      for(LayoutIterator lit = m_fineGrids.layoutIterator();
          lit.ok(); ++lit)
        {
          localIrreg -= m_fineGrids.get(lit());
        }

    }
  //now create the coarse lo hi and center intvectsets--
  //these sets determine where you have to do one sided differences
  //because of the configuration of the coarse layout.
  for(DataIterator dit = m_coarsenedFineGrids.dataIterator();
      dit.ok(); ++dit)
    {
      const Box& fineBox = m_fineGrids.get(dit());
      Box coarsenedFineBox = coarsen(grow(fineBox, m_radius), m_refRat);
      coarsenedFineBox &= m_coarDomain;

      IntVectSet coarsenedFineInterp(coarsenedFineBox);

      // Iterate over boxes in coarsened fine domain, and subtract off
      // from the set of coarse cells from which the fine ghost cells
      // will be interpolated.

      for(LayoutIterator litCF = m_coarsenedFineGrids.layoutIterator();
          litCF.ok(); ++litCF)
        {
          const Box& otherCoarsenedBox = m_coarsenedFineGrids.get(litCF());
          coarsenedFineInterp -= otherCoarsenedBox;
        }

      // Now that we have the coarsened cells required for interpolation,
      // construct IntvectSets specifying the one-sided and centered
      // stencil locations.
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          IntVectSet& coarseCeInterp = m_coarCeInterp[idir][dit()];
          IntVectSet& coarseLoInterp = m_coarLoInterp[idir][dit()];
          IntVectSet& coarseHiInterp = m_coarHiInterp[idir][dit()];

          coarseCeInterp = coarsenedFineInterp;

          coarseLoInterp = coarseCeInterp;
          coarseLoInterp.shift(BASISV(idir));

          coarseHiInterp = coarseCeInterp;
          coarseHiInterp.shift(-BASISV(idir));

          // We iterate over the coarse grids and subtract them off of the
          // one-sided stencils.
          for(LayoutIterator litC = m_coarGrids.layoutIterator(); litC.ok(); ++litC)
            {
              const Box& bx = m_coarGrids.get(litC());
              coarseLoInterp -= bx;
              coarseHiInterp -= bx;
            }

          coarseLoInterp.shift(-BASISV(idir));
          coarseHiInterp.shift(BASISV(idir));

          coarseCeInterp -= coarseLoInterp;
          coarseCeInterp -= coarseHiInterp;
        }
    }

  /////////////////////////
    for(int derivDir = 0; derivDir < SpaceDim; derivDir++)
      {
        m_hiStencils[derivDir].define(m_coarsenedFineGrids);
        m_loStencils[derivDir].define(m_coarsenedFineGrids);

        LayoutData<BaseIVFAB<VoFStencil > >&
          loStenLDF= m_loStencils[derivDir];
        LayoutData<BaseIVFAB<VoFStencil > >&
          hiStenLDF= m_hiStencils[derivDir];

        for(DataIterator dit = m_coarsenedFineGrids.dataIterator();
            dit.ok(); ++dit)
          {
            const EBISBox& ebisBox  =  m_coarsenedFineEBISL[dit()];
            const IntVectSet& fineIrreg  = m_irregRegionsFine[dit()];
            IntVectSet coarIrreg = coarsen(fineIrreg, m_refRat);
            BaseIVFAB<VoFStencil >& hiStenFAB = hiStenLDF[dit()];
            BaseIVFAB<VoFStencil >& loStenFAB = loStenLDF[dit()];

            hiStenFAB.define(coarIrreg, ebisBox.getEBGraph(), 1);
            loStenFAB.define(coarIrreg, ebisBox.getEBGraph(), 1);

            for(VoFIterator vofit(coarIrreg, ebisBox.getEBGraph());
                vofit.ok(); ++vofit)
              {
                const VolIndex& vof = vofit();
                Vector<FaceIndex> loFaces=
                  ebisBox.getFaces(vof, derivDir, Side::Lo);
                Vector<FaceIndex> hiFaces=
                  ebisBox.getFaces(vof, derivDir, Side::Hi);
                Vector<VolIndex> loVoFs;
                Vector<VolIndex> hiVoFs;
                for(int iface = 0; iface < loFaces.size(); iface++)
                  {
                    //use one-sided diff at edge of domain boundary
                    if(!loFaces[iface].isBoundary())
                      loVoFs.push_back(loFaces[iface].getVoF(Side::Lo));
                  }
                for(int iface = 0; iface < hiFaces.size(); iface++)
                  {
                    //use one-sided diff at edge of domain boundary
                    if(!hiFaces[iface].isBoundary())
                      hiVoFs.push_back(hiFaces[iface].getVoF(Side::Hi));
                  }
                if(loVoFs.size() > 0)
                  {
                    // have vofs on the low side.
                    //one sided diff with low side
                    Vector<VolIndex> stenVoFs;
                    Vector<Real>     stenWeig;

                    Real rfacesLo = Real(loVoFs.size());
                    Vector<Real> loWeig(loVoFs.size(), -1.0/rfacesLo);
                    stenVoFs = loVoFs;
                    stenWeig = loWeig;
                    stenVoFs.push_back(vof);
                    stenWeig.push_back(1.0);

                    //finally put the stencil into its container
                    VoFStencil& vofsten = loStenFAB(vof, 0);
                    vofsten.clear();
                    CH_assert(stenVoFs.size() == stenWeig.size());
                    for(int isten = 0; isten < stenVoFs.size(); isten++)
                      {
                        vofsten.add(stenVoFs[isten], stenWeig[isten]);
                      }
                  }
                if(hiVoFs.size() > 0)
                  {
                    // have vofs on the high side.
                    //one sided diff with high side
                    Vector<VolIndex> stenVoFs;
                    Vector<Real>     stenWeig;

                    Real rfacesHi = Real(hiVoFs.size());
                    Vector<Real> hiWeig(hiVoFs.size(),  1.0/rfacesHi);
                    stenVoFs = hiVoFs;
                    stenWeig = hiWeig;
                    stenVoFs.push_back(vof);
                    stenWeig.push_back(-1.0);

                    //finally put the stencil into its container
                    VoFStencil& vofsten = hiStenFAB(vof, 0);
                    vofsten.clear();
                    CH_assert(stenVoFs.size() == stenWeig.size());
                    for(int isten = 0; isten < stenVoFs.size(); isten++)
                      {
                        vofsten.add(stenVoFs[isten], stenWeig[isten]);
                      }
                  }
              } //end loop over vofs
          } //end loop over grids
      } //end loop over derivative directions
}
/************************************/
/************************************/
bool
EBPWQuadFillPatch::isDefined() const
{
  return m_isDefined;
}
/************************************/
/************************************/
void
EBPWQuadFillPatch::interpolate(LevelData<EBCellFAB>& a_fineData,
                               const LevelData<EBCellFAB>& a_coarDataOld,
                               const LevelData<EBCellFAB>& a_coarDataNew,
                               const Real& a_coarTimeOld,
                               const Real& a_coarTimeNew,
                               const Real& a_fineTime,
                               const Interval& a_variables)
{
  CH_assert(isDefined());
  CH_assert(a_coarTimeNew >= a_coarTimeOld);
  CH_assert(a_coarTimeNew >= a_fineTime);
  CH_assert(a_fineTime >= a_coarTimeOld);

  //level data copy fills coarse ghost cells as well as interior
  a_coarDataOld.copyTo(a_variables, m_coarOnFDataOld, a_variables);
  a_coarDataNew.copyTo(a_variables, m_coarOnFDataNew, a_variables);

  m_coarOnFDataOld.exchange(a_variables);
  m_coarOnFDataNew.exchange(a_variables);
#if CH_SPACEDIM==2
  //  Box bdebug(IntVect(16, 16),IntVect(23, 23));
#endif

  for(DataIterator fineit = m_coarsenedFineGrids.dataIterator();
      fineit.ok(); ++fineit)
    {
      Box fineBox = m_fineGrids.get(fineit());
      Box coarsenedFineBox = m_coarsenedFineGrids.get(fineit());
      EBCellFAB&  fineData = a_fineData[fineit()];
      const EBCellFAB&  coarDataOld = m_coarOnFDataOld[fineit()];
      const EBCellFAB&  coarDataNew = m_coarOnFDataNew[fineit()];
      // interpolateFAB interpolates from an entire coarse grid onto an
      // entire fine grids coarse-fine regions.
      interpolateFAB(fineData,
                     coarDataOld,
                     coarDataNew,
                     a_coarTimeOld,
                     a_coarTimeNew,
                     a_fineTime, fineit(),
                     a_variables);
    }
  //exchange
  a_fineData.exchange(a_variables);
  //exchange the corners
  a_fineData.exchange(a_variables, m_cornerCopier);
}
/************************************/
/************************************/
void
EBPWQuadFillPatch::interpolateFAB(EBCellFAB& a_fine,
                                  const EBCellFAB& a_coarDataOld,
                                  const EBCellFAB& a_coarDataNew,
                                  const Real& a_coarTimeOld,
                                  const Real& a_coarTimeNew,
                                  const Real& a_fineTime,
                                  const DataIndex& a_datInd,
                                  const Interval& a_variables) const
{
  //
  //interpolation factor
  Real factor = 0.0;
  if((a_coarTimeNew - a_coarTimeOld) > 1.0e-8)
    factor = (a_fineTime - a_coarTimeOld)/(a_coarTimeNew - a_coarTimeOld);

  Box refBox(IntVect::Zero,
             (m_refRat-1)*IntVect::Unit);

  const EBISBox& fineEBISBox = m_fineEBISL[a_datInd];
  const IntVectSet& ivsIrreg = m_irregRegionsFine[a_datInd];
  for(VoFIterator vofit(ivsIrreg, fineEBISBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& fineVoF = vofit();
      VolIndex coarVoF = m_fineEBISL.coarsen(fineVoF, m_refRat, a_datInd);
      const IntVect& coarIV = coarVoF.gridIndex();
      const IntVect& fineIV = fineVoF.gridIndex();

      Box stenBox = Box(coarIV,coarIV);
      stenBox.grow(1);
      bool doQuad = false;
      if(m_coarDomain.contains(stenBox))
        {
          if(m_coarsenedFineEBISL[a_datInd].isRegular(stenBox))
            {
              if(Abs(factor) < 1.e-8)
                {
                  doQuad = true;
                }
              else
                {
                  //                  MayDay::Warning("Dropped order in EBPWQuadFillPatch because of time difference");
                }
            }
        }

      const BaseFab<Real>& coarRegFAB = a_coarDataOld.getSingleValuedFAB();
      //make a temp fine data holder (because in fortran we do all the fine vofs at once)
      Box fineTempBox = Box(coarIV,coarIV);
      fineTempBox.refine(m_refRat);
      BaseFab<Real> fineRegTemp(fineTempBox,a_fine.nComp());


      for(int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++)
        {
          Real interpVal;
          if(doQuad)
            {//do bi/tri-quadratic
              IntVect ishift = IntVect::Zero;
              Box ivBox = Box(coarIV,coarIV);
              FORT_EBINTERPQUADSHIFT(CHF_FRA1(fineRegTemp, ivar),
                                     CHF_CONST_FRA1(coarRegFAB,ivar),
                                     CHF_BOX(ivBox),
                                     CHF_CONST_INT(m_refRat),
                                     CHF_CONST_INTVECT(ishift),
                                     CHF_BOX(refBox));
              //now copy the fine value over
              interpVal = fineRegTemp(fineIV,ivar);
            }//end quadratic
          else
            {//do linear
              //first set the fine stuff to the coar stuff
              Real fineValOld = a_coarDataOld(coarVoF, ivar);
              Real fineValNew = a_coarDataNew(coarVoF, ivar);

              //now add first derivative terms
              for(int derivDir = 0; derivDir  < SpaceDim; derivDir++)
                {
                  const LayoutData<BaseIVFAB<VoFStencil > >&
                    loDirSten = m_loStencils[derivDir];
                  const LayoutData<BaseIVFAB<VoFStencil > >&
                    hiDirSten = m_hiStencils[derivDir];
                  const LayoutData<IntVectSet>& loCoarInterp = m_coarLoInterp[derivDir];
                  const LayoutData<IntVectSet>& hiCoarInterp = m_coarHiInterp[derivDir];
                  const LayoutData<IntVectSet>& ceCoarInterp = m_coarCeInterp[derivDir];

                  const BaseIVFAB<VoFStencil >& loStenBF = loDirSten[a_datInd];
                  const BaseIVFAB<VoFStencil >& hiStenBF = hiDirSten[a_datInd];

                  const IntVectSet& loInterpSet = loCoarInterp[a_datInd];
                  const IntVectSet& hiInterpSet = hiCoarInterp[a_datInd];
                  const IntVectSet& ceInterpSet = ceCoarInterp[a_datInd];

                  Real deltaOld = computeDMinMod(loStenBF, hiStenBF, coarVoF,
                                                 a_coarDataOld, ivar,
                                                 loInterpSet, hiInterpSet, ceInterpSet);
                  Real deltaNew = computeDMinMod(loStenBF, hiStenBF, coarVoF,
                                                 a_coarDataNew, ivar,
                                                 loInterpSet, hiInterpSet, ceInterpSet);

                  //these are the locations in space of the data
                  //the factor of 1/nref is because we normalized
                  //for dxcoar == 1
                  Real coarLoc = Real(coarIV[derivDir]) + 0.5;
                  Real fineLoc  = (Real(fineIV[derivDir]) + 0.5)/Real(m_refRat);

                  fineValOld += deltaOld*(fineLoc-coarLoc);
                  fineValNew += deltaNew*(fineLoc-coarLoc);
                  //add in contribution of deriv*dist
                } //end loop over dirs
              //finally interpolate in time
              interpVal =  fineValOld + factor*(fineValNew-fineValOld);
            }//end linear

          //the first vof gets the interpolated value.
          a_fine(fineVoF,ivar) = interpVal;
        }//end loop over variables
    }//end loop over fine vofs at the coarse/fine interface
}
/************************************/
/************************************/
bool
EBPWQuadFillPatch::extractFromSten(Real& a_value,
                                   const VoFStencil& a_vofsten,
                                   const EBCellFAB& a_coarData,
                                   const int& a_ivar) const
{
  Real retval = 0.0;
  for(int isten = 0; isten < a_vofsten.size(); isten++)
    {
      const VolIndex&  stenvof = a_vofsten.vof(isten);
      const Real& stenwgt = a_vofsten.weight(isten);
      retval += stenwgt*a_coarData(stenvof, a_ivar);
    }
  a_value = retval;
  return(a_vofsten.size() > 0);
}
/************************************/
/************************************/
Real
EBPWQuadFillPatch::computeDMinMod(const BaseIVFAB<VoFStencil>& a_loStenBF,
                                  const BaseIVFAB<VoFStencil>& a_hiStenBF,
                                  const VolIndex& a_coarVoF,
                                  const EBCellFAB& a_coarData,
                                  const int& a_ivar,
                                  const IntVectSet& loInterpSet,
                                  const IntVectSet& hiInterpSet,
                                  const IntVectSet& ceInterpSet) const
{
  //compute delta derivDir direction (normalized s.t. dxcoar == 1)
  CH_assert(isDefined());
  const IntVect& ivCoar = a_coarVoF.gridIndex();
  Real deltaminmod;
  if(ceInterpSet.contains(ivCoar))
    {
      //centered diffs allowed on coar union of rectangles
      //but might be disallowed because of EB
      Real deltahi = 0.0;
      bool hasHiSten = extractFromSten(deltahi, a_hiStenBF(a_coarVoF,0),
                                       a_coarData, a_ivar);
      Real deltalo = 0.0;
      bool hasLoSten = extractFromSten(deltalo, a_loStenBF(a_coarVoF,0),
                                       a_coarData, a_ivar);
      if(hasHiSten && hasLoSten)
        {
          //have deltas in both directions, do minmod thing.
          //should do this most of the time.
          Real mono = deltahi*deltalo;
          if(mono > 0.0)
            {
              Real rsign = 1.0;
              if((deltahi + deltalo) < 0.0)
                rsign = -1.0;
              deltaminmod = rsign*Min(Abs(deltalo), Abs(deltahi));
            }
          else
            {
              deltaminmod = 0.0;
            }
        }
      else if(hasHiSten)
        {
          deltaminmod = deltahi;
        }
      else if(hasLoSten)
        {
          deltaminmod = deltalo;
        }
      else
        {
          //no derivs exist because of EB
          deltaminmod = 0.0;
        }
    }//end coar union of rectangles allows centered diff
  else if(loInterpSet.contains(ivCoar))
    {
      Real deltalo = 0.0;
      bool hasLoSten = extractFromSten(deltalo, a_loStenBF(a_coarVoF,0),
                                       a_coarData, a_ivar);
      if(hasLoSten)
        {
          deltaminmod = deltalo;
        }
      else
        {
          //no derivs exist because of EB
          deltaminmod = 0.0;
        }
    }//end if doing lo one sided diffs
  else if(hiInterpSet.contains(ivCoar))
    {
      Real deltahi = 0.0;
      bool hasHiSten = extractFromSten(deltahi, a_hiStenBF(a_coarVoF,0),
                                       a_coarData, a_ivar);
      if(hasHiSten)
        {
          deltaminmod = deltahi;
        }
      else
        {
          //no derivs exist because of EB
          deltaminmod = 0.0;
        }
    }//end if doing lo one sided diffs
  else
    {
      //has to be in one of the sets
      MayDay::Error("sets in ebpwl inconsistant");
    }
  return deltaminmod;
}

/************************************/
/************************************/
#include "NamespaceFooter.H"
