#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// dtgraves fri aug 17, 2001

#include "EBCoarseAverage.H"
#include "EBAverageF_F.H"
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "DisjointBoxLayout.H"
#include "EBCellFactory.H"
#include "EBFluxFactory.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

/************************************/
void
EBCoarseAverage::setDefaultValues()
{
  m_isDefined = false;
  m_refRat = -1;
  m_nComp = -1;
}
/************************************/
EBCoarseAverage::EBCoarseAverage()
{
  setDefaultValues();
}
/************************************/
EBCoarseAverage::~EBCoarseAverage()
{
}
/************************************/
EBCoarseAverage::EBCoarseAverage(const DisjointBoxLayout& a_dblFine,
                                 const DisjointBoxLayout& a_dblCoar,
                                 const EBISLayout& a_ebislFine,
                                 const EBISLayout& a_ebislCoar,
                                 const ProblemDomain& a_domainCoar,
                                 const int& a_nref,
                                 const int& a_nvar,
                                 const EBIndexSpace* ebisPtr)
{
  setDefaultValues();

  define(a_dblFine, a_dblCoar, a_ebislFine, a_ebislCoar,
         a_domainCoar, a_nref, a_nvar, ebisPtr);
}
/************************************/
EBCoarseAverage::EBCoarseAverage(const EBLevelGrid& a_eblgFine,
                                 const EBLevelGrid& a_eblgCoar,
                                 const EBLevelGrid& a_eblgCoarsenedFine,
                                 const int& a_nref,
                                 const int& a_nvar)
{
  setDefaultValues();

  define(a_eblgFine, a_eblgCoar, a_eblgCoarsenedFine, a_nref, a_nvar);
}
/************************************/
void
EBCoarseAverage::define(const DisjointBoxLayout& a_dblFine,
                        const DisjointBoxLayout& a_dblCoar,
                        const EBISLayout& a_ebislFine,
                        const EBISLayout& a_ebislCoar,
                        const ProblemDomain& a_domainCoar,
                        const int& a_nref,
                        const int& a_nvar,
                        
                        const EBIndexSpace* ebisPtr)
{
  CH_TIME("EBCoarseAverage::define");
  CH_assert(ebisPtr->isDefined());

  ProblemDomain domainFine = a_domainCoar;
  domainFine.refine(a_nref);
  EBLevelGrid eblgFine = EBLevelGrid(a_dblFine,  domainFine,2,Chombo_EBIS::instance());
  EBLevelGrid eblgCoar = EBLevelGrid(a_dblCoar,a_domainCoar,2,Chombo_EBIS::instance());
  EBLevelGrid eblgCoarsenedFine;
  coarsen(eblgCoarsenedFine,eblgFine,a_nref);
  if(a_nref >= 2)
    {
      eblgCoarsenedFine.getEBISL().setMaxRefinementRatio(a_nref,ebisPtr);
    }

  define(eblgFine, eblgCoar, eblgCoarsenedFine,a_nref, a_nvar);
}
/************************************/
void
EBCoarseAverage::define(const EBLevelGrid& a_eblgFine,
                        const EBLevelGrid& a_eblgCoar,
                        const EBLevelGrid& a_eblgCoarsenedFine,
                        const int& a_nref,
                        const int& a_nvar)
{
  CH_TIME("EBCoarseAverage::define(EBLG)");
  CH_assert(a_nref > 0);
  CH_assert(a_nvar > 0);
  CH_assert(a_eblgFine.getDBL().coarsenable(a_nref));
  CH_assert(a_eblgCoarsenedFine.getEBISL().getGhost()>0);
  CH_assert(a_eblgCoarsenedFine.getEBISL().getMaxRefinementRatio()>=a_nref);

  m_isDefined = true;
  m_refRat = a_nref;
  m_nComp = a_nvar;
  m_coarEBLG          = a_eblgCoar;
  m_fineEBLG          = a_eblgFine;
  m_coarsenedFineEBLG = a_eblgCoarsenedFine;
  IntVect ghostZero = IntVect::Zero;
  IntVect ghostUnit = IntVect::Unit;
  LayoutData<IntVectSet> irregSets(m_coarsenedFineEBLG.getDBL());
  for(DataIterator dit = m_coarsenedFineEBLG.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      irregSets[dit()] = m_coarsenedFineEBLG.getEBISL()[dit()].getIrregIVS(m_coarsenedFineEBLG.getDBL().get(dit()));
    }
  EBCellFactory        ebcellfact(m_coarsenedFineEBLG.getEBISL());
  EBFluxFactory        ebfluxfact(m_coarsenedFineEBLG.getEBISL());
  BaseIVFactory<Real>  baseivfact(m_coarsenedFineEBLG.getEBISL(), irregSets);
  m_coarsenedFineData.define(    m_coarsenedFineEBLG.getDBL(), m_nComp,ghostZero, ebcellfact);
  m_coarsenedFineFluxData.define(m_coarsenedFineEBLG.getDBL(), m_nComp,ghostUnit, ebfluxfact);
  m_coarsenedFineBIVFData.define(m_coarsenedFineEBLG.getDBL(), m_nComp,ghostZero, baseivfact);
}
/************************************/
bool
EBCoarseAverage::isDefined() const
{
  return m_isDefined;
}
/************************************/
void
EBCoarseAverage::average(LevelData<EBCellFAB>& a_coarData,
                         const LevelData<EBCellFAB>& a_fineData,
                         const Interval& a_variables)
{
  CH_TIME("EBCoarseAverage::average(LD<EBCellFAB>)");
  CH_assert(isDefined());
  //first average the data onto the coarsenedFine data
  //then copy over
  for(DataIterator fineit = m_fineEBLG.getDBL().dataIterator();
      fineit.ok(); ++fineit)
    {
      averageFAB(m_coarsenedFineData[fineit()],
                 a_fineData[fineit()],
                 fineit(),
                 a_variables);
    }

  //this is the part that is blocking
  m_coarsenedFineData.copyTo(a_variables, a_coarData, a_variables);
}
/************************************/
void
EBCoarseAverage::average(LevelData<EBFluxFAB>&       a_coarData,
                         const LevelData<EBFluxFAB>& a_fineData,
                         const Interval& a_variables)
{
  CH_TIME("EBCoarseAverage::average(LD<EBFluxFAB>)");
  CH_assert(isDefined());
  //first average the data onto the coarsenedFine data
  //then copy over
  for(DataIterator fineit = m_fineEBLG.getDBL().dataIterator();
      fineit.ok(); ++fineit)
    {
      EBFluxFAB& coarsenedFineFluxFAB = m_coarsenedFineFluxData[fineit()];
      const EBFluxFAB&    fineFluxFAB = a_fineData[fineit()];
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          averageFAB(coarsenedFineFluxFAB[idir],
                     fineFluxFAB[idir],
                     fineit(),
                     a_variables,
                     idir);
        }
    }

  //this is the part that is blocking
  m_coarsenedFineFluxData.copyTo(a_variables, a_coarData, a_variables);
}

/************************************/
void
EBCoarseAverage::average(LevelData<BaseIVFAB<Real> >&        a_coarData,
                         const LevelData<BaseIVFAB<Real> >&  a_fineData,
                         const Interval& a_variables)
{
  CH_TIME("EBCoarseAverage::average(LD<BaseIVFAB>)");
  CH_assert(isDefined());
  //first average the data onto the coarsenedFine data
  //then copy over
  for(DataIterator fineit = m_fineEBLG.getDBL().dataIterator();
      fineit.ok(); ++fineit)
    {
      BaseIVFAB<Real>& coarsenedFineFAB = m_coarsenedFineBIVFData[fineit()];
      const BaseIVFAB<Real>&    fineFAB = a_fineData[fineit()];
      averageFAB(coarsenedFineFAB,
                 fineFAB,
                 fineit(),
                 a_variables);
    }

  //this is the part that is blocking
  m_coarsenedFineBIVFData.copyTo(a_variables, a_coarData, a_variables);
}
/************************************/
void
EBCoarseAverage::averaRZ(LevelData<EBCellFAB>&       a_coarData,
                         const LevelData<EBCellFAB>& a_fineData,
                         const Interval&             a_variables,
                         const Real&                 a_dx)
{
  CH_assert(isDefined());
  //first average the data onto the coarsenedFine data
  //then copy over
  for(DataIterator fineit = m_fineEBLG.getDBL().dataIterator();
      fineit.ok(); ++fineit)
    {
      averageFRZ(m_coarsenedFineData[fineit()],
                 a_fineData[fineit()],
                 fineit(),
                 a_variables,
                 a_dx);
    }

  //this is the part that is blocking
  m_coarsenedFineData.copyTo(a_variables, a_coarData, a_variables);
}
/************************************/
void
EBCoarseAverage::averageFAB(EBCellFAB&       a_coar,
                            const EBCellFAB& a_fine,
                            const DataIndex& a_datInd,
                            const Interval&  a_variables) const
{
  CH_TIME("EBCoarseAverage::averageFAB(EBCellFAB)");
  CH_assert(isDefined());
  //do all cells as if they were regular
  BaseFab<Real>& coarRegFAB =  a_coar.getSingleValuedFAB();
  const BaseFab<Real>& fineRegFAB = a_fine.getSingleValuedFAB();
  Box refbox(IntVect::Zero,
             (m_refRat-1)*IntVect::Unit);
  const Box& coarBox = m_coarsenedFineEBLG.getDBL().get(a_datInd);

#ifndef NDEBUG
  Box fineBox = refine(coarBox, m_refRat);
  CH_assert(coarRegFAB.box().contains(coarBox));
  CH_assert(fineRegFAB.box().contains(fineBox));
#endif

  for(int ivar = a_variables.begin();
      ivar <= a_variables.end(); ivar++)
    {
      FORT_EBAVERAGE(CHF_FRA1(coarRegFAB,ivar),
                     CHF_CONST_FRA1(fineRegFAB,ivar),
                     CHF_BOX(coarBox),
                     CHF_CONST_INT(m_refRat),
                     CHF_BOX(refbox));
    }

  //overwrite irregular cells

  //recall that datInd is from the fine layout.
  const EBISBox& ebisBoxCoar = m_coarsenedFineEBLG.getEBISL()[a_datInd];
  const EBISBox& ebisBoxFine = m_fineEBLG.getEBISL()[a_datInd];
  IntVectSet coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

  //fine cell volume is normalized to one.
  //compute coarse volume
  Real dxCoar = Real(m_refRat);
  Real cellVolCoar = 1.0;
  for(int idir = 0; idir < SpaceDim; idir++)
    cellVolCoar *= dxCoar;

  for(VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph());
      vofitCoar.ok(); ++vofitCoar)
    {
      const VolIndex& coarVoF = vofitCoar();
      Vector<VolIndex> fineVoFs =
        m_coarsenedFineEBLG.getEBISL().refine(coarVoF, m_refRat, a_datInd);

      Real volCoar = cellVolCoar*ebisBoxCoar.volFrac(coarVoF);
      for(int ivar = a_variables.begin();
          ivar <= a_variables.end(); ivar++)
        {
          Real dataVal = 0.0;
          if(volCoar > 0.)
            {
              for(int ifine = 0; ifine < fineVoFs.size(); ifine++)
                {
                  const VolIndex& fineVoF = fineVoFs[ifine];
                  //fine cell volume is normalized to one...
                  Real volFine  = ebisBoxFine.volFrac(fineVoF);
                  Real fineVal =  a_fine(fineVoF, ivar);
                  if(volFine > 0.)
                    {
                      dataVal += fineVal*volFine;
                    }
                }
              dataVal /= volCoar;
            }
          else
            {
              //if there is no real volume, just take the ave
              //of fine values
              for(int ifine = 0; ifine < fineVoFs.size(); ifine++)
                {
                  const VolIndex& fineVoF = fineVoFs[ifine];
                  Real fineVal =  a_fine(fineVoF, ivar);
                  dataVal += fineVal;
                }
              if(fineVoFs.size() > 0)
                {
                  dataVal /= Real(fineVoFs.size());
                }
            }
          a_coar(coarVoF, ivar) = dataVal;
        }
    }
}

/**
   bndry area weighted averaging of irregular data
**********************************/

void
EBCoarseAverage::averageFAB(BaseIVFAB<Real>&       a_coar,
                            const BaseIVFAB<Real>& a_fine,
                            const DataIndex&       a_datInd,
                            const Interval&        a_variables) const
{
  CH_assert(isDefined());
  //recall that datInd is from the fine layout.
  const EBISBox& ebisBoxCoar = m_coarsenedFineEBLG.getEBISL()[a_datInd];
  const EBISBox& ebisBoxFine = m_fineEBLG.getEBISL()[a_datInd];
  const IntVectSet& coarIrregIVS = a_coar.getIVS();
  const IntVectSet& fineIrregIVS = a_fine.getIVS();

  for(VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph());
      vofitCoar.ok(); ++vofitCoar)
    {
      const VolIndex& coarVoF = vofitCoar();
      Vector<VolIndex> fineVoFs =
        m_coarsenedFineEBLG.getEBISL().refine(coarVoF, m_refRat, a_datInd);

      for(int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++)
        {
          int  numVoFs = 0;
          Real areaTot = 0;
          Real dataVal = 0;
          for(int ifine = 0; ifine < fineVoFs.size(); ifine++)
            {
              const VolIndex& fineVoF = fineVoFs[ifine];
              if(fineIrregIVS.contains(fineVoF.gridIndex()))
                {
                  Real bndryArea = ebisBoxFine.bndryArea(fineVoF);
                  if(bndryArea > 0)
                    {
                      areaTot += bndryArea;
                      numVoFs++;
                      dataVal += a_fine(fineVoF, ivar);
                    }
                }
            }
          if(numVoFs > 1)
            {
              dataVal /= Real(numVoFs);
            }
          a_coar(coarVoF, ivar) = dataVal;
        }
    }
}
/************************************/
void
EBCoarseAverage::averageFAB(EBFaceFAB&       a_coar,
                            const EBFaceFAB& a_fine,
                            const DataIndex& a_datInd,
                            const Interval&  a_variables,
                            const int&       a_dir) const
{
  CH_TIME("EBCoarseAverage::averageFAB(EBFaceFAB)");
  CH_assert(isDefined());

  //do all cells as if they were regular
  BaseFab<Real>& coarRegFAB =  a_coar.getSingleValuedFAB();
  const BaseFab<Real>& fineRegFAB = a_fine.getSingleValuedFAB();
  Box refbox(IntVect::Zero,
             (m_refRat-1)*IntVect::Unit);
  refbox.surroundingNodes(a_dir);
  const Box& coarDBLBox = m_coarsenedFineEBLG.getDBL().get(a_datInd);
  Box coarFaceBox = coarDBLBox;
  coarFaceBox.surroundingNodes(a_dir);

#ifndef NDEBUG
  Box fineBox = refine(coarFaceBox, m_refRat);
  CH_assert(coarRegFAB.box().contains(coarFaceBox));
  CH_assert(fineRegFAB.box().contains(fineBox));
#endif

  for(int ivar = a_variables.begin();
      ivar <= a_variables.end(); ivar++)
    {
      FORT_EBAVERAGEFACE(CHF_FRA1(coarRegFAB,ivar),
                         CHF_CONST_FRA1(fineRegFAB,ivar),
                         CHF_BOX(coarFaceBox),
                         CHF_CONST_INT(m_refRat),
                         CHF_BOX(refbox),
                         CHF_CONST_INT(a_dir));
    }

  //overwrite irregular cells

  //recall that datInd is from the fine layout.
  const EBISBox& ebisBoxCoar = m_coarsenedFineEBLG.getEBISL()[a_datInd];
  const EBISBox& ebisBoxFine = m_fineEBLG.getEBISL()[a_datInd];
  IntVectSet coarIrregIVS = ebisBoxCoar.getIrregIVS(coarDBLBox);

  //fine face area is normalized to one.
  //compute coarse volume
  Real dxCoar = Real(m_refRat);
  Real faceAreaCoar = 1.0;
  for(int idir = 0; idir < (SpaceDim - 1); idir++)
    faceAreaCoar *= dxCoar;

  for(FaceIterator faceitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph(),a_dir,
                              FaceStop::SurroundingWithBoundary);
      faceitCoar.ok(); ++faceitCoar)
    {
      const FaceIndex& coarFace = faceitCoar();
      Vector<FaceIndex> fineFaces = m_coarsenedFineEBLG.getEBISL().refine(coarFace,m_refRat,a_datInd);

      Real areaCoar = faceAreaCoar*ebisBoxCoar.areaFrac(coarFace);
      for(int ivar = a_variables.begin();
          ivar <= a_variables.end(); ivar++)
        {
          Real dataVal = 0.0;
          if(areaCoar > 0.)
            {
              for(int ifine = 0; ifine < fineFaces.size(); ifine++)
                {
                  const FaceIndex& fineFace = fineFaces[ifine];
                  //fine face area is normalized to one...
                  Real areaFine  = ebisBoxFine.areaFrac(fineFace);
                  Real fineVal =  a_fine(fineFace, ivar);
                  if(areaFine > 0.)
                    {
                      dataVal += fineVal*areaFine;
                    }
                }
              dataVal /= areaCoar;
            }
          else
            {
              //if there is no real area, just take the ave
              //of fine values
              for(int ifine = 0; ifine < fineFaces.size(); ifine++)
                {
                  const FaceIndex& fineFace = fineFaces[ifine];
                  Real fineVal =  a_fine(fineFace, ivar);
                  dataVal += fineVal;
                }
              if(fineFaces.size() > 0)
                {
                  dataVal /= Real(fineFaces.size());
                }
            }

          a_coar(coarFace, ivar) = dataVal;
        }
    }
}
/************************************/
void
EBCoarseAverage::averageFRZ(EBCellFAB&       a_coar,
                            const EBCellFAB& a_fine,
                            const DataIndex& a_datInd,
                            const Interval&  a_variables,
                            const Real&      a_dx) const
{
  CH_assert(isDefined());
  //do all cells as if they were regular
  BaseFab<Real>& coarRegFAB =  a_coar.getSingleValuedFAB();
  const BaseFab<Real>& fineRegFAB = a_fine.getSingleValuedFAB();
  Box refbox(IntVect::Zero,
             (m_refRat-1)*IntVect::Unit);
  const Box& coarBox = m_coarsenedFineEBLG.getDBL().get(a_datInd);

#ifndef NDEBUG
  Box fineBox = refine(coarBox, m_refRat);
  CH_assert(coarRegFAB.box().contains(coarBox));
  CH_assert(fineRegFAB.box().contains(fineBox));
#endif
  Real dxFine = a_dx;
  Real dxCoar = Real(m_refRat)*a_dx;

  for(int ivar = a_variables.begin();
      ivar <= a_variables.end(); ivar++)
    {
      FORT_EBAVERARZ(CHF_FRA1(coarRegFAB,ivar),
                     CHF_CONST_FRA1(fineRegFAB,ivar),
                     CHF_BOX(coarBox),
                     CHF_CONST_INT(m_refRat),
                     CHF_BOX(refbox),
                     CHF_REAL(dxCoar),
                     CHF_REAL(dxFine));
    }

  //overwrite irregular cells

  //recall that datInd is from the fine layout.
  const EBISBox& ebisBoxCoar = m_coarsenedFineEBLG.getEBISL()[a_datInd];
  const EBISBox& ebisBoxFine = m_fineEBLG.getEBISL()[a_datInd];
  IntVectSet coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

  //fine cell volume is normalized to one.

  for(VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph());
      vofitCoar.ok(); ++vofitCoar)
    {
      const VolIndex& vofCoar = vofitCoar();
      Real volCoar;
      EBArith::getCompVolRZ(volCoar, ebisBoxCoar, dxCoar, vofCoar);

      Vector<VolIndex> vofsFine =
        m_coarsenedFineEBLG.getEBISL().refine(vofCoar, m_refRat, a_datInd);

      for(int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++)
        {
          Real dataVal = 0.0;
          for(int ifine = 0; ifine < vofsFine.size(); ifine++)
            {
              const VolIndex& vofFine = vofsFine[ifine];
              Real  volFine;
              EBArith::getCompVolRZ(volFine, ebisBoxFine, dxFine, vofFine);
              Real fineVal =  a_fine(vofFine, ivar);
              if(volFine > 0.)
                {
                  dataVal += fineVal*volFine;
                }
            }
          if(volCoar > 0.0)
            {
              dataVal /= volCoar;
            }
          a_coar(vofCoar, ivar) = dataVal;
        }
    }
}
/************************************/
#include "NamespaceFooter.H"
