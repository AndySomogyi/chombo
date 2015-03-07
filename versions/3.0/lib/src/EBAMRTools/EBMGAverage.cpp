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

#include "EBMGAverage.H"
#include "EBMGAverageF_F.H"
#include "VoFIterator.H"
#include "DisjointBoxLayout.H"
#include "EBCellFactory.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"
/************************************/
void
EBMGAverage::setDefaultValues()
{
  m_isDefined = false;
  m_refRat = -1;
  m_nComp = -1;
}
/************************************/
EBMGAverage::EBMGAverage()
{
  setDefaultValues();
}
/************************************/
EBMGAverage::~EBMGAverage()
{
}
 /************************************/
EBMGAverage::EBMGAverage(const DisjointBoxLayout& a_dblFine,
                         const DisjointBoxLayout& a_dblCoar,
                         const EBISLayout& a_ebislFine,
                         const EBISLayout& a_ebislCoar,
                         const ProblemDomain& a_domainCoar,
                         const int& a_nref,
                         const int& a_nvar,
                         const EBIndexSpace* ebisPtr,
                         const IntVect&              a_ghostCellsRHS)
{
  setDefaultValues();

  define(a_dblFine, a_dblCoar, a_ebislFine, a_ebislCoar,
         a_domainCoar, a_nref, a_nvar, ebisPtr, a_ghostCellsRHS);

}
/************************************/
void
EBMGAverage::define(const DisjointBoxLayout&    a_dblFine,
                    const DisjointBoxLayout&    a_dblCoar,
                    const EBISLayout&           a_ebislFine,
                    const EBISLayout&           a_ebislCoar,
                    const ProblemDomain&        a_domainCoar,
                    const int&                  a_nref,
                    const int&                  a_nvar,
                    const EBIndexSpace*         ebisPtr,
                    const IntVect&              a_ghostCellsRHS)
{
  CH_TIMERS("EBMGAverage::define");
  CH_TIMER("fillEBISLayout", t1);
  CH_TIMER("defines", t2);
  CH_TIMER("setVal", t3);
  CH_TIMER("defineStencils", t4);
  CH_assert(a_nref > 0);
  CH_assert(a_nvar > 0);
  m_ghost = a_ghostCellsRHS;
  m_isDefined = true;
  m_refRat = a_nref;
  m_nComp = a_nvar;
  m_coarGrids = a_dblCoar;
  m_fineGrids = a_dblFine;
  m_coarEBISL = a_ebislCoar;
  m_fineEBISL = a_ebislFine;
  m_coarDomain = a_domainCoar;
  m_fineDomain = refine(m_coarDomain, m_refRat);


  //const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  CH_assert(ebisPtr->isDefined());
  int nghost = 4;
  m_refinedCoarseGrids = DisjointBoxLayout();
  refine(m_refinedCoarseGrids, m_coarGrids, m_refRat);
  m_copier.define(a_dblFine, m_refinedCoarseGrids, a_ghostCellsRHS);

  CH_START(t1);
  ebisPtr->fillEBISLayout(m_refinedCoarseEBISL,
                          m_refinedCoarseGrids,
                          m_fineDomain, nghost);
  CH_STOP(t1);
  if(m_refRat > 2)
    {
      m_coarEBISL.setMaxRefinementRatio(m_refRat, ebisPtr);
    }
  EBCellFactory ebcellfact(m_refinedCoarseEBISL);
  CH_START(t2);
  m_refinedCoarseData.define(m_refinedCoarseGrids, m_nComp,
                             m_ghost, ebcellfact);
  CH_STOP(t2);
  CH_START(t3);
  for(DataIterator dit = m_refinedCoarseGrids.dataIterator(); dit.ok(); ++dit)
    {
      m_refinedCoarseData[dit()].setVal(0.0);
    }
  CH_STOP(t3);

  CH_START(t4);
  LayoutData<Vector<VoFStencil> > averageStencil;
  LayoutData<VoFIterator > vofItIrregCoar;
  LayoutData<VoFIterator > vofItIrregRefCoar;
  m_averageEBStencil.define(m_coarGrids);
  averageStencil.define(    m_coarGrids);
  vofItIrregRefCoar.define(m_coarGrids);
  vofItIrregCoar.define(m_coarGrids);

  for(DataIterator dit = m_coarGrids.dataIterator(); dit.ok(); ++dit)
    {
      Box refBox(IntVect::Zero, IntVect::Zero);
      refBox.refine(m_refRat);
      int numFinePerCoar = refBox.numPts();

      const Box& boxRefCoar =          m_refinedCoarseGrids.get(dit());
      const EBISBox& ebisBoxRefCoar =  m_refinedCoarseEBISL[dit()];
      const EBGraph& ebGraphRefCoar = ebisBoxRefCoar.getEBGraph();

      const Box& boxCoar = a_dblCoar.get(dit());
      const EBISBox& ebisBoxCoar = a_ebislCoar[dit()];
      const EBGraph& ebGraphCoar = ebisBoxCoar.getEBGraph();

      IntVectSet notRegularCoar = ebisBoxCoar.getIrregIVS(boxCoar);
      vofItIrregCoar[dit()].define(notRegularCoar, ebGraphCoar);

      IntVectSet ivsRefCoar = refine(notRegularCoar, m_refRat);
      vofItIrregRefCoar[dit()].define(ivsRefCoar, ebGraphRefCoar);

      const Vector<VolIndex>& allCoarVofs = vofItIrregCoar[dit()].getVector();

      Vector<VoFStencil>& averageStencils = averageStencil[dit()];

      averageStencils.resize(allCoarVofs.size());

      for(int icoar = 0; icoar < allCoarVofs.size(); icoar++)
        {
          Vector<VolIndex> fineVofs;
          if(m_refRat > 2)
            {
              fineVofs = a_ebislCoar.refine(allCoarVofs[icoar], m_refRat, dit());
            }
          else
            {
              fineVofs = a_ebislCoar[dit()].refine(allCoarVofs[icoar]);
            }

          VoFStencil& averageStencil = averageStencils[icoar];

          for (int ifine = 0; ifine < fineVofs.size(); ifine++)
            {
              averageStencil.add(fineVofs[ifine], 1./numFinePerCoar);
            }
        }

      m_averageEBStencil[dit()] = RefCountedPtr<EBStencil>(new EBStencil(allCoarVofs, averageStencil[dit()], boxCoar,  boxRefCoar, ebisBoxCoar,  ebisBoxRefCoar,  a_ghostCellsRHS, a_ghostCellsRHS));

    }
  CH_STOP(t4);
}
/************************************/
bool
EBMGAverage::isDefined() const
{
  return m_isDefined;
}
/************************************/
void
EBMGAverage::average(LevelData<EBCellFAB>&       a_coarData,
                     const LevelData<EBCellFAB>& a_fineData,
                     const Interval&             a_variables)
{
  CH_assert(a_fineData.ghostVect() == m_ghost);
  CH_assert(a_coarData.ghostVect() == m_ghost);
  CH_assert(isDefined());

  a_fineData.copyTo(a_variables, m_refinedCoarseData, a_variables, m_copier);

  for(DataIterator dit = m_coarGrids.dataIterator();
      dit.ok(); ++dit)
    {
      averageFAB(a_coarData[dit()],
                 m_refinedCoarseData[dit()],
                 dit(),
                 a_variables);
    }
}
/************************************/
void
EBMGAverage::averageMG(LevelData<EBCellFAB>&       a_coarData,
                       const LevelData<EBCellFAB>& a_fineData,
                       const Interval&             a_variables)
{
  CH_assert(a_fineData.ghostVect() == m_ghost);
  CH_assert(a_coarData.ghostVect() == m_ghost);
  CH_assert(isDefined());

  for(DataIterator dit = m_coarGrids.dataIterator();
      dit.ok(); ++dit)
    {
      averageFAB(a_coarData[dit()],
                 a_fineData[dit()],
                 dit(),
                 a_variables);
    }
}
/************************************/
void
EBMGAverage::averageFAB(EBCellFAB&       a_coar,
                        const EBCellFAB& a_refCoar,
                        const DataIndex& a_datInd,
                        const Interval&  a_variables) const
{
  CH_TIMERS("EBMGAverage::average");
  CH_TIMER("regular_average", t1);
  CH_TIMER("irregular_average", t2);
  CH_assert(isDefined());

  const Box& coarBox = m_coarGrids.get(a_datInd);

  //do all cells as if they were regular
  Box refBox(IntVect::Zero, IntVect::Zero);
  refBox.refine(m_refRat);
  int numFinePerCoar = refBox.numPts();

  BaseFab<Real>& coarRegFAB =             a_coar.getSingleValuedFAB();
  const BaseFab<Real>& refCoarRegFAB = a_refCoar.getSingleValuedFAB();

  //set to zero because the fortran is a bit simpleminded
  //and does stuff additively
  a_coar.setVal(0.);
  CH_START(t1);
  for(int comp = a_variables.begin();  comp <= a_variables.end(); comp++)
    {
      FORT_REGAVERAGE(CHF_FRA1(coarRegFAB,comp),
                      CHF_CONST_FRA1(refCoarRegFAB,comp),
                      CHF_BOX(coarBox),
                      CHF_BOX(refBox),
                      CHF_CONST_INT(numFinePerCoar),
                      CHF_CONST_INT(m_refRat));
    }
  CH_STOP(t1);

  //this is really volume-weighted averaging even though it does
  //not look that way.

  //so (in the traditional sense) we want to preserve
  //rhoc * volc = sum(rhof * volf)
  //this translates to
  //volfrac_C * rhoC = (1/numFinePerCoar)(sum(volFrac_F * rhoF))
  //but the data input to this routine is all kappa weigthed so
  //the volumefractions have already been multiplied
  //which means
  // rhoC = (1/numFinePerCoar)(sum(rhoF))
  //which is what this does

  CH_START(t2);
  m_averageEBStencil[a_datInd]->apply(a_coar, a_refCoar, false);
  CH_STOP(t2);

}
/************************************/
#include "NamespaceFooter.H"
