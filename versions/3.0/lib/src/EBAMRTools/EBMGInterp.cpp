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

#include "EBMGInterp.H"
#include "EBMGInterpF_F.H"
#include "VoFIterator.H"
#include "DisjointBoxLayout.H"
#include "EBCellFactory.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "CH_Timer.H"
#include "EBLevelGrid.H"
#include "NamespaceHeader.H"

/************************************/
void
EBMGInterp::setDefaultValues()
{
  m_isDefined = false;
  m_refRat = -1;
  m_nComp = -1;
}
/************************************/
EBMGInterp::EBMGInterp()
{
  setDefaultValues();
}
/************************************/
EBMGInterp::~EBMGInterp()
{
}
// /************************************/
EBMGInterp::EBMGInterp(const DisjointBoxLayout& a_dblFine,
                       const DisjointBoxLayout& a_dblCoar,
                       const EBISLayout& a_ebislFine,
                       const EBISLayout& a_ebislCoar,
                       const ProblemDomain& a_domainCoar,
                       const int& a_nref,
                       const int& a_nvar,
                       const EBIndexSpace*  ebisPtr,
                       const IntVect& a_ghost)
{
  setDefaultValues();

  define(a_dblFine, a_dblCoar, a_ebislFine, a_ebislCoar,
         a_domainCoar, a_nref, a_nvar, ebisPtr, a_ghost);

}
/************************************/
void
EBMGInterp::define(const DisjointBoxLayout&    a_dblFine,
                   const DisjointBoxLayout&    a_dblCoar,
                   const EBISLayout&           a_ebislFine,
                   const EBISLayout&           a_ebislCoar,
                   const ProblemDomain&        a_domainCoar,
                   const int&                  a_nref,
                   const int&                  a_nvar,
                   const EBIndexSpace*         ebisPtr,
                   const IntVect&              a_ghostCellsPhi)
{
  CH_TIMERS("EBMGInterp::define");
  CH_TIMER("fillEBISLayout", t1);
  CH_TIMER("defines", t2);
  CH_TIMER("setVal", t3);
  CH_TIMER("defineStencils", t4);
  CH_assert(a_nref > 0);
  m_ghost = a_ghostCellsPhi;

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
  m_copierRCtoF.define(m_refinedCoarseGrids, a_dblFine, a_ghostCellsPhi);
  m_copierFtoRC.define(a_dblFine, m_refinedCoarseGrids, a_ghostCellsPhi);

  CH_START(t1);
  ebisPtr->fillEBISLayout(m_refinedCoarseEBISL,
                          m_refinedCoarseGrids,
                          m_fineDomain, nghost);
  if(m_refRat > 2)
    {
      m_refinedCoarseEBISL.setMaxCoarseningRatio(m_refRat);
    }
  CH_STOP(t1);

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
  LayoutData<Vector<VoFStencil> > interpStencil;
  LayoutData<VoFIterator > vofItIrregRefCoar;

  m_interpEBStencil.define(m_coarGrids);
  interpStencil.define(    m_coarGrids);
  vofItIrregRefCoar.define(m_coarGrids);


  for(DataIterator dit = m_coarGrids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& boxRefCoar =          m_refinedCoarseGrids.get(dit());
      const EBISBox& ebisBoxRefCoar =  m_refinedCoarseEBISL[dit()];
      const EBGraph& ebGraphRefCoar = ebisBoxRefCoar.getEBGraph();

      const Box& boxCoar = a_dblCoar.get(dit());
      const EBISBox& ebisBoxCoar = a_ebislCoar[dit()];
      IntVectSet notRegularCoar = ebisBoxCoar.getIrregIVS(boxCoar);

      IntVectSet ivsRefCoar = refine(notRegularCoar, m_refRat);
      vofItIrregRefCoar[dit()].define(ivsRefCoar, ebGraphRefCoar);

      const Vector<VolIndex>& allFineVofs = vofItIrregRefCoar[dit()].getVector();

      Vector<VoFStencil>& interpStencils = interpStencil[dit()];
      interpStencils.resize(allFineVofs.size());

      for(int ifine = 0; ifine < allFineVofs.size(); ifine++)
        {
          VolIndex coarseVof;
          if(m_refRat > 2)
            {
              coarseVof = m_refinedCoarseEBISL.coarsen(allFineVofs[ifine], m_refRat, dit());
            }
          else
            {
              coarseVof = m_refinedCoarseEBISL[dit()].coarsen(allFineVofs[ifine]);
            }

          VoFStencil& interpStencil = interpStencils[ifine];

          //weight is one because there is only one vof to the stencil
          //and we are doing piecewise constant interpolation
          //for the copy stencil the weight is one because it is a copy
          interpStencil.add(coarseVof, 1.);
        }

      m_interpEBStencil[dit()] = RefCountedPtr<EBStencil>(new EBStencil(allFineVofs, interpStencil[dit()],  boxRefCoar, boxCoar,    ebisBoxRefCoar, ebisBoxCoar,       a_ghostCellsPhi, a_ghostCellsPhi));

    }
  CH_STOP(t4);
}
/************************************/
bool
EBMGInterp::isDefined() const
{
  return m_isDefined;
}

/************************************/
class EBAddOp : public LDOperator<EBCellFAB>
{
public:

  virtual void linearIn(EBCellFAB& arg,  void* buf, const Box& R,
                        const Interval& comps) const
  {
    EBCellFAB tmp;
    tmp.clone(arg);
    tmp.linearIn(buf, R, comps);
    arg.plus(tmp, R, comps.begin(), comps.begin(), comps.size());
  }

  void op(EBCellFAB& dest,
          const Box& RegionFrom,
          const Interval& Cdest,
          const Box& RegionTo,
          const EBCellFAB& src,
          const Interval& Csrc) const
  {
    dest.plus(src, RegionFrom, Csrc.begin(), Cdest.begin(), Cdest.size());
  }

};
void
EBMGInterp::pwcInterp(LevelData<EBCellFAB>&       a_fineData,
                      const LevelData<EBCellFAB>& a_coarData,
                      const Interval&             a_variables)
{
  CH_assert(a_fineData.ghostVect() == m_ghost);
  CH_assert(a_coarData.ghostVect() == m_ghost);

//   a_fineData.copyTo(a_variables, m_refinedCoarseData, a_variables, m_copierFtoRC);

  for(DataIterator dit = m_coarGrids.dataIterator();
      dit.ok(); ++dit)
    {
      m_refinedCoarseData[dit()].setVal(0.);
      pwcInterpFAB(m_refinedCoarseData[dit()],
                   a_coarData[dit()],
                   dit(),
                   a_variables);
    }

  EBAddOp op;
  m_refinedCoarseData.copyTo(a_variables, a_fineData, a_variables, m_copierRCtoF, op);
}
/************************************/
void
EBMGInterp::pwcInterpMG(LevelData<EBCellFAB>&       a_fineData,
                        const LevelData<EBCellFAB>& a_coarData,
                        const Interval&             a_variables)
{
  CH_assert(a_fineData.ghostVect() == m_ghost);
  CH_assert(a_coarData.ghostVect() == m_ghost);

  for(DataIterator dit = m_coarGrids.dataIterator();
      dit.ok(); ++dit)
    {
      pwcInterpFAB(a_fineData[dit()],
                   a_coarData[dit()],
                   dit(),
                   a_variables);
    }
}
/************************************/
void
EBMGInterp::pwcInterpFAB(EBCellFAB&       a_refCoar,
                         const EBCellFAB& a_coar,
                         const DataIndex& a_datInd,
                         const Interval&  a_variables) const
{
  CH_TIMERS("EBMGInterp::interp");
  CH_TIMER("regular_interp", t1);
  CH_TIMER("irregular_interp", t2);
  CH_assert(isDefined());

  const Box& coarBox = m_coarGrids.get(a_datInd);

  m_interpEBStencil[a_datInd]->cache(a_refCoar);

  //do all cells as if they were regular
  Box refBox(IntVect::Zero, IntVect::Zero);
  refBox.refine(m_refRat);

  const BaseFab<Real>& coarRegFAB =    a_coar.getSingleValuedFAB();
  BaseFab<Real>& refCoarRegFAB    = a_refCoar.getSingleValuedFAB();

  CH_START(t1);
  for(int ivar = a_variables.begin();  ivar <= a_variables.end(); ivar++)
    {

      FORT_REGPROLONG(CHF_FRA1(refCoarRegFAB,ivar),
                      CHF_CONST_FRA1(coarRegFAB,ivar),
                      CHF_BOX(coarBox),
                      CHF_BOX(refBox),
                      CHF_CONST_INT(m_refRat));
    }
  CH_STOP(t1);

  m_interpEBStencil[a_datInd]->uncache(a_refCoar);

  CH_START(t2);
  m_interpEBStencil[a_datInd]->apply(a_refCoar, a_coar, true);
  CH_STOP(t2);
}
/************************************/
#include "NamespaceFooter.H"
