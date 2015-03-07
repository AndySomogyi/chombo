#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EBTensorCFInterp.H"
#include "LayoutIterator.H"
#include "DataIterator.H"
#include "EBAlias.H"
#include "NamespaceHeader.H"

/***********************/
// default constructor
/***********************/
EBTensorCFInterp::
EBTensorCFInterp(const DisjointBoxLayout& a_fineBoxes,
                 const DisjointBoxLayout& a_coarBoxes,
                 const EBISLayout&        a_ebislFine,
                 const EBISLayout&        a_ebislCoar,
                 const ProblemDomain&     a_domainCoar,
                 const int&               a_nref,
                 const int&               a_nvar,
                 Real                     a_dxFine)
{
  CH_assert(a_nvar > 0);
  CH_assert(a_nref > 0);
  CH_assert(a_fineBoxes.isClosed());
  CH_assert(a_coarBoxes.isClosed());
  CH_assert(!a_domainCoar.isEmpty());

  m_nComp  = a_nvar;
  m_refRat = a_nref;
  m_domainCoar = a_domainCoar;
  m_domainFine = refine(m_domainCoar, m_refRat);
  m_gridsFine = a_fineBoxes;
  m_gridsCoar = a_coarBoxes;
  m_ebislFine = a_ebislFine;
  m_ebislCoar = a_ebislCoar;
  if(m_ebislFine.getMaxCoarseningRatio() < m_refRat)
    {
      m_ebislFine.setMaxCoarseningRatio(m_refRat);
    }
  Real dxFine = a_dxFine;
  CH_assert(a_nvar > 0);
  CH_assert(a_nref > 0);
  CH_assert(a_fineBoxes.isClosed());
  CH_assert(a_coarBoxes.isClosed());
  CH_assert(!a_domainCoar.isEmpty());

  m_nComp  = a_nvar;
  m_refRat = a_nref;
  m_domainCoar = a_domainCoar;
  m_domainFine = refine(m_domainCoar, m_refRat);
  m_gridsFine = a_fineBoxes;
  m_gridsCoar = a_coarBoxes;
  m_ebislFine = a_ebislFine;
  m_ebislCoar = a_ebislCoar;
  if(m_ebislFine.getMaxCoarseningRatio() < m_refRat)
    {
      m_ebislFine.setMaxCoarseningRatio(m_refRat);
    }
  m_tensorCFInterp.define(a_fineBoxes, &a_coarBoxes, dxFine, a_nref, a_nvar, m_domainFine);
}
/***********************/
void
EBTensorCFInterp::coarseFineInterp(LevelData<EBCellFAB>&       a_fineData,
                                   LevelData<EBCellFAB>&       a_tanGradF,
                                   const LevelData<EBCellFAB>& a_coarData)
{
  CH_assert(a_fineData.nComp() == m_nComp);
  CH_assert(a_coarData.nComp() == m_nComp);
  CH_assert(a_tanGradF.nComp() == SpaceDim*m_nComp);
  LevelData<FArrayBox> fineDataLDFAB, coarDataLDFAB, tanGradFLDFAB;

  aliasEB(fineDataLDFAB, a_fineData);
  aliasEB(tanGradFLDFAB, a_tanGradF);
  aliasEB(coarDataLDFAB, (LevelData<EBCellFAB>&)a_coarData);
  m_tensorCFInterp.coarseFineInterp(fineDataLDFAB, tanGradFLDFAB, coarDataLDFAB);
}
void
EBTensorCFInterp::coarseFineInterpH(LevelData<EBCellFAB>& a_fineData,
                                    LevelData<EBCellFAB>& a_tanGradF)
{
  CH_assert(a_fineData.nComp() == m_nComp);
  CH_assert(a_tanGradF.nComp() == SpaceDim*m_nComp);
  LevelData<FArrayBox> fineDataLDFAB, tanGradFLDFAB;

  aliasEB(fineDataLDFAB, a_fineData);
  aliasEB(tanGradFLDFAB, a_tanGradF);
  m_tensorCFInterp.coarseFineInterpH(fineDataLDFAB, tanGradFLDFAB);
}
/***********************/
EBTensorCFInterp::~EBTensorCFInterp()
{
}











#include "NamespaceFooter.H"
