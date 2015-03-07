#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif


#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "AMRMultiGrid.H"
#include "EBBackwardEuler.H"
#include "EBLevelDataOps.H"
#include "NamespaceHeader.H"

void EBBackwardEuler::
resetAlphaAndBeta(const Real& a_alpha,
                  const Real& a_beta)
{
  Vector<MGLevelOp< LevelData<EBCellFAB> >* > ops = m_solver->getAllOperators();
  for(int iop = 0; iop < ops.size(); iop++)
    {

      TGAHelmOp< LevelData<EBCellFAB> >* helmop = (TGAHelmOp< LevelData<EBCellFAB> >*) ops[iop];
      helmop->setAlphaAndBeta(a_alpha, a_beta);
    }
}

/*****/
TGAHelmOp<LevelData<EBCellFAB> >* 
EBBackwardEuler::
newOp(const ProblemDomain&                             a_indexSpace,
      const AMRLevelOpFactory<LevelData<EBCellFAB> >&  a_opFact)
{
  AMRLevelOpFactory<LevelData<EBCellFAB> >& opFact = (AMRLevelOpFactory<LevelData<EBCellFAB> >&) a_opFact;
  TGAHelmOp<LevelData<EBCellFAB> >* retval = (TGAHelmOp<LevelData<EBCellFAB> >*) opFact.AMRnewOp(a_indexSpace);
  return retval;
}
/*****/
EBBackwardEuler::
EBBackwardEuler(const RefCountedPtr<AMRMultiGrid< LevelData<EBCellFAB> > >& a_solver,
                const AMRLevelOpFactory<LevelData<EBCellFAB> > &            a_opFact,
                const ProblemDomain&                                        a_level0Domain,
                const Vector<int>&                                          a_refRat,
                int a_numLevels, int a_verbosity)
{
  m_verbosity = a_verbosity;
  m_level0Domain = a_level0Domain;
  m_refRat = a_refRat;
  m_solver  = a_solver;
  m_numLevels = a_numLevels;
  if(m_numLevels < 0)
    {
      m_numLevels = a_refRat.size();
    }

  m_ops.resize(m_numLevels);

  AMRLevelOpFactory<LevelData<EBCellFAB> >& opFact  = 
    (AMRLevelOpFactory<LevelData<EBCellFAB> >&) a_opFact;
  ProblemDomain curDom = m_level0Domain;
  for(int ilev = 0; ilev < m_numLevels; ilev++)
    {
      m_ops[ilev] = RefCountedPtr<TGAHelmOp<LevelData<EBCellFAB> > >(newOp(curDom, opFact));
      curDom.refine(a_refRat[ilev]);
    }

  m_dataCreated = false;
}

/*****/
EBBackwardEuler::~EBBackwardEuler()
{
  for(int ilev = 0; ilev < m_rhst.size(); ilev++)
    {
      if(m_rhst[ilev] != NULL)
        {
          delete m_rhst[ilev];
          m_rhst[ilev] = NULL;
        }
    }
}
/***/
void 
EBBackwardEuler::
createData(Vector<LevelData<EBCellFAB>* >&       a_source,
           int a_lbase,int a_lmax)
{
  m_dataCreated = true;
  m_rhst.resize(a_source.size(), NULL);
  for(int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      m_rhst[ilev] = new LevelData<EBCellFAB>();
      m_ops[ilev]->create(*m_rhst[ilev], *a_source[ilev]);
    }
}
/***/
void EBBackwardEuler::
oneStep(Vector<LevelData<EBCellFAB>* >&       a_phiNew,
        Vector<LevelData<EBCellFAB>* >&       a_phiOld,
        Vector<LevelData<EBCellFAB>* >&       a_source,
        const Real&        a_dt,
        int                a_lbase,
        int                a_lmax,
        bool               a_zeroPhi)
{
  if(!m_dataCreated)
    {
      createData(a_source, a_lbase, a_lmax);
    }

  if(m_verbosity >= 3)
    {
      pout() << "  EBBackwardEuler:: making rhs" << std::endl;
    }
  createEulerRHS(m_rhst, a_source, a_phiOld, a_lbase, a_lmax, a_dt);

  if(m_verbosity >= 3)
    { 
      pout() << "  EBBackwardEuler:: solving" << std::endl;
    }
  solveHelm(a_phiNew, m_rhst, a_lbase, a_lmax,a_dt, a_zeroPhi);
}
/*******/
//fills a_ans = dt*kappa*a_source + kappa*acoef*phiOld
void EBBackwardEuler::
createEulerRHS(Vector<LevelData<EBCellFAB>* >&   a_ans,
               Vector<LevelData<EBCellFAB>* >&   a_rho,
               Vector<LevelData<EBCellFAB>* >&   a_phiOld,
               int                               a_lbase,
               int                               a_lmax,
               Real                              a_dt)

{
  for(int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      DisjointBoxLayout grids = a_phiOld[ilev]->disjointBoxLayout();
      for(DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
        {
          //this makes rhs  = 0
          (*a_ans[ilev])[dit()].setVal(0.);
          //this makes rhs = phiOld
          (*a_ans[ilev])[dit()] +=  ((*a_phiOld[ilev])[dit()]);
        }
      //this makes rhs = kappa*acoef*(phi^n)
      m_ops[ilev]->diagonalScale(*a_ans[ilev]);      
      
      for(DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
        {
          const EBCellFAB& rho  = (*a_rho[ilev])[dit()];
          EBCellFAB scaleRho(rho.getEBISBox(), rho.box(), rho.nComp());
          scaleRho.setVal(0.);

          scaleRho += rho;
          //this makes scalerho = kappa*a_rho
          EBLevelDataOps::kappaWeight(scaleRho);

          //this makes scalerho = dt*kappa*a_rho          
          scaleRho *= a_dt;
          
          //this makes rhs = kappa*acoef*(phi^n) + dt*kappa*a_rho
          (*a_ans[ilev])[dit()] += scaleRho;
        }
    }
}
/*******/
void EBBackwardEuler::
solveHelm(Vector<LevelData<EBCellFAB>*>&       a_ans,
          Vector<LevelData<EBCellFAB>*>&       a_rhs,
          int               a_lbase,
          int               a_lmax,
          Real              a_dt,
          bool              a_zeroPhi)
{
  Real factor  = -a_dt;
  resetAlphaAndBeta(1.0, factor);

  m_solver->solveNoInit(a_ans, a_rhs, a_lmax, a_lbase, true);

  if(m_solver->m_exitStatus==1 && m_verbosity>3)
    {
      pout() << "EBBackwardEuler:: WARNING: solver exitStatus == 1" << std::endl;
    }
}

#include "NamespaceFooter.H"

