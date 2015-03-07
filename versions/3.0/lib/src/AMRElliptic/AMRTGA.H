#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _AMRTGA_H_
#define _AMRTGA_H_

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "AMRMultiGrid.H"
#include "NamespaceHeader.H"

///
/**
   This operator is meant to represent the general helmholtz operator that
   AMRTGA will be solving.  This operator is of the form
   alpha A(x) phi + beta div B(x) grad phi = rho.
   AMRTGA needs to reset the constants alpha and beta often.
 */
template <class T>
class TGAHelmOp : public AMRLevelOp<T>
{
public:
  TGAHelmOp()
  {;}
  virtual ~TGAHelmOp()
  {;}

  ///set the constants in the equation
  virtual void setAlphaAndBeta(const Real& a_alpha,
                               const Real& a_beta) = 0;

  ///
  /**
     Set the diagonal scaling of the operator.  If you are solving
     rho(x) dphi/dt = L(phi),
     this would mean multiply by rho.
     In EB applications, even for constant coefficients, it means to
     multiply by kappa.
   */
  virtual void diagonalScale(T& a_rhs) = 0;

  /// apply operator without any boundary or coarse-fine boundary conditions and no finer level
  virtual void applyOpNoBoundary(T& a_ans, const T& a_phi) = 0;


};

///
/**
 */
template <class T, class TFlux>
class LevelTGAHelmOp : public TGAHelmOp< T >
{
public:
   ///  These functions are part of the LevelTGA interface......

  ///
  virtual void fillGrad(const T& a_phi) = 0;

  ///
  virtual void getFlux(TFlux&           a_flux,
                       const T&         a_data,
                       const Box&       a_grid,
                       const DataIndex& a_dit,
                       Real             a_scale) = 0;

  virtual ~LevelTGAHelmOp()
  {;}
};



///
/**
   Template implementation of the TGA algorithm
   to solve the heat equation.
 **/
template <class T>
class AMRTGA
{
public:

  ///
  TGAHelmOp<T>* newOp(const ProblemDomain&                                  a_indexSpace,
                      const AMRLevelOpFactory<T>&                           a_opFact)
  {
    AMRLevelOpFactory<T>& opFact = (AMRLevelOpFactory<T>&) a_opFact;
    TGAHelmOp<T>* retval = (TGAHelmOp<T>*) opFact.AMRnewOp(a_indexSpace);
    return retval;
  }

  ///
  ~AMRTGA();

  ///
  /**
   **/
  AMRTGA(const RefCountedPtr<AMRMultiGrid<T> >     &           a_solver,
         const AMRLevelOpFactory<T>&                           a_factory,
         const ProblemDomain&                                  a_level0Domain,
         const Vector<int>&                                    a_refRat,
         int a_numLevels = -1,
         int a_verbosity = 0);


  ///
  /**
     This advances a parabolic pde from a_phiOld to a_phiNew using TGA on a non-moving domain with source term a_source
  **/
  void oneStep(Vector<T*>&             a_phiNew,
               Vector<T*>&             a_phiOld,
               Vector<T*>&             a_source,
               const Real&             a_dt,
               int                     a_lbase,
               int                     a_lmax);

  ///
  void resetAlphaAndBeta(const Real& a_alpha,
                         const Real& a_beta);

protected:
  void solveHelm(Vector<T*>&       a_ans,
                 Vector<T*>&       a_rhs,
                 int               a_lbase,
                 int               a_lmax,
                 Real              a_mu,
                 Real              a_dt);

  void applyHelm(Vector<T*>&       a_ans,
                 Vector<T*>&       a_source,
                 int               a_lbase,
                 int               a_lmax,
                 Real              a_mu,
                 Real              a_dt,
                 bool              a_homogeneousBC);

  void setMu();

  void createData(Vector<T* >&       a_source,
                  int                a_lbase,
                  int                a_lmax);

private:
  ProblemDomain                                  m_level0Domain;
  Vector<int>                                    m_refRat;
  RefCountedPtr<AMRMultiGrid<T> >                m_solver;
  Vector<RefCountedPtr<TGAHelmOp<T> > >          m_ops;
  Real m_mu1, m_mu2, m_mu3, m_mu4;
  int m_verbosity, m_numLevels;
  bool                                           m_dataCreated;
  Vector<T*>                                     m_rhst;

  //copy constructor and operator= disallowed for all the usual reasons
  AMRTGA(const AMRTGA<T>& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  void operator=(const AMRTGA<T>& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  /// weak construction is bad.   Ref Counted pointers are your friends.
  AMRTGA()
  {
    MayDay::Error("invalid operator");
  }


};

/*****/
template <class T>
void AMRTGA<T>::
resetAlphaAndBeta(const Real& a_alpha,
                  const Real& a_beta)
{
  Vector<MGLevelOp<T>* > ops = m_solver->getAllOperators();
  for(int iop = 0; iop < ops.size(); iop++)
    {

      TGAHelmOp<T>* helmop = (TGAHelmOp<T>*) ops[iop];
      helmop->setAlphaAndBeta(a_alpha, a_beta);
    }
}

/*****/
template <class T>
AMRTGA<T>::
AMRTGA(const RefCountedPtr<AMRMultiGrid<T> >&                a_solver,
       const AMRLevelOpFactory<T> &                          a_opFact,
       const ProblemDomain&                                  a_level0Domain,
       const Vector<int>&                                    a_refRat,
       int a_numLevels,
       int a_verbosity)
{
  //cast to remove const because base class definition is weird
  AMRLevelOpFactory<T>& opFact  = (AMRLevelOpFactory<T>&) a_opFact;

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

  ProblemDomain curDom = m_level0Domain;
  for(int ilev = 0; ilev < m_numLevels; ilev++)
    {
      m_ops[ilev] = RefCountedPtr<TGAHelmOp<T> >(newOp(curDom, opFact));
      curDom.refine(a_refRat[ilev]);
    }
  setMu();
  m_dataCreated = false;
}
template <class T>
void AMRTGA<T>::
setMu()
{
  Real tgaEpsilon = 1.e-12;
  Real a = 2.0 - sqrt(2.0) - tgaEpsilon;
  m_mu1 = (a - sqrt(pow(a,2) - 4.0*a + 2.0))/2.0 ;
  m_mu2 = (a + sqrt(pow(a,2) - 4.0*a + 2.0))/2.0 ;
  m_mu3 = (1.0 - a);
  m_mu4 = 0.5 - a;
  if(m_verbosity > 4)
    {
      pout() << "   AMRTGA:: epsilon = " << tgaEpsilon << std::endl;
    }
}

template <class T>
void AMRTGA<T>::
createData(Vector<T* >&       a_source,
           int                a_lbase,
           int                a_lmax)
{
  m_dataCreated = true;
  m_rhst.resize(a_source.size(), NULL);
  for(int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      m_rhst[ilev] = new T();
      m_ops[ilev]->create(*m_rhst[ilev], *a_source[ilev]);
    }
}
/*****/
template <class T>
AMRTGA<T>::~AMRTGA()
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
template <class T>
void AMRTGA<T>::
oneStep(Vector<T* >&       a_phiNew,
        Vector<T* >&       a_phiOld,
        Vector<T* >&       a_source,
        const Real&        a_dt,
        int                a_lbase,
        int                a_lmax)
{
  if(!m_dataCreated)
    {
      createData(a_source, a_lbase, a_lmax);
    }

  if(m_verbosity > 3)
    {
      pout() << "  AMRTGA:: starting mu4 operation" << std::endl;
    }
  applyHelm(m_rhst,a_source, a_lbase, a_lmax, m_mu4, a_dt, true);
  for(int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      m_ops[ilev]->scale(*m_rhst[ilev], a_dt);
    }

  if(m_verbosity > 3)
    {
      pout() << "  AMRTGA:: starting mu3 operation" << std::endl;
    }
  applyHelm(a_phiNew,a_phiOld, a_lbase, a_lmax, m_mu3, a_dt, false);

  for(int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      m_ops[ilev]->incr(*m_rhst[ilev], *a_phiNew[ilev], 1.0);
    }

  if(m_verbosity > 2)
    {
      pout() << "  AMRTGA:: starting mu2 operation" << std::endl;
    }
  solveHelm(a_phiNew,m_rhst, a_lbase, a_lmax, m_mu2, a_dt);
  for(int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      m_ops[ilev]->assign(*m_rhst[ilev], *a_phiNew[ilev]);
    }

  if(m_verbosity > 2)
    {
      pout() << "  AMRTGA:: starting mu1 operation" << std::endl;
    }
  for(int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      m_ops[ilev]->diagonalScale(*m_rhst[ilev]);
    }
  solveHelm(a_phiNew,m_rhst, a_lbase, a_lmax, m_mu1, a_dt);
}
/*******/
template <class T>
void AMRTGA<T>::
applyHelm(Vector<T*>&      a_ans,
          Vector<T*>&      a_phi,
          int              a_lbase,
          int              a_lmax,
          Real             a_mu,
          Real             a_dt,
          bool             a_homogeneousBC)
{
  Real factor  = a_mu*a_dt;

  resetAlphaAndBeta(1.0, factor);

  m_solver->computeAMROperator(a_ans,
                               a_phi,
                               a_lmax,
                               a_lbase,
                               a_homogeneousBC);

  for(int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      m_ops[ilev]->scale(*a_ans[ilev], -1);
    }
}
/*******/
template <class T>
void AMRTGA<T>::
solveHelm(Vector<T*>&       a_ans,
          Vector<T*>&       a_rhs,
          int               a_lbase,
          int               a_lmax,
          Real              a_mu,
          Real              a_dt)
{
  Real factor  = -a_mu*a_dt;

  resetAlphaAndBeta(1.0, factor);

  m_solver->solveNoInit(a_ans,
                        a_rhs,
                        a_lmax, a_lbase);
  if(m_solver->m_exitStatus==1 && m_verbosity>3)
    {
      pout() << "AMRTGA:: WARNING: solver exitStatus == 1" << std::endl;
    }
}

#include "NamespaceFooter.H"
#endif
