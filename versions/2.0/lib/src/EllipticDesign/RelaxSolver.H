#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// BVS, Feb 7, 2005

#ifndef _RELAXSOLVER_H_
#define _RELAXSOLVER_H_

#include "LinearSolver.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

///
/**
   Iterative solver which only uses the operator's preconditioner.
   Probably only useful as a bottom solver.
 */
template <class T>
class RelaxSolver : public LinearSolver<T>
{
public:
  RelaxSolver();
  virtual ~RelaxSolver();

  ///
  /**
     Set whether the solver uses only homogeneous boundary conditions
   */
  virtual void setHomogeneous(bool a_homogeneous)
     {m_homogeneous = a_homogeneous;}

  ///
  /**
     Define the solver.   a_op is the linear operator.
   */
  virtual void define(LinearOp<T>* a_op, bool a_homogeneous);

  ///
  /**
     Solve L(phi) = rho
   */
  virtual void solve(T& a_phi, const T& a_rhs);

  ///
  /**
     public member data: whether or not to use inhomogeneous boundary conditions.
   */
  bool m_homogeneous;

  ///
  /**
     public member data:  linear operator.
   */
  LinearOp<T>* m_op;

  ///
  /**
     public member data: maximum number of iterations
   */
  int m_imax;

  ///
  /**
     public member data: how much screen output user wants
   */
  int m_verbosity;

  ///
  /**
     public member data:  solver tolerance
   */
  Real m_eps;

};

// *******************************************************
// RelaxSolver Implementation
// *******************************************************

template <class T>
RelaxSolver<T>::RelaxSolver()
  :m_homogeneous(false), m_op(NULL), m_imax(40), m_verbosity(0), m_eps(1.0E-6) {;}

template <class T>
RelaxSolver<T>::~RelaxSolver(){;}

template <class T>
void RelaxSolver<T>::define(LinearOp<T>* a_operator, bool a_homogeneous)
{
  m_homogeneous = a_homogeneous;
  m_op = a_operator;
}

template <class T>
void RelaxSolver<T>::solve(T& a_phi, const T& a_rhs)
{
  CH_assert(m_op != NULL);
  CH_TIME("RelaxSolver::solve");
  T r, e;

  m_op->create(r, a_rhs);
  m_op->create(e, a_phi);
  


  m_op->residual(r, a_phi, a_rhs, m_homogeneous);

  Real norm = m_op->norm(r, 2);
  Real initialNorm = norm;
  int iter =0;
  while(iter <  m_imax){
    m_op->setToZero(e);
    m_op->preCond(e, r);
    m_op->incr(a_phi, e, 1.0);
    m_op->residual(r, a_phi, a_rhs, m_homogeneous);
    norm = m_op->norm(r, 2);
    if(m_verbosity > 3)
      {
        pout() << "iter# " << iter << " L2 Norm of resid: " << norm << endl;
      }
    iter++;
    if(norm <  m_eps*initialNorm) break;
  }

  m_op->clear(r);
  m_op->clear(e);
}

#include "NamespaceFooter.H"
#endif /*_RELAXSOLVER_H_*/
