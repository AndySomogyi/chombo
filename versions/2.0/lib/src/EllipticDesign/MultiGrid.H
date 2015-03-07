#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// BVS, June 18, 2003

// We can assume that template class T has null construction.

#ifndef _MULTIGRID_H_
#define _MULTIGRID_H_

#include "CH_Timer.H"
#include "LinearSolver.H"
#include "ProblemDomain.H"
#include "REAL.H"
#include "Box.H"
#include "parstream.H"
#include "NamespaceHeader.H"

///
/**
   class MGLevelOp handles the additional tasks of coordinating operations
   between this level, and the level one coarser than this 'level'
*/
template <class T>
class MGLevelOp : public LinearOp<T>
{
public:
  ///
  /**
   */
  virtual ~MGLevelOp(){;}

  ///
  /**
     Create a coarsened  (by two) version of the input data.  This does not include averaging 
     the data.   So if a_fine is over a Box of (0, 0, 0) (63, 63, 63), a_fine should
     be over a Box (0, 0, 0) (31, 31, 31).   
   */
  virtual void createCoarser(T& a_coarse, const T& a_fine, bool ghosted) = 0;

  ///
  /**
     Use your relaxtion operator to remove the high frequency wave numbers from
     the correction so that it may be averaged to a coarser refinement.   
     A point relaxtion scheme, for example takes the form  
     a_correction -= lambda*(L(a_correction) - a_residual).   
   */
  virtual void relax(T& a_correction, const T& a_residual, int a_iterations) = 0 ;


  ///
  /**
     calculate restricted residual
     a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
  */
  virtual void restrictResidual(T& a_resCoarse, T& a_phiFine, const T& a_rhsFine) = 0;

  ///
  /**
     correct the fine solution based on coarse correction
     a_phiThisLevel += I[2h->h](a_correctCoarse)
  */
  virtual void prolongIncrement(T& a_phiThisLevel, const T& a_correctCoarse) = 0;

};

///
/**
   Factory class for generating MGLevelOps
 */
template <class T>
class MGLevelOpFactory 
{
public:
  virtual ~MGLevelOpFactory(){;}

  ///
  /**
     Create an operator at an  index space = coarsen(a_fineIndexSpace, 2^a_depth)
     Return NULL if no such Multigrid level can be created at this a_depth.
     If a_homoOnly = true, then only homogeneous boundary conditions will be needed.  
  */
  virtual MGLevelOp<T>* MGnewOp(const ProblemDomain& a_FineindexSpace,
                                int                  a_depth,
                                bool                 a_homoOnly = true) = 0;
};


///
/**
   Class to execute geometric multigrid.   This class is not meant to be particularly
   user-friendly and a good option for people who want something a tad less raw is to
   use AMRMultiGrid instead.
 */
template <class T>
class MultiGrid
{
public:
  MultiGrid();
  virtual ~MultiGrid();

  /// Function to define a MultiGrid object.
  /**
     a_factory is the factory for generating operators.   
     a_bottomSolver is called at the bottom of v-cycle.
     a_domain is the problem domain at the top of the vcycle.
     maxDepth defines the location of the bottom of the v-cycle.
     The vycle will terminate (hit bottom) when the factory returns NULL for a paticular
     depth if maxdepth = -1.  Otherwise the vycle terminates at maxdepth.
   */
  virtual void define(MGLevelOpFactory<T>& a_factory,
                      LinearSolver<T>* a_bottomSolver,
                      const ProblemDomain& a_domain,
                      int  a_maxDepth = -1);

  ///
  /**
     solve L(a_phi) = a_rhs.   Tolerance is how much you want the norm of the error reduced.
     verbosity is how chatty you want the function to be.   maxIterations is the maximum number
     of v-cycles.   This does the whole residual correction switcharoo and calls oneCycle up to
     maxIterations times, evaluating the residual as it goes.
   */
  virtual void solve(T& a_phi, const T& a_rhs, Real a_tolerance, int a_maxIterations, int verbosity= 0);

  /// Execute ONE v-cycle of multigrid.
  /** 
      If you want the solution to converge, you need to iterate this.
      See solve() or AMRMultiGrid::solve for a more automatic solve() function.
      This operates residual-correction form of solution
      so all boundary conditions are assumed to be homogeneous.
      L(a_e) = a_res
  */
  virtual void oneCycle(T& a_e, const T& a_res);

  
  // used by AMRMultiGrid--not a part of the public API
  void init(const T& a_correction, const T& a_residual); 

  //internal function
  void cycle(int a_depth, T& a_correction, const T& a_residual);

  //internal function
  void clear();

  // used by AMRMultiGrid--not a part of the public API
  void setBottomSolver(LinearSolver<T>* a_bottomSolver);

  ///
  /**
     Public solver parameters.   m_pre and m_post are the ones 
     that usually get set and are the number of relaxations performed
     before and after multigrid recursion.    See AMRMultiGrid for a
     more user-friendly interface.
  */
  int  m_depth, m_pre, m_post, m_cycle,m_numMG;
  bool m_homogeneous;
  LinearSolver<T>*  m_bottomSolver;

  ///
  /**
     for changing coefficients --- not for the faint of heart.
  */
  Vector< MGLevelOp<T>* > getAllOperators();

protected:
  bool m_defined;

  ProblemDomain m_topLevelDomain;
  Vector<MGLevelOp<T>*> m_op;
  Vector< T* >   m_residual;
  Vector< T* >   m_correction;
};

// *******************************************************
// MultiGrid Implementation
// *******************************************************

template <class T>
Vector< MGLevelOp<T>* > 
MultiGrid<T>::getAllOperators()
{
  Vector<MGLevelOp<T>*> retval = m_op;
  return retval;
}

template <class T>
MultiGrid<T>::MultiGrid():
   m_pre(3), m_post(3), m_cycle(1), m_numMG(1), m_homogeneous(true),
  m_defined(false){;}

template <class T>
MultiGrid<T>::~MultiGrid()
{
  clear();
}

template <class T>
void MultiGrid<T>::define(MGLevelOpFactory<T>& a_factory,
                          LinearSolver<T>*     a_bottomSolver,
                          const ProblemDomain& a_domain,
                          int                  a_maxDepth)
{

  CH_TIME("MultiGrid::define");
  m_depth = 0;
  m_bottomSolver = a_bottomSolver;
  m_residual.resize(m_depth);   // the zero'th entry is left null constructed
  m_correction.resize(m_depth); // the zero'th entry is left null constructed
  m_op.resize(m_depth);
  m_topLevelDomain = a_domain;

  MGLevelOp<T>* nextOp =  a_factory.MGnewOp(a_domain, 0, true);
  while(nextOp != NULL)
    {
      m_residual.push_back(new T());
      m_correction.push_back(new T());
      m_op.push_back(nextOp);
      m_depth++;
      if((m_depth < a_maxDepth) || (a_maxDepth < 0))
        {
          nextOp = a_factory.MGnewOp(a_domain , m_depth, true);
        }
      else
        {
          nextOp = NULL;
        }
    }
  m_bottomSolver->define(m_op[m_depth-1], true);
  m_defined = true;
}

template <class T>
void MultiGrid<T>::setBottomSolver(LinearSolver<T>* a_bottomSolver)
{
  m_bottomSolver = a_bottomSolver;
  m_bottomSolver->define(m_op[m_depth-1], true);
}

template <class T>
void MultiGrid<T>::init(const T& a_e, const T& a_residual)
{
  CH_TIME("MutliGrid::init");
  if (m_depth > 1)
    {
      m_op[0]->createCoarser(*(m_residual[1]), a_residual, false);
      m_op[0]->createCoarser(*(m_correction[1]), a_e, true);
    }

  for (int i = 2; i < m_depth; i++)
    {
      m_op[i-1]->createCoarser(*(m_residual[i]), *(m_residual[i-1]), false);
      m_op[i-1]->createCoarser(*(m_correction[i]), *(m_correction[i-1]), true);
    }
}

template <class T>
void MultiGrid<T>::solve(T& a_phi, const T& a_rhs, Real a_tolerance, int a_maxIterations, int verbosity)
{
  CH_TIME("MultiGrid::solve");
  this->init(a_phi, a_rhs);

  T correction, residual;
  m_op[0]->create(correction, a_phi);
  m_op[0]->create(residual, a_rhs);
  m_op[0]->setToZero(a_phi);
  m_op[0]->residual(residual, a_phi, a_rhs, false);

  Real errorno = m_op[0]->norm(residual, 0);
  if(verbosity > 2)
    {
      pout() << "multigrid::solve initial residual = " << errorno << std::endl;
    }
  Real compval = a_tolerance*errorno;
  Real epsilon = 1.0e-16;
  compval = Max(compval, epsilon);
  Real error = errorno;
  int iter = 0;
  while((error > compval) && (error > a_tolerance*errorno) && (iter < a_maxIterations))
    {
      m_op[0]->setToZero(correction);
      m_op[0]->residual(residual, a_phi, a_rhs, false);
      error = m_op[0]->norm(residual, 0);
      if(verbosity > 3)
        {
          pout() << "multigrid::solve iter = " << iter <<  ",  residual = " << error << std::endl;
        }

      this->cycle(0, correction, residual);
      m_op[0]->incr(a_phi, correction, 1.0);

      iter++;
    }
  if(verbosity > 2)
    {
      pout() << "multigrid::solve final residual = " << error << std::endl;
    }

  m_op[0]->clear(correction);
  m_op[0]->clear(residual);
}

template <class T>
void MultiGrid<T>::oneCycle(T& a_e, const T& a_residual)
{
  this->init(a_e, a_residual);

  if (m_homogeneous)
    {
      cycle(0, a_e, a_residual);
    }
  else
    {
      T correction, residual;
      m_op[0]->create(correction, a_e);
      m_op[0]->create(residual, a_residual);
      m_op[0]->residual(residual, a_e, a_residual);
      m_op[0]->setToZero(correction);
      this->cycle(0, correction, residual);
      m_op[0]->incr(a_e, correction, 1.0);

      m_op[0]->clear(correction);
      m_op[0]->clear(residual);
    }
}

template <class T>
void MultiGrid<T>::cycle(int depth, T& correction, const T& residual)
{
  // currently I can't drop a timer in this function because it is recursive.  doh

  if(depth == m_depth-1)
    {
      m_bottomSolver->solve(correction, residual);
    }
  else
    {

      int cycles = m_cycle;
      if(depth == m_depth-1) cycles=1;
      for(int c=0; c<cycles; ++c)
        {

          m_op[depth  ]->relax(correction, residual, m_pre);
          m_op[depth  ]->restrictResidual(*(m_residual[depth+1]), correction, residual);
          m_op[depth  ]->setToZero(*(m_correction[depth+1]));

          for(int img = 0; img < m_numMG; img++)
            {
              cycle(depth+1, *(m_correction[depth+1]), *(m_residual[depth+1]));
            }

          m_op[depth  ]->prolongIncrement(correction, *(m_correction[depth+1]));
          m_op[depth  ]->relax(correction, residual, m_post);

        }
    }

}

template <class T>
void MultiGrid<T>::clear()
{
  for(int i=0; i<m_op.size(); ++i)
    {
      delete m_op[i];
      delete m_correction[i];
      delete m_residual[i];
    }
  m_op.resize(0);
  m_residual.resize(0);
  m_correction.resize(0);
  m_defined = false;
}

#include "NamespaceFooter.H"
#endif /*_MULTIGRID_H_*/
