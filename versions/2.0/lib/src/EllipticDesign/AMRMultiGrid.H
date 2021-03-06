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

#ifndef _AMRMULTIGRID_H_
#define _AMRMULTIGRID_H_

#include "MultiGrid.H"
#include "REAL.H"
#include "Box.H"
#include "NoOpSolver.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "Copier.H"
#include "SPMD.H"
#include "Misc.H"

#include "NamespaceHeader.H"

///
/**
   Operator class for AMR Multigrid
 */
template <class T>
class AMRLevelOp : public MGLevelOp<T>
{
public:

  virtual ~AMRLevelOp(){;}

  virtual Real AMRNorm(const T&   a_coarResid,
               const T&   a_fineResid,
               const int& a_refRat,
               const int& a_ord)
  { return this->norm(a_coarResid, 0);} 

  ///
  /**
     return the refinement ratio to next coarser level.
     return 1 when there are no coarser AMRLevelOp objects.
   */
  virtual int refToCoarser() = 0;

  ///
  /**
      a_residual = a_rhs - L(a_phiFine, a_phi, a_phiCoarse)
  */
  virtual void AMRResidual(T& a_residual, const T& a_phiFine, const T& a_phi,
                           const T& a_phiCoarse, const T& a_rhs,
                           bool a_homogeneousDomBC,
                           AMRLevelOp<T>*  a_finerOp) = 0;

  ///
  /** 
      a_residual = a_rhs - L^nf(a_phi, a_phiCoarse)  
      assume no finer AMR level
  */
  virtual void AMRResidualNF(T& a_residual, const T& a_phi, const T& a_phiCoarse,
                             const T& a_rhs, bool a_homogeneousBC) = 0;

  ///
  /** 
      a_residual = a_rhs - L(a_phiFine, a_phi)   
      assume no coarser AMR level
  */
  virtual void AMRResidualNC(T& a_residual, const T& a_phiFine, const T& a_phi,
                             const T& a_rhs, bool a_homogeneousBC,
                             AMRLevelOp<T>* a_finerOp) = 0;



  ///
  /**
     Apply the AMR operator, including coarse-fine matching
  */
  virtual void AMROperator(T& a_LofPhi, 
                           const T& a_phiFine, const T& a_phi,
                           const T& a_phiCoarse, 
                           bool a_homogeneousDomBC,
                           AMRLevelOp<T>*  a_finerOp) = 0;

  ///
  /** 
      Apply the AMR operator, including coarse-fine matching.
      assume no finer AMR level
  */
  virtual void AMROperatorNF(T& a_LofPhi, 
                             const T& a_phi, 
                             const T& a_phiCoarse,
                             bool a_homogeneousBC) = 0;

  ///
  /** 
      Apply the AMR operator, including coarse-fine matching
      assume no coarser AMR level
  */
  virtual void AMROperatorNC(T& a_LofPhi, 
                             const T& a_phiFine, 
                             const T& a_phi,
                             bool a_homogeneousBC,
                             AMRLevelOp<T>* a_finerOp) = 0;


  ///
  /** a_resCoarse = I[h-2h]( a_residual - L(a_correction, a_coarseCorrection)) */
  virtual void AMRRestrict(T& a_resCoarse, const T& a_residual, const T& a_correction,
                           const T& a_coarseCorrection) = 0;



  ///
  /** a_correction += I[2h->h](a_coarseCorrection) */
  virtual void AMRProlong(T& a_correction, const T& a_coarseCorrection) = 0;



  ///
  /** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
  virtual void AMRUpdateResidual(T& a_residual, const T& a_correction,
                                 const T& a_coarseCorrection) = 0;
 

  ///
  /**
   */
  virtual void createCoarsened(T&       a_lhs,
                               const T& a_rhs,
                               const int&     a_refRat) = 0;



  //===================================================================
  // optional optimizations for an AMRLevelOp.  These are not pure virtual
  // functions, since we can build the equivalent algorithmic components from
  // pure virtual functions.  The AMRMultiGrid algorithm actually calls *these*
  // functions, which a smart operator can perform in faster ways.
  //===================================================================


  
  virtual void buildCopier(Copier& a_copier, const T& a_lhs, const T& a_rhs)
  {;}

  virtual void assignCopier(T& a_lhs, const T& a_rhs, const Copier& a_copier)
  {
    this->assign(a_lhs, a_rhs);
  }

  virtual void zeroCovered(T& a_lhs, T& a_rhs, const Copier& a_copier)
  {
    this->setToZero(a_rhs);
    this->assignCopier(a_lhs, a_rhs, a_copier);
  }
  virtual Real localMaxNorm(const T& a_phi) { return this->norm(a_phi, 0);}

    /** optimization of AMRProlong that sends in the existing temporary and copier */
  virtual void AMRProlongS(T& a_correction, const T& a_coarseCorrection,
                           T& a_temp, const Copier& a_copier)
  {
    AMRProlong(a_correction, a_coarseCorrection);
  }
  virtual void AMRRestrictS(T& a_resCoarse, const T& a_residual, const T& a_correction,
                            const T& a_coarseCorrection, T& scratch)
  {
    AMRRestrict(a_resCoarse, a_residual, a_correction, a_coarseCorrection);
  }
};

///
/**
   Factory to create AMRLevelOps
 */
template <class T>
class AMRLevelOpFactory : public MGLevelOpFactory<T>
{
public:
  virtual ~AMRLevelOpFactory(){;}

  ///
  /**
     return a new operator.  this is done with a new call.
     caller is responsible for deletion
   */
  virtual AMRLevelOp<T>* AMRnewOp(const ProblemDomain& a_indexSpace)=0;

  ///
  /**
     return refinement ratio  to next finer level.
   */
  virtual int refToFiner(const ProblemDomain& a_indexSpace) const =0;

};

///
/**
   Class to solve elliptic equations using the Martin and Cartwright algorithm.
 */
template <class T>
class AMRMultiGrid
{
public:

  AMRMultiGrid();
  
  virtual ~AMRMultiGrid();

  ///
  /**
     Define the solver.   
     a_coarseDomain is the index space on the coarsest AMR level.
     a_factory is the operator factory through which all special information is conveyed.
     a_bottomSolver is the solver to be used at the termination of multigrid coarsening.
     a_numLevels is the number of AMR levels.
   */
  virtual void define(const ProblemDomain& a_coarseDomain,
                      AMRLevelOpFactory<T>& a_factory,
                      LinearSolver<T>* a_bottomSolver,
                      int a_numLevels);

  ///
  /**
     Solve L(phi) = rho from l_base to l_max.  To solve over all levels,
     l_base = 0 and l_max = max_level = numLevels-1.
   */
  virtual void solve(Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                     int l_max, int l_base);

  void relaxOnlyHomogeneous(Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                            int l_max, int l_base);
  
  virtual void AMRVCycle(Vector<T*>& a_correction,
                         Vector<T*>& a_residual,
                         int l, int l_max, int l_base);

  void setMGCycle(int a_numMG);

  void init(const Vector<T*>& a_phi, const Vector<T*>& a_rhs,
            int l_max, int l_base);

  Real m_eps, m_hang, m_normThresh;
  int m_imin, m_iterMax, m_verbosity, m_exitStatus;
  int m_pre, m_post, m_bottom, m_numMG;

  AMRLevelOp<T>& levelOp(int level);

  Real computeAMRResidual(Vector<T*>&       a_resid,
                          Vector<T*>&       a_phi,
                          const Vector<T*>& a_rhs,
                          int               l_max, 
                          int               l_base,
                          bool              a_homogeneousBC=false,
                          bool              a_computeNorm=true);

  ///
  /**
     For changing coefficients.  Use at thy own peril.
  */
  Vector< MGLevelOp<T>* > getAllOperators();

  ///
  /**
     Set parameters of the solve. 
     a_pre is the number of smoothings before averaging.
     a_post is the number of smoothings after averaging.
     a_bottom is the number of smoothings at the bottom level.
     a_numMG = 1 for vcycle, =2 for wcycle (use 1).
     a_itermax is the max number of v cycles.
     a_hang is the minimum amount of change per vcycle.
     a_eps is the solution tolerance.
     a_normThresh is how close to zero eps*resid is allowed to get.
   */
  void setSolverParameters(const int&   a_pre,
                           const int&   a_post,
                           const int&   a_bottom,
                           const int&   a_numMG,
                           const int&   a_iterMax,
                           const Real&  a_eps,
                           const Real&  a_hang,
                           const Real&  a_normThresh);

protected:

  void relax(T& phi, T& R, int depth);

  void computeAMRResidualLevel(Vector<T*>&       a_resid,
                               Vector<T*>&       a_phi,
                               const Vector<T*>& a_rhs,
                               int l_max, int l_base, int ilev, 
                               bool a_homogeneousBC);

  Vector<AMRLevelOp<T>*>          m_op;
  Vector<MultiGrid<T> >           m_mg;

  Vector<T*>  m_correction;
  Vector<T*>  m_residual;
  Vector<T*>  m_resC;
  Vector<Copier> m_resCopier;
  Vector<Copier> m_reverseCopier;

  NoOpSolver<T>    m_nosolve;

  LinearSolver<T>* m_bottomSolver;

  void clear();

};

//*******************************************************
// AMRMultigrid Implementation
//*******************************************************

//===================================================================

template <class T>
void
AMRMultiGrid<T>::setSolverParameters(const int&   a_pre,        
                                     const int&   a_post,       
                                     const int&   a_bottom,     
                                     const int&   a_numMG,     
                                     const int&   a_iterMax,
                                     const Real&  a_eps,        
                                     const Real&  a_hang,        
                                     const Real&  a_normThresh)
{
  m_pre        =    a_pre;       
  m_post       =    a_post;      
  m_bottom     =    a_bottom;    
  m_eps        =    a_eps;      
  m_hang       =    a_hang;     
  m_normThresh =    a_normThresh;
  m_iterMax    =    a_iterMax;
  for(int img = 0; img < m_mg.size(); img++)
    {
      m_mg[img].m_pre    = a_pre;
      m_mg[img].m_post   = a_post;
    }
  setMGCycle(a_numMG);
}
template <class T>
Vector< MGLevelOp<T> * > 
AMRMultiGrid<T>::getAllOperators()
{
  Vector< MGLevelOp<T>* > retval;
  for(int iop = 0;  iop < m_op.size(); iop++)
    {
      MGLevelOp<T>* operPtr = (MGLevelOp<T>*) m_op[iop];
      retval.push_back(operPtr);
    }
  
  for(int img = 0; img < m_mg.size(); img++)
    {
      Vector< MGLevelOp<T>* > mgOps = m_mg[img].getAllOperators();
      retval.append(mgOps);
    }
  return retval;
}

template <class T>
AMRMultiGrid<T>::AMRMultiGrid():m_eps(1E-6), m_hang(1E-15), m_normThresh(1E-30),
                                m_imin(5), m_iterMax(20), m_verbosity(3),
                                m_pre(2), m_post(2), m_bottom(2),
                                m_numMG(1),
                                m_bottomSolver(NULL){}

template <class T>
void AMRMultiGrid<T>::setMGCycle(int a_numMG)
{
  for(int ilev=0; ilev<m_op.size(); ilev++)
    {
      m_mg[ilev].m_numMG = a_numMG;
    }
  m_numMG = a_numMG;
}
template <class T>
void AMRMultiGrid<T>::relax(T& a_correction, T& a_residual, int depth)
{
  CH_TIME("AMRMultiGrid::relax");
  m_op[depth]->relax(a_correction, a_residual, m_pre);  //numSmoothDown

  if(m_op[depth]->refToCoarser() > 2)
    {
      int intermediateDepth = 0;
      int r = m_op[depth]->refToCoarser();
      while(r>2)
        {
          r/=2;
          intermediateDepth++;
        }
      int tmp = m_mg[depth].m_depth;
      m_mg[depth].m_depth = intermediateDepth;
      //ok, there is an intermediate multigrid level that is
      // not an AMR level.  We use regular MultiGrid smoothing algorithm for this.
      m_mg[depth].cycle(0, a_correction, a_residual);
      m_mg[depth].m_depth = tmp;
    }
}
/************/
template <class T>
AMRMultiGrid<T>::~AMRMultiGrid()
{
  CH_TIME("~AMRMultiGrid");
  clear();
}

template <class T>
AMRLevelOp<T>&  AMRMultiGrid<T>::levelOp(int level)
{
  return *(m_op[level]);
}

/************/
template <class T>
Real AMRMultiGrid<T>::computeAMRResidual(Vector<T*>&       a_resid,
                                         Vector<T*>&       a_phi,
                                         const Vector<T*>& a_rhs,
                                         int               l_max, 
                                         int               l_base,
                                         bool              a_homogeneousBC,
                                         bool              a_computeNorm)
{
  CH_TIME("AMRMultiGrid::computeAMRResidual");
  Real rnorm = 0;
  Real localNorm = 0;
  for(int ilev=l_base; ilev<= l_max; ilev++)
    {
      //always used at top level where bcs are inhomogeneous
      computeAMRResidualLevel(a_resid,
                              a_phi,
                              a_rhs,
                              l_max, l_base, ilev, a_homogeneousBC);
      if(a_computeNorm){
        if(ilev == l_max)
          {
            localNorm =m_op[ilev]->localMaxNorm(*a_resid[ilev]);
          }
        else
          {
            //m_op[ilev+1]->setToZero(*m_resC[ilev+1]);
            //m_op[ilev]->assignCopier(*a_resid[ilev], *m_resC[ilev+1], m_resCopier[ilev+1]);
            m_op[ilev]->zeroCovered(*a_resid[ilev], *m_resC[ilev+1], m_resCopier[ilev+1]);
            localNorm = m_op[ilev]->localMaxNorm(*a_resid[ilev]);
          }
      }
      rnorm = Max(localNorm, rnorm);
    }
#ifdef CH_MPI
  if(a_computeNorm){
    Real recv;
    int result = MPI_Allreduce(&rnorm, &recv, 1, MPI_CH_REAL,
                               MPI_MAX, Chombo_MPI::comm);
    if(result != MPI_SUCCESS){ //bark!!!
      MayDay::Error("sorry, but I had a communcation error on norm");
    }
    rnorm = recv;
  }
#endif

  return rnorm; // if a_computeNorm is false, then this just returns zero.
}
/************/
template <class T>
void AMRMultiGrid<T>::computeAMRResidualLevel(Vector<T*>&       a_resid,
                                              Vector<T*>&       a_phi,
                                              const Vector<T*>& a_rhs,
                                              int l_max, int l_base, int ilev,
                                              bool a_homogeneousBC)
{
  if(l_max != l_base)
    {
      if(ilev == l_max)
        {
          m_op[l_max]->AMRResidualNF(*(a_resid[l_max]), *(a_phi[l_max]),
                                     *(a_phi[l_max-1]), *(a_rhs[l_max]),
                                     a_homogeneousBC);
        }
      else if(ilev == l_base)
        {
          if(l_base == 0)
            {
              m_op[l_base]->AMRResidualNC(*(a_resid[l_base]), *(a_phi[l_base+1]),
                                          *(a_phi[l_base]),  *(a_rhs[l_base]),
                                           a_homogeneousBC, m_op[l_base+1]);
            }
          else
            {
              m_op[l_base]->AMRResidual(*a_resid[l_base], *a_phi[l_base+1], *a_phi[l_base],
                                        *a_phi[l_base-1], *a_rhs[l_base],
                                         a_homogeneousBC, m_op[l_base+1]);
            }
        }
      else
        {
          m_op[ilev]->AMRResidual(*a_resid[ilev], *a_phi[ilev+1], *a_phi[ilev],
                                  *a_phi[ilev-1], *a_rhs[ilev],
                                   a_homogeneousBC, m_op[ilev+1]);
        }
    }
  else
    {
      CH_assert(ilev == l_base);
      if(l_base == 0)
        {
          m_op[l_max]->residual(*a_resid[l_max], *a_phi[l_max], *a_rhs[l_max],a_homogeneousBC);
        }
      else
        {
          m_op[l_max]->AMRResidualNF(*(a_resid[l_max]), *(a_phi[l_max]),
                                     *(a_phi[l_max-1]), *(a_rhs[l_max]),
                                      a_homogeneousBC);
        }
    }

}

template<class T>
void AMRMultiGrid<T>::solve(Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                            int l_max, int l_base)
{
  CH_TIMERS("AMRMultiGrid::solve");
  CH_TIMER("AMRMultiGrid::AMRVcycle", vtimer);
  init(a_phi, a_rhs, l_max, l_base);

  CH_assert(l_base <= l_max);
  CH_assert(a_rhs.size() == a_phi.size());

  //these correspond to the residual and correction
  //that live in AMRSolver
  Vector<T*> uberCorrection(a_rhs.size());
  Vector<T*> uberResidual(a_rhs.size());

  int lowlim = l_base;
  if( l_base > 0 )
    lowlim--;
  //  for(int ilev=l_base; ilev<= l_max; ilev++)
  for(int ilev=lowlim; ilev<= l_max; ilev++)
    {
      uberCorrection[ilev] = new T();
      uberResidual[ilev] = new T();
      m_op[ilev]->create(*uberCorrection[ilev], *a_phi[ilev]);
      if(ilev >= l_base)
        {
          m_op[ilev]->create(*uberResidual[ilev], *a_rhs[ilev]);
        }
      m_op[ilev]->setToZero(*(uberCorrection[ilev]));
    }

  //compute initial residual and initialize internal residual to it
  
  Real initial_rnorm = 0;
  {
    CH_TIME("Initial AMR Residual");
    initial_rnorm=computeAMRResidual(uberResidual, a_phi, a_rhs, l_max, l_base);
  }
  Real rnorm = initial_rnorm;
  Real norm_last = 2*initial_rnorm;

  /// set bottom solver convergence norm and solver tolerance
  m_bottomSolver->setConvergenceMetrics(initial_rnorm, m_eps);

  int iter=0;
  if(m_verbosity>=4)
    {
      pout() << "    AMRMultiGrid:: iteration = " << iter << ", error norm = " << rnorm << std::endl;
    }
  bool goNorm = rnorm > m_normThresh;                        //iterate if norm is not small enough
  bool goRedu = rnorm > m_eps*initial_rnorm;                 //iterate if initial norm is not reduced enough
  bool goIter = iter < m_iterMax;                            //iterate if iter < max iteration count
  bool goHang = iter < m_imin || rnorm <(1-m_hang)*norm_last;//iterate if we didn't hang
  while (goIter && goRedu && goHang && goNorm)
    {
      norm_last = rnorm;

      //this generates a correction from the current residual
      CH_START(vtimer);
      AMRVCycle(uberCorrection, uberResidual, l_max, l_max, l_base);
      CH_STOP(vtimer);

      //increment phi by correction and reset correction to zero
      for(int ilev=l_base; ilev<= l_max; ilev++)
        {
          m_op[ilev]->incr(*(a_phi[ilev]), *(uberCorrection[ilev]), 1.0);
          m_op[ilev]->setToZero(*(uberCorrection[ilev]));
        }

      //recompute residual
      rnorm = computeAMRResidual(uberResidual, a_phi, a_rhs, l_max, l_base);
      iter++;
      if(m_verbosity>=4)
        {
          pout() << "    AMRMultiGrid:: iteration = " << iter << ", error norm = " << rnorm 
                 << ", rate = " << norm_last/rnorm << std::endl;
        }
      
      goNorm = rnorm > m_normThresh;                        //keep iterating if norm is not small enough
      goRedu = rnorm > m_eps*initial_rnorm;                 //keep iterating if initial norm is not reduced enough
      goIter = iter < m_iterMax;                            //keep iterating if iter < max iteration count
      goHang = iter < m_imin || rnorm <(1-m_hang)*norm_last;//keep iterating if we didn't hang
    }
  m_exitStatus = int(!goRedu) + int(!goIter)*2 + int(!goHang)*4 + int(!goNorm)*8;
  if(m_verbosity==3)
    {
      pout() << "    AMRMultiGrid:: iteration = " << iter << ", error norm = " << rnorm << std::endl;
    }
  if(m_verbosity>1)
    {
      if(!goIter && goRedu && goNorm)
        {
          pout() << "    AMRMultiGrid:: WARNING: Exit because max iteration count exceeded" << std::endl;
        }
      if(!goHang && goRedu && goNorm)
        {
          pout() << "    AMRMultiGrid:: WARNING: Exit because of solver hang" << std::endl;
        }
      if(m_verbosity>4)
        {
          pout() << "    AMRMultiGrid:: exitStatus = " << m_exitStatus << std::endl;
        }
    }
  for(int i=l_base; i <= l_max; i++)
    {
      m_op[i]->clear(*uberCorrection[i]);
      m_op[i]->clear(*uberResidual[i]);
      delete uberCorrection[i];
      delete uberResidual[i];
    }
}


template<class T>
void AMRMultiGrid<T>::relaxOnlyHomogeneous(Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                                           int l_max, int l_base)
{
  CH_TIME("AMRMultiGrid::relaxOnly");
  CH_assert(l_max == l_base);
  CH_assert(l_max == 0);
  init(a_phi, a_rhs, l_max, l_base);

  CH_assert(a_rhs.size() == a_phi.size());

  Vector<T*> uberResidual(a_rhs.size());

  for(int ilev=0; ilev<m_op.size(); ilev++)
    {
      uberResidual[ilev] = new T();
      m_op[ilev]->create(*uberResidual[ilev], *a_rhs[ilev]);
    }

  //compute initial residual and initialize internal residual to it
  m_op[0]->residual(*uberResidual[0], *a_phi[0], *a_rhs[0], true);  
  Real initial_rnorm = m_op[0]->norm(*uberResidual[0], 0);
  Real rnorm = initial_rnorm;
  Real norm_last = 2*initial_rnorm;


  int iter=0;
  if(m_verbosity>=3)
    {
      pout() << "    AMRMultiGrid::relaxOnly iteration = " << iter << ", error norm = " << rnorm  << std::endl;
    }
  for(iter = 1; iter < 40; iter++)
    {

      norm_last = rnorm;
      m_op[0]->relax(*a_phi[0], *a_rhs[0], 100);  

      //recompute residual
      m_op[0]->residual(*uberResidual[0], *a_phi[0], *a_rhs[0], true);  
      rnorm = m_op[0]->norm(*uberResidual[0], 2);
      if(m_verbosity>=4)
        {
          pout() << "    AMRMultiGrid::relaxOnly iteration = " << iter << ", error norm = " << rnorm 
                 << ", rate = " << norm_last/rnorm << std::endl;
        }
      
    }
  m_exitStatus = 0;
  for(int i=0; i<m_op.size(); i++)
    {
      m_op[i]->clear(*uberResidual[i]);
      delete uberResidual[i];
    }
}

template <class T>
void AMRMultiGrid<T>::clear()
{
  for(int i=0; i<m_op.size(); i++)
    {
      m_op[i]->clear(*m_correction[i]);
      m_op[i]->clear(*m_residual[i]);
      m_op[i]->clear(*m_resC[i]);

      delete m_correction[i];
      delete m_residual[i];
      delete m_resC[i];
      delete m_op[i];
    }
}

template <class T>
void AMRMultiGrid<T>::init(const Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                           int l_max, int l_base)
{
  CH_TIME("AMRMultiGrid::init");
  CH_assert(a_phi.size()<=m_op.size());
  for(int i=l_base; i<= l_max; i++)
    {
      AMRLevelOp<T>& op = *(m_op[i]);
  
      op.create(*m_correction[i], *a_phi[i]);
      op.create(*m_residual[i],   *a_rhs[i]);

      m_mg[i].init(*a_phi[i],     *a_rhs[i]);

      int r = op.refToCoarser();
      op.createCoarsened(*m_resC[i], *a_rhs[i], r);
      if(i!=l_base) {
        op.buildCopier(m_resCopier[i],  *a_rhs[i-1], *m_resC[i]);
        m_reverseCopier[i] = m_resCopier[i];
        m_reverseCopier[i].reverse();
      }
    }

  //  for(int ilev=0; ilev<=m_mg.size(); ilev++)
  for(int ilev=l_base; ilev<=l_max; ilev++)
    {
      m_mg[ilev].m_pre=1;
      m_mg[ilev].m_post=1;
      m_mg[ilev].m_bottomSolver = &m_nosolve;
    }
  
  m_mg[l_base].m_pre  = m_pre;
  m_mg[l_base].m_post = m_post;
  m_mg[l_base].setBottomSolver(m_bottomSolver);

}

template <class T>
void AMRMultiGrid<T>::define(const ProblemDomain& a_coarseDomain,
                             AMRLevelOpFactory<T>& a_factory,
                             LinearSolver<T>* a_bottomSolver,
                             int a_maxAMRLevels)
{
  CH_TIME("AMRMultiGrid::define");
  this->clear();
  m_op.        resize( a_maxAMRLevels, NULL);
  m_mg.        resize( a_maxAMRLevels);

 
  m_correction.resize( a_maxAMRLevels, NULL);
  m_residual.  resize( a_maxAMRLevels, NULL);
  m_resC.      resize( a_maxAMRLevels, NULL);
  m_resCopier. resize( a_maxAMRLevels);
  m_reverseCopier.resize(  a_maxAMRLevels);
  m_bottomSolver = a_bottomSolver;

  ProblemDomain current = a_coarseDomain;
  for(int i=0; i<a_maxAMRLevels; i++)
    {
 
      m_correction[i] = new T();
      m_residual[i]   = new T();
      m_resC[i]       = new T();
      m_op[i] = a_factory.AMRnewOp(current);
      int maxDepth = -1; //coarsen all the way for now
      m_mg[i].define(a_factory, &m_nosolve, current, maxDepth);
      current.refine(a_factory.refToFiner(current));
    }
}

template<class T>
void AMRMultiGrid<T>::AMRVCycle(Vector<T*>& a_uberCorrection,
                                Vector<T*>& a_uberResidual,
                                int ilev, int l_max, int l_base)
{
  if (ilev == l_max)
    {
      for(int level=l_base; level<= l_max; level++)
        {
          m_op[level]->assignLocal(*m_residual[level], *a_uberResidual[level]);
          m_op[level]->setToZero(*m_correction[level]);
        }
    }
  
  if(l_max == l_base)
    {
      CH_assert(ilev == l_base);
      m_mg[l_base].oneCycle(*(a_uberCorrection[ilev]), *(a_uberResidual[ilev]));
    }
  else if(ilev == l_base)
    {
      m_mg[l_base].oneCycle(*(m_correction[ilev]), *(m_residual[ilev]));

      m_op[ilev]->incr(*(a_uberCorrection[ilev]), *(m_correction[ilev]), 1.0);
    }
  else
    {

      //downsweep
      // in AMRLevelMG, it is a coarsened version of this level
      //but here we need to make it the same as the coarser
      //level so I will use increment instead of copyTo

      this->relax(*(m_correction[ilev]), *(m_residual[ilev]), ilev);
      m_op[ilev]->incr(*(a_uberCorrection[ilev]), *(m_correction[ilev]), 1.0);

      //set next level correction to zero
      m_op[ilev-1]->setToZero(*(m_correction[ilev-1]));

      //==== Compute residual for the next coarser level.=====

      //recompute residual for next coarser level
      //need homogeneous BCs because in residual-correction form
      computeAMRResidualLevel(m_residual,
                              a_uberCorrection,
                              a_uberResidual,
                              l_max, l_base, ilev-1,
                              true);

      // compute the restriction of the residual to the coarser level resC.
      m_op[ilev]->AMRRestrictS(*(m_resC[ilev]),
                               *(m_residual[ilev]),
                               *(m_correction[ilev]),
                               *(m_correction[ilev-1]),
                               *(a_uberCorrection[ilev]));

 
      //overwrite next coarser's residual with coarsened residual from
      //this level
      m_op[ilev-1]->assignCopier(*m_residual[ilev-1], *(m_resC[ilev]), m_resCopier[ilev]);

      //============finish Compute residual fo rthe next coarser level======

      for(int img = 0; img < m_numMG; img++)
        {
          AMRVCycle(a_uberCorrection, a_uberResidual, ilev-1, l_max, l_base);
        }

      //upsweep
      //increment the correction with coarser version
      //m_op[ilev]->AMRProlong(*(m_correction[ilev]), *(m_correction[ilev-1]));
      m_op[ilev]->AMRProlongS(*(m_correction[ilev]), *(m_correction[ilev-1]), 
                              *m_resC[ilev], m_reverseCopier[ilev]);
      //recompute residual
      m_op[ilev]->AMRUpdateResidual(*(m_residual[ilev]), *(m_correction[ilev]), *(m_correction[ilev-1]));
      

      //compute correction to the correction
      T& dCorr = *(a_uberCorrection[ilev]); // user uberCorrection as holder for correction to correction
      m_op[ilev]->setToZero(dCorr);
      this->relax(dCorr, *(m_residual[ilev]), ilev);

      //correct the correction with the correction to the correction
      m_op[ilev]->incr(*(m_correction[ilev]), dCorr, 1.0);

 
      m_op[ilev]->assignLocal(*(a_uberCorrection[ilev]), *(m_correction[ilev]));

    }
}

#include "NamespaceFooter.H"
#endif
