#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// LevelNodeSolver.cpp
// adapted from LevelSolver by DTGraves, Fri, July 23, 1999
// petermc, Fri, Jan 26, 2001

#include <cmath>
#include <iostream>
using std::cout;
using std::endl;

#include "LevelNodeSolver.H"
#include "IntVectSet.H"
#include "NodeSetOperations.H"
#include "LayoutIterator.H"
#include "Norms.H"

// ---------------------------------------------------------
// default constructor
LevelNodeSolver::LevelNodeSolver()
{
  setDefaultValues();
}


// ---------------------------------------------------------
void
LevelNodeSolver::define(const DisjointBoxLayout& a_grids,
                        const DisjointBoxLayout* a_gridsCoarsePtr,
                        const Box& a_domain,
                        Real a_dx,
                        int a_refToCoarse,
                        const NodeLevelOp* const a_opin,
                        int a_minLength)
{
  ProblemDomain probdomain(a_domain);
  define(a_grids, a_gridsCoarsePtr, probdomain,
         a_dx, a_refToCoarse, a_opin, a_minLength);
}


// ---------------------------------------------------------
void
LevelNodeSolver::define(const DisjointBoxLayout& a_grids,
                        const DisjointBoxLayout* a_gridsCoarsePtr,
                        const ProblemDomain& a_domain,
                        Real a_dx,
                        int a_refToCoarse,
                        const NodeLevelOp* const a_opin,
                        int a_minLength)
{
  clearMemory();
  m_grids = a_grids;
  m_domain = a_domain;
  m_refToCoarse = a_refToCoarse;
  m_dx = a_dx;
  m_resid.define(m_grids, 1, IntVect::Zero);
  m_scratch.define(m_grids, 1, IntVect::Zero);
  m_corr.define(m_grids, 1, IntVect::Unit);
  m_isDefined = true;

  // added 27 March 2003, for use with norms
  interiorBoundaryNodes(m_IVSV, m_grids, m_domain);
  exteriorBoundaryNodes(m_IVSVext, m_IVSV, m_grids);

  m_levelOpPtr = a_opin->new_levelop();
  bool homogeneousOnly = false;
  int ncomp = 1;
  m_levelOpPtr->define(m_grids, a_gridsCoarsePtr, m_dx,
                       m_refToCoarse, m_domain,
                       homogeneousOnly, ncomp, m_verbose);

 CH_assert(a_minLength >= 0);

  int nCoarserLevels = 0;
  // a_minLength == 0 means do not coarsen.

  if (a_minLength > 0)
    nCoarserLevels = countCoarserLevels(a_minLength);

  m_levelMG.define(m_grids, a_gridsCoarsePtr, m_domain, m_dx,
                   m_refToCoarse, a_opin, nCoarserLevels);
  m_levelMG.setVerbose(m_verbose);
}


// ---------------------------------------------------------
int
LevelNodeSolver::countCoarserLevels(int a_minLength)
{
  // Coarsen boxes once to start, because coarsest level boxes
  // must also be coarsenable (for proper coarse-fine interpolation).
  // However, if grid merely covers entire physical domain, then
  // want to coarsen all the way down.

  int numCoarser = -1;

  // If m_grids cover all of m_domain, then set numCoarser = 0.
  // Else retain numCoarser = -1.
  IntVectSet ivsCheck(m_domain.domainBox());
  LayoutIterator lit = m_grids.layoutIterator();
  for (lit.reset(); lit.ok(); ++lit)
    ivsCheck -= m_grids.get(lit());
  if (ivsCheck.isEmpty()) numCoarser++;

  // Find the largest power of 2 that divides the lengths of EVERY grid
  // at this level.
  // Increment numCoarser by the exponent.

  int refRatioTotal = 2 * a_minLength;
  while (true)
    {
      bool addOne = true;
      for (lit.reset(); lit.ok(); ++lit)
        {
          Box grid = m_grids.get(lit());
          // check whether grid is refinable by refRatioTotal
          if (refine(coarsen(grid, refRatioTotal), refRatioTotal) != grid)
            addOne = false;
          if (!addOne) break;
        }
      if (!addOne) break;
      // succeeded; then try refining even more.
      numCoarser++;
      refRatioTotal *= 2;
    }
  numCoarser = Max(numCoarser, 0);
  return numCoarser;
}

// ---------------------------------------------------------
// complete constructor
LevelNodeSolver::LevelNodeSolver(const DisjointBoxLayout& a_grids,
                                 const DisjointBoxLayout* a_gridsCoarsePtr,
                                 const Box& a_domain,
                                 Real a_dx,
                                 int a_refToCoarse,
                                 const NodeLevelOp* const a_opin,
                                 int a_minLength)
{
  ProblemDomain probdomain(a_domain);
  setDefaultValues();
  define(a_grids, a_gridsCoarsePtr, probdomain, a_dx,
         a_refToCoarse, a_opin, a_minLength);
}


// ---------------------------------------------------------
// complete constructor
LevelNodeSolver::LevelNodeSolver(const DisjointBoxLayout& a_grids,
                                 const DisjointBoxLayout* a_gridsCoarsePtr,
                                 const ProblemDomain& a_domain,
                                 Real a_dx,
                                 int a_refToCoarse,
                                 const NodeLevelOp* const a_opin,
                                 int a_minLength)
{
  setDefaultValues();
  define(a_grids, a_gridsCoarsePtr, a_domain, a_dx,
         a_refToCoarse, a_opin, a_minLength);
}


// ---------------------------------------------------------
// destructor
LevelNodeSolver::~LevelNodeSolver()
{
  clearMemory();
}


// ---------------------------------------------------------
bool
LevelNodeSolver::isDefined() const
{
  return m_isDefined;
}


// ---------------------------------------------------------
void
LevelNodeSolver::setDefaultValues()
{
  m_levelOpPtr = NULL;
  m_isDefined = false;
  m_refToCoarse = -1;
  m_dx = -1.0;
  m_maxIter = 33;
  m_minIter = 4;
  m_tolerance = 1.0e-10;
  m_operatorTolerance = 1.0e-4;
  m_bottomSolveFlag = true;
#ifdef CH_USE_FLOAT
  m_tolerance = sqrt(m_tolerance);
#endif
  m_verbose = false;
}


// ---------------------------------------------------------
void
LevelNodeSolver::setnumSmoothUp(int a_numSmoothUp)
{
  m_levelMG.setnumSmoothUp(a_numSmoothUp);
}


// ---------------------------------------------------------
void
LevelNodeSolver::setnumSmoothDown(int a_numSmoothDown)
{
  m_levelMG.setnumSmoothDown(a_numSmoothDown);
}


// ---------------------------------------------------------
void
LevelNodeSolver::setBottomSmoothing(bool a_bottomSolveFlag)
{
  m_bottomSolveFlag = a_bottomSolveFlag;
}


// ---------------------------------------------------------
void
LevelNodeSolver::setTolerance(Real a_tolerance)
{
  m_tolerance = a_tolerance;
}


// ---------------------------------------------------------
void
LevelNodeSolver::setOperatorTolerance(Real a_tolerance)
{
 CH_assert (a_tolerance >= 0.0);
  m_operatorTolerance = a_tolerance;
}


// ---------------------------------------------------------
void
LevelNodeSolver::setVerbose(bool a_verbose)
{
  m_verbose = a_verbose;
  if (isDefined())
    m_levelMG.setVerbose(a_verbose);
  if (m_levelOpPtr != NULL)
    m_levelOpPtr->setVerbose(a_verbose);
}


// ---------------------------------------------------------
void
LevelNodeSolver::setMaxIter(int a_maxIter)
{
  m_maxIter = a_maxIter;
}


// ---------------------------------------------------------
void
LevelNodeSolver::setMinIter(int a_minIter)
{
  m_minIter = a_minIter;
}


// ---------------------------------------------------------
// returns LevelNodeSolver to a basically undefined state
void
LevelNodeSolver::clearMemory()
{
  if(m_levelOpPtr != NULL)
    delete m_levelOpPtr;
  m_levelOpPtr = NULL;
}


// ---------------------------------------------------------
// solve on just this level, inhomogeneous bcs
void
LevelNodeSolver::levelSolve(LevelData<NodeFArrayBox>& a_phi,
                            const LevelData<NodeFArrayBox>* a_phiCoarse,
                            const LevelData<NodeFArrayBox>& a_rhs,
                            bool a_initializePhiToZero)
{
 CH_assert(isDefined());
 CH_assert(a_phi.nComp() == 1);
 CH_assert(a_rhs.nComp() == 1);

  /* compute initial residual
     -- can't use residual() fn because don't want to zero out higher levels.
     -- also note that we use no-fine operator, because this is a LEVEL solve.
  */

 CH_assert(m_levelOpPtr->isDefined());
  DataIterator dit = a_phi.dataIterator();

  if (a_initializePhiToZero)
    {
      //initial guess of phi is zero
      for (dit.reset(); dit.ok(); ++dit)
        a_phi[dit()].getFab().setVal(0.0);
    }

  // set m_resid = rhs - operator(a_phi)
  // with the two-level inhomogeneous operator
  m_levelOpPtr->residualI(m_resid, a_phi, a_phiCoarse, a_rhs);

  // set m_corr to zero.
  for (dit.reset(); dit.ok(); ++dit)
    m_corr[dit()].getFab().setVal(0.0);

  // Real initRes = norm(m_resid, m_dx, 0, m_resid.interval(), m_verbose);
  Real initRes = norm(m_resid, m_domain, m_IVSVext,
                      m_dx, m_resid.interval(), 0, m_verbose);

  // if initRes = 0, don't bother solving
  if (initRes < 0.01*m_tolerance)
    return;

  bool done = false;
  Real oldRes = initRes;
  int iter = 0;
  while (!done)
    {
      m_levelMG.mgRelax(m_corr, m_resid, m_bottomSolveFlag);

      // m_scratch = m_resid - operator(m_corr)
      // with the homogeneous two-level operator
      m_levelOpPtr->residualH(m_scratch, m_corr, m_resid);

      //      Real currentRes = norm(m_scratch, m_dx, 0, m_scratch.interval(), m_verbose);
      Real currentRes = norm(m_scratch, m_domain, m_IVSVext,
                             m_dx, m_scratch.interval(), 0, m_verbose);

      iter++;

      if (m_verbose)
        cout << "LevelNodeSolver::levelSolve iteration #= " << iter
             << ": Max(res) = " << currentRes
             << endl;

      // done = (currentRes <= m_tolerance*initRes || iter > m_maxIter);
      // petermc, 4 Oct 2002, added third test is to catch case where
      // convergence hangs due to machine precision (or
      // solvability) issues.  This follows the approach
      // used in AMRSolver
      done = ((iter > m_maxIter)
              || (currentRes <= m_tolerance * initRes)
              || ((iter > m_minIter)
                  && ((currentRes/oldRes) > 1.0-m_operatorTolerance)));

      oldRes = currentRes;

    }  // end solver iterations

  // update phi
  for(dit.reset(); dit.ok(); ++dit)
    a_phi[dit()].getFab() += m_corr[dit()].getFab();
}


// ---------------------------------------------------------
// solve on just this level, homogeneous bcs
void
LevelNodeSolver::levelSolveH(LevelData<NodeFArrayBox>& a_phi,
                             const LevelData<NodeFArrayBox>& a_rhs,
                             bool a_initializePhiToZero)
{
 CH_assert(isDefined());
 CH_assert(a_phi.nComp() == 1);
 CH_assert(a_rhs.nComp() == 1);

  /* compute initial residual
     -- can't use residual() fn because don't want to zero out higher levels.
     -- also note that we use no-fine operator, because this is a LEVEL solve.
  */

 CH_assert(m_levelOpPtr->isDefined());

  DataIterator dit = a_phi.dataIterator();

  if (a_initializePhiToZero) {
    //initial guess of phi is zero
    for (dit.reset(); dit.ok(); ++dit)
      a_phi[dit()].getFab().setVal(0.0);
  }

  // m_resid = a_rhs - operator(a_phi)
  // with the two-level homogeneous operator
  m_levelOpPtr->residualH(m_resid, a_phi, a_rhs);

  for (dit.reset(); dit.ok(); ++dit)
    m_corr[dit()].getFab().setVal(0.0);

  //   Real initRes = norm(m_resid, m_dx, 0, m_resid.interval(), m_verbose);
  Real initRes = norm(m_resid, m_domain, m_IVSVext,
                      m_dx, m_resid.interval(), 0, m_verbose);

  // if initRes = 0, don't bother solving
  if (initRes < 0.01*m_tolerance)
    return;

  bool done = false;
  Real oldRes = initRes;
  int iter = 0;
  while (!done)
    {
      // compute modified residual LofPhi = res - LNf(corr)
      // use homogeneous form of ApplyOp because want to use
      // homogeneous CF BC's for correction
      m_levelMG.mgRelax(m_corr, m_resid, m_bottomSolveFlag);

      // set m_scratch = m_resid - operator(m_corr)
      // with the two-level homogeneous operator
      m_levelOpPtr->residualH(m_scratch, m_corr, m_resid);

      // Real currentRes = norm(m_scratch, m_dx, 0, m_scratch.interval(), m_verbose);
      Real currentRes = norm(m_scratch, m_domain, m_IVSVext,
                             m_dx, m_scratch.interval(), 0, m_verbose);

      iter++;

      if (m_verbose)
        cout << "LevelNodeSolver::levelSolveH iteration #= " << iter
             << ": Max(res) = " << currentRes
             << endl;

      // done = (currentRes <= m_tolerance*initRes || iter > m_maxIter);
      // petermc, 4 Oct 2002, added third test is to catch case where
      // convergence hangs due to machine precision (or
      // solvability) issues.  This follows the approach
      // used in AMRSolver
      done = ((iter > m_maxIter)
              || (currentRes <= m_tolerance * initRes)
              || ((iter > m_minIter)
                  && ((currentRes/oldRes) > 1.0-m_operatorTolerance)));

      oldRes = currentRes;

    }  // end solver iterations

  // update phi
  for(dit.reset(); dit.ok(); ++dit)
    a_phi[dit()].getFab() += m_corr[dit()].getFab();
}
