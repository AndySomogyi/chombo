#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// AMRNodeSolver.cpp
// adapted from AMRSolver by DTGraves, Tues, July 6, 1999
// petermc, Wed, Nov 22, 2000

#include <cmath>
#include <iostream>
using std::cout;
using std::endl;

#include "AMRNodeSolver.H"
#include "LayoutIterator.H"
// #include "DatasetClient.H"


// ---------------------------------------------------------
void
AMRNodeSolver::setDefaultValues()
{
  m_numLevels = -1;
  m_finestLevel = -1;
  m_refRatio.resize(0);
  m_tolerance = 1.0e-10;
#ifdef CH_USE_FLOAT
  m_tolerance = 10.*sqrt(m_tolerance);
#endif
  m_operatorTolerance = 1.0e-4;
  m_maxIter = 42;
  m_minIter = 5;
  m_numSmoothUp = 4;
  m_numSmoothDown = 4;
  m_isDefined = false;
  m_verbose = true;
}


// ---------------------------------------------------------
AMRNodeSolver::AMRNodeSolver() : m_amrmgLevel()
{
  setDefaultValues();
}


// ---------------------------------------------------------
void
AMRNodeSolver::clear()
{
  for (int ilev = 0; ilev < m_amrmgLevel.size(); ilev++)
    {
      if (m_amrmgLevel[ilev] != NULL)
        delete m_amrmgLevel[ilev];
    }
}


// ---------------------------------------------------------
AMRNodeSolver::~AMRNodeSolver()
{
  clear();
}


// ---------------------------------------------------------
AMRNodeSolver::AMRNodeSolver(const Vector<DisjointBoxLayout>& a_gridsLevel,
                             const Vector<Box>& a_domainLevel,
                             const Vector<Real>& a_dxLevel,
                             const Vector<int>& a_refRatio,
                             int a_numLevels,
                             int a_lBase,
                             const NodeLevelOp* const a_opin,
                             int a_minLength) : m_amrmgLevel()
{
  setDefaultValues();
  Vector<ProblemDomain> probdomainLevel(a_numLevels);
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    probdomainLevel[ilev] = ProblemDomain(a_domainLevel[ilev]);

  define(a_gridsLevel, probdomainLevel, a_dxLevel,
         a_refRatio, a_numLevels, a_lBase, a_opin, a_minLength);
}


// ---------------------------------------------------------
AMRNodeSolver::AMRNodeSolver(const Vector<DisjointBoxLayout>& a_gridsLevel,
                             const Vector<ProblemDomain>& a_domainLevel,
                             const Vector<Real>& a_dxLevel,
                             const Vector<int>& a_refRatio,
                             int a_numLevels,
                             int a_lBase,
                             const NodeLevelOp* const a_opin,
                             int a_minLength) : m_amrmgLevel()
{
  setDefaultValues();
  define(a_gridsLevel, a_domainLevel, a_dxLevel,
         a_refRatio, a_numLevels, a_lBase, a_opin, a_minLength);
}


// ---------------------------------------------------------
void
AMRNodeSolver::define(const Vector<DisjointBoxLayout>& a_gridsLevel,
                      const Vector<Box>& a_domainLevel,
                      const Vector<Real>& a_dxLevel,
                      const Vector<int>& a_refRatio,
                      int a_numLevels,
                      int a_lBase,
                      const NodeLevelOp* const a_opin,
                      int a_minLength)
{
  Vector<ProblemDomain> probdomainLevel(a_numLevels);
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    probdomainLevel[ilev] = ProblemDomain(a_domainLevel[ilev]);

  define(a_gridsLevel, probdomainLevel, a_dxLevel,
         a_refRatio, a_numLevels, a_lBase, a_opin, a_minLength);
}


// ---------------------------------------------------------
void
AMRNodeSolver::define(const Vector<DisjointBoxLayout>& a_gridsLevel,
                      const Vector<ProblemDomain>& a_domainLevel,
                      const Vector<Real>& a_dxLevel,
                      const Vector<int>& a_refRatio,
                      int a_numLevels,
                      int a_lBase,
                      const NodeLevelOp* const a_opin,
                      int a_minLength)
{
  clear();
  m_isDefined = true;

  m_gridsLevel = a_gridsLevel;
  m_domainLevel = a_domainLevel;
  m_numLevels = a_numLevels;
  m_finestLevel = a_numLevels - 1;
  m_refRatio = a_refRatio;
  m_dxLevel = a_dxLevel;
  m_lBase = a_lBase;

  // for each level, check that all grids are in the physical domain.
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      const Box& domainbox = m_domainLevel[ilev].domainBox();
      const DisjointBoxLayout& grids = m_gridsLevel[ilev];
      for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
        CH_assert(domainbox.contains(grids[dit()]));
    }

  // define an AMRNodeLevelMG object for each level.
  m_amrmgLevel.resize(m_numLevels, NULL);
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      m_amrmgLevel[ilev] = new AMRNodeLevelMG;
      m_amrmgLevel[ilev]->setVerbose(m_verbose);
      m_amrmgLevel[ilev]->define(this, ilev, a_opin);
      m_amrmgLevel[ilev]->setnumSmoothUp(m_numSmoothUp);
      m_amrmgLevel[ilev]->setnumSmoothDown(m_numSmoothDown);
    }
  CH_assert(m_lBase >= 0);

  // define a LevelNodeSolver object for the base level.
  // the grids remain cell-centered, so this part should stay the same.
  const DisjointBoxLayout& levsolvGrids= m_gridsLevel[m_lBase];
  const ProblemDomain&     levsolvDom  = m_domainLevel[m_lBase];
  const Real&              levsolvDx   = m_dxLevel[m_lBase];
  const DisjointBoxLayout* levsolvBase = NULL;
  int levsolvRef = 2;
  if (m_lBase > 0)
    {
      levsolvBase = &m_gridsLevel[m_lBase-1];
      levsolvRef = m_refRatio[m_lBase-1];
    }
  m_levelSolver.define(levsolvGrids, levsolvBase,
                       levsolvDom,   levsolvDx,
                       levsolvRef, a_opin, a_minLength);
  m_levelSolver.setVerbose(m_verbose);
  m_levelSolver.setnumSmoothUp(m_numSmoothUp);
  m_levelSolver.setnumSmoothDown(m_numSmoothDown);
}


// ---------------------------------------------------------
bool
AMRNodeSolver::isDefined() const
{
  return m_isDefined;
}


// ---------------------------------------------------------
void
AMRNodeSolver::setNumSmoothUp(int a_numSmoothUp)
{
  CH_assert(isDefined());
  CH_assert(a_numSmoothUp >= 0);
  m_numSmoothUp = a_numSmoothUp;
  m_levelSolver.setnumSmoothUp(m_numSmoothUp);
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    m_amrmgLevel[ilev]->setnumSmoothUp(m_numSmoothUp);
}


// ---------------------------------------------------------
void
AMRNodeSolver::setNumSmoothDown(int a_numSmoothDown)
{
  CH_assert(isDefined());
  CH_assert(a_numSmoothDown >= 0);
  m_numSmoothDown = a_numSmoothDown;
  m_levelSolver.setnumSmoothDown(m_numSmoothDown);
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    m_amrmgLevel[ilev]->setnumSmoothDown(m_numSmoothDown);
}


// ---------------------------------------------------------
void
AMRNodeSolver::setTolerance(Real a_tolerance)
{
  CH_assert(a_tolerance >= 0);
  m_tolerance = a_tolerance;
}


// ---------------------------------------------------------
void
AMRNodeSolver::setBottomSmoothing(bool a_doBottomSmooth)
{
  m_levelSolver.setBottomSmoothing(a_doBottomSmooth);
}


// ---------------------------------------------------------
void
AMRNodeSolver::setBottomTolerance(Real a_tolerance)
{
  CH_assert(a_tolerance >= 0);
  m_levelSolver.setTolerance(a_tolerance);
}


// ---------------------------------------------------------
void
AMRNodeSolver::setOperatorTolerance(Real a_tolerance)
{
  CH_assert(a_tolerance >= 0);
  m_operatorTolerance = a_tolerance;
}


// ---------------------------------------------------------
void
AMRNodeSolver::setMaxIter(int a_maxIter)
{
  CH_assert(a_maxIter >= 0);
  m_maxIter = a_maxIter;
}


// ---------------------------------------------------------
void
AMRNodeSolver::setMinIter(int a_minIter)
{
  CH_assert(a_minIter >= 0);
  m_minIter = a_minIter;
}


// ---------------------------------------------------------
void
AMRNodeSolver::setBottomMaxIter(int a_maxIter)
{
  CH_assert(isDefined());
  CH_assert(a_maxIter >= 0);
  m_levelSolver.setMaxIter(a_maxIter);
}


// ---------------------------------------------------------
void
AMRNodeSolver::setVerbose(bool a_verbose)
{
  m_verbose = a_verbose;
  if (isDefined())
    {
      m_levelSolver.setVerbose(a_verbose);
      for (int ilev = 0; ilev < m_numLevels; ilev++)
        m_amrmgLevel[ilev]->setVerbose(a_verbose);
    }
}


// ---------------------------------------------------------
void
AMRNodeSolver::solveAMR(Vector<LevelData<NodeFArrayBox> *>& a_phiLevel,
                        const Vector<LevelData<NodeFArrayBox> *>& a_rhsLevel)
{
  CH_assert(isDefined());
  CH_assert(a_phiLevel.size() > m_finestLevel);
  CH_assert(a_rhsLevel.size() > m_finestLevel);

  // As initial guess, set a_phiLevel == 0 at all levels.
  for (int ilev = m_lBase; ilev <= m_finestLevel; ilev++)
    {
      LevelData<NodeFArrayBox>& phiLev = *a_phiLevel[ilev];
      for (DataIterator dit = phiLev.dataIterator(); dit.ok(); ++dit)
        phiLev[dit()].getFab().setVal(0.);
    }

  // First compute the residual itself.  Added by petermc, 11 Jun 2002
  for (int ilev = m_finestLevel; ilev >= m_lBase; ilev--)
    m_amrmgLevel[ilev]->computeAMRResidual(a_phiLevel, a_rhsLevel);

  // Removed by petermc, 2 Aug 2002, because you don't have
  // the residual calculated yet.
  // Real initRes = computeResidualNorm(a_phiLevel, a_rhsLevel, 0);

  // petermc, 10 apr 2001; commented out 5 Jun 2003
  // to deal with huge initial residuals, as from inhomo bcs
  // (where boundary is nonzero and interior is zero)
  // Real initRHSmax = 0.;
  // for (int ilev = m_finestLevel; ilev >= m_lBase; ilev--)
  // {
  // Real rhsMaxLevel = m_amrmgLevel[ilev]->computeNorm(*a_rhsLevel[ilev], 0);
  // if (rhsMaxLevel > initRHSmax) initRHSmax = rhsMaxLevel;
  // }
  // if (initRHSmax > 0.)
  // initRes = Min(initRes, initRHSmax);

  Real initRes = computeResidualNorm(0); //initRHSmax;
  Real currentRes = initRes;
  Real oldRes = currentRes;

  bool done = (initRes == 0.0);
  // If initRes == 0, don't bother solving...
  if (done)
    {
      if (m_verbose)
        {
          cout << "AMRNodeSolver::solveAMR(): initial residual == 0" << endl;
          cout << "returning..." << endl;
        }
      return;
    }

  if (m_verbose)
    {
      cout << "AMRNodeSolver, levels " << m_lBase
           << " to " << m_finestLevel
           << " initial max(residual) = " << initRes << endl;
    }
  oldRes = initRes; //  * 2.; // kludge to get through test below

  // BEGIN MAIN LOOP

  int iter = 0;
  while (!done)
    {
      // residuals at all levels get calculated in AMRVCycleMG
      AMRVCycleMG(a_phiLevel, a_rhsLevel);

      currentRes = computeResidualNorm(0);

      iter++;

      // DFM (3/24/2000) second test is in case where convergence
      // hangs due to machine precision (or solvability) issues...
      done = ((iter >= m_maxIter)
              || ( (iter >= m_minIter)
                   && (currentRes/oldRes) > 1.0-m_operatorTolerance)
              || (currentRes <= m_tolerance*initRes) );

      if (m_verbose)
        cout << "AMRNodeSolver iteration #  " << iter
             << ": NORM Max(res) = " << currentRes << endl;

      oldRes = currentRes;
    } // end while (!done)

  // END MAIN LOOP

  if (currentRes <= m_tolerance*initRes) // convergence succeeded
    {
      if (m_verbose)
        cout << "AMRNodeSolver: " << iter << " iterations, converged final NORM max(res) = "
             << currentRes << endl;
    }
  else if (iter >= m_maxIter)
    {
      if (m_verbose)
        cout << "AMRNodeSolver: reached maximum number of "
             << iter << " iterations, final NORM max(res) = "
             << currentRes << endl;
    }
  else if ((currentRes/oldRes) > 1.0-m_operatorTolerance) // hang
    {
      if (m_verbose)
        cout << "AMRNodeSolver: reached solver hang point; " << iter
             << " iterations, final NORM max(res) = " << currentRes << endl;
    }
  else // convergence failed
    {
      if (m_verbose)
        cout << "AMRNodeSolver NOT CONVERGED! - final NORM max(res) = "
             << currentRes << " after " << iter << " iterations" << endl;
    }
}


// ---------------------------------------------------------
void
AMRNodeSolver::AMRVCycleMG(Vector<LevelData<NodeFArrayBox> *>& a_phiLevel,
                           const Vector<LevelData<NodeFArrayBox> *>& a_rhsLevel)
{
  CH_assert(isDefined());
  CH_assert(a_phiLevel.size() > m_finestLevel);
  CH_assert(a_rhsLevel.size() > m_finestLevel);

  // Compute residual on level m_finestLevel.
  m_amrmgLevel[m_finestLevel]->computeAMRResidual(a_phiLevel, a_rhsLevel);

  // SWEEP DOWN V-CYCLE, computing residual at each level in turn.
  for (int ilev = m_finestLevel; ilev > m_lBase; ilev--)
    m_amrmgLevel[ilev]->downSweep(a_phiLevel, a_rhsLevel);

  // Solve at level m_lBase:
  // rhs in bottomRes, initial guess 0 and final estimate in bottomCorr.
  LevelData<NodeFArrayBox>& bottomCorr = m_amrmgLevel[m_lBase]->m_corr;
  const LevelData<NodeFArrayBox>& bottomRes = m_amrmgLevel[m_lBase]->m_resid;
  m_levelSolver.levelSolveH(bottomCorr, bottomRes);

  // Add correction to phi at m_lBase:
  // bottomPhi += bottomCorr.
  LevelData<NodeFArrayBox> & bottomPhi = *a_phiLevel[m_lBase];
  for (DataIterator dit = bottomPhi.dataIterator(); dit.ok(); ++dit)
    bottomPhi[dit()].getFab() += bottomCorr[dit()].getFab();

  // SWEEP BACK UP V-CYCLE, adding correction to phi at each level.
  for (int ilev = m_lBase+1; ilev <= m_finestLevel; ilev++)
    m_amrmgLevel[ilev]->upSweep(a_phiLevel, a_rhsLevel);
}


// ---------------------------------------------------------
Real
AMRNodeSolver::computeResidualNorm(int a_normType)
{
  CH_assert(isDefined());
  CH_assert((a_normType >= 0) && (a_normType <= 2));
  Real normTotal = 0.;

  for (int ilev = m_finestLevel; ilev >= m_lBase; ilev--)
    {
      // Store residual in  m_amrmgLevel[ilev]->m_resid.
      // Removed by petermc, 2 Aug 2002.
      // Avoid this redundant computation, and assume you already have
      // the residual at every level.
      // m_amrmgLevel[ilev]->computeAMRResidual(a_phiLevel, a_rhsLevel);
      Real normLevel = m_amrmgLevel[ilev]->computeResidualNorm(a_normType);
      if (a_normType == 0)
        {
          if (m_verbose)
            cout << "level " << ilev << " res = " << normLevel << endl;
          if (normLevel > normTotal)
            normTotal = normLevel;
        }
      else if (a_normType == 1)
        normTotal += normLevel;
      else if (a_normType == 2)
        normTotal += normLevel * normLevel;
      else
        normTotal += pow(normLevel, Real(a_normType));
    }

  if (! (a_normType == 0 || a_normType == 1))
    if (a_normType == 2)
      normTotal = sqrt(normTotal);
    else
      normTotal = pow(normTotal, (Real)1.0/Real(a_normType));

  return normTotal;
}


// ---------------------------------------------------------
void
AMRNodeSolver::computeAMRResidual(LevelData<NodeFArrayBox> & a_res,
                                  Vector<LevelData<NodeFArrayBox> *>& a_phiLevel,
                                  const Vector<LevelData<NodeFArrayBox> *>& a_rhsLevel,
                                  int a_ilev)
{
  CH_assert(isDefined());
  CH_assert(a_ilev <= m_finestLevel);
  CH_assert(a_ilev >= 0);
  CH_assert(a_phiLevel.size() == m_numLevels);
  CH_assert(a_rhsLevel.size() == m_numLevels);
  CH_assert(a_rhsLevel[a_ilev]->getBoxes() == m_gridsLevel[a_ilev]);
  CH_assert(a_phiLevel[a_ilev]->getBoxes() == m_gridsLevel[a_ilev]);
  CH_assert(a_res.getBoxes() == m_gridsLevel[a_ilev]);

  //compute residual internal to amrmgLevel
  m_amrmgLevel[a_ilev]->computeAMRResidual(a_phiLevel, a_rhsLevel);
  const LevelData<NodeFArrayBox> & amrmgRes = m_amrmgLevel[a_ilev]->m_resid;
  CH_assert(amrmgRes.getBoxes() == m_gridsLevel[a_ilev]);

  //copy residual into output array
  for (DataIterator dit = amrmgRes.dataIterator(); dit.ok(); ++dit)
    a_res[dit()].copy(amrmgRes[dit()]);
}


// ---------------------------------------------------------
void
AMRNodeSolver::applyAMROperator(LevelData<NodeFArrayBox> & a_lofPhi,
                                Vector<LevelData<NodeFArrayBox> *>& a_phiLevel,
                                int a_ilev)
{
  CH_assert(isDefined());
  CH_assert(a_ilev <= m_finestLevel);
  CH_assert(a_ilev >= 0);
  CH_assert(a_phiLevel.size() > m_finestLevel);
  m_amrmgLevel[a_ilev]->applyAMROperator(a_lofPhi, a_phiLevel);
}


// ---------------------------------------------------------
void
AMRNodeSolver::applyAMRGradient(LevelData<NodeFArrayBox> & a_gradPhi,
                                Vector<LevelData<NodeFArrayBox> *>& a_phiLevel,
                                int a_ilev)
{
  CH_assert(isDefined());
  CH_assert(a_ilev <= m_finestLevel);
  CH_assert(a_ilev >= 0);
  CH_assert(a_phiLevel.size() > m_finestLevel);
  m_amrmgLevel[a_ilev]->applyAMRGradient(a_gradPhi, a_phiLevel);
}
