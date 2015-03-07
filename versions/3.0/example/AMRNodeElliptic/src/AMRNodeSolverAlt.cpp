#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// AMRNodeSolverAlt.cpp
// adapted from AMRSolver by DTGraves, Tues, July 6, 1999
// AMRNodeSolver by petermc, 28 Nov 2000
// AMRNodeSolverAlt by petermc, 22 May 2002
// petermc, 23 Aug 2002, added a_interpolationDegree
// petermc, 10 Oct 2002, removed a_interpolationDegree

#include <cmath>
#include <iostream>
using std::cout;
using std::endl;

#include "AMRNodeSolverAlt.H"
#include "LayoutIterator.H"
// #include "DatasetClient.H"


// ---------------------------------------------------------
void
AMRNodeSolverAlt::setDefaultValues()
{
  m_numLevels = -1;
  m_finestLevel = -1;
  m_refRatio.resize(0);
  m_tolerance = 1.0e-10;
#ifdef CH_USE_FLOAT
  m_tolerance = 10.*sqrt(m_tolerance);
#endif
  m_numSmoothUp = 4;
  m_numSmoothDown = 4;
  m_verbose = true;
  m_isDefined = false;
}


// ---------------------------------------------------------
AMRNodeSolverAlt::AMRNodeSolverAlt() : m_levelSolver()
{
  setDefaultValues();
}


// ---------------------------------------------------------
void
AMRNodeSolverAlt::clear()
{
  for (int ilev = 0; ilev < m_levelSolver.size(); ilev++)
    {
      if (m_levelSolver[ilev] != NULL)
        delete m_levelSolver[ilev];
    }
}


// ---------------------------------------------------------
AMRNodeSolverAlt::~AMRNodeSolverAlt()
{
  clear();
}


// ---------------------------------------------------------
AMRNodeSolverAlt::AMRNodeSolverAlt(const Vector<DisjointBoxLayout>& a_gridsLevel,
                                   const Vector<Box>& a_domainLevel,
                                   const Vector<Real>& a_dxLevel,
                                   const Vector<int>& a_refRatio,
                                   int a_numLevels,
                                   int a_lBase,
                                   const NodeLevelOp* const a_opin,
                                   int a_minLength) : m_levelSolver()
{
  setDefaultValues();
  Vector<ProblemDomain> probdomainLevel(a_numLevels);
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    probdomainLevel[ilev] = ProblemDomain(a_domainLevel[ilev]);

  define(a_gridsLevel, probdomainLevel, a_dxLevel,
         a_refRatio, a_numLevels, a_lBase, a_opin, a_minLength);
}


// ---------------------------------------------------------
AMRNodeSolverAlt::AMRNodeSolverAlt(const Vector<DisjointBoxLayout>& a_gridsLevel,
                                   const Vector<ProblemDomain>& a_domainLevel,
                                   const Vector<Real>& a_dxLevel,
                                   const Vector<int>& a_refRatio,
                                   int a_numLevels,
                                   int a_lBase,
                                   const NodeLevelOp* const a_opin,
                                   int a_minLength) : m_levelSolver()
{
  setDefaultValues();
  define(a_gridsLevel, a_domainLevel, a_dxLevel,
         a_refRatio, a_numLevels, a_lBase, a_opin, a_minLength);
}


// ---------------------------------------------------------
void
AMRNodeSolverAlt::define(const Vector<DisjointBoxLayout>& a_gridsLevel,
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
AMRNodeSolverAlt::define(const Vector<DisjointBoxLayout>& a_gridsLevel,
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
  CH_assert(m_lBase >= 0);

  // for each level, check that all grids are in the physical domain.
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      const Box& domainbox = m_domainLevel[ilev].domainBox();
      const DisjointBoxLayout& grids = m_gridsLevel[ilev];
      for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
        CH_assert(domainbox.contains(grids[dit()]));
    }

  // define a LevelNodeSolver object for each level.
  m_levelSolver.resize(m_numLevels, NULL);
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      const DisjointBoxLayout& levsolvGrids = m_gridsLevel[ilev];
      const ProblemDomain&     levsolvDom   = m_domainLevel[ilev];
      const Real&              levsolvDx    = m_dxLevel[ilev];
      const DisjointBoxLayout* levsolvCoarser = NULL;
      int levsolvRef = 2;
      if (ilev > 0)
        {
          levsolvCoarser = &m_gridsLevel[ilev-1];
          levsolvRef = m_refRatio[ilev-1];
        }
      m_levelSolver[ilev] = new LevelNodeSolver;
      m_levelSolver[ilev]->setVerbose(m_verbose);
      m_levelSolver[ilev]->define(levsolvGrids, levsolvCoarser,
                                  levsolvDom,   levsolvDx,
                                  levsolvRef, a_opin,
                                  a_minLength);
      m_levelSolver[ilev]->setnumSmoothUp(m_numSmoothUp);
      m_levelSolver[ilev]->setnumSmoothDown(m_numSmoothDown);
    }
}


// ---------------------------------------------------------
bool
AMRNodeSolverAlt::isDefined() const
{
  return m_isDefined;
}


// ---------------------------------------------------------
void
AMRNodeSolverAlt::setNumSmoothUp(int a_numSmoothUp)
{
  CH_assert(a_numSmoothUp >= 0);
  m_numSmoothUp = a_numSmoothUp;
  if (isDefined())
    for (int ilev = m_lBase; ilev <= m_finestLevel; ilev++)
      m_levelSolver[ilev]->setnumSmoothUp(m_numSmoothUp);
}


// ---------------------------------------------------------
void
AMRNodeSolverAlt::setNumSmoothDown(int a_numSmoothDown)
{
  CH_assert(a_numSmoothDown >= 0);
  m_numSmoothDown = a_numSmoothDown;
  if (isDefined())
    for (int ilev = m_lBase; ilev <= m_finestLevel; ilev++)
      m_levelSolver[ilev]->setnumSmoothDown(m_numSmoothDown);
}


// ---------------------------------------------------------
void
AMRNodeSolverAlt::setTolerance(Real a_tolerance)
{
  m_tolerance = a_tolerance;
  if (isDefined())
    for (int ilev = m_lBase; ilev <= m_finestLevel; ilev++)
      m_levelSolver[ilev]->setTolerance(a_tolerance);
}


// ---------------------------------------------------------
void
AMRNodeSolverAlt::setBottomSmoothing(bool a_doBottomSmooth)
{
  CH_assert(isDefined());
  for (int ilev = m_lBase; ilev <= m_finestLevel; ilev++)
    m_levelSolver[ilev]->setBottomSmoothing(a_doBottomSmooth);
}


// ---------------------------------------------------------
void
AMRNodeSolverAlt::setMaxIter(int a_maxIter)
{
  CH_assert(isDefined());
  CH_assert(a_maxIter >= 0);
  for (int ilev = m_lBase; ilev <= m_finestLevel; ilev++)
    m_levelSolver[ilev]->setMaxIter(a_maxIter);
}


// ---------------------------------------------------------
void
AMRNodeSolverAlt::setVerbose(bool a_verbose)
{
  m_verbose = a_verbose;
  if (isDefined())
    for (int ilev = m_lBase; ilev <= m_finestLevel; ilev++)
      m_levelSolver[ilev]->setVerbose(a_verbose);
}


// ---------------------------------------------------------
void
AMRNodeSolverAlt::solveAMR(Vector<LevelData<NodeFArrayBox> *>& a_phiLevel,
                        const Vector<LevelData<NodeFArrayBox> *>& a_rhsLevel)
{
  CH_assert(isDefined());
  CH_assert(a_phiLevel.size() > m_finestLevel);
  CH_assert(a_rhsLevel.size() > m_finestLevel);

  if (m_verbose)
    cout << "AMRNodeSolverAlt, levels " << m_lBase
         << " to " << m_finestLevel
         << endl;

  // Solve at level m_lBase.

  if (m_verbose)
    cout << "AMRNodeSolverAlt::solveAMR() solving at level "
         << m_lBase << " ************* " << endl;

  if (m_lBase > 0)
    // interpolate phi from m_lBase-1 to m_lBase on c/f boundary.
    m_levelSolver[m_lBase]->levelSolve(*a_phiLevel[m_lBase],
                                       a_phiLevel[m_lBase-1],
                                       *a_rhsLevel[m_lBase],
                                       true); // initialize phi to 0
  else
    m_levelSolver[m_lBase]->levelSolve(*a_phiLevel[m_lBase],
                                       NULL, // no coarser level
                                       *a_rhsLevel[m_lBase],
                                       true); // initialize phi to 0

  for (int ilev = m_lBase+1; ilev <= m_finestLevel; ilev++)
    {
      if (m_verbose)
        cout << "AMRNodeSolverAlt::solveAMR() solving at level "
             << ilev << " ************* " << endl;

      // Solve at level ilev.
      // levelSolve() initializes a_phiLevel[ilev] with zeroes.
      // levelSolve() interpolates phi from ilev-1 to ilev on c/f boundary.
      m_levelSolver[ilev]->levelSolve(*a_phiLevel[ilev],
                                      a_phiLevel[ilev-1],
                                      *a_rhsLevel[ilev],
                                      true); // initialize phi to 0
    }
}
