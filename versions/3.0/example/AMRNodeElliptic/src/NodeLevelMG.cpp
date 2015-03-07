#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// NodeLevelMG.cpp
// adapted from LevelMG by DTGraves, Sat, July 18, 1999
// petermc, 30 May 2001

#include "NodeLevelMG.H"
#include "NodeLevelMGF_F.H"

void
NodeLevelMG::setDefaultValues()
{
  // set defaults
  // m_numBottomGSRB = 16;
  m_numSmoothUp = 4;
  m_numSmoothDown = 4;
  m_verbose = false;
  m_lCoarsePtr = NULL;
  m_levelopPtr = NULL;
  m_dx = -1;
  m_isDefined = false;
}

bool
NodeLevelMG::isDefined() const
{
  return m_isDefined;
}

// petermc removed this, 28 Jul 2003
//  void
//  NodeLevelMG::setnumBottomGSRB(int a_numBottomGSRB)
//  {
//    m_numBottomGSRB = a_numBottomGSRB;
//    if (m_lCoarsePtr != NULL)
//      m_lCoarsePtr->setnumBottomGSRB(m_numBottomGSRB);
//  }


void
NodeLevelMG::setnumSmoothUp(int a_numSmoothUp)
{
  m_numSmoothUp = a_numSmoothUp;
  if (m_lCoarsePtr != NULL)
    m_lCoarsePtr->setnumSmoothUp(m_numSmoothUp);
}


void
NodeLevelMG::setnumSmoothDown(int a_numSmoothDown)
{
  m_numSmoothDown = a_numSmoothDown;
  if (m_lCoarsePtr != NULL)
    m_lCoarsePtr->setnumSmoothDown(m_numSmoothDown);
}


void
NodeLevelMG::setVerbose(bool a_verbose)
{
  m_verbose = a_verbose;
  if (m_lCoarsePtr != NULL)
    m_lCoarsePtr->setVerbose(a_verbose);
  if (m_levelopPtr != NULL)
    m_levelopPtr->setVerbose(a_verbose);
}


void
NodeLevelMG::clearMemory()
{
  if (m_lCoarsePtr != NULL)
    {
      delete m_lCoarsePtr;
      m_lCoarsePtr = NULL;
    }
  if(m_levelopPtr != NULL)
    {
      delete m_levelopPtr;
      m_levelopPtr = NULL;
    }
  m_isDefined = false;
}


// default constructor
NodeLevelMG::NodeLevelMG()
{
  setDefaultValues();
}


/// Constructor with DisjointBoxLayout, number of coarser levels
//this calls the levelmg big define function
NodeLevelMG::NodeLevelMG(const DisjointBoxLayout& a_grids,
                         const DisjointBoxLayout* a_gridsCoarsePtr,
                         const Box& a_domain,
                         Real a_dx,
                         int a_refToCoarse,
                         const NodeLevelOp* const a_opin,
                         int a_nCoarserLevels)
{
  setDefaultValues();
  ProblemDomain probdomain(a_domain);
  define(a_grids, a_gridsCoarsePtr, probdomain, a_dx,
         a_refToCoarse, a_opin, a_nCoarserLevels);

}


/// Constructor with DisjointBoxLayout, number of coarser levels
//this calls the levelmg big define function
NodeLevelMG::NodeLevelMG(const DisjointBoxLayout& a_grids,
                         const DisjointBoxLayout* a_gridsCoarsePtr,
                         const ProblemDomain& a_domain,
                         Real a_dx,
                         int a_refToCoarse,
                         const NodeLevelOp* const a_opin,
                         int a_nCoarserLevels)
{
  setDefaultValues();
  define(a_grids, a_gridsCoarsePtr, a_domain, a_dx,
         a_refToCoarse, a_opin, a_nCoarserLevels);

}


// corresponding define function
void
NodeLevelMG::define(const DisjointBoxLayout& a_grids,
                    const DisjointBoxLayout* a_gridsCoarsePtr,
                    const Box& a_domain,
                    Real a_dx,
                    int a_refToCoarse,
                    const NodeLevelOp* const a_opin,
                    int a_nCoarserLevels)
{
  ProblemDomain probdomain(a_domain);
  define(a_grids, a_gridsCoarsePtr, probdomain, a_dx,
         a_refToCoarse, a_opin, a_nCoarserLevels);
}


// corresponding define function
void
NodeLevelMG::define(const DisjointBoxLayout& a_grids,
                    const DisjointBoxLayout* a_gridsCoarsePtr,
                    const ProblemDomain& a_domain,
                    Real a_dx,
                    int a_refToCoarse,
                    const NodeLevelOp* const a_opin,
                    int a_nCoarserLevels)
{
  CH_assert(a_opin != NULL);
  CH_assert(a_nCoarserLevels >= 0);
  CH_assert(!a_domain.isEmpty());

  //either there are no coarser levels
  //or we have a legitimate ref ratio
  CH_assert((a_refToCoarse > 0) || (a_gridsCoarsePtr == NULL));
  CH_assert(a_dx > 0);

  clearMemory();

  m_grids = a_grids;
  m_gridsCoarsePtr = a_gridsCoarsePtr;

  m_dx = a_dx;
  m_domain  = a_domain;
  m_numBottomPoints = m_domain.domainBox().numPts();

  m_nCoarserLevels = a_nCoarserLevels;

  m_isDefined = true;
  int ncomp = 1;
  // need a ghost cell layer for averaging.
  m_resid.define(m_grids, ncomp, (a_refToCoarse/2) * IntVect::Unit);
  m_levelopPtr = a_opin->new_levelop();
  //   int new_ncomp = a_opin->m_ncomp;

  //no need to do inhomogeneous CFInterp here
  //all NodeLevelMG operations (smooth, etc) use
  //only homogeneous c/f interpolation
  bool homogeneousOnly = true;
  m_levelopPtr->define(m_grids,  m_gridsCoarsePtr,
                       m_dx,  a_refToCoarse,
                       m_domain, homogeneousOnly, ncomp, m_verbose);
  // recursively define MG hierarchy
  if (m_nCoarserLevels != 0)
    {

      //this constructor is called at the finest mg level.
      //this means that it might need an averaging operator
      //but it does not need an interpolation operator
      m_refToCoarsened = 2;

      // These are used in FORT_NODEINTERPMG.
      m_boxRef = Box(IntVect::Zero, (m_refToCoarsened-1)*IntVect::Unit);
      Box corners(IntVect::Zero, IntVect::Unit);
      m_weights.define(corners, m_boxRef.numPts());
      FORT_NODEINTERPMG_GETWEIGHTS(CHF_CONST_INT(m_refToCoarsened),
                                   CHF_BOX(m_boxRef),
                                   CHF_FRA(m_weights));

      // OLD stuff; the assertion fails.  I think gridsCoarsePtr isn't needed.
      // CH_assert(m_gridsCoarsePtr != NULL);
      // DisjointBoxLayout gridsCoarse = *m_gridsCoarsePtr;
      // m_averageOp.define(m_grids, gridsCoarse, 1, m_refToCoarsened, m_domain);
      if (m_coarsenedGrids.isClosed())
        {
          m_coarsenedGrids = DisjointBoxLayout();
        }
      coarsen(m_coarsenedGrids, m_grids, m_refToCoarsened);

      //      m_averageOp.define(m_grids, m_coarsenedGrids, 1, m_refToCoarsened, m_domain);
      // petermc, 11 Apr 2003:
      // m_coarsenedGrids is just coarsened version of grids, so use
      // the alternative define
      m_averageOp.define(m_coarsenedGrids, 1, m_refToCoarsened, m_domain);
      m_crseCorr.define(m_coarsenedGrids, 1, IntVect::Unit);
      m_crseResid.define(m_coarsenedGrids, 1, IntVect::Zero);
      m_lCoarsePtr = new NodeLevelMG(*this, m_refToCoarsened, a_opin);
      if (m_lCoarsePtr == NULL)
        {
          MayDay::Error("Out of Memory in NodeLevelMG::define");
        }
    }
}


// --------------------------------------------------------------
//  constructor for coarsened version of NodeLevelMG
// refcoarse is the ratio from coarsened thing to this thing
// --------------------------------------------------------------
NodeLevelMG::NodeLevelMG(NodeLevelMG& a_L,
                         int a_refToCoarse,
                         const NodeLevelOp* a_opin)
{
  //set pointers to NULL so define does not try to delete them
  setDefaultValues();
  define(a_L, a_refToCoarse, a_opin);
}


// --------------------------------------------------------------
//  constructor for coarsened version of NodeLevelMG
// a_refToCoarse is the ratio from coarsened thing to this thing
// --------------------------------------------------------------
void
NodeLevelMG::define(const NodeLevelMG& a_L,
                    int a_refToCoarse,
                    const NodeLevelOp* a_opin)
{
  CH_assert(a_refToCoarse == 2);
  CH_assert(a_opin != NULL);
  clearMemory();

  m_isDefined = true;
  m_verbose = a_L.m_verbose;
  m_grids = a_L.m_coarsenedGrids;
  Real dxCoarse = a_refToCoarse*a_L.m_dx;

  ProblemDomain coarsenedDomain = coarsen(a_L.m_domain, a_refToCoarse);
  m_nCoarserLevels = a_L.m_nCoarserLevels - 1;
  CH_assert(m_nCoarserLevels >= 0);

  m_gridsCoarsePtr = a_L.m_gridsCoarsePtr;
  m_dx = dxCoarse;
  m_domain = coarsenedDomain;
  m_numBottomPoints = m_domain.domainBox().numPts();

  // need a ghost layer for averaging.
  m_resid.define(m_grids, 1, (a_refToCoarse/2) * IntVect::Unit);
  m_levelopPtr = a_L.m_levelopPtr->new_levelop();

  m_levelopPtr->define(a_L.m_levelopPtr, a_refToCoarse);

  // recursively define MG hierarchy
  if (m_nCoarserLevels != 0)
    {
      m_refToCoarsened = 2;

      // These are used in FORT_NODEINTERPMG.
      m_boxRef = Box(IntVect::Zero, (m_refToCoarsened-1)*IntVect::Unit);
      Box corners(IntVect::Zero, IntVect::Unit);
      m_weights.define(corners, m_boxRef.numPts());
      FORT_NODEINTERPMG_GETWEIGHTS(CHF_CONST_INT(m_refToCoarsened),
                                   CHF_BOX(m_boxRef),
                                   CHF_FRA(m_weights));

      // OLD causes segfault because m_gridsCoarsePtr == NULL.
      // m_averageOp.define(m_grids, *m_gridsCoarsePtr, 1, m_refToCoarsened, m_domain);
      coarsen(m_coarsenedGrids, m_grids, m_refToCoarsened);
      // m_averageOp.define(m_grids, m_coarsenedGrids, 1, m_refToCoarsened, m_domain);
      // petermc, 11 Apr 2003:
      // m_coarsenedGrids is just coarsened version of m_grids, so use
      // the alternative define
      m_averageOp.define(m_coarsenedGrids, 1, m_refToCoarsened, m_domain);

      m_crseCorr.define(m_coarsenedGrids, 1, IntVect::Unit);
      m_crseResid.define(m_coarsenedGrids, 1, IntVect::Zero);
      m_lCoarsePtr = new NodeLevelMG(*this, m_refToCoarsened, a_opin);
      if(m_lCoarsePtr == NULL)
        {
          MayDay::Error("Out of Memory in NodeLevelMG::define");
        }
    }
  else
    {
      m_lCoarsePtr = NULL;
    }
}

// ------------------------------------------------------------
// destructor
// ------------------------------------------------------------
NodeLevelMG::~NodeLevelMG()
{
  clearMemory();
}


// -------------------------------------------------------------
void
NodeLevelMG::mgRelax(LevelData<NodeFArrayBox>& a_phi,
                     const LevelData<NodeFArrayBox>& a_rhs,
                     bool a_bottomSolveFlag)
{
  CH_assert(isDefined());
  int ncomp = a_rhs.nComp();
  CH_assert(ncomp == 1);
  CH_assert(a_phi.nComp() == 1);
  CH_assert(a_phi.ghostVect() == IntVect::Unit);

  // now, either coarsen or call bottom solver
  if (m_lCoarsePtr != NULL)  // if level > 0
    {
      // coarser level exists

      for (DataIterator ditc = m_grids.dataIterator(); ditc.ok(); ++ditc)
        {
          m_crseCorr[ditc()].getFab().setVal(0.0);
          m_crseResid[ditc()].getFab().setVal(0.0);
        }

      // first, smooth on this level, m_numSmoothDown iterations
      for (int iter = 0; iter < m_numSmoothDown; iter++)
        m_levelopPtr->smooth(a_phi, a_rhs);

      // define residual on this grid using homogeneous c/f bcs.
      // Set m_resid = rhs - operator(a_phi)
      // with the two-level homogeneous operator.
      m_levelopPtr->residualH(m_resid, a_phi, a_rhs);

      { // begin recursive part:
        // relax averaged residual, and interpolate correction.

        // coarsen residual and physical boundary condition
        // -- what physical boundary condition?

        // Set  m_crseResid = average(m_resid)  at interior points of
        // the coarsened fine-level domain.
        CH_assert(m_averageOp.isDefined());
        m_averageOp.averageToCoarse(m_crseResid, m_resid);

        // relax recursively on coarsened level
        m_lCoarsePtr->mgRelax(m_crseCorr, m_crseResid, a_bottomSolveFlag);

        // interpolate correction m_crseCorr back up and modify solution
        crseCorrect(a_phi);

      } // end recursive part

      // smooth again on this level, m_numSmoothUp iterations
      for (int iter = 0; iter < m_numSmoothUp; iter++)
        m_levelopPtr->smooth(a_phi, a_rhs);

    }
  // there is no coarser level
  // petermc, 1 Oct 2002
  // PREVIOUSLY in the no-coarser-level case:
  // if (m_domain.numPts() == 1), do m_numBottomGSRB smoothings.
  // else if (a_bottomSolveFlag), do m_numBottomGSRB smoothings,
  //                              and also call bottomSmoother.
  // else, do m_numSmoothDown smoothings.

  // So all I have done is in the last case, change the number
  // of iterations from m_numSmoothDown to m_numBottomGSRB.

  // petermc, 28 Jul 2003
  // I used to have m_numBottomGSRB smoothings, and then
  // if (a_bottomSolveFlag && m_domain.domainBox().numPts() > 1)
  //   m_levelopPtr->bottomSmoother(a_phi, a_rhs);
  else if (a_bottomSolveFlag)
    {
      // really at the bottom

      // for (int iter = 0; iter < m_numBottomPoints; iter++)
      for (int iter = 0; iter < m_numSmoothUp + m_numSmoothDown; iter++)
        m_levelopPtr->smooth(a_phi, a_rhs);

      if (m_numBottomPoints > 1)
        m_levelopPtr->bottomSmoother(a_phi, a_rhs);
    }
  else
    {
      // bottom of mini V-cycle

      for (int iter = 0; iter < m_numSmoothUp + m_numSmoothDown; iter++)
        m_levelopPtr->smooth(a_phi, a_rhs);
    }
}


// ---------------------------------------------------------
/** Correct a_phi on intersection with m_crseCorr.
    Called from mgRelax() after the recursive call to the coarser level.
 */
void
NodeLevelMG::crseCorrect(LevelData<NodeFArrayBox>& a_phi)
{
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      Box crseBox(m_coarsenedGrids.get(dit()));

      const FArrayBox& crseFab = m_crseCorr[dit()].getFab();
      FArrayBox& fineFab = a_phi[dit()].getFab();

      FORT_NODEINTERPMG(CHF_FRA(fineFab),
                        CHF_CONST_FRA(crseFab),
                        CHF_BOX(crseBox),
                        CHF_CONST_INT(m_refToCoarsened),
                        CHF_BOX(m_boxRef),
                        CHF_FRA(m_weights));
    }
}


NodeLevelOp*
NodeLevelMG::levelOpPtr()
{
  return m_levelopPtr;
}

NodeLevelMG*
NodeLevelMG::lCoarsePtr()
{
  return m_lCoarsePtr;
}
