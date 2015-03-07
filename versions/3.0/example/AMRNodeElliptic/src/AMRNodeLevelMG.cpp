#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// AMRNodeLevelMG.cpp
// from AMRLevelMG by DTGraves, Thurs, July 15, 1999
// petermc, Tues, Nov 28, 2000
// petermc, 4 Apr 2002, added using std::cout and using std::endl

#include "SPMD.H"
#include "AMRNodeSolver.H"
#include "AMRNodeLevelMG.H"
#include "NodeSetOperations.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
// #include "DatasetClient.H"
#include "Norms.H"
using std::cout;
using std::endl;

// ---------------------------------------------------------
/// constructor
AMRNodeLevelMG::AMRNodeLevelMG()
{
  setDefaultValues();
}


// ---------------------------------------------------------
/// constructor
AMRNodeLevelMG::AMRNodeLevelMG(const AMRNodeSolver* const a_parent,
                               int a_level,
                               const NodeLevelOp* const a_opin)
{
  setDefaultValues();
  define(a_parent, a_level, a_opin);
}


// ---------------------------------------------------------
/// destructor
AMRNodeLevelMG::~AMRNodeLevelMG()
{
  clearMemory();
}


// ---------------------------------------------------------
/// define level
void
AMRNodeLevelMG::define(const AMRNodeSolver* const a_parent,
                       int a_level,
                       const NodeLevelOp* const a_opin)
{
  clearMemory();
 CH_assert(a_parent != NULL);
 CH_assert(a_parent->isDefined());
 CH_assert(a_level >= 0);
 CH_assert(a_level <= a_parent->m_finestLevel);
  m_isDefined = true;
  m_parent = a_parent;
  m_level = a_level;
  int nCoarserLevels;
  const DisjointBoxLayout& grids = a_parent->m_gridsLevel[a_level];
  const DisjointBoxLayout* gridsCoarsePtr = NULL;
  //put reasonable value into refratio
  int refToCoarse = 2;
  m_domain = a_parent->m_domainLevel[a_level];
  if (m_level > 0)
    {
      const DisjointBoxLayout& baseGrids =
        a_parent->m_gridsLevel[a_level-1];
      gridsCoarsePtr = &baseGrids;
      refToCoarse = a_parent->m_refRatio[a_level-1];
      // nCoarserLevels = log2(reftoCoarse) - 1
      // bug fixed by petermc, 15 Apr 2002
      nCoarserLevels = -1;
      for (int interRatio = refToCoarse; interRatio >= 2; interRatio /= 2)
        nCoarserLevels++;
    }
  else
    {
      nCoarserLevels = 0;
    }

  m_dx = m_parent->m_dxLevel[m_level];

  m_levelMG.define(grids, gridsCoarsePtr, m_domain, m_dx,
                   refToCoarse, a_opin, nCoarserLevels);
  m_levelMG.setVerbose(m_verbose);

  int ncomps = 1;
  m_mginterp.define(grids, ncomps, refToCoarse, m_domain);

  //residual,lofphi has no ghost cells, corr and phi have 1 ghost cell
  m_LofPhi.define(grids, 1, (refToCoarse/2) * IntVect::Unit);
  m_resid.define(grids, 1, IntVect::Zero);
  // m_residInterior.define(grids, 1, IntVect::Zero);
  m_corr.define(grids, 1, IntVect::Unit);
  m_dcorr.define(grids, 1, IntVect::Unit);
  m_phiSave.define(grids, 1,IntVect::Unit);

  // copyInteriorNodes(m_resid, m_residInterior, m_domain);
  if (m_verbose)
    cout << "IBN/EBN AMRNodeLevelMG on " << m_domain.domainBox().size()
         << ", from "
         << grids.size() << " grids to selves on " << m_domain.domainBox().size()
         << endl;

  interiorBoundaryNodes(m_IVSV, grids, m_domain);

  exteriorBoundaryNodes(m_IVSVext, m_IVSV, grids);

  if (m_level > 0)
    {
      coarsen(m_coarsenedGrids, grids, refToCoarse);

      ProblemDomain coarsenedBase = coarsen(m_domain, refToCoarse);

      m_residCoarsened.define(m_coarsenedGrids, 1,IntVect::Zero);

      // copyInteriorNodes(nextCoarser.m_resid, m_residCoarsened, m_coarsenedBase);

      if (m_verbose)
        cout << "IBN AMRNodeLevelMG on " << m_domain.domainBox().size()
             << ", from "
             << m_coarsenedGrids.size() << " coarsened grids to "
             << (a_parent->m_gridsLevel[a_level-1]).size()
             << " coarser-level grids on " << coarsenedBase.domainBox().size()
             << endl;

      interiorBoundaryNodes(m_IVSVcoarsened, a_parent->m_gridsLevel[a_level-1], m_coarsenedGrids, coarsenedBase);

      // m_averageOp.define(grids, *gridsCoarsePtr, 1, refToCoarse, m_domain);
      // m_averageOp.define(grids, m_coarsenedGrids, 1, refToCoarse, m_domain);
      // petermc, 11 Apr 2003:
      // m_coarsenedGrids is just coarsened version of grids, so use
      // the alternative define
      m_averageOp.define(m_coarsenedGrids, 1, refToCoarse, m_domain);
      m_averageOp.setVerbose(m_verbose);
    }

  m_levelOpPtr = a_opin->new_levelop();
  bool homogeneousOnly = false;
  m_levelOpPtr->define(grids, gridsCoarsePtr,
                       m_dx, refToCoarse,
                       m_domain, homogeneousOnly, ncomps, m_verbose);

  if (m_level < m_parent->m_finestLevel)
    {
      m_finerGrids.deepCopy(m_parent->m_gridsLevel[m_level+1]);
      m_finerGrids.close();

      int refToFine =  m_parent->m_refRatio[m_level];
      coarsen(m_coarsenedFineGrids, m_finerGrids, refToFine);

      if (m_verbose)
        cout << "IBN AMRNodeLevelMG on " << m_domain.domainBox().size()
             << ", from "
             << m_coarsenedFineGrids.size() << " coarsened finer-level grids to "
             << grids.size() << " grids on " << m_domain.domainBox().size()
             << endl;

      interiorBoundaryNodes(m_IVSVcoarsenedFine, grids, m_coarsenedFineGrids, m_domain);
    }
}


// ---------------------------------------------------------
/// has define function been called?  if not, most fcns won't work
bool
AMRNodeLevelMG::isDefined() const
{
  return m_isDefined;
}


// ---------------------------------------------------------
/// complete, 3-level operator
// petermc, 23 Jan 2001
void
AMRNodeLevelMG::applyAMROperator(LevelData<NodeFArrayBox>& a_LofPhi,
                                 Vector<LevelData<NodeFArrayBox>*>& a_phiLevel)
{
 CH_assert(isDefined());
  LevelData<NodeFArrayBox>& phi = *a_phiLevel[m_level];

  // Initialize a_LofPhi with zero (on boxes in this proc).
  for (DataIterator dit = a_LofPhi.dataIterator(); dit.ok(); ++dit)
    a_LofPhi[dit()].getFab().setVal(0.0);

  // If there is a finer level, then project a_phiLevel[m_level+1]
  // on interior coarse nodes down to a_phiLevel[m_level].

  if (m_level < m_parent->m_finestLevel)
    projectFineInterior(phi, *a_phiLevel[m_level+1]);

  LevelData<NodeFArrayBox>* crsePhi = NULL;
  if (m_level > 0) crsePhi = a_phiLevel[m_level-1];

  // Evaluate operator on this level including interpolated bc's from Level-1.
  m_levelOpPtr->applyOpI(a_LofPhi, phi, crsePhi);

  // In the CELL-centered version, now you reflux to enforce
  // flux-matching from finer levels.
}



// ---------------------------------------------------------
/// complete, 3-level gradient
// petermc, 23 Jan 2001
void
AMRNodeLevelMG::applyAMRGradient(LevelData<NodeFArrayBox>& a_gradPhi,
                                 Vector<LevelData<NodeFArrayBox>*>& a_phiLevel)
{
 CH_assert(isDefined());
  LevelData<NodeFArrayBox>& phi = *a_phiLevel[m_level];

  // Initialize a_gradPhi with zero (on boxes in this proc).
  for (DataIterator dit = a_gradPhi.dataIterator(); dit.ok(); ++dit)
    a_gradPhi[dit()].getFab().setVal(0.0);

  // If there is a finer level, then project a_phiLevel[m_level+1]
  // on interior coarse nodes down to a_phiLevel[m_level].

  if (m_level < m_parent->m_finestLevel)
    projectFineInterior(phi, *a_phiLevel[m_level+1]);

  LevelData<NodeFArrayBox>* crsePhi = NULL;
  if (m_level > 0) crsePhi = a_phiLevel[m_level-1];

  // Evaluate operator on this level including interpolated bc's from Level-1.
  m_levelOpPtr->gradient(a_gradPhi, phi, crsePhi);

  // In the CELL-centered version, now you reflux to enforce
  // flux-matching from finer levels.
}



// ---------------------------------------------------------
// set m_resid = *a_rhsLevel[m_level] - operator(*a_phiLevel[m_level])
// with the complete, 3-level operator
void
AMRNodeLevelMG::computeAMRResidual(Vector<LevelData<NodeFArrayBox>* >& a_phiLevel,
                                   const Vector<LevelData<NodeFArrayBox>* >& a_rhsLevel)
{
 CH_assert(isDefined());
  const LevelData<NodeFArrayBox>& rhs = *a_rhsLevel[m_level];

  // First set  m_resid = operator(a_phiLevel).
  applyAMROperator(m_resid, a_phiLevel);
  // Then set  m_resid = rhs - m_resid.
  for (DataIterator dit = rhs.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& residFab = m_resid[dit()].getFab();
      const FArrayBox& rhsFab = rhs[dit()].getFab();

      residFab -= rhsFab;
      residFab.negate();
    }
  // We need to set m_resid == 0 on phys boundary and c/f boundary.
  // Otherwise m_resid == rhs here, and this will affect norm(m_resid).
  zeroBoundaryNodes(m_resid, m_IVSVext);
}


// ---------------------------------------------------------
// Project phiFine on coarse interior nodes down to phi.
void
AMRNodeLevelMG::projectFineInterior(LevelData<NodeFArrayBox>& a_phi,
                                    const LevelData<NodeFArrayBox>& a_phiFine)
{
  int refToFine =  m_parent->m_refRatio[m_level];
  int numComps = a_phi.nComp();

  // coarsen a_phiFine, call it phiFineCoarsened
  LevelData<NodeFArrayBox> phiFineCoarsened(m_coarsenedFineGrids, numComps, IntVect::Zero);

  for (DataIterator dit = m_coarsenedFineGrids.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& phiFineCoarsenedFab = phiFineCoarsened[dit()].getFab();
      const FArrayBox& phiFineFab = a_phiFine[dit()].getFab();

      const Box coarsenedNodes = surroundingNodes(m_coarsenedFineGrids.get(dit()));

      for (BoxIterator bit(coarsenedNodes); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          for (int ivar = 0; ivar < numComps; ivar++)
            phiFineCoarsenedFab(iv, ivar) = phiFineFab(refToFine*iv, ivar);
        }
    }

  // Copy phiFineCoarsened to a_phi on interior nodes of coarsenedFineGrids.
  // Code is in NodeFArrayBox.

  copyInteriorNodes(a_phi, phiFineCoarsened, m_IVSVcoarsenedFine);
}


// ---------------------------------------------------------
//      Sweep up V-cycle
void
AMRNodeLevelMG::upSweep(Vector<LevelData<NodeFArrayBox> *>& a_phiLevel,
                        const Vector<LevelData<NodeFArrayBox> *>& a_rhsLevel)
{
 CH_assert(isDefined());
 CH_assert(m_level > 0);

  LevelData<NodeFArrayBox>& phi = *a_phiLevel[m_level];
  const LevelData<NodeFArrayBox>& rhs = *a_rhsLevel[m_level];
 CH_assert(phi.nComp() == 1);
 CH_assert(rhs.nComp() == 1);

  // Interpolate corrections  m_corr  from next coarser level to
  // the correction field  m_corr  at this level.

  // Initialize correction at the next coarser level.
  //this uses the two-level only operator (Lnf)
  AMRNodeLevelMG& nextCoarser = *m_parent->m_amrmgLevel[m_level-1];

  m_mginterp.interpToFine(m_corr, nextCoarser.m_corr);

  DataIterator dit = m_resid.dataIterator();

  // two-level operator on the correction:
  // modify m_corr on boundary with coarse level by interpolation
  // from coarser m_corr; then apply homogeneous physical bcs.
  m_levelOpPtr->applyOpIcfHphys(m_LofPhi, m_corr, &nextCoarser.m_corr);
  // Finish resid calculation, m_resid -= op(m_corr)
  for (dit.reset(); dit.ok(); ++dit)
    m_resid[dit()].getFab() -= m_LofPhi[dit()].getFab();

  // initialize the correction to the correction, and apply smoother.
  for (dit.reset(); dit.ok(); ++dit)
    m_dcorr[dit()].getFab().setVal(0.);

  smooth(m_dcorr, m_resid);

  // Increment correction:  m_corr += m_dcorr
  // Increment saved value of phi with corr:  m_phiSave += m_corr
  // Overwrite phi:  phi = m_phiSave.
  for (dit.reset(); dit.ok(); ++dit)
    {
      m_corr[dit()].getFab() += m_dcorr[dit()].getFab();
      m_phiSave[dit()].getFab() += m_corr[dit()].getFab();
      phi[dit()].copy(m_phiSave[dit()]);
    }
}


// ---------------------------------------------------------
//    Sweep down V-cycle
void
AMRNodeLevelMG::downSweep(Vector<LevelData<NodeFArrayBox>* >& a_phiLevel,
                          const Vector<LevelData<NodeFArrayBox>* >& a_rhsLevel)
{
 CH_assert(isDefined());
  LevelData<NodeFArrayBox>& phi = *a_phiLevel[m_level];
  LevelData<NodeFArrayBox>& rhs = *a_rhsLevel[m_level];
 CH_assert(phi.nComp() == 1);
 CH_assert(rhs.nComp() == 1);

  // You should already have m_resid at this level.

  // set m_phiSave = phi
  // set m_corr = 0
  // smooth(m_corr, m_resid)
  // set phi += m_corr

  // If coarser level exists:
  // set m_corr[coarser] = 0
  // set m_LofPhi = m_resid - op(m_corr, m_corr[coarser])
  // set m_residCoarsened = average(m_LofPhi)
  // set m_resid[coarser] = rhs[coarser] - op(phi[coarser])
  // set m_resid[coarser] = m_residCoarsened  at valid nodes of this level

  DataIterator dit = phi.dataIterator();

  // Copy phi into m_phiSave, and set m_corr to zero.
  for (dit.reset(); dit.ok(); ++dit)
    {
      m_phiSave[dit()].copy(phi[dit()]);
      m_corr[dit()].getFab().setVal(0.);
    }

  // Apply smoother on m_resid to get new m_corr.
  smooth(m_corr, m_resid);

  // Update phi:  phi += m_corr.
  for (dit.reset(); dit.ok(); ++dit)
    phi[dit()].getFab() += m_corr[dit()].getFab();

  // Form residual on next coarser level.
  if (m_level > 0)
    {
      // Initialize correction at the next coarser level.
      // this has to be done before applyOpI
      //int refToCoarse = m_parent->m_refRatio[m_level-1];
      AMRNodeLevelMG& nextCoarser = *m_parent->m_amrmgLevel[m_level-1];

      for (DataIterator ditc = nextCoarser.m_corr.dataIterator(); ditc.ok(); ++ditc)
        nextCoarser.m_corr[ditc()].getFab().setVal(0.);

      // Set m_LofPhi = m_resid - operator(m_corr)
      // with the two-level IcfHphys operator.

      // m_levelOpPtr->residualIcfHphys(m_LofPhi, m_corr, &nextCoarser.m_corr, m_resid);
      // nextCoarser.m_corr is zero anyway
      m_levelOpPtr->residualH(m_LofPhi, m_corr, m_resid);

      // Set m_residCoarsened to average of m_LofPhi.
      m_averageOp.averageToCoarse(m_residCoarsened, m_LofPhi);

      nextCoarser.computeAMRResidual(a_phiLevel, a_rhsLevel);

      // overwrite residual on next coarser level (nextCoarser.m_resid) with
      // coarsened residual from this level (this.m_residCoarsened),
      // on the interior coarse nodes of this level.

      copyInteriorNodes(nextCoarser.m_resid, m_residCoarsened, m_IVSVcoarsened);
    }
}


// ---------------------------------------------------------
// compute and return norm of internal data m_resid.
Real
AMRNodeLevelMG::computeResidualNorm(int a_normType) const
{
 CH_assert(isDefined());
  Real normPeterson = computeNorm(m_resid, a_normType);
  return normPeterson;
}


// ---------------------------------------------------------
void
AMRNodeLevelMG::setnumSmoothUp(int a_numSmoothUp)
{
 CH_assert(isDefined());
  m_levelMG.setnumSmoothUp(a_numSmoothUp);
}


// ---------------------------------------------------------
void
AMRNodeLevelMG::setnumSmoothDown(int a_numSmoothDown)
{
 CH_assert(isDefined());
  m_levelMG.setnumSmoothDown(a_numSmoothDown);
}


// ---------------------------------------------------------
void
AMRNodeLevelMG::setVerbose(bool a_verbose)
{
  m_verbose = a_verbose;
  if ( isDefined() )
    {
      m_levelMG.setVerbose(a_verbose);
      if (m_levelOpPtr != NULL)
        m_levelOpPtr->setVerbose(a_verbose);
    }
}


// ---------------------------------------------------------
void
AMRNodeLevelMG::setDefaultValues()
{
  m_isDefined = false;
  m_parent = NULL;
  m_level = -1;
  m_isDefined = false;
  m_levelOpPtr = NULL;
  m_verbose = true;
}


// ---------------------------------------------------------
void
AMRNodeLevelMG::clearMemory()
{
  if (m_levelOpPtr != NULL)
    delete m_levelOpPtr;
  m_levelOpPtr = NULL;
}


// ---------------------------------------------------------
//  smooth a_phi
void
AMRNodeLevelMG::smooth(LevelData<NodeFArrayBox>& a_phi,
                       const LevelData<NodeFArrayBox>& a_rhs)
{
 CH_assert(isDefined());
 CH_assert(a_phi.ghostVect() >= IntVect::Unit);
 CH_assert(a_phi.getBoxes() == m_parent->m_gridsLevel[m_level]);
 CH_assert(a_rhs.getBoxes() == m_parent->m_gridsLevel[m_level]);
  //no bottom solve here
  m_levelMG.mgRelax(a_phi, a_rhs, false);
}


// ---------------------------------------------------------
// returns normType norm of mfab
Real
AMRNodeLevelMG::computeNorm(const LevelData<NodeFArrayBox>& a_nfinput,
                            int a_normType) const
{
 CH_assert(isDefined());
 CH_assert((a_normType >= 0) && (a_normType <= 2));
 CH_assert(a_nfinput.nComp() == 1);

  const Interval& intvl = a_nfinput.interval();
  Real normLevel = 0.;
  if (m_level < m_parent->m_finestLevel)
    {
      int refToFine =  m_parent->m_refRatio[m_level];
      normLevel = norm(a_nfinput, m_domain, m_coarsenedFineGrids,
                       m_IVSVext, m_IVSVcoarsenedFine,
                       refToFine, m_dx, intvl, a_normType, m_verbose);
//        if (a_normType == 0)
//          // We have saved m_IVSVcoarsenedFine so that it can be reused here.
//          return maxnorm(a_nfinput, m_domain, m_coarsenedFineGrids,
//                         m_IVSVcoarsenedFine, refToFine, intvl, m_verbose);
//        else
//          return norm(a_nfinput, m_domain, &m_finerGrids,
//                      refToFine, m_dx, intvl, a_normType, m_verbose);
    }
  else
    {
      // no finer level
      // return norm(a_nfinput, m_dx, a_normType, intvl, m_verbose);
      normLevel = norm(a_nfinput, m_domain,
                       m_IVSVext,
                       m_dx, intvl, a_normType, m_verbose);
    }
  return normLevel;
}
