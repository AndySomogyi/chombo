#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// NodePoissonOp.cpp
// adapted from PoissonOp.cpp by DTGraves, Weds, July 21, 1999
// petermc, Tue, Dec 12, 2000
// petermc, 7 Apr 2002, added using std::cout and using std::endl

#include "NodePoissonOp.H"
#include "NodePoissonOpF_F.H"
#include "NodeSetOperations.H"
#include "LayoutIterator.H"
// #include "DatasetClient.H"
#include "BoxIterator.H"
#include "LoHiSide.H"
#include "Tuple.H"
#include "NodeCGSmoother.H"
#include "NodeIntegrals.H"
using std::cout;
using std::cerr;
using std::endl;

// ---------------------------------------------------------
// default constructor
NodePoissonOp::NodePoissonOp() : NodeLevelOp()
{
  setDefaultValues();
}


// ---------------------------------------------------------
NodeLevelOp*
NodePoissonOp::new_levelop() const
{
  //only check to see if the boundary conditions are
  //defined.  not necessary to have whole thing defined.
  //the boundary conditions need to be defined because
  //the interface does not know about the boundary
  //conditions.  the solvers call the other define functions
  //in our inimitable two-stage construction
  CH_assert(m_isBCDefined);
  NodePoissonOp* newop = new NodePoissonOp();
  if(newop == NULL)
    {
      MayDay::Error("Out of Memory in NodePoissonOp::new_levelop");
    }
  newop->setDomainNodeBC(m_dombc);
  // petermc, 10 Oct 2002:
  newop->setInterpolationDegree(m_interpolationDegree);
  newop->setVerbose(m_verbose);
  if (m_bottomSmootherPtr != NULL)
    newop->setBottomSmoother(*m_bottomSmootherPtr);
  return static_cast<NodeLevelOp*>(newop);
}


// ---------------------------------------------------------
NodePoissonOp::~NodePoissonOp()
{
  clearMemory();
  if (m_bottomSmootherPtr != NULL)
    {
      delete m_bottomSmootherPtr;
      m_bottomSmootherPtr = NULL;
    }
}


// ---------------------------------------------------------
void
NodePoissonOp::define(
                      const DisjointBoxLayout& a_grids,
                      const DisjointBoxLayout* a_gridsCoarsePtr,
                      Real  a_dx,
                      int a_refToCoarse,
                      const Box& a_domain,
                      bool a_homogeneousOnly,
                      int a_ncomp,
                      bool a_verbose)
{
  ProblemDomain probdomain(a_domain);
  define(a_grids, a_gridsCoarsePtr, a_dx, a_refToCoarse, probdomain,
         a_homogeneousOnly, a_ncomp, a_verbose);
}


// ---------------------------------------------------------
void
NodePoissonOp::define(
                      const DisjointBoxLayout& a_grids,
                      const DisjointBoxLayout* a_gridsCoarsePtr,
                      Real  a_dx,
                      int a_refToCoarse,
                      const ProblemDomain& a_domain,
                      bool a_homogeneousOnly,
                      int a_ncomp,
                      bool a_verbose)
{
  clearMemory();
  m_isDefined = true;
  m_verbose = a_verbose;
  CH_assert(!a_domain.isEmpty());

  m_dx = a_dx;
  m_domain = a_domain;
  m_domainInteriorNodes = surroundingNodes(m_domain.domainBox());
  for (int idir = 0; idir < SpaceDim; idir++)
    if (m_domain.isPeriodic(idir))
      m_domainInteriorNodes.setPeriodic(idir, true);
    else
      {
        m_domainInteriorNodes.setPeriodic(idir, false);
        m_domainInteriorNodes.grow(idir, -1);
      }
  m_grids = a_grids;
  // petermc, 13 May 2003:  added m_domain to support periodic.
  m_exchangeCopier.define(m_grids, m_grids, m_domain, IntVect::Unit);
  m_refToCoarse = a_refToCoarse;
  m_dxCoarse = m_refToCoarse*m_dx;
  m_ihcfiEnabled = (!a_homogeneousOnly);
  m_ncomp = a_ncomp;
  if (a_gridsCoarsePtr != NULL)
    {
      m_gridsCoarse = *a_gridsCoarsePtr;
    }
  setCFIVS();

  if (m_ihcfiEnabled && m_refToCoarse > 1)
    {
      m_qcfi.define(m_grids, m_dx, m_domain,
                    m_loCFIVS, m_hiCFIVS,
                    m_refToCoarse, m_interpolationDegree, m_ncomp, m_verbose);
      if (m_isBCDefined) m_qcfi.setDomainNodeBC(m_dombc);
    }
}


// ---------------------------------------------------------
void
NodePoissonOp::define(const NodeLevelOp* a_opfine,
                      int a_refToFine)
{
  clearMemory();
  const NodePoissonOp* opfineptr =
    dynamic_cast<const NodePoissonOp*>(a_opfine);
  if(opfineptr == NULL)
    MayDay::Error("NodePoissonOp::define: casting failed");
  const NodePoissonOp& opfine = *opfineptr;
  CH_assert(opfine.isDefined());
  CH_assert(a_refToFine > 0);
  m_isDefined = true;
  m_ihcfiEnabled = false;
  setDomainNodeBC(opfine.m_dombc);
  // petermc, 10 Oct 2002:
  setInterpolationDegree(opfine.m_interpolationDegree);

  m_dx = (opfine.m_dx)*a_refToFine;
  m_domain = coarsen(opfine.m_domain, a_refToFine);
  m_domainInteriorNodes = surroundingNodes(m_domain.domainBox());
  for (int idir = 0; idir < SpaceDim; idir++)
    if (m_domain.isPeriodic(idir))
      m_domainInteriorNodes.setPeriodic(idir, true);
    else
      {
        m_domainInteriorNodes.setPeriodic(idir, false);
        m_domainInteriorNodes.grow(idir, -1);
      }
  coarsen(m_grids, opfine.m_grids, a_refToFine);
  // petermc, 13 May 2003:  added m_domain to support periodic.
  m_exchangeCopier.define(m_grids, m_grids, m_domain, IntVect::Unit);
  m_dxCoarse = opfine.m_dxCoarse;
  m_ncomp = opfine.m_ncomp;

  //this refratio is the ratio between this level and the
  //next coarser level in the amr hierarchy
  m_refToCoarse = -1;
  //don't need to define the cfinterp since inhomogeneous cfinterpolaton
  //is never done here.  We shall leave m_qcfi and m_gridsCoarse undefined
  //and just leave the refratio -1.  We DO need fineinterp_ivs for doing
  //homogeneous cfinterpolation
  setCFIVS();
  setBottomSmoother(*(opfine.m_bottomSmootherPtr));
}


// ---------------------------------------------------------
void
NodePoissonOp::setCFIVS()
{
  if (m_verbose)
    cout << "setCFIVS in NodePoissonOp on "
         << m_grids.size() << " grids on "
         << m_domain.domainBox().size()
         << endl;

  DataIterator dit = m_grids.dataIterator();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_loCFIVS[idir].define(m_grids);
      m_hiCFIVS[idir].define(m_grids);
      for (dit.begin(); dit.ok(); ++dit)
        {
          Box bx(m_grids.get(dit()));
          m_loCFIVS[idir][dit()].define(m_domain, bx, m_grids, idir, Side::Lo);
          m_hiCFIVS[idir][dit()].define(m_domain, bx, m_grids, idir, Side::Hi);
        }
    }

  if (m_verbose)
    cout << "IBN/EBN NodePoissonOp on "
         << (m_domain.domainBox().size())
         << ", from "
         << m_grids.size() << " grids to selves on "
         << (m_domain.domainBox().size())
         << endl;

  interiorBoundaryNodes(m_IVSV, m_grids, m_domain);

  // m_IVSVfull added by petermc, 21 Jul 2003
  fullIntVectSets(m_IVSVfull, m_IVSV);

  exteriorBoundaryNodes(m_IVSVext, m_IVSV, m_grids);
}


// ---------------------------------------------------------
bool
NodePoissonOp::isDefined() const
{
  return (m_isDefined && m_isBCDefined);
}


// ---------------------------------------------------------
void
NodePoissonOp::setVerbose(bool a_verbose)
{
  m_verbose = a_verbose;
  if (m_bottomSmootherPtr != NULL)
    m_bottomSmootherPtr->setVerbose(a_verbose);
}


// ---------------------------------------------------------
void
NodePoissonOp::setBottomSmoother(const NodeBaseBottomSmoother& a_bottomSmoother)
{
  if (m_bottomSmootherPtr != NULL) delete m_bottomSmootherPtr;

  m_bottomSmootherPtr = a_bottomSmoother.new_bottomSmoother();
  m_bottomSmootherPtr->setVerbose(m_verbose);
}


// ---------------------------------------------------------
// interpolate a_phif from a_phiCoarse on c/f boundary.
// petermc, 31 Jan 2001
void
NodePoissonOp::CFInterp(LevelData<NodeFArrayBox>& a_phif,
                        const LevelData<NodeFArrayBox>& a_phiCoarse,
                        bool a_inhomogeneous)
{
  CH_assert(isDefined());
  CH_assert(m_ihcfiEnabled);
  CH_assert(a_phif.ghostVect() >= IntVect::Unit);
  CH_assert(m_refToCoarse > 1);
  CH_assert(m_qcfi.isDefined());
  m_qcfi.coarseFineInterp(a_phif, a_phiCoarse, a_inhomogeneous);
}


// ---------------------------------------------------------
// does homogeneous coarse/fine interpolation
// petermc, 22 Jan 2001, for nodes.
void
NodePoissonOp::homogeneousCFInterp(LevelData<NodeFArrayBox>& a_phi)
{
  CH_assert(isDefined());
  // CH_assert(a_phi.ghostVect() >= IntVect::Unit);

  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      const DataIndex& datInd = dit();
      FArrayBox& phiFab = a_phi[datInd].getFab();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              // fill phiFab on side sit(), direction idir, box m_grids[datInd]
              homogeneousCFInterp(phiFab, datInd, idir, sit());
            }
        }
    }
}


// ---------------------------------------------------------
// does homogeneous coarse/fine interpolation.
// petermc, 22 Jan 2001
void
NodePoissonOp::homogeneousCFInterp(FArrayBox& a_phiFab,
                                   const DataIndex& a_datInd,
                                   const int a_idir,
                                   const Side::LoHiSide a_hiorlo)
{
  CH_assert(isDefined());
  CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
  CH_assert (m_ncomp == a_phiFab.nComp());

  const NodeCFIVS* cfivs_ptr = NULL;
  switch (a_hiorlo)
    {
    case Side::Lo:
      {
        cfivs_ptr = &m_loCFIVS[a_idir][a_datInd];
        break;
      }
    case Side::Hi:
      {
        cfivs_ptr = &m_hiCFIVS[a_idir][a_datInd];
        break;
      }
    default:
      {
        cerr << "NodePoissonOp::homogeneousCFInterp(): bogus side" << endl;
        abort();
      }
    }

  // set a_phif to 0 on coarse-fine boundary nodes.
  if (! cfivs_ptr->isEmpty() )
    if (cfivs_ptr->isPacked())
      {
        Box cfivsBox(cfivs_ptr->packedBox());
        // shift from CELL centering to NODE centering
        // because the indices of bx actually represent NODEs.
        cfivsBox.shiftHalf(-IntVect::Unit);
        cfivsBox &= a_phiFab.box();
        // zero out all m_ncomp components starting with 0
        a_phiFab.setVal(0., cfivsBox, 0, m_ncomp);
      }
    else
      {
        IntVectSet interp_ivs = cfivs_ptr->getFineIVS();
        IVSIterator ivsit(interp_ivs);
        for (ivsit.begin(); ivsit.ok(); ++ivsit)
          {
            IntVect iv = ivsit();
            for (int ivar = 0; ivar < m_ncomp; ivar++)
              {
                a_phiFab(iv, ivar) = 0.;
              }
          }
      }
}


// ---------------------------------------------------------
//this smoother assumes problem has
//already been put into residual-correction
//  form, so that crse/fine BC's are homogeneous
void
NodePoissonOp::smooth(LevelData<NodeFArrayBox>& a_phi,
                      const LevelData<NodeFArrayBox>& a_rhs)
{
  CH_assert(isDefined());
  levelGSRB(a_phi, a_rhs);
}


// ---------------------------------------------------------
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
void
NodePoissonOp::applyOpI(LevelData<NodeFArrayBox>& a_LofPhi,
                        LevelData<NodeFArrayBox>& a_phi,
                        const LevelData<NodeFArrayBox>* a_phiCoarsePtr)
{
  CH_assert(isDefined());
  CH_assert(m_ihcfiEnabled);
  CH_assert (a_phi.nComp() == a_LofPhi.nComp());
  CH_assert (a_phi.nComp() == m_ncomp);

  // If there is a coarser level, interpolate from it on C/F boundary.
  if (a_phiCoarsePtr != NULL)
    CFInterp(a_phi, *a_phiCoarsePtr, true); // inhomogeneous phys bc


  // fill in intersection of ghostcells and a_phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  // inhomogeneous physical boundary conditions
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    m_dombc.applyInhomogeneousBCs(a_phi[dit()], m_domain, m_dx);

  applyOpOnly(a_LofPhi, a_phi);
}


// ---------------------------------------------------------
// evaluate gradient, inhomogeneous C/F boundary conditions
// note that gradient applies C/F interpolation and copy bc's
void
NodePoissonOp::gradient(LevelData<NodeFArrayBox>& a_gradPhi,
                        LevelData<NodeFArrayBox>& a_phi,
                        const LevelData<NodeFArrayBox>* a_phiCoarsePtr)
{
  CH_assert(isDefined());
  CH_assert(m_ihcfiEnabled);
  CH_assert(SpaceDim * a_phi.nComp() == a_gradPhi.nComp());
  CH_assert(a_phi.nComp() == m_ncomp);

  // If there is a coarser level, interpolate from it on C/F boundary.
  if (a_phiCoarsePtr != NULL)
    CFInterp(a_phi, *a_phiCoarsePtr, true); // inhomogeneous phys bc

  // fill in intersection of ghostcells and a_phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  // inhomogeneous physical boundary conditions
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    m_dombc.applyInhomogeneousBCs(a_phi[dit()], m_domain, m_dx);

  gradientOnly(a_gradPhi, a_phi);
}


// ---------------------------------------------------------
// evaluate operator, homogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
// homogeneous physical bcs
void
NodePoissonOp::applyOpH(LevelData<NodeFArrayBox>& a_LofPhi,
                        LevelData<NodeFArrayBox>& a_phi)
{
  CH_assert(isDefined());
  CH_assert (a_phi.nComp() == a_LofPhi.nComp());
  CH_assert (a_phi.nComp() == m_ncomp);

  // fill in zeros on coarse-fine interfaces of a_phi.
  homogeneousCFInterp(a_phi);

  // homogeneous physical boundary conditions
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    m_dombc.applyHomogeneousBCs(a_phi[dit()], m_domain, m_dx);

  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  applyOpOnly(a_LofPhi, a_phi);
}


// ---------------------------------------------------------
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
// homogeneous physical bcs
void
NodePoissonOp::applyOpIcfHphys(
                               LevelData<NodeFArrayBox>& a_LofPhi,
                               LevelData<NodeFArrayBox>& a_phi,
                               const LevelData<NodeFArrayBox>* a_phiCoarsePtr)
{
  CH_assert(isDefined());
  CH_assert(m_ihcfiEnabled);
  CH_assert (a_phi.nComp() == a_LofPhi.nComp());
  CH_assert (a_phi.nComp() == m_ncomp);

  // apply C/F boundary conditions...
  if (a_phiCoarsePtr != NULL)
    CFInterp(a_phi, *a_phiCoarsePtr, false); // not inhomogeneous phys bc

  //fill in intersection of ghostcells and phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  // homogeneous physical boundary conditions
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    m_dombc.applyHomogeneousBCs(a_phi[dit()], m_domain, m_dx);

  // removed by petermc, 29 Jul 2003
  // a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  applyOpOnly(a_LofPhi, a_phi);
}


// ---------------------------------------------------------
// evaluate operator, homogeneous C/F boundary conditions
// note -- grids for phi need not correspond to bx
// inhomogeneous physical boundary conditions
void
NodePoissonOp::applyOpHcfIphys(LevelData<NodeFArrayBox>& a_LofPhi,
                               LevelData<NodeFArrayBox>& a_phi)
{
  CH_assert(isDefined());
  CH_assert (a_phi.nComp() == a_LofPhi.nComp());
  CH_assert (a_phi.nComp() == m_ncomp);

  homogeneousCFInterp(a_phi);

  //fill in intersection of ghostcells and phi's boxes
  // petermc removed this, 29 Jul 2003
  // a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    m_dombc.applyInhomogeneousBCs(a_phi[dit()], m_domain, m_dx);

  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  applyOpOnly(a_LofPhi, a_phi);
}


// ---------------------------------------------------------
// This does a GSRB Pre/Conditioned CG on a level
// for the bottom solver.
void
NodePoissonOp::bottomSmoother(LevelData<NodeFArrayBox>& a_phi,
                              const LevelData<NodeFArrayBox>& a_rhs)
{
  CH_assert(isDefined());
  CH_assert(m_bottomSmootherPtr != NULL);
  m_bottomSmootherPtr->doBottomSmooth(a_phi, a_rhs, this);
  // added by petermc, 13 Mar 2003:
  // if periodic in ALL directions, then modify a_phi so that
  // it averages to zero.
  if (m_domain.isPeriodic())
    {
      bool allPeriodic = true;
      for (int idir = 0; idir < SpaceDim; idir++)
        if (! m_domain.isPeriodic(idir) ) allPeriodic = false;

      if (allPeriodic)
        {
          Real dxOne = 1.;  // fake dx to make things simpler
          Real phiIntegral = integral(a_phi, dxOne, a_phi.interval(), false);
          Real volume = Real(m_domain.domainBox().numPts());
          Real phiOffset = phiIntegral / volume;
          DataIterator dit =  a_phi.dataIterator();
          for (dit.reset(); dit.ok(); ++dit)
            {
              FArrayBox& phiFab = a_phi[dit()].getFab();
              phiFab -= phiOffset;
            }
        }
    }
}


//---------------------------------------------------------
// this preconditioner first initializes phihat to (IA)phihat = rhshat
// (diagonization of L -- A is the matrix version of L)
// then smooths with a couple of passes of levelGSRB
void
NodePoissonOp::levelPreconditioner(LevelData<NodeFArrayBox>& a_phihat,
                                   const LevelData<NodeFArrayBox>& a_rhshat)
{
  // diagonal term of this operator is 4/h/h in 2D, 6/h/h in 3D,
  // so inverse of this is our initial multiplier
  Real mult = - (m_dx*m_dx) / (2.0*SpaceDim);
  Interval comps = a_phihat.interval();
  CH_assert(a_phihat.nComp() == a_rhshat.nComp());

  // Use copy() now, more efficient than copyTo().
  for (DataIterator dit = a_phihat.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& phiFab = a_phihat[dit()].getFab();
      const FArrayBox& rhsFab = a_rhshat[dit()].getFab();
      Box bx(rhsFab.box());

      phiFab.copy(bx, comps, bx, rhsFab, comps);

      phiFab *= mult;
    }

  // There are problems with copyTo() on LevelData<NodeFArrayBox>
  // when src and dest are shaped differently.
  // But here the underlying layouts should be the same.
  // a_rhshat.copyTo(comps, a_phihat, comps);
  //
  //  for (dit.begin(); dit.ok(); ++dit)
  //    {
  //      a_phihat[dit()].getFab() *= mult;
  //    }

  levelGSRB(a_phihat, a_rhshat);
  levelGSRB(a_phihat, a_rhshat);
}


//---------------------------------------------------------
void
NodePoissonOp::setDomainNodeBC(const DomainNodeBC& a_dombcIn)
{
  m_dombc = a_dombcIn;
  if (m_ihcfiEnabled && m_refToCoarse > 2) m_qcfi.setDomainNodeBC(m_dombc);
  m_isBCDefined = true;
}


//---------------------------------------------------------
void
NodePoissonOp::setInterpolationDegree(int a_interpolationDegree)
{
  m_interpolationDegree = a_interpolationDegree;
}


// ---------------------------------------------------------
// return to undefined state
void
NodePoissonOp::setDefaultValues()
{
  m_dx = -1.0;
  m_refToCoarse = -1;
  m_isDefined = false;
  m_ihcfiEnabled = false;
  m_isBCDefined = false;
  m_verbose = false;
  // default is to use CG as bottom smoother.
  m_bottomSmootherPtr = new NodeCGSmoother;
  // petermc, 10 Oct 2002:
  // default is to use quadratic/biquadratic interpolation in m_qcfi.
  m_interpolationDegree = 2;
}


// ---------------------------------------------------------
// return to undefined state
void
NodePoissonOp::clearMemory()
{
  // delete storage
  // now reset all other variables
  m_dx = -1.0;
  m_refToCoarse = -1;
  m_isDefined = false;
  m_ihcfiEnabled = false;
  // don't touch bottom smoother here, since we want it to be able
  // to pass through the define statement unaltered
}


// ---------------------------------------------------------
//  levelGSRB does GSRB on a level.  it has no knowledge of overlying
//  finer grids, and assumes coarse/fine BC's are homogeneous
void
NodePoissonOp::levelGSRB(LevelData<NodeFArrayBox>& a_phi,
                         const LevelData<NodeFArrayBox>& a_rhs)
{
  CH_assert(isDefined());
  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == m_ncomp);
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  DataIterator dit = a_phi.dataIterator();

  // do red pass, then black pass
  for (int whichPass = 0; whichPass <= 1; whichPass++)
    {
      // do coarse/fine and copy bc's.
      // should be done patch by patch.

      homogeneousCFInterp(a_phi);

      //fill in intersection of ghostcells and a_phi's boxes
      a_phi.exchange(a_phi.interval(), m_exchangeCopier);

      // invoke physical BC's where necessary
      for (dit.begin(); dit.ok(); ++dit)
        m_dombc.applyHomogeneousBCs(a_phi[dit()], m_domain, m_dx);

      //fill in intersection of ghostcells and a_phi's boxes
      // dfm -- for a 5 point stencil, this should not be necessary
      // added back by petermc, 10 Mar 2003
      // deleted by petermc, 29 Jul 2003
      // a_phi.exchange(a_phi.interval(), m_exchangeCopier);
      // this includes wraparound for periodic

      for (dit.begin(); dit.ok(); ++dit)
        {
          Box thisBox(surroundingNodes(m_grids.get(dit())));
          // removed by petermc, 10 Mar 2003.  ghosts will be used.
          // brought back by petermc, 13 Mar 2003, with new definition
          // of m_domainInteriorNodes depending on boundary conditions.
          thisBox &= m_domainInteriorNodes;
          if (! thisBox.isEmpty() )
            {
              FArrayBox& phiFab = a_phi[dit()].getFab();
              const FArrayBox& rhsFab = a_rhs[dit()].getFab();
              CH_assert(phiFab.box().contains(grow(thisBox, 1)));
              CH_assert(rhsFab.box().contains(thisBox));
              FORT_NODEGSRBLEVELLAP(CHF_FRA(phiFab),
                                    CHF_CONST_FRA(rhsFab),
                                    CHF_BOX(thisBox),
                                    CHF_CONST_REAL(m_dx),
                                    CHF_CONST_INT(whichPass));
            }
        } // end loop through grids

      // deleted by petermc, 29 Jul 2003:  don't need this now
      // a_phi.exchange(a_phi.interval(), m_exchangeCopier);

    } // end loop through red-black
}


// ---------------------------------------------------------
// apply operator; called by applyOpI and applyOpH after
// coarse-fine interpolation.
void
NodePoissonOp::applyOpOnly(LevelData<NodeFArrayBox>& a_LofPhi,
                           const LevelData<NodeFArrayBox>& a_phi)
{
  // Apply operator to a_phi _at interior nodes_.
  // See NodePoissonOp.ChF for NODEOPLAP().

  CH_assert(a_phi.ghostVect() >= IntVect::Unit); // need this for stencil
  Box alldomainbox = m_domain.domainBox();
  // added by petermc, 21 Jul 2003, used in C++ method below.
  // Real dxinv2 = 1. / (m_dx * m_dx);
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      const Box& thisBox = m_grids.get(dit());

      const BaseFab<Real>& thisPhiFab = a_phi[dit()].getFab();

      BaseFab<Real>& thisLPhiFab = a_LofPhi[dit()].getFab();
      // initialize with 0, because otherwise it can contain garbage
      thisLPhiFab.setVal(0.);

      // (1)
      // set thisLPhi at inner NODES of box
      // (thisLPhiFab.box() is NODE-centered)

      // const Box& inner = grow(thisLPhiFab.box(), -1);
      // ... but remember thisLPhiFab has ghost cells now.
      Box boxNodes = surroundingNodes(thisBox);
      // get the enclosed nodes
      boxNodes.grow(-1);

      // Box inner(thisBox);
      // for (int dir = 0; dir < CH_SPACEDIM; dir++)
      // inner.growLo(dir, -1);

      if (! boxNodes.isEmpty() )
        {
          // In the periodic case, you must have copied the edge data
          // to the ghost nodes on the other side.
          FORT_NODEOPLAP(CHF_FRA(thisLPhiFab),
                         CHF_CONST_FRA(thisPhiFab),
                         CHF_BOX(boxNodes),
                         CHF_CONST_REAL(m_dx));
        }

      // (2)
      // remove inner nodes from interior, leaving only grid-boundary nodes.
      // cell indices of innerCells are node indices of interior nodes.

      // new:  actually innerCells are nodes.
      // const Box& innerCells = grow(surroundingNodes(thisBox), -1);

      // now do averaging for each interior boundary node.
      Vector<IntVectSet>& IVSvec = m_IVSV[dit()];
      BitSet& fullvec = m_IVSVfull[dit()];

      for (int lcomp = 0; lcomp < IVSvec.size(); lcomp++)
        {
          IntVectSet& IVS = IVSvec[lcomp];
          if (fullvec[lcomp])
            {
              Box container(IVS.minBox());
              FORT_NODEOPLAP(CHF_FRA(thisLPhiFab),
                             CHF_CONST_FRA(thisPhiFab),
                             CHF_BOX(container),
                             CHF_CONST_REAL(m_dx));
            }
          else
            for (IVSIterator it(IVS); it.ok(); ++it)
              {
                IntVect iv = it();
                // Set thisLPhiFab(iv, 0) to operator on thisPhiFab.
                // Use either Fortran or C++.
                // begin Fortran method
                FORT_NODEOPLAPPOINT(CHF_FRA(thisLPhiFab),
                                    CHF_CONST_FRA(thisPhiFab),
                                    CHF_CONST_INTVECT(iv),
                                    CHF_CONST_REAL(m_dx));
                // end Fortran method
                // begin C++ method
//                  // petermc, 21 Jul 2003
//                  Real lphi = 0.;
//                  for (int idir = 0; idir < SpaceDim; idir++)
//                    {
//                      IntVect e = BASISV(idir);
//                      lphi += dxinv2 *
//                        ((thisPhiFab(iv + e, 0) - thisPhiFab(iv, 0)) -
//                         (thisPhiFab(iv, 0) - thisPhiFab(iv - e, 0)));
//                    }
//                  thisLPhiFab(iv, 0) = lphi;
                // end C++ method
              }
        }
    }
}


// ---------------------------------------------------------
// apply operator; called by applyOpI and applyOpH after
// coarse-fine interpolation.
void
NodePoissonOp::gradientOnly(LevelData<NodeFArrayBox>& a_gradPhi,
                            const LevelData<NodeFArrayBox>& a_phi)
{
  // Apply operator to a_phi _at interior nodes_.
  // See NodePoissonOp.ChF for NODEOPLAP().

  CH_assert(a_phi.ghostVect() >= IntVect::Unit); // need this for stencil
  Box alldomainbox = m_domain.domainBox();

  // added by petermc, 21 Jul 2003, used in C++ method below.
  // Real dxinvh = 0.5 / m_dx;

  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      const Box& thisBox = m_grids.get(dit());

      const BaseFab<Real>& thisPhiFab = a_phi[dit()].getFab();

      BaseFab<Real>& thisGradFab = a_gradPhi[dit()].getFab();
      // initialize with 0, because otherwise it can contain garbage
      thisGradFab.setVal(0.);

      // (1)
      // set thisGrad at inner NODES of box
      // (thisGradFab.box() is NODE-centered)

      // const Box& inner = grow(thisGradFab.box(), -1);
      // ... but remember thisGradFab has ghost cells now.
      Box boxNodes = surroundingNodes(thisBox);
      // get the enclosed nodes
      boxNodes.grow(-1);

      // Box inner(thisBox);
      // for (int dir = 0; dir < CH_SPACEDIM; dir++)
      // inner.growLo(dir, -1);

      if (! boxNodes.isEmpty() )
        {
          // In the periodic case, you must have copied the edge data
          // to the ghost nodes on the other side.
          FORT_NODEGRAD(CHF_FRA(thisGradFab),
                        CHF_CONST_FRA(thisPhiFab),
                        CHF_BOX(boxNodes),
                        CHF_CONST_REAL(m_dx));
        }

      // (2)
      // remove inner nodes from interior, leaving only grid-boundary nodes.
      // cell indices of innerCells are node indices of interior nodes.

      // new:  actually innerCells are nodes.
      // const Box& innerCells = grow(surroundingNodes(thisBox), -1);

      // now do averaging for each interior boundary node.

      Vector<IntVectSet>& IVSvec = m_IVSV[dit()];
      BitSet& fullvec = m_IVSVfull[dit()];

      for (int lcomp = 0; lcomp < IVSvec.size(); lcomp++)
        {
          IntVectSet& IVS = IVSvec[lcomp];
          if (fullvec[lcomp])
            {
              Box container(IVS.minBox());
              FORT_NODEGRAD(CHF_FRA(thisGradFab),
                            CHF_CONST_FRA(thisPhiFab),
                            CHF_BOX(container),
                            CHF_CONST_REAL(m_dx));
            }
          else
            for (IVSIterator it(IVS); it.ok(); ++it)
              {
                IntVect iv = it();
                // Set thisGradFab(iv, 0:SpaceDim-1) to gradient of thisPhiFab.
                // Use either Fortran or C++.
                // begin Fortran method
                FORT_NODEGRADPOINT(CHF_FRA(thisGradFab),
                                   CHF_CONST_FRA(thisPhiFab),
                                   CHF_CONST_INTVECT(iv),
                                   CHF_CONST_REAL(m_dx));
                // end Fortran method
                // begin C++ method
//                  // petermc, 21 Jul 2003
//                  for (int idir = 0; idir < SpaceDim; idir++)
//                    {
//                      IntVect e = BASISV(idir);
//                      thisGradFab(iv, idir) = dxinvh *
//                        (thisPhiFab(iv + e, 0) - thisPhiFab(iv - e, 0));
//                    }
                // end C++ method
              }
        }
    }
}


// ---------------------------------------------------------
void
NodePoissonOp::residualI(LevelData<NodeFArrayBox>& a_resid,
                         LevelData<NodeFArrayBox>& a_phi,
                         const LevelData<NodeFArrayBox>* a_phiCoarsePtr,
                         const LevelData<NodeFArrayBox>& a_rhs)
{
  applyOpI(a_resid, a_phi, a_phiCoarsePtr);
  residualDiff(a_resid, a_rhs);

  // Zero out the difference at the boundary nodes, because
  // the data there is not valid.
  // added by petermc, 24 May 2002
  zeroBoundaryNodes(a_resid, m_IVSVext);
}


// ---------------------------------------------------------
void
NodePoissonOp::residualH(LevelData<NodeFArrayBox>& a_resid,
                         LevelData<NodeFArrayBox>& a_phi,
                         const LevelData<NodeFArrayBox>& a_rhs)
{
  applyOpH(a_resid, a_phi);
  residualDiff(a_resid, a_rhs);
}


// ---------------------------------------------------------
void
NodePoissonOp::residualHcfIphys(LevelData<NodeFArrayBox>& a_resid,
                                LevelData<NodeFArrayBox>& a_phi,
                                const LevelData<NodeFArrayBox>& a_rhs)
{
  applyOpHcfIphys(a_resid, a_phi);
  residualDiff(a_resid, a_rhs);
}


// ---------------------------------------------------------
void
NodePoissonOp::residualIcfHphys(LevelData<NodeFArrayBox>& a_resid,
                                LevelData<NodeFArrayBox>& a_phi,
                                const LevelData<NodeFArrayBox>* a_phiCoarsePtr,
                                const LevelData<NodeFArrayBox>& a_rhs)
{
  applyOpIcfHphys(a_resid, a_phi, a_phiCoarsePtr);
  residualDiff(a_resid, a_rhs);

  // Zero out the difference at the boundary nodes, because
  // the data there is not valid.
  // added by petermc, 8 Oct 2002
  zeroBoundaryNodes(a_resid, m_IVSVext);
}


// ---------------------------------------------------------
void
NodePoissonOp::residualDiff(LevelData<NodeFArrayBox>& a_resid,
                            const LevelData<NodeFArrayBox>& a_rhs)
{
  for (DataIterator dit = a_resid.dataIterator(); dit.ok(); ++dit)
    {
      NodeFArrayBox& resid = a_resid[dit()];
      FArrayBox& residFab = resid.getFab();
      const FArrayBox& rhsFab = a_rhs[dit()].getFab();

      residFab -= rhsFab;
      residFab.negate();
    }
}
