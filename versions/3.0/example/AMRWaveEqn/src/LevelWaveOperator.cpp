#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <cstdio>

#include "LevelWaveOperator.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelFluxRegister.H"
#include "QuadCFInterp.H"
#include "timeInterp.H"
#include "AMRIO.H"

#include "LoHiSide.H"

// Constructor - set up some defaults
LevelWaveOperator::LevelWaveOperator()
{
  m_defined = false;
  m_dx = 0.0;
  m_refineCoarse = 0;
}

// Destructor - free up storage
LevelWaveOperator::~LevelWaveOperator()
{
}

// Define the object so that time stepping can begin
void LevelWaveOperator::define(
                          const DisjointBoxLayout&  a_thisDisjointBoxLayout,
                          const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
                          const ProblemDomain&      a_domain,
                          const int&                a_refineCoarse,
                          const int&                a_numFields,
                          const Real&               a_dx,
                          const bool&               a_hasCoarser,
                          const bool&               a_hasFiner)
{
  // Sanity checks
 CH_assert(a_refineCoarse > 0);
 CH_assert(a_dx > 0.0);

  // Make a copy of the current grids
  m_grids  = a_thisDisjointBoxLayout;

  // Cache data
  m_dx = a_dx;
  m_domain = a_domain;
  m_refineCoarse = a_refineCoarse;
  m_hasCoarser = a_hasCoarser;
  m_hasFiner = a_hasFiner;
  m_numFields = a_numFields;

  // Set up the interpolator if there is a coarser level
  if (m_hasCoarser)
    {
      m_patcher.define(a_thisDisjointBoxLayout,
                       &a_coarserDisjointBoxLayout,
                       m_dx,
                       a_refineCoarse,
                       m_numFields,
                       a_domain);
    }

  // Everything is defined
  m_defined = true;
}

// Evaluate the wave operator at time a_time.
// "a_finerFluxRegister" is the flux register with the next finer level.
// "a_coarseFluxRegister" is flux register with the next coarser level.
// The flux registers are incremented with the normal derivatives of a_phi
// at the boundary, times a fraction of the time step corresponding to
// which stage of an explicit time-stepping scheme (e.g. Runge-Kutta)
// is being invoked.
void LevelWaveOperator::eval(
                        LevelData<FArrayBox>&       a_phi,
                        LevelData<FArrayBox>&       a_LOfPhi,
                        LevelFluxRegister&          a_finerFluxRegister,
                        LevelFluxRegister&          a_coarserFluxRegister,
                        const LevelData<FArrayBox>& a_phiCoarseOld,
                        const Real&                 a_TCoarseOld,
                        const LevelData<FArrayBox>& a_phiCoarseNew,
                        const Real&                 a_TCoarseNew,
                        Real                        a_time,
                        Real                        a_weight)
{
  // Make sure everything is defined
 CH_assert(m_defined);

  // Create temporary storage with a layer of "m_numGhost" ghost cells
  IntVect ivGhost = m_numGhost*IntVect::Unit;

  if (m_hasCoarser)
  {
    // Check that current fine-level time "a_time" falls between the old and new coarse times
    Real alpha = (a_time - a_TCoarseOld) / (a_TCoarseNew - a_TCoarseOld);
    Real dtCoarse = (a_TCoarseNew - a_TCoarseOld);

    // Truncate the fraction to the range [0,1] to remove floating-point
    // subtraction roundoff effects
    Real eps = 0.04 * dtCoarse / m_refineCoarse;

    if (Abs(alpha) < eps)
    {
      alpha = 0.0;
    }

    if (Abs(1.0-alpha) < eps)
    {
      alpha = 1.0;
    }

    // Current time before old coarse time
    if (alpha < 0.0)
    {
      MayDay::Error( "LevelWaveOperator::eval: alpha < 0.0");
    }

    // Current time after new coarse time
    if (alpha > 1.0)
    {
      MayDay::Error( "LevelWaveOperator::eval: alpha > 1.0");
    }

    // Interpolate ghost cells from next coarser level using both space
    // and time interpolation.
    {
      // first, interpolate coarse data to the current time
      LevelData<FArrayBox> phiCoarse; phiCoarse.define(a_phiCoarseOld);
      timeInterp(phiCoarse       ,a_time,
                 a_phiCoarseOld  ,a_TCoarseOld,
                 a_phiCoarseNew  ,a_TCoarseNew,
                 a_phi.interval());
      // use currrent-time coarse data to fill fine ghost cells
      m_patcher.coarseFineInterp(a_phi,phiCoarse);
    }
  }

  // Exchange all the ghost cell data between grids at this level
  a_phi.exchange(a_phi.interval());

  // Beginning of loop through patches/grids.
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    // The current box
    Box curBox = m_grids.get(dit());

    // The current grid of conserved variables
    FArrayBox& curPhi = a_phi[dit()];
    FArrayBox& curLOfPhi = a_LOfPhi[dit()];

    // The current source terms if they exist

    // The fluxes computed for this grid - used for refluxing and returning
    // other face centered quantities
    FArrayBox flux[SpaceDim];

    curLOfPhi.setVal(0.0);

    // Evaluate the single-grid operator.

    for (int idir = 0; idir < SpaceDim;idir++)
    {
        Box fluxBox = surroundingNodes(curBox,idir);
        flux[idir].resize(fluxBox,curLOfPhi.nComp());
        // curLOfPhi is incremented by d^2 phi/dx_idir^2
        // flux[idir] gets d phi/dx_idir for the low-side face
        FORT_APPLYLAP(CHF_FRA(curLOfPhi),
                      CHF_FRA(flux[idir]),
                      CHF_CONST_FRA(curPhi),
                      CHF_BOX(curBox),
                      CHF_BOX(fluxBox),
                      CHF_CONST_REAL(m_dx),
                      CHF_CONST_INT(idir));
    }

    // Do flux register updates
    //
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      // Increment coarse flux register between this level and the next
      // finer level - this level is the next coarser level with respect
      // to the next finer level
      if (m_hasFiner && a_finerFluxRegister.isDefined())
      {
        a_finerFluxRegister.incrementCoarse(flux[idir],a_weight,
                                            dit(),a_phi.interval(),
                                            a_phi.interval(),idir);
      }

      // Increment fine flux registers between this level and the next
      // coarser level - this level is the next finer level with respect
      // to the next coarser level
      if (m_hasCoarser && a_coarserFluxRegister.isDefined())
      {
        a_coarserFluxRegister.incrementFine(flux[idir],a_weight,
                                            dit(),a_phi.interval(),
                                            a_phi.interval(),idir,Side::Lo);
        a_coarserFluxRegister.incrementFine(flux[idir],a_weight,
                                            dit(),a_phi.interval(),
                                            a_phi.interval(),idir,Side::Hi);
      }
    }
  }

}
void LevelWaveOperator::avgdown(LevelData<FArrayBox>& a_phi,
                                LevelData<FArrayBox>& a_phiCoarse)
{
  // Make sure everything is defined

 CH_assert(m_defined);

  // Make sure there is a coarser level.

 CH_assert(m_hasCoarser);

  // Interpolate fine data from next coarser level.
  // We need to use this only if the number of cells over which we are
  // computing the average of l(phi) in FORT_HOAVGDOWN is equal to
  // the size of the refined region. This is always true if nrefine = 2,
  // but need not be otherwise.

  m_patcher.coarseFineInterp(a_phi,a_phiCoarse);

  // Build the stencils for averaging phi and Laplacian(phi).

  Box unitbox(IntVect::Zero,IntVect::Zero);
  Box avStencil = refine(unitbox,m_refineCoarse);
  int halfRefine = m_refineCoarse/2;
  Box lapStencil = refine(unitbox,2) + IntVect::Unit*(halfRefine - 1);
  int sizeLapStencil = lapStencil.volume();
  int sizeAvStencil = avStencil.volume();

  // Exchange all the ghost cell data between grids at this level
  a_phi.exchange(a_phi.interval());

  // Allocate LevelData for corsened version of fine data.

  DisjointBoxLayout coarsenedGrids;
  coarsen(coarsenedGrids,m_grids,m_refineCoarse);

  LevelData<FArrayBox> coarsenedPhi(coarsenedGrids,a_phi.nComp());

  // Beginning of loop through patches/grids.

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    // The current box
    Box curBox = coarsenedGrids.get(dit());

    // The current grid of conserved variables
    FArrayBox& curPhi = a_phi[dit()];
    FArrayBox& curCoarsenedPhi = coarsenedPhi[dit()];

    FORT_HOAVGDOWN(CHF_CONST_FRA(curPhi),
                       CHF_FRA(curCoarsenedPhi),
                   CHF_CONST_INT(m_refineCoarse),
                   CHF_BOX(curBox),
                   CHF_BOX(avStencil),
                   CHF_CONST_INT(sizeAvStencil),
                   CHF_BOX(lapStencil),
                   CHF_CONST_INT(sizeLapStencil)
                   );

  }
  // copy coarsened data onto coarse data.

    coarsenedPhi.copyTo(a_phi.interval(),a_phiCoarse,a_phi.interval());
}
