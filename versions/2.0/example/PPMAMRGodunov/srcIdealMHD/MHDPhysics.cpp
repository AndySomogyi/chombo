#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "LoHiSide.H"
#include "LoHiCenter.H"

#include "MHDPhysics.H"

#include "MHDPhysicsF_F.H"
#include "LGintegrator.H"

#include "CH_HDF5.H"

MHDPhysics::MHDPhysics(const Real& a_smallPressure)
  :
  GodunovPhysics()
{
  m_smallPressure = a_smallPressure;
}

MHDPhysics::~MHDPhysics()
{
}

// Compute the maximum wave speed
Real MHDPhysics::getMaxWaveSpeed(const FArrayBox& a_U,
                                 const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.contains(a_box));

  Real speed = 0.0;

  FORT_MAXWAVESPEEDF(CHF_REAL(speed),
                     CHF_CONST_FRA(a_U),
                     CHF_BOX(a_box));

  return speed;
}

GodunovPhysics* MHDPhysics::new_godunovPhysics() const
{
  CH_assert(m_isBCSet);

  GodunovPhysics* retval = static_cast<GodunovPhysics*>
                           (new MHDPhysics(m_smallPressure));

  retval->setPhysIBC(m_bc);
  return retval;
}

// Number of conserved variables
int MHDPhysics::numConserved()
{
  CH_assert(isDefined());

  return UNUM;
}

// Names of the conserved variables
Vector<string> MHDPhysics::stateNames()
{
  CH_assert(isDefined());

  Vector<string> retval;

  retval.push_back("density");

  retval.push_back("x-momentum");
  retval.push_back("y-momentum");
  retval.push_back("z-momentum");

  retval.push_back("x-magnetic");
  retval.push_back("y-magnetic");
  retval.push_back("z-magnetic");

  retval.push_back("energy-density");

  return retval;
}

// Names of the plot variables
Vector<string> MHDPhysics::plotNames()
{
  CH_assert(isDefined());

  Vector<string> retval;

  retval.push_back("density");

  retval.push_back("x-velocity");
  retval.push_back("y-velocity");
  retval.push_back("z-velocity");

  retval.push_back("x-magnetic");
  retval.push_back("y-magnetic");
  retval.push_back("z-magnetic");

  retval.push_back("pressure");

  return retval;
}
// Number of flux variables
int MHDPhysics::numFluxes()
{
  CH_assert(isDefined());

  // In some computations there may be more fluxes than conserved variables
  return UNUM;
}

void MHDPhysics::getFlux(FArrayBox&       a_flux,
                         const FArrayBox& a_whalf,
                         const int&       a_dir,
                         const Box&       a_box)
{
  CH_assert(isDefined());

  FORT_GETFLUXF(CHF_FRA(a_flux),
                CHF_CONST_FRA(a_whalf),
                CHF_CONST_INT(a_dir),
                CHF_BOX(a_box));
}

// Number of primitive variables
int MHDPhysics::numPrimitives()
{
  CH_assert(isDefined());

  // This doesn't equal the number of conserved variables because
  // auxiliary/redundant variable may be computed and stored
  return WNUM;
}

// Compute expansion amplitudes of dW in right eigenvectors.
void MHDPhysics::charAnalysis(FArrayBox&       a_dW,
                              const FArrayBox& a_W,
                              const int&       a_dir,
                              const Box&       a_box)
{
  CH_assert(isDefined());

  FORT_CHARANALYSISF(CHF_FRA(a_dW),
                     CHF_CONST_FRA(a_W),
                     CHF_CONST_INT(a_dir),
                     CHF_BOX(a_box));
}

void MHDPhysics::charSynthesis(FArrayBox&       a_dW,
                               const FArrayBox& a_W,
                               const int&       a_dir,
                               const Box&       a_box)
{
  CH_assert(isDefined());

  FORT_CHARSYNTHESISF(CHF_FRA(a_dW),
                      CHF_CONST_FRA(a_W),
                      CHF_CONST_INT(a_dir),
                      CHF_BOX(a_box));
}

void MHDPhysics::charValues(FArrayBox&       a_lambda,
                            const FArrayBox& a_W,
                            const int&       a_dir,
                            const Box&       a_box)
{
  CH_assert(isDefined());

  FORT_CHARVALUESF(CHF_FRA(a_lambda),
                   CHF_CONST_FRA(a_W),
                   CHF_CONST_INT(a_dir),
                   CHF_BOX(a_box));
}

// Add to (increment) the source terms given the current state
// Empty implementation - does nothing.
void MHDPhysics::incrementSource(FArrayBox&       a_S,
                                 const FArrayBox& a_W,
                                 const Box&       a_box)
{
}

// Compute a Riemann problem and generate fluxes at the faces.
void MHDPhysics::riemann(FArrayBox&       a_WGdnv,
                         const FArrayBox& a_WLeft,
                         const FArrayBox& a_WRight,
                         const FArrayBox& a_W,
                         const Real&      a_time,
                         const int&       a_dir,
                         const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert(a_WGdnv.box().contains(a_box));

  // Get the numbers of relevant variables
  int numPrim = numPrimitives();

  CH_assert(a_WGdnv .nComp() == numPrim);
  CH_assert(a_WLeft .nComp() == numPrim);
  CH_assert(a_WRight.nComp() == numPrim);

  // Cast away "const" inputs so their boxes can be shifted left or right
  // 1/2 cell and then back again (no net change is made!)
  FArrayBox& shiftWLeft  = (FArrayBox&)a_WLeft;
  FArrayBox& shiftWRight = (FArrayBox&)a_WRight;

  // Solution to the Riemann problem

  // Shift the left and right primitive variable boxes 1/2 cell so they are
  // face centered
  shiftWLeft .shiftHalf(a_dir, 1);
  shiftWRight.shiftHalf(a_dir,-1);

  CH_assert(shiftWLeft .box().contains(a_box));
  CH_assert(shiftWRight.box().contains(a_box));

  // Riemann solver computes Wgdnv all edges that are not on the physical
  // boundary.
  FORT_RIEMANNF(CHF_FRA(a_WGdnv),
                CHF_CONST_FRA(shiftWLeft),
                CHF_CONST_FRA(shiftWRight),
                CHF_CONST_INT(a_dir),
                CHF_BOX(a_box));

  // Call boundary Riemann solver (note: periodic BC's are handled there).
  m_bc->primBC(a_WGdnv,shiftWLeft ,a_W,a_dir,Side::Hi,a_time);
  m_bc->primBC(a_WGdnv,shiftWRight,a_W,a_dir,Side::Lo,a_time);

  // Shift the left and right primitive variable boxes back to their original
  // position.
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
}

void MHDPhysics::postNormalPred(FArrayBox&       a_dWMinus,
                                FArrayBox&       a_dWPlus,
                                const FArrayBox& a_W,
                                const Real&      a_dt,
                                const Real&      a_dx,
                                const int&       a_dir,
                                const Box&       a_box)

{
  Box inbox(a_box);
  inbox.grow(a_dir,1);

  Box loBox,hiBox,centerBox,domain;
  int hasLo,hasHi;

  loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,domain,
             inbox,m_domain,a_dir);

#if 0
  pout() << "postNormalPred Direction" << a_dir << endl;

  if (hasLo)
  {
    pout() << "postNormalPred loBox" << loBox << endl;
  }

  if (hasHi)
  {
    pout() << "postNormalPred hiBox" << hiBox << endl;
  }

  pout() << "postNormalPred centerBox" << centerBox << endl;
  pout() << "postNormalPred a_box" << a_box << endl;
  pout() << "postNormalPred domain" << domain << endl;
#endif

  // Stone correction
  Real dtbydx = a_dt / a_dx;
  FORT_STONECORRECTION(CHF_FRA(a_dWMinus),
                       CHF_CONST_FRA(a_W),
                       CHF_BOX(loBox),
                       CHF_CONST_INT(hasLo),
                       CHF_BOX(hiBox),
                       CHF_CONST_INT(hasHi),
                       CHF_BOX(centerBox),
                       CHF_CONST_INT(a_dir),
                       CHF_CONST_REAL(dtbydx));

  FORT_STONECORRECTION(CHF_FRA(a_dWPlus),
                       CHF_CONST_FRA(a_W),
                       CHF_BOX(loBox),
                       CHF_CONST_INT(hasLo),
                       CHF_BOX(hiBox),
                       CHF_CONST_INT(hasHi),
                       CHF_BOX(centerBox),
                       CHF_CONST_INT(a_dir),
                       CHF_CONST_REAL(dtbydx));
}

// Compute increment of primitive variables.
void MHDPhysics::quasilinearUpdate(FArrayBox&       a_dWdx,
                                   const FArrayBox& a_WHalf,
                                   const FArrayBox& a_W,
                                   const Real&      a_scale,
                                   const int&       a_dir,
                                   const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_dWdx.box().contains(a_box));

  FORT_GETADWDXF(CHF_FRA(a_dWdx),
                 CHF_CONST_FRA(a_WHalf),
                 CHF_CONST_FRA(a_W),
                 CHF_CONST_REAL(a_scale),
                 CHF_CONST_INT(a_dir),
                 CHF_BOX(a_box));
}

// Compute the primitive variables from the conserved variables
void MHDPhysics::consToPrim(FArrayBox&       a_W,
                            const FArrayBox& a_U,
                            const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  FORT_CONSTOPRIMF(CHF_FRA(a_W),
                   CHF_CONST_FRA(a_U),
                   CHF_BOX(a_box));
}

// Interval within the primitive variables corresponding to the velocities
Interval MHDPhysics::velocityInterval()
{
  CH_assert(isDefined());

#if CH_SPACEDIM==2
  Interval retval(WVELX,WVELY);
#elif CH_SPACEDIM==3
  Interval retval(WVELX,WVELZ);
#else
  bogus_spacedim();
#endif

  return retval;
}

// Component index within the primitive variables of the pressure
int MHDPhysics::pressureIndex()
{
  CH_assert(isDefined());

  return WPRES;
}

// Used to limit the absolute value of a "pressure" difference
Real MHDPhysics::smallPressure()
{
  CH_assert(isDefined());

  return m_smallPressure;
}

// Component index within the primitive variables of the bulk modulus
int MHDPhysics::bulkModulusIndex()
{
  CH_assert(isDefined());

  return WPRES;
}
