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
#include <string>
using std::string;

#include "LoHiSide.H"

#include "GodunovUtilities.H"
#include "LoHiCenter.H"
#include "GodunovPhysics.H"
#include "GodunovUtilitiesF_F.H"
#include "NamespaceHeader.H"

// Flag everything as not defined or set
GodunovUtilities::GodunovUtilities()
{
  m_isDefined = false;
}

GodunovUtilities::~GodunovUtilities()
{
}

// Define this object and the boundary condition object
void GodunovUtilities::define(ProblemDomain& a_domain,
                              const Real&    a_dx)
{
  // Store the domain and grid spacing
  m_domain    = a_domain;
  m_dx        = a_dx;
  m_isDefined = true;
}

// Compute the flattening coefficients from the primitive variables
void GodunovUtilities::computeFlattening(FArrayBox&       a_flattening,
                                         const FArrayBox& a_W,
                                         const Interval&  a_velInt,
                                         const int&       a_presInd,
                                         const Real&      a_smallPres,
                                         const int&       a_bulkModulusInd,
                                         const Box&       a_box)
{
  CH_assert(m_isDefined);
  CH_assert(a_W.box().contains(a_box));

  // The current direction
  int idir;

  // The directional flattening coefficients
  FArrayBox zetaDir(a_box,SpaceDim);

  // The divergence of the velocity
  FArrayBox dVel(a_box,SpaceDim);

  // The interval of the primitive variables corresponding to the velocity
  Interval velInterval= a_velInt;
  int v0index = velInterval.begin();

  // Get the directional flattening coefficients in each direction
  for (idir = 0; idir < SpaceDim; idir++)
  {
    // A box one larger (in direction "idir") than the final result box
    Box box1 = a_box;
    box1.grow(idir,1);

    // A box two larger (in direction "idir") than the final result box
    Box box2 = a_box;
    box2.grow(idir,2);

    // A box three larger (in direction "idir") than the final result box
    Box box3 = a_box;
    box3.grow(idir,3);

    // The primitive variables need to be defined over the largest box
    // CH_assert(a_W.box().contains(box3));

    // Compute where centered differences can be used and where one sided
    // differences need to be used.  The data used for the differences is
    // defined on "box3"
    Box loBox,hiBox,centerBox,entireBox;
    int hasLo,hasHi;

    loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
               box3,m_domain,idir);

    // Compute the first differences in "pressure"
    FArrayBox delta1p(entireBox,1);

    FORT_GETGRADF(CHF_FRA1(delta1p,0),
                  CHF_CONST_FRA1(a_W,a_presInd),
                  CHF_CONST_INT(idir),
                  CHF_BOX(loBox),
                  CHF_CONST_INT(hasLo),
                  CHF_BOX(hiBox),
                  CHF_CONST_INT(hasHi),
                  CHF_BOX(centerBox));

    // Compute where centered differences can be used and where one sided
    // differences need to be used.  The data used for the differences is
    // defined on "box2"
    loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
               box2,m_domain,idir);

    // Compute the second differences in "pressure"
    FArrayBox delta2p(entireBox,1);

    FORT_GETDPTWOF(CHF_FRA1(delta2p,0),
                   CHF_CONST_FRA1(delta1p,0),
                   CHF_CONST_INT(idir),
                   CHF_BOX(loBox),
                   CHF_CONST_INT(hasLo),
                   CHF_BOX(hiBox),
                   CHF_CONST_INT(hasHi),
                   CHF_BOX(centerBox));

    // Compute a 3-way minimum of the "bulk modulus"
    FArrayBox bulkMin(entireBox,1);

    int bulkIndex = a_bulkModulusInd;

    FORT_MIN3PTSF(CHF_FRA1(bulkMin,0),
                  CHF_CONST_FRA1(a_W,bulkIndex),
                  CHF_CONST_INT(idir),
                  CHF_BOX(loBox),
                  CHF_CONST_INT(hasLo),
                  CHF_BOX(hiBox),
                  CHF_CONST_INT(hasHi),
                  CHF_BOX(centerBox));

    // Use the first and second differences normalized by the 3-way minimum
    // computed above to generate directional flattening coefficients
    FArrayBox zetaTwiddleDir(entireBox,1);

    FORT_GETFLATF(CHF_FRA1(zetaTwiddleDir,0),
                  CHF_CONST_FRA1(delta1p,0),
                  CHF_CONST_FRA1(delta2p,0),
                  CHF_CONST_REAL(a_smallPres),
                  CHF_CONST_FRA1(bulkMin,0),
                  CHF_BOX(entireBox));

    // Compute where centered differences can be used and where one sided
    // differences need to be used.  The data used for the differences is
    // defined on "box1"
    loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
               box1,m_domain,idir);

    // Take a 3-way minimum of the directional flattening coefficients
    FORT_MIN3PTSF(CHF_FRA1(zetaDir,idir),
                  CHF_CONST_FRA1(zetaTwiddleDir,0),
                  CHF_CONST_INT(idir),
                  CHF_BOX(loBox),
                  CHF_CONST_INT(hasLo),
                  CHF_BOX(hiBox),
                  CHF_CONST_INT(hasHi),
                  CHF_BOX(centerBox));

    // Compute each component of the divergence of the velocity
    FORT_GETGRADF(CHF_FRA1(dVel,idir),
                  CHF_CONST_FRA1(a_W,v0index+idir),
                  CHF_CONST_INT(idir),
                  CHF_BOX(loBox),
                  CHF_CONST_INT(hasLo),
                  CHF_BOX(hiBox),
                  CHF_CONST_INT(hasHi),
                  CHF_BOX(centerBox));
  }

  // At each point, set the flattening coefficient to the minimum of all
  // the directional flattening coefficients if the divergence of the velocity
  // is negative, otherwise set it to 1 (no flattening).
  FORT_MINFLATF(CHF_FRA1(a_flattening,0),
                CHF_CONST_FRA(zetaDir),
                CHF_CONST_FRA(dVel),
                CHF_BOX(a_box));
}

void GodunovUtilities::applyFlattening(      FArrayBox& a_dW,
                                       const FArrayBox& a_flat,
                                       const Box&       a_box)
{
  CH_assert(a_dW.box().contains(a_box));
  CH_assert(a_flat.box().contains(a_box));

  int numSlopes = a_dW.nComp();

  for (int islope = 0;islope < numSlopes;islope++)
  {
    a_dW.mult(a_flat,a_box,0,islope);
  }
}

void GodunovUtilities::vanLeerSlopes(FArrayBox&       a_dW,
                                     const FArrayBox& a_W,
                                     const int&       a_numSlopes,
                                     const bool&      a_useLimiting,
                                     const int&       a_dir,
                                     const Box&       a_box)
{
  // A box one larger (in direction "a_dir") than the final result box
  Box box1 = a_box;
  box1.grow(a_dir,1);

  // Compute where centered differences can be used and where one sided
  // differences need to be used.
  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
             box1,m_domain,a_dir);

  // Compute 2nd order slopes - including one sided differences
  FArrayBox dWMinus(entireBox,a_numSlopes);
  FArrayBox dWPlus (entireBox,a_numSlopes);

  slopes(a_dW, dWMinus, dWPlus, a_W, a_numSlopes, a_dir,
         loBox, hiBox, centerBox, entireBox, hasLo, hasHi);

  // Apply the slope limiter if requested
  if (a_useLimiting)
  {
    slopeLimiter(a_dW,dWMinus,dWPlus,a_numSlopes,centerBox);
  }
}

void GodunovUtilities::fourthOrderSlopes(FArrayBox&       a_dW,
                                         const FArrayBox& a_W,
                                         const FArrayBox& a_dWvL,
                                         const int&       a_numSlopes,
                                         const int&       a_dir,
                                         const Box&       a_box)
{
  // Number of slopes to compute
  int numSlope = a_numSlopes;

  CH_assert(a_dW.nComp() == numSlope);
  CH_assert(a_W.nComp() >= numSlope);

  // A box one larger (in direction "a_dir") than the final result box
  Box box1 = a_box;
  box1.grow(a_dir,1);

  // Compute where centered differences can be used and where one sided
  // differences need to be used.
  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
             box1,m_domain,a_dir);

  CH_assert(a_dW.box().contains(entireBox));

  FORT_FOURTHSLOPEDIFFSF(CHF_FRA(a_dW),
                         CHF_CONST_FRA(a_W),
                         CHF_CONST_FRA(a_dWvL),
                         CHF_CONST_INT(numSlope),
                         CHF_CONST_INT(a_dir),
                         CHF_BOX(loBox),
                         CHF_CONST_INT(hasLo),
                         CHF_BOX(hiBox),
                         CHF_CONST_INT(hasHi),
                         CHF_BOX(centerBox));
}

void GodunovUtilities::oneSidedDifferences(FArrayBox&       a_dWMinus,
                                           FArrayBox&       a_dWPlus,
                                           const FArrayBox& a_W,
                                           const int&       a_dir,
                                           const Box&       a_box)
{
  int numSlopes = a_dWMinus.nComp();

  // A box one larger (in direction "a_dir") than the final result box
  Box box1 = a_box;
  box1.grow(a_dir,1);

  // Compute where centered differences can be used and where one sided
  // differences need to be used.
  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
             box1,m_domain,a_dir);

  // Compute 2nd order slopes - including one sided differences
  FArrayBox deltaWC(entireBox, numSlopes);

  slopes(deltaWC, a_dWMinus, a_dWPlus, a_W, numSlopes, a_dir,
         loBox, hiBox, centerBox, entireBox, hasLo, hasHi);
}

void GodunovUtilities::slopes(FArrayBox&       a_dWCent,
                              FArrayBox&       a_dWMinus,
                              FArrayBox&       a_dWPlus,
                              const FArrayBox& a_W,
                              const int&       a_numSlopes,
                              const int&       a_dir,
                              const Box&       a_loBox,
                              const Box&       a_hiBox,
                              const Box&       a_centerBox,
                              const Box&       a_entireBox,
                              const int&       a_hasLo,
                              const int&       a_hasHi)
{
  CH_assert(a_dWCent .nComp() == a_numSlopes);
  CH_assert(a_dWMinus.nComp() == a_numSlopes);
  CH_assert(a_dWPlus .nComp() == a_numSlopes);
  CH_assert(a_W.nComp() >= a_numSlopes);

  CH_assert(a_dWCent .box().contains(a_entireBox));
  CH_assert(a_dWMinus.box().contains(a_entireBox));
  CH_assert(a_dWPlus .box().contains(a_entireBox));
  CH_assert(a_W.box().contains( a_entireBox));

  CH_assert((a_dir >= 0 )&& (a_dir < SpaceDim));

  FORT_SECONDSLOPEDIFFSF(CHF_FRA(a_dWCent),
                         CHF_FRA(a_dWMinus),
                         CHF_FRA(a_dWPlus),
                         CHF_CONST_FRA(a_W),
                         CHF_CONST_INT(a_numSlopes),
                         CHF_CONST_INT(a_dir),
                         CHF_BOX(a_loBox),
                         CHF_CONST_INT(a_hasLo),
                         CHF_BOX(a_hiBox),
                         CHF_CONST_INT(a_hasHi),
                         CHF_BOX(a_centerBox));
}

// Apply a van Leer limiter directly to the slopes.
void GodunovUtilities::slopeLimiter(FArrayBox&       a_dW,
                                    const FArrayBox& a_dWLeft,
                                    const FArrayBox& a_dWRigh,
                                    const int&       a_numSlopes,
                                    const Box&       a_box)
{
  CH_assert(m_isDefined);
  CH_assert(a_dW.nComp()     == a_numSlopes);
  CH_assert(a_dWLeft.nComp() == a_numSlopes);
  CH_assert(a_dWRigh.nComp() == a_numSlopes);
  CH_assert(a_dW.box().contains(a_box));
  CH_assert(a_dWLeft.box().contains(a_box));
  CH_assert(a_dWRigh.box().contains(a_box));

  FORT_VANLEERLIMITERF(CHF_FRA(a_dW),
                       CHF_CONST_FRA(a_dWLeft),
                       CHF_CONST_FRA(a_dWRigh),
                       CHF_CONST_INT(a_numSlopes),
                       CHF_BOX(a_box));
}

// Piecewise linear normal predictor. Computes increments in the
// characteristic amplitudes.
void GodunovUtilities::PLMNormalPred(FArrayBox&       a_dWCharLo,
                                     FArrayBox&       a_dWCharHi,
                                     const FArrayBox& a_dWChar,
                                     const FArrayBox& a_Lambda,
                                     const Real&      a_dtbydx,
                                     const Box&       a_box)
{
  int numPrim = a_dWChar.nComp();
  CH_assert(a_dWCharLo.nComp() == numPrim);
  CH_assert(a_dWCharHi.nComp() == numPrim);
  CH_assert(a_Lambda.nComp() == numPrim);
  CH_assert(a_dWCharLo.box().contains(a_box));
  CH_assert(a_dWCharHi.box().contains(a_box));
  CH_assert(a_dWChar.box().contains(a_box));
  CH_assert(a_Lambda.box().contains(a_box));

  FORT_PLMNORMALPREDF(CHF_FRA(a_dWCharLo),
                      CHF_FRA(a_dWCharHi),
                      CHF_CONST_FRA(a_dWChar),
                      CHF_CONST_FRA(a_Lambda),
                      CHF_CONST_REAL(a_dtbydx),
                      CHF_CONST_INT(numPrim),
                      CHF_BOX(a_box));
}

void GodunovUtilities::PPMFaceValues(FArrayBox&            a_WFace,
                                     const FArrayBox&      a_W,
                                     const int&            a_numSlopes,
                                     const bool&           a_useLimiting,
                                     const int&            a_dir,
                                     const Box&            a_box,
                                     const Real&           a_time,
                                     const GodunovPhysics* a_physPtr)
{
  // Input box is the face-centered box on which a_WFace is computed.

  // No limiter needs to be applied to these values, since the van Leer
  // limiter already takes care of such issues.

  CH_assert(a_WFace.box().contains(a_box));

  Box vlBox = a_box;
  vlBox.enclosedCells();
  vlBox.grow(a_dir,1);

  FArrayBox dW(vlBox,a_numSlopes);
  vanLeerSlopes(dW,a_W,a_numSlopes,a_useLimiting,a_dir,vlBox);

  if (a_physPtr != NULL)
  {
    Real a_time = 0.0;
    a_physPtr->getPhysIBC()->setBdrySlopes(dW,a_W,a_dir,a_time);
  }

  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  loHiCenterFace(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
                 vlBox,m_domain,a_dir);

  FORT_PPMFACEVALUESF(CHF_FRA(a_WFace),
                      CHF_CONST_FRA(a_W),
                      CHF_CONST_FRA(dW),
                      CHF_CONST_INT(a_numSlopes),
                      CHF_CONST_INT(a_dir),
                      CHF_BOX(loBox),
                      CHF_CONST_INT(hasLo),
                      CHF_BOX(hiBox),
                      CHF_CONST_INT(hasHi),
                      CHF_BOX(centerBox));
}

void GodunovUtilities::PPMLimiter(FArrayBox& a_dWMinus,
                                  FArrayBox& a_dWPlus,
                                  const int& a_numSlopes,
                                  const Box& a_box)
{
  FORT_PPMLIMITERF(CHF_FRA(a_dWMinus),
                   CHF_FRA(a_dWPlus),
                   CHF_CONST_INT(a_numSlopes),
                   CHF_BOX(a_box));
}

void GodunovUtilities::PPMNormalPred(FArrayBox&       a_dWMinus,
                                     FArrayBox&       a_dWPlus,
                                     const FArrayBox& a_Lambda,
                                     const Real&      a_dtbydx,
                                     const int&       a_numSlopes,
                                     const Box&       a_box)
{
  FORT_PPMNORMALPREDF(CHF_FRA(a_dWMinus),
                      CHF_FRA(a_dWPlus),
                      CHF_CONST_FRA(a_Lambda),
                      CHF_CONST_REAL(a_dtbydx),
                      CHF_CONST_INT(a_numSlopes),
                      CHF_BOX(a_box));
}

// Compute a face centered divergence of the velocity
void GodunovUtilities::divVel(FArrayBox&       a_divVel,
                              const FArrayBox& a_W,
                              const Interval&  a_velInt,
                              const int&       a_dir,
                              const Box&       a_box)
{
  CH_assert((a_dir >= 0 )&& (a_dir < SpaceDim));
  CH_assert(a_divVel.box().contains(a_box));

  Box dveltanBox = a_box;
  dveltanBox.enclosedCells(a_dir);
  dveltanBox.grow(1);
  dveltanBox &= m_domain;

  // First, we need to calculate the directional derivatives of
  // the tangential components of velocity at the cell-centers.
  // want to make sure numTanComps is > 0 even in 1D, since
  // we pass it in to the divergence (although it isn't used there)
  int numTanComps = std::max(1, SpaceDim-1);
  FArrayBox dveltan(dveltanBox,numTanComps);

  // Get the interval of the primitive variables corresponding to the velocity
  int v0index = a_velInt.begin();

  // Go through the tangential directions
  for (int i = 0, dir = 0; dir < SpaceDim; ++dir)
  {
    if (dir != a_dir)
    {
      // This velocity component is tangential.  Build the box in which
      // d(v[dir])/d(x[dir]) is to be computed.
      Box primBox = a_box;
      primBox.enclosedCells(a_dir);
      primBox.grow(dir,1).grow(a_dir,1);

      // Compute where centered differences can be used and where one sided
      // differences need to be used.  The data used for the differences is
      // defined on "dveltanBox"
      Box gradBox,hiBox,loBox,centerBox;
      int hasLo,hasHi;

      loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,gradBox,
                 primBox,m_domain,dir);

      // Compute d(v[dir])/d(x[dir]).
      FORT_GETGRADF(CHF_FRA1(dveltan,i),
                    CHF_CONST_FRA1(a_W,v0index+dir),
                    CHF_CONST_INT(dir),
                    CHF_BOX(loBox),
                    CHF_CONST_INT(hasLo),
                    CHF_BOX(hiBox),
                    CHF_CONST_INT(hasHi),
                    CHF_BOX(centerBox));
      i++;
    }
  }

  // Now, we can calculate the divergence of the normal velocity
  // at the center normal-direction faces. To do this, we determine
  // the faces at which we have sufficient data to compute centered
  // estimates of h*(div(u)). At the remaining faces. i.e. those
  // corresponding to the physical boundaries, we use zeroth-order
  // extrapolation.

  Box divBox = a_box;
  divBox.enclosedCells(a_dir);
  divBox.grow(a_dir,1);

  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  loHiCenterFace(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
                 divBox,m_domain,a_dir);

  // All of the boxes computed above are shifted so as to be cell-centered,
  // with the index of the cell center being identified with the low face.
  // We then shift a_divVel to be compatible with that convention on input,
  // and undo the shift on output.  Basically, everything is made compatible
  // a_W (which is cell-centered).

  loBox.shiftHalf(a_dir,1);
  centerBox.shiftHalf(a_dir,1);
  hiBox.shiftHalf(a_dir,1);

  a_divVel.shiftHalf(a_dir,1);

  FORT_DIVUEDGEF(CHF_FRA1(a_divVel,0),
                 CHF_CONST_FRA1(a_W,v0index+a_dir),
                 CHF_CONST_FRA(dveltan),
                 CHF_CONST_INT(a_dir),
                 CHF_BOX(loBox),
                 CHF_CONST_INT(hasLo),
                 CHF_BOX(hiBox),
                 CHF_CONST_INT(hasHi),
                 CHF_BOX(centerBox));

  a_divVel.shiftHalf(a_dir,-1);
}

// Increment fluxes with artificial viscosity.
void GodunovUtilities::artificialViscosity(FArrayBox&       a_F,
                                           const FArrayBox& a_U,
                                           const FArrayBox& a_divVel,
                                           const Real&      a_coeff,
                                           const int&       a_dir,
                                           const Box&       a_box)
{
  CH_assert((a_dir >= 0 )&& (a_dir < SpaceDim));
  CH_assert(a_F.box().contains(a_box));
  CH_assert(a_divVel.box().contains(a_box));
  CH_assert(a_coeff >= 0.0);

  // Note: a_box is face centered in the a_dir direction and is set up for
  // updating a_F

  // Remove any boundary faces in direction a_dir from a_box
  Box noBoundaryBox = a_box;
  noBoundaryBox.enclosedCells(a_dir);
  noBoundaryBox.grow(a_dir,1);
  noBoundaryBox &= m_domain;
  noBoundaryBox.grow(a_dir,-1);
  noBoundaryBox.surroundingNodes(a_dir);

  FORT_ARTVISCF(CHF_FRA(a_F),
                CHF_CONST_FRA(a_U),
                CHF_CONST_FRA1(a_divVel,0),
                CHF_CONST_REAL(a_coeff),
                CHF_CONST_INT(a_dir),
                CHF_BOX(noBoundaryBox));
}
#include "NamespaceFooter.H"
