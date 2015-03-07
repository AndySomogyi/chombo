#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _GODUNOVUTILITIES_H_
#define _GODUNOVUTILITIES_H_

#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "REAL.H"
#include "FArrayBox.H"
#include "NamespaceHeader.H"

class GodunovPhysics;

///
/**
   Utility class for higher-order Godunov methods: slopes, parabolic
   interpolants, limiters. Contains no physics-dependent methods, but
   one of the member functions (PPMFaceLimiter) may require a pointer
   to a physics class analysis class in order to perform limiting in
   characteristic variables.
*/
class GodunovUtilities
{
public:
  /// Constructor
  /**
   */
  GodunovUtilities();

  /// Destructor
  /**
   */
  ~GodunovUtilities();

  /// Define the object
  /**
   */
  void define(ProblemDomain& a_domain,
              const Real&    a_dx);

  /// Compute the slope flattening coefficients
  /**
     Compute the slope flattening coefficients, a_flattening, using the
     primitive variables, a_W, within a_box.
  */
  void computeFlattening(FArrayBox&       a_flattening,
                         const FArrayBox& a_W,
                         const Interval&  a_velInt,
                         const int&       a_presInd,
                         const Real&      a_smallPres,
                         const int&       a_bulkModulusInd,
                         const Box&       a_box);

  /// Apply the flattening to slopes.
  /**
   */
  void applyFlattening(FArrayBox&       a_dW,
                       const FArrayBox& a_flat,
                       const Box&       a_box);

  /// Compute componentwise van Leer slopes.
  /**
     Given cell averages W, compute van Leer slopes dW on a
     component-by-component basis.
   */
  void vanLeerSlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_numSlopes,
                     const bool&      a_useLimiting,
                     const int&       a_dir,
                     const Box&       a_box);

  /// Compute fourth-order slopes.
  /**
     Given cell averages W and van Leer slopes dWvL, compute fourth-order
     slopes dW4. Limiting is performed in a separate pass.
  */
  void fourthOrderSlopes(FArrayBox&       a_dW4,
                         const FArrayBox& a_W,
                         const FArrayBox& a_dWvL,
                         const int&       a_numSlopes,
                         const int&       a_dir,
                         const Box&       a_box);

  /// Compute slopes (dW- and dW+) using one sided differences
  /**
   */
  void oneSidedDifferences(FArrayBox&       a_dWMinus,
                           FArrayBox&       a_dWPlus,
                           const FArrayBox& a_W,
                           const int&       a_dir,
                           const Box&       a_box);

  /// Compute slopes (dW (center), dW-, and dW+)
  /**
   */
  void slopes(FArrayBox&       a_dWCent,
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
              const int&       a_hasHi);

  /// van Leer slope limiter.
  /**
     On input, dW contains the centered, unlimited slopes, and
     dW{Minus,Plus} contain the one-sided slopes from the minus, plus sides.
     On output, dW contains the limited slopes.
     slopes dW4. Limiting is performed in a separate pass.
  */
  void slopeLimiter(FArrayBox&       a_dW,
                    const FArrayBox& a_dWMinus,
                    const FArrayBox& a_dWPlus,
                    const int&       a_numSlopes,
                    const Box&       a_box);

  /// PLM normal predictor.
  /**
     Compute the increments in the characteristic amplitudes
  */
  void PLMNormalPred(FArrayBox&       a_dWCharMinus,
                     FArrayBox&       a_dWCharPlus,
                     const FArrayBox& a_dWChar,
                     const FArrayBox& a_Lambda,
                     const Real&      a_dtbydx,
                     const Box&       a_box);

  /// PPM face-centered interpolant.
  /**
     Given the cell average W and the van Leer slope dWvL, compute fourth-order
     accurate face-centered values WFace given by differentiating the indefinite
     integral. Limiting is performed in a separate pass.
  */
  void PPMFaceValues(FArrayBox&            a_WFace,
                     const FArrayBox&      a_W,
                     const int&            a_numSlopes,
                     const bool&           a_useLimiting,
                     const int&            a_dir,
                     const Box&            a_box,
                     const Real&           a_time,
                     const GodunovPhysics* a_physPtr = NULL);

  /// PPM Limiter.
  /**
     On input, dWMinus and dWPlus are the differences between the face
     values on the minus and plus sides of cells and the average in the cell.
     The PPM limiter is applied to these values to obtain a monotone
     interpolant in the cell.
  */
  void PPMLimiter(FArrayBox& a_dWMinus,
                  FArrayBox& a_dWPlus,
                  const int& a_numSlopes,
                  const Box& a_box);

  /// PPM normal predictor.
  /**
     On input, dW{Minus,Plus}, contain the characteristic
     expansions of the differences between the {minus, plus} face values
     and the cell average. On output, dW{Minus,Plus} contain the
     characteristic amplitudes of the corrections required to compute
     the normal predictor.
  */
  void PPMNormalPred(FArrayBox&       a_dWMinus,
                     FArrayBox&       a_dWPlus,
                     const FArrayBox& a_Lambda,
                     const Real&      a_dtbydx,
                     const int&       a_numSlopes,
                     const Box&       a_box);

  /// Compute face-centered velocity divergence.
  /**
     Returns face-centered velocity divergence on for faces in the
     direction a_dir. The velocities are the components a_velInterval
     of a_W.
  */
  void divVel(FArrayBox&       a_divVel,
              const FArrayBox& a_W,
              const Interval&  a_velInt,
              const int&       a_dir,
              const Box&       a_box);

  /// Compute face-centered artificial viscosity flux.
  /**
     Increments face-centered flux in the a_dir direction with quadratic
     artificial viscosity.
  */
  void artificialViscosity(FArrayBox&       a_F,
                           const FArrayBox& a_U,
                           const FArrayBox& a_divVel,
                           const Real&      a_scale,
                           const int&       a_dir,
                           const Box&       a_box);


protected:
  // Problem domain and grid spacing
  ProblemDomain m_domain;
  Real          m_dx;

  // Has this object been defined
  bool          m_isDefined;


private:
  // We may allow copying and assignment later.
  // Disallowed for all the usual reasons
  void operator=(const GodunovUtilities& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  GodunovUtilities(const GodunovUtilities& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
