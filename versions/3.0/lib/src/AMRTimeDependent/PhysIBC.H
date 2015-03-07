#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//
//  PhysIBC.H
// ============
//
// Virtual base class through which a user specifies the initial and boundary
// conditions for a hyperbolic system of PDEs.
//

#ifndef _PHYSIBC_H_
#define _PHYSIBC_H_

#include "FArrayBox.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

/// Physical/domain initial and boundary conditions
/**
   Virtual base class through which a user specifies the initial and
   boundary conditions for a hyperbolic system of PDEs.
 */
class PhysIBC
{
public:
  /// Constructor
  /**
   */
  PhysIBC();

  /// Destructor
  /**
   */
  virtual ~PhysIBC();

  /// Define the object
  /**
     Set the problem domain index space and the grid spacing for this
     initial and boundary condition object.
   */
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx);

  /// Factory method - this object is its own factory
  /**
     Return a point to a new PhysIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used).
   */
  virtual PhysIBC* new_physIBC() = 0;

  /// Set up initial conditions
  /**
   */
  virtual void initialize(LevelData<FArrayBox>& a_U) = 0;

  /// Set boundary fluxes
  /**
   */
  virtual void primBC(FArrayBox&            a_WGdnv,
                      const FArrayBox&      a_Wextrap,
                      const FArrayBox&      a_W,
                      const int&            a_dir,
                      const Side::LoHiSide& a_side,
                      const Real&           a_time) = 0;

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one sided difference
     approximations.  If this function doesn't change them they will be
     used for the slopes at the boundaries.
   */
  virtual
  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
                     const Real&      a_time) = 0;

  /// Adjust boundary fluxes to account for artificial viscosity
  /**
   */
  virtual
  void artViscBC(FArrayBox&       a_F,
                 const FArrayBox& a_U,
                 const FArrayBox& a_divVel,
                 const int&       a_dir,
                 const Real&      a_time) = 0;


protected:

  /// This function is called by primBC to get boundary faces.
  virtual
  void getBoundaryFaces(/// face-centered box of boundary faces to fill in
                        Box&                   a_boundaryBox,
                        /// face-centered box
                        const Box&             a_dataFaceBox,
                        /// which direction
                        const int&             a_dir,
                        /// which side
                        const Side::LoHiSide&  a_side);

  // Has this object been defined
  bool          m_isDefined;

  // The current level's problem domain
  ProblemDomain m_domain;

  // The current level's grid spacing
  Real          m_dx;


private:
  // Disallowed for all the usual reasons
  void operator=(const PhysIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  PhysIBC(const PhysIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
