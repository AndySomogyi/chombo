#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//  PhysIBC.H
// ============
//
// Virtual base class through which a user specifies the initial and boundary
// conditions for a hyperbolic system of PDEs.
//

#ifndef _PHYSIBC_H
#define _PHYSIBC_H

#include "FArrayBox.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "UsingNamespace.H"

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
  PhysIBC(){};

  /// Destructor
  /**
   */
  virtual ~PhysIBC() {};

  /// Define the object
/**
     Set the problem domain index space and the grid spacing for this
     initial and boundary condition object.
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx);

   */
  /// Factory method - this object is its own factory
  /**
     Return a point to a new PhysIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used).
   */
  virtual PhysIBC* new_physIBC() = 0;

  /// Set up initial conditions
  /**
   */
  virtual void initialize(LevelData<FArrayBox>& a_phi,
                          LevelData<FArrayBox>& a_pi,Real a_dx) = 0;

protected:
  // define() has been called
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

#endif
