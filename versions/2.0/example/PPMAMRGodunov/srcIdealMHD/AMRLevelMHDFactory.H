#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _AMRLEVELMHDFACTORY_H_
#define _AMRLEVELMHDFACTORY_H_

#include "AMRLevelFactory.H"

#include "AMRLevelMHD.H"

/// factory object for creating AMRLevelMHD objects
/**
 */
class AMRLevelMHDFactory : public AMRLevelFactory
{
public:
  /// Null constructor
  /**
   */
  AMRLevelMHDFactory();

  /// Destructor
  /**
   */
  virtual ~AMRLevelMHDFactory();

  /// Define the parameters of the factory
  /**
   */
  virtual void define(const Real&                 a_cfl,
                      const Real&                 a_domainLength,
                      const int&                  a_verbosity,
                      const Real&                 a_refineThresh,
                      const int&                  a_tagBufferSize,
                      const Real&                 a_initialDtMultiplier,
                      const GodunovPhysics* const a_godunovPhysics,
                      const int&                  a_normalPredOrder,
                      const bool&                 a_useFourthOrderSlopes,
                      const bool&                 a_usePrimLimiting,
                      const bool&                 a_useCharLimiting,
                      const bool&                 a_useFlattening,
                      const bool&                 a_useArtificialViscosity,
                      const Real&                 a_artificialViscosity,
                      const bool&                 a_doFilterBField,
                      const bool&                 a_projectB);

  /// Create a new AMRLevel (for MHD)
  /**
   */
  virtual AMRLevel* new_amrlevel() const;

  /// Check that this object has been defined
  /**
   */
  bool isDefined() const;


protected:
  // CFL number
  Real m_cfl;

  // Physical dimension of the longest side of the domain
  Real m_domainLength;

  // Verbosity of the object
  int m_verbosity;

  // Refinement threshold for gradient
  Real m_refineThresh;

  // Tag buffer size
  int m_tagBufferSize;

  // Initial dt multiplier
  Real m_initialDtMultiplier;

  // Pointer to the object that supplies the physics needed by the integrator
  // (used as a factory here)
  const GodunovPhysics* m_godunovPhysics;

  // Order of the normal predictor (1 -> PLM, 2 -> PPM)
  int m_normalPredOrder;

  // Use 4th other slope computations (otherwise, use 2nd order)
  bool m_useFourthOrderSlopes;

  // Do slope limiting in the primitive or characteristic variables, respect.
  bool m_usePrimLimiting;
  bool m_useCharLimiting;

  // Do slope flattening - MUST BE USING 4th order slopes
  bool m_useFlattening;

  // Apply artificial viscosity of a set value
  bool m_useArtificialViscosity;
  Real m_artificialViscosity;

  // Apply filter to magnetic field
  bool m_doFilterBField;

  // project magnetic field
  bool m_projectB;

  // Has this object been defined
  bool m_isDefined;


private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelMHDFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelMHDFactory(const AMRLevelMHDFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
