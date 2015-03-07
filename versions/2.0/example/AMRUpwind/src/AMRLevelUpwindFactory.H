#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _AMRLEVELUPWINDFACTORY_H_
#define _AMRLEVELUPWINDFACTORY_H_

#include "AMRLevelFactory.H"
#include "AMRLevelUpwind.H"
#include "UsingNamespace.H"

/// AMR Godunov factory
/**
 */
class AMRLevelUpwindFactory : public AMRLevelFactory
{
public:
  /// Null constructor
  /**
   */
  AMRLevelUpwindFactory();

  /// Virtual constructor
  /**
   */
  virtual AMRLevel* new_amrlevel() const;

  /// Destructor
  /**
   */
  virtual ~AMRLevelUpwindFactory();

  /// CFL number
  /**
   */
  virtual void CFL(Real a_cfl);

  /// advection velocity
  /**
   */
  virtual void advectionVel(const RealVect& a_advVel);
  
  /// Physical dimension of the longest side of the domain
  /**
   */
  virtual void domainLength(Real a_domainLength);

  ///
  /**
   */
  virtual void verbosity(const int& verbosity);

  /// Refinement threshold
  /**
   */
  virtual void refinementThreshold(Real a_refineThresh);

  /// Tag buffer size
  /**
   */
  void tagBufferSize(int a_tagBufferSize);

  /// Initial dt multiplier
  /**
   */
  void initialDtMultiplier(Real a_initialDtMultiplier);


  /// Check that everything is defined
  /**
   */
  bool isDefined() const;


protected:
  // Some default values
  void setDefaultValues();

  int m_verbosity;

  // CFL number
  Real m_cfl;
  bool m_cflSet;

  // Physical dimension of the longest side of the domain
  Real m_domainLength;
  bool m_domainLengthSet;

  // Refinement threshold for gradient
  Real m_refineThresh;
  bool m_refineThreshSet;

  // Tag buffer size
  int  m_tagBufferSize;
  bool m_tagBufferSizeSet;

  // Initial dt multiplier
  Real m_initialDtMultiplier;
  bool m_initialDtMultiplierSet;

  // advection velocity
  RealVect m_advectionVel;
  bool m_advectionVelSet;

private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelUpwindFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelUpwindFactory(const AMRLevelUpwindFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
