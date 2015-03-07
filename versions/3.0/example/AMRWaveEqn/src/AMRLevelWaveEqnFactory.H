#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _AMR_LEVEL_WAVE_EQN_FACTORY_H_
#define _AMR_LEVEL_WAVE_EQN_FACTORY_H_

#include "AMRLevelFactory.H"
#include "WaveIBC.H"
#include "AMRLevelWaveEquation.H"
#include "UsingNamespace.H"

/// AMR Wave Equation factory
/**
 */
class AMRLevelWaveEqnFactory : public AMRLevelFactory
{
public:
  /// Null constructor
  /**
   */
  AMRLevelWaveEqnFactory();

  /// Virtual constructor
  /**
   */
  virtual AMRLevel* new_amrlevel() const;

  /// Destructor
  /**
   */
  virtual ~AMRLevelWaveEqnFactory();

  /// CFL number
  /**
   */
  virtual void CFL(Real a_cfl);

  /// Physical dimension of the longest side of the domain
  /**
   */
  virtual void domainLength(Real a_domainLength);

  void x0(const Real& a_x0);
  void IBC(PhysIBC* a_w0){m_wave_ibc = a_w0->new_physIBC();};
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
  Real m_r0;
  Real m_cfl;
  bool m_cflSet;
  PhysIBC* m_wave_ibc;

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

private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelWaveEqnFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelWaveEqnFactory(const AMRLevelWaveEqnFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
