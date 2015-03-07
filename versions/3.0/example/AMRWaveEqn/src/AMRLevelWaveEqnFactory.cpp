#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "AMRLevel.H"
#include "AMRLevelWaveEquation.H"
#include "AMRLevelWaveEqnFactory.H"

AMRLevelWaveEqnFactory::AMRLevelWaveEqnFactory()
{
  setDefaultValues();
}

// Virtual constructor
AMRLevel* AMRLevelWaveEqnFactory::new_amrlevel() const
{
  // Make sure everything is defined
 CH_assert(isDefined());

  // Create a new AMRLevelWaveEquation
  AMRLevelWaveEquation* amrWavePtr = new AMRLevelWaveEquation();

  // Set up new object
  amrWavePtr->CFL(m_cfl);
  amrWavePtr->domainLength(m_domainLength);
  amrWavePtr->refinementThreshold(m_refineThresh);
  amrWavePtr->tagBufferSize(m_tagBufferSize);
  amrWavePtr->initialDtMultiplier(m_initialDtMultiplier);
  amrWavePtr->verbosity(m_verbosity);
  amrWavePtr->IBC(m_wave_ibc);
  // Return it
  return (static_cast <AMRLevel*> (amrWavePtr));
}

AMRLevelWaveEqnFactory::~AMRLevelWaveEqnFactory()
{
}

// CFL number
void AMRLevelWaveEqnFactory::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
  m_cflSet = true;
}

void AMRLevelWaveEqnFactory::verbosity(const int& a_verbosity)
{
  m_verbosity = a_verbosity;
}

// Physical dimension of the longest side of the domain
void AMRLevelWaveEqnFactory::domainLength(Real a_domainLength)
{
  m_domainLength = a_domainLength;
  m_domainLengthSet = true;
}

// Refinement threshold
void AMRLevelWaveEqnFactory::refinementThreshold(Real a_refineThresh)
{
  m_refineThresh = a_refineThresh;
  m_refineThreshSet = true;
}

// Tag buffer size
void AMRLevelWaveEqnFactory::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
  m_tagBufferSizeSet = true;
}

// Initial dt multiplier
void AMRLevelWaveEqnFactory::initialDtMultiplier(Real a_initialDtMultiplier)
{
  m_initialDtMultiplier = a_initialDtMultiplier;
  m_initialDtMultiplierSet = true;
}

// Check that everything is defined
bool AMRLevelWaveEqnFactory::isDefined() const
{
  return (m_cflSet &&
          m_domainLengthSet &&
          m_refineThreshSet &&
          m_tagBufferSizeSet &&
          m_initialDtMultiplierSet);
}

// Some default values
void AMRLevelWaveEqnFactory::setDefaultValues()
{
  CFL(0.8);
  domainLength(1.0);
  refinementThreshold(0.2);
  tagBufferSize(2);
  initialDtMultiplier(0.1);
  m_verbosity = 0;
}
