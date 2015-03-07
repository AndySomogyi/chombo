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
#include "AMRLevelUpwind.H"
#include "AMRLevelUpwindFactory.H"

AMRLevelUpwindFactory::AMRLevelUpwindFactory()
{
  setDefaultValues();
}

// Virtual constructor
AMRLevel* AMRLevelUpwindFactory::new_amrlevel() const
{
  // Make sure everything is defined
 CH_assert(isDefined());

  // Create a new AMRLevelUpwind
  AMRLevelUpwind* amrLevelPTr = new AMRLevelUpwind();

  // Set up new object
  amrLevelPTr->CFL(m_cfl);
  amrLevelPTr->advectionVel(m_advectionVel);
  amrLevelPTr->domainLength(m_domainLength);
  amrLevelPTr->refinementThreshold(m_refineThresh);
  amrLevelPTr->tagBufferSize(m_tagBufferSize);
  amrLevelPTr->initialDtMultiplier(m_initialDtMultiplier);
  amrLevelPTr->verbosity(m_verbosity);
  // Return it
  return (static_cast <AMRLevel*> (amrLevelPTr));
}

AMRLevelUpwindFactory::~AMRLevelUpwindFactory()
{
}

// CFL number
void AMRLevelUpwindFactory::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
  m_cflSet = true;
}

void AMRLevelUpwindFactory::advectionVel(const RealVect& a_advVel)
{
  m_advectionVel = a_advVel;
  m_advectionVelSet = true;
}

void AMRLevelUpwindFactory::verbosity(const int& a_verbosity)
{
  m_verbosity = a_verbosity;
}

// Physical dimension of the longest side of the domain
void AMRLevelUpwindFactory::domainLength(Real a_domainLength)
{
  m_domainLength = a_domainLength;
  m_domainLengthSet = true;
}

// Refinement threshold
void AMRLevelUpwindFactory::refinementThreshold(Real a_refineThresh)
{
  m_refineThresh = a_refineThresh;
  m_refineThreshSet = true;
}

// Tag buffer size
void AMRLevelUpwindFactory::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
  m_tagBufferSizeSet = true;
}

// Initial dt multiplier
void AMRLevelUpwindFactory::initialDtMultiplier(Real a_initialDtMultiplier)
{
  m_initialDtMultiplier = a_initialDtMultiplier;
  m_initialDtMultiplierSet = true;
}


// Check that everything is defined
bool AMRLevelUpwindFactory::isDefined() const
{
  return (m_cflSet &&
          m_domainLengthSet &&
          m_refineThreshSet &&
          m_tagBufferSizeSet &&
          m_initialDtMultiplierSet &&
          m_advectionVelSet);
}

// Some default values
void AMRLevelUpwindFactory::setDefaultValues()
{
  CFL(0.8);
  domainLength(1.0);
  refinementThreshold(0.2);
  tagBufferSize(2);
  initialDtMultiplier(0.1);
  m_advectionVelSet = false;
  m_verbosity = 0;
}
