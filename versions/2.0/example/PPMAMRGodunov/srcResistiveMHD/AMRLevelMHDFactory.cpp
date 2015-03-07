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

#include "AMRLevelMHDFactory.H"
#include "AMRLevelMHD.H"

AMRLevelMHDFactory::AMRLevelMHDFactory()
{
  m_godunovPhysics = NULL;
  m_isDefined = false;
}

AMRLevelMHDFactory::~AMRLevelMHDFactory()
{
  if (m_godunovPhysics != NULL)
  {
    delete m_godunovPhysics;
    m_godunovPhysics = NULL;
  }

  m_isDefined = false;
}

void AMRLevelMHDFactory::define(const Real&                 a_cfl,
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
				const BaseHeatSolver*       a_diffusionSolverPtr,
				const DomainGhostBC&        a_dombcB,
				const DomainGhostBC&        a_dombcV,
				const DomainGhostBC&        a_dombcT,
				const Real&                 a_gamma,
				const Real&                 a_mu,
				const Real&                 a_eta,
				const Real&                 a_kappa)
{
  // Store the CFL number
  m_cfl = a_cfl;

  // Store the physical dimension of the longest side of the domain
  m_domainLength = a_domainLength;

  // Store the verbosity of the object
  m_verbosity = a_verbosity;

  // Store the refinement threshold for gradient
  m_refineThresh = a_refineThresh;

  // Store the tag buffer size
  m_tagBufferSize = a_tagBufferSize;

  // Store the initial dt multiplier
  m_initialDtMultiplier = a_initialDtMultiplier;

  // Delete any existing physics object
  if (m_godunovPhysics != NULL)
  {
    delete m_godunovPhysics;
    m_godunovPhysics = NULL;
  }

  // Store the object that supplies the physics needed by the integrator
  // (used as a factory)
  m_godunovPhysics = a_godunovPhysics->new_godunovPhysics();

  // Store the order of the normal predictor (1 -> PLM, 2 -> PPM)
  m_normalPredOrder = a_normalPredOrder;

  // Store the slope computation parameters
  m_useFourthOrderSlopes = a_useFourthOrderSlopes;
  m_usePrimLimiting      = a_usePrimLimiting;
  m_useCharLimiting      = a_useCharLimiting;
  m_useFlattening        = a_useFlattening;

  // Artificial viscosity coefficient must be greater than zero
  CH_assert(!a_useArtificialViscosity || (a_artificialViscosity >= 0.0));

  // Store the artificial viscosity flag and coefficient
  m_useArtificialViscosity = a_useArtificialViscosity;
  m_artificialViscosity    = a_artificialViscosity;

  //Filter magnetic field
  m_doFilterBField = a_doFilterBField;

  // Setting the base heat solver for diffusion terms
  m_diffusionSolverPtr=a_diffusionSolverPtr->new_heatSolver();

  // Domain ghost BCs for B, V, and T
  m_dombcB=a_dombcB;
  m_dombcV=a_dombcV;
  m_dombcT=a_dombcT;

  // Plasma properties
  m_gamma=a_gamma;
  m_mu = a_mu;
  m_kappa = a_kappa; 
  m_eta = a_eta; 
  // The object is defined
  m_isDefined = true;
}

// Virtual constructor
AMRLevel* AMRLevelMHDFactory::new_amrlevel() const
{
  // Make sure everything is defined
  CH_assert(isDefined());

  // Create a new AMRLevelMHD
  AMRLevelMHD* amrGodPtr = new AMRLevelMHD();

  // Define the new object
  amrGodPtr->defineParams(m_cfl,
                          m_domainLength,
                          m_verbosity,
                          m_refineThresh,
                          m_tagBufferSize,
                          m_initialDtMultiplier,
                          m_godunovPhysics,
                          m_normalPredOrder,
                          m_useFourthOrderSlopes,
                          m_usePrimLimiting,
                          m_useCharLimiting,
                          m_useFlattening,
                          m_useArtificialViscosity,
                          m_artificialViscosity,
                          m_doFilterBField,
			  m_diffusionSolverPtr,
			  m_dombcB,
			  m_dombcV,
			  m_dombcT,
			  m_gamma,
			  m_mu,
			  m_eta,
			  m_kappa);

  // Return it
  return (static_cast <AMRLevel*> (amrGodPtr));
}

// Check that everything is defined
bool AMRLevelMHDFactory::isDefined() const
{
  return m_isDefined;
}
