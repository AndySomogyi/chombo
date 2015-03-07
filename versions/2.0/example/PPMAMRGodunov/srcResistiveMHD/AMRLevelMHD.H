#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _AMRLEVELMHD_H_
#define _AMRLEVELMHD_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "AMRLevel.H"
#include "CoarseAverage.H"
#include "FineInterp.H"
#include "LevelFluxRegister.H"
#include "LevelSolver.H"
#include "PoissonOp.H"
#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"

#include "LevelGodunov.H"
#include "PhysIBC.H"
#include "VCHelmholtzOp.H"
#include "BaseHeatSolver.H"

/// AMR Godunov
/**
 */
class AMRLevelMHD: public AMRLevel
{
public:
  /// Constructor
  /**
   */
  AMRLevelMHD();

  /// Destructor
  /**
   */
  virtual ~AMRLevelMHD();

  /// Define the parameters the object needs
  /**
   */
  void defineParams(const Real&                 a_cfl,
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
		    const Real&                 a_kappa);

  /// This instance should never get called - historical
  /**
   */
  virtual void define(AMRLevel*  a_coarserLevelPtr,
                      const Box& a_problemDomain,
                      int        a_level,
                      int        a_refRatio);

  /// Define new AMR level
  /**
   */
  virtual void define(AMRLevel*            a_coarserLevelPtr,
                      const ProblemDomain& a_problemDomain,
                      int                  a_level,
                      int                  a_refRatio);

  /// Advance by one timestep
  /**
   */
  virtual Real advance();

  /// Things to do after a timestep
  /**
   */
  virtual void postTimeStep();

  /// Create tags for regridding
  /**
   */
  virtual void tagCells(IntVectSet& a_tags);

  /// Create tags at initialization
  /**
   */
  virtual void tagCellsInit(IntVectSet& a_tags);

  /// Set up data on this level after regridding
  /**
   */
  virtual void regrid(const Vector<Box>& a_newGrids);

  /// Initialize grids
  /**
   */
  virtual void initialGrid(const Vector<Box>& a_newGrids);

  /// Initialize data
  /**
   */
  virtual void initialData();

  /// Things to do after initialization
  /**
   */
  virtual void postInitialize();

#ifdef CH_USE_HDF5
  /// Write checkpoint header
  /**
   */
  virtual void writeCheckpointHeader(HDF5Handle& a_handle) const;

  /// Write checkpoint data for this level
  /**
   */
  virtual void writeCheckpointLevel(HDF5Handle& a_handle) const;

  /// Read checkpoint header
  /**
   */
  virtual void readCheckpointHeader(HDF5Handle& a_handle);

  /// Read checkpoint data for this level
  /**
   */
  virtual void readCheckpointLevel(HDF5Handle& a_handle);

  /// Write plotfile header
  /**
   */
  virtual void writePlotHeader(HDF5Handle& a_handle) const;

  /// Write plotfile data for this level
  /**
   */
  virtual void writePlotLevel(HDF5Handle& a_handle) const;
#endif

  /// Returns the dt computed earlier for this level
  /**
   */
  virtual Real computeDt();

  /// Compute dt using initial data
  /**
   */
  virtual Real computeInitialDt();

  ///
  const LevelData<FArrayBox>& getStateNew() const;

  ///
  const LevelData<FArrayBox>& getStateOld() const;

  ///
  bool allDefined() const;


protected:
  // Create a load-balanced DisjointBoxLayout from a collection of Boxes
  DisjointBoxLayout loadBalance(const Vector<Box>& a_grids);

  // Setup menagerie of data structures
  void levelSetup();

  // Compute Divergence
  void DivergenceB(LevelData<FArrayBox>&      a_divb,
                   const LayoutData<FluxBox>& a_WHalf);

  // Project the magnetic field
  void ProjectBField(LayoutData<FluxBox>&        a_WHalf,
                     const LevelData<FArrayBox>& a_phi);

  //Filter the magnetic field - Old function
  void FilterBFieldOld(LevelData<FArrayBox> &a_U);

  //Filter the magnetic field
  void FilterBField(LevelData<FArrayBox>& a_U,
                    LevelFluxRegister&    a_finerFluxRegister,
                    LevelFluxRegister&    a_coarserFluxRegister);

  // source is div(diffusive flux)
  void parabolicFluxSource(LevelData<FArrayBox>&       a_U,
			   LevelData<FArrayBox>&       a_S,
			   const LevelData<FArrayBox>& a_UCoarseOld,
			   const Real&                 a_TCoarseOld,
			   const LevelData<FArrayBox>& a_UCoarseNew,
			   const Real&                 a_TCoarseNew,
			   const DisjointBoxLayout & a_coarseGrids,
			   const int&                a_nrefCoarse,
			   const Real&                 a_time,
			   const Real&                 a_dt);

  //Solution interpolation in time
 void InterpolateSolnInTime(LevelData<FArrayBox>& a_U,
			    const LevelData<FArrayBox> &a_UOld,
			    const LevelData<FArrayBox> &a_UNew,
			    const DisjointBoxLayout & a_grids,
			    const Real& alpha);

  // Velocity from Momemtum before calling diffusion solve
  void computeVelocityFromMomentum(LevelData<FArrayBox> &a_Mom,
				   const LevelData<FArrayBox> &a_UCons,
				   DisjointBoxLayout & a_grids);
  // updateRho: updates density and computes time average of density field
  void updateRho(LevelData<FArrayBox> &a_rhoHalf,
		 LevelData<FArrayBox> &a_U,
		 const LevelData<FArrayBox> &a_dU);

  // RHS for velocity diffusion solve
  void computeRHSVelocity(LevelData<FArrayBox> &a_srcMomentum,
			  const LevelData<FArrayBox> &a_dU,
			  const LayoutData<FluxBox> &a_WHalf);

  //divE_v term in energy equation for temperature diffusion solve
  void diffusiveEnergyTerm(LevelData<FArrayBox>&       a_divEv,
			   LevelData<FArrayBox>&       a_U,
			   LevelFluxRegister&          a_finerFluxRegister,
			   LevelFluxRegister&          a_coarserFluxRegister,
			   const LevelData<FArrayBox>& a_UCoarseOld,
			   const Real&                 a_TCoarseOld,
			   const LevelData<FArrayBox>& a_UCoarseNew,
			   const Real&                 a_TCoarseNew,
			   const DisjointBoxLayout & a_coarseGrids,
			   const int&                a_nrefCoarse,
			   const Real&                 a_time);


// Temperature field at cell centers from conserved vars
  void computeTemperatureFromConserved(LevelData<FArrayBox> &a_Temp,
				     const LevelData<FArrayBox> &a_UCons,
				     DisjointBoxLayout & a_grids);

  // RHS for diffusion solve for Temperature
  void computeRHSTemperature(LevelData<FArrayBox>& a_srcTemperature,
			     const LevelData<FArrayBox>&  a_UNew,
			     const LevelData<FArrayBox>&  a_UOld,
			     const LevelData<FArrayBox>&  a_divEv,
			     const LevelData<FArrayBox> &a_dU);

  //Update of conserved variable with hyperbolic fluxes and diffusion terms
  //for momentum and magnetic field equations
  void updateMomentumAndMagneticFields(LevelData<FArrayBox>&a_UNew,
				   const LevelData<FArrayBox> &a_dU,
				   const LevelData<FArrayBox> &a_etaLaplacianB,
				   const LevelData<FArrayBox> &a_muLaplacianV);

  //Update of energy variable with hyperbolic flux and diffusion term
  void updateEnergy(LevelData<FArrayBox>&a_UNew,
		    const LevelData<FArrayBox> &a_dU,
		    const LevelData<FArrayBox> &a_kappaLaplacianT,
		    const LevelData<FArrayBox> &a_divEv);

// Get the next coarser level
  AMRLevelMHD* getCoarserLevel() const;

  // Get the next finer level
  AMRLevelMHD* getFinerLevel() const;

  // Conserved state, U, at old and new time
  LevelData<FArrayBox> m_UOld,m_UNew;

  // Phi for projection
  LevelData<FArrayBox> m_Phi;

  //Whalf: Time centered face values of primitive variables
   LayoutData<FluxBox> m_whalf;

 //Source term which uses the parabolic fluxes
  LevelData<FArrayBox> m_source;

 // rhs=(divB) for projection
  LevelData<FArrayBox> m_rhs;
  // CFL number
  Real m_cfl;

  // Grid spacing
  Real m_dx;

  // Interpolation from coarse to fine level
  FineInterp m_fineInterp;

  // Averaging from fine to coarse level
  CoarseAverage m_coarseAverage;

  // New time step
  Real m_dtNew;

  // Number of conserved states
  int m_numStates;

  // Number of primitive states
  int m_numPrims;

  // Number of fluxes
  int m_numFluxes;

  // Names of conserved states
  Vector<string> m_stateNames;

  // Names of plot variables
  Vector<string> m_plotNames;

  // Number of ghost cells (in each direction)
  int m_numGhost;

  // Physical dimension of the longest side of the domain
  Real m_domainLength;

  // Level integrator
  LevelGodunov m_levelGodunov;

  // Level solver for Poisson equation for Phi
  LevelSolver m_levelSolver;

  // Flux register
  LevelFluxRegister m_fluxRegister;

  // velocity Flux register
  LevelFluxRegister m_momentumFluxRegister;
  // Magnetic Flux register
  LevelFluxRegister m_magneticFluxRegister;
  // Energy/temperature Flux register
  LevelFluxRegister m_temperatureFluxRegister;

  // Pointer to the class defining the physics of the problem
  GodunovPhysics* m_gdnvPhysics;

  //Poisson operator
  PoissonOp m_poissonOp;

  // Order of the normal predictor (1 -> PLM, 2-> PPM)
  int m_normalPredOrder;

  // Use 4th order slope computations (otherwise, use 2nd order)
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


  // Heat-type equation solver
  BaseHeatSolver* m_bhsolver;

   // Domain ghost BCs
  DomainGhostBC        m_dombcB;
  DomainGhostBC        m_dombcV;
  DomainGhostBC        m_dombcT;

 //Helmholtz operators
  VCHelmholtzOp m_bhelmop;
  VCHelmholtzOp m_vhelmop;
  VCHelmholtzOp m_thelmop;

  // Refinement threshold for gradient
  Real m_refineThresh;

  // Tag buffer size
  int m_tagBufferSize;

  // Flag coarser and finer levels
  bool m_hasCoarser;
  bool m_hasFiner;

  // Grid layout for this level
  DisjointBoxLayout m_grids;

  //Plasma properties
  Real m_mu;
  Real m_eta;
  Real m_kappa;
  Real m_gamma;
  // True if all the parameters for this object are defined
  bool m_paramsDefined;


private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelMHD& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelMHD(const AMRLevelMHD& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
