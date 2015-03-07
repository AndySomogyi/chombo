#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _AMR_LEVEL_WAVE_EQUATION_H_
#define _AMR_LEVEL_WAVE_EQUATION_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "AMRLevel.H"
#include "CoarseAverage.H"
#include "FineInterp.H"
#include "LevelFluxRegister.H"
#include "PatchGodunovF_F.H"

#include "LevelWaveOperator.H"
#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "PhysIBC.H"
#include "UsingNamespace.H"

/// AMR Level for second-order wave equation.
/**
 */
class AMRLevelWaveEquation : public AMRLevel
{
  friend class AMRLevelWaveEqnFactory ;
public:
  /// Constructor
  /**
   */
  AMRLevelWaveEquation();

  /// Destructor
  /**
   */
  virtual ~AMRLevelWaveEquation();

  /// Define new AMR level
  /**
   */
  virtual void define(AMRLevel*            a_coarserLevelPtr,
                      const ProblemDomain& a_problemDomain,
                      int                  a_level,
                      int                  a_refRatio);

  /// Advance by one timestep -- RK4
  /**
   */
  virtual Real advance();

  /// Things to do after a timestep -- reflux
  /**
   */
  virtual void postTimeStep();

  /// Create tags for regridding
  /**
   */
  virtual void tagCells(IntVectSet& a_tags) ;

  /// Create tags at initialization
  /**
   */
  virtual void tagCellsInit(IntVectSet& a_tags) ;

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

  /// Compute integral of energy over this and all finer levels
  /**
   */
  Real computeEnergy();

  /// Set the CFL number
  /**
   */
  virtual void CFL(Real a_cfl);

  /// Set the physical dimension of the longest side of the domain
  /**
   */
  virtual void domainLength(Real a_domainLength);

  /// Set the refinement threshold
  /**
   */
  virtual void refinementThreshold(Real a_refineThresh);

  /// Set the tag buffer size
  /**
   */
  void tagBufferSize(int a_tagBufferSize);

  void IBC(PhysIBC* a_ibc) {m_physIBCPtr = a_ibc;};

  /// Accessors
  /**
   */
  LevelData<FArrayBox> * phi(){ return &m_phiNew ; } ;
  LevelData<FArrayBox> * pi() { return &m_piNew  ; } ;

protected:
  // Create a load-balanced DisjointBoxLayout from a collection of Boxes
  DisjointBoxLayout loadBalance(const Vector<Box>& a_grids);

  // Setup menagerie of data structures
  void levelSetup();

  // Get the next coarser level
  AMRLevelWaveEquation* getCoarserLevel() const;

  // Get the next finer level
  AMRLevelWaveEquation* getFinerLevel() const;


private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelWaveEquation& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelWaveEquation(const AMRLevelWaveEquation& a_input)
  {
    MayDay::Error("invalid operator");
  }

protected:

  // member variables


  // lower corner of problem domain in physical space
  Real m_r0 ;


  // Fields phi, pi = d(phi)/dt at old and new time
  LevelData<FArrayBox> m_phiOld,m_phiNew,m_piOld,m_piNew;

  // CFL number
  Real m_cfl;

  // Grid spacing
  Real m_dx;

  // Interpolation from fine to coarse level
  FineInterp m_fineInterp;

  // Averaging from coarse to fine level
  CoarseAverage m_coarseAverage;

  // New time step
  Real m_dtNew;

  // Number of field variables (phi) = number of conserved states (pi).
  int m_numStates;

  // Names of conserved states
  Vector<std::string> m_stateNames;

  // Number of ghost cells (in each direction)
  int m_numGhost;

  // Physical dimension of the longest side of the domain
  Real m_domainLength;

  // Level integrator
  LevelWaveOperator m_levelWaveOperator;

  // Flux register
  LevelFluxRegister m_fluxRegister;

  // Refinement threshold for gradient
  Real m_refineThresh;

  // Tag buffer size
  int m_tagBufferSize;

  // Flag coarser and finer levels
  bool m_hasCoarser;
  bool m_hasFiner;

  DisjointBoxLayout m_grids;

  // Initial / boundary condition interface class.

  PhysIBC* m_physIBCPtr;

  void initializeODE(LevelData<FArrayBox>& a_piOld,LevelData<FArrayBox>& a_pi);
  void updateODE(LevelData<FArrayBox>& a_LOfPhi,LevelData<FArrayBox>& a_pi,
                 Real a_dt);

};

#endif
