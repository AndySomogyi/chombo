#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _AMRLEVELUPWIND_H_
#define _AMRLEVELUPWIND_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "AMRLevel.H"
#include "CoarseAverage.H"
#include "FineInterp.H"
#include "LevelFluxRegister.H"

#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "UsingNamespace.H"

/// AMR simple Upwind 
/**
 */
class AMRLevelUpwind : public AMRLevel
{
public:
  /// Constructor
  /**
   */
  AMRLevelUpwind();

  /// Destructor
  /**
   */
  virtual ~AMRLevelUpwind();

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

  /// Set the advection velocity
  /**
   */
  void advectionVel(const RealVect& a_vel);

  /// accessor for the advection velocity
  /**
   */
  const RealVect& advectionVel() const;

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

  ///
  const LevelData<FArrayBox>& getStateNew() const;

  ///
  const LevelData<FArrayBox>& getStateOld() const;


protected:
  // Create a load-balanced DisjointBoxLayout from a collection of Boxes
  DisjointBoxLayout loadBalance(const Vector<Box>& a_grids);

  // Setup menagerie of data structures
  void levelSetup();

  // Get the next coarser level
  AMRLevelUpwind* getCoarserLevel() const;

  // Get the next finer level
  AMRLevelUpwind* getFinerLevel() const;

  // Conserved state, U, at old and new time
  LevelData<FArrayBox> m_UOld,m_UNew;

  /// Advection velocity (constant in space and time
  RealVect m_advectionVel;

  // CFL number
  Real m_cfl;

  // Grid spacing
  Real m_dx;

  /// Averaging from fine to coarse level
  CoarseAverage m_coarseAverage;

  /// Interpolation from coarse to fine level. 
  FineInterp m_fineInterp;

  // New time step
  Real m_dtNew;

  // Number of conserved states
  int m_numStates;

  // Names of conserved states
  Vector<std::string> m_stateNames;

  // Number of ghost cells (in each direction)
  int m_numGhost;

  // Physical dimension of the longest side of the domain
  Real m_domainLength;

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


private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelUpwind& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelUpwind(const AMRLevelUpwind& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
