#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef LEVEL_WAVEOPERATOR_H
#define LEVEL_WAVEOPERATOR_H

// Original design:                     Phil Colella
// Original implementation:
// Implementation date:

#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "QuadCFInterp.H"
#include "ExtrapFillPatch.H"
#include "LevelFluxRegister.H"
#include "ProblemDomain.H"
#include "WaveOpF_F.H"
#include "UsingNamespace.H"

/// Level Wave Operator - a class for evaluating the Laplacian.

/**
 */
class LevelWaveOperator
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  LevelWaveOperator();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  ~LevelWaveOperator();

  /// Actual constructor.
  /**
     Inside the routine, we cast away const-ness on the data members
     for the assignment.  The arguments passed in are maintained const
     (coding standards).  a_refine is the refinement ratio between
     this level and the next coarser level.  For the coarsest level, an
     empty DisjointBoxLayout is passed in for coarserDisjointBoxLayout.
   */
  void define(const DisjointBoxLayout&  a_thisDisjointBoxLayout,
              const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
              const ProblemDomain&      a_domain,
              const int&                a_refineCoarse,
              const int&                a_numStates,
              const Real&               a_dx,
              const bool&               a_hasCoarser,
              const bool&               a_hasFiner);

  /// Evaluate the operator at a given time.
  /**
     For the coarsest level, phiCoarseOld and phiCoarseNew are empty
     LevelData<FArrayBox> objects.  Also, either phiCoarseOld or
     phiCoarseNew might be empty to indicate that t(nf) = t(nc) the
     one grid is at the current time and interpolation is not required
     for boundary condition generation.  phi must be defined on the same
     DisjointBoxLayouts as were used in define().  Coarse flux register
     is flux register with the next coarser level.  Fine flux register is
     the flux register with the next finer level.  To the finer level FR,
     this level is the coarse level.  To the coarser level FR, this level
     is the fine level.
   */
  void eval(LevelData<FArrayBox>&       a_Phi,
            LevelData<FArrayBox>&       a_LOfPhi,
            LevelFluxRegister&          a_finerFluxRegister,
            LevelFluxRegister&          a_coarserFluxRegister,
            const LevelData<FArrayBox>& a_phiCoarseOld,
            const Real&                 a_TCoarseOld,
            const LevelData<FArrayBox>& a_phiCoarseNew,
            const Real&                 a_TCoarseNew,
            Real                        a_time,
            Real                        a_fluxweight);
  /// Average down from the fine grid to the coarse grid.
  /**
     This operator averages the solution from the fine grid onto the
     grid at the next coarser AMR level. Our averaging operator given as

     phiCoarse = <phi> - h^2 <L(phi)>_2 / 24 .

     Here <.> denotes the arithmetic average over the fine cells covered
     by each coarse cell, L(phi) is the (2*SpaceDim + 1)-point
     discretization of the Laplacian, with ghost cells set by QuadCFInterp,
     and <.>_2 is the average over the fine subblock of size 2 located in the
     center of the coarse cell.

  */

  void avgdown(
            LevelData<FArrayBox>&       a_Phi,
            LevelData<FArrayBox>&       a_phiCoarse);

protected:
  // define() has been called
  bool m_defined;

  // Box layout for this level
  DisjointBoxLayout m_grids;

  // Interpolator for filling in ghost cells from the next coarser level
  QuadCFInterp m_patcher;

  // Grid spacing
  Real m_dx;

  // Problem domain - index space for this level
  ProblemDomain m_domain;

  // Refinement ratio between this level and the next coarser
  int m_refineCoarse;

  // Number of conservative fluxes and other face centered quantities
  int m_numFluxes;

  int m_numGhost;
  int m_numFields;

  // Flags for the existence of coarser and finer levels, respectively
  bool m_hasCoarser;
  bool m_hasFiner;

private:
  // Disallowed for all the usual reasons
  void operator=(const LevelWaveOperator& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  LevelWaveOperator(const LevelWaveOperator& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
