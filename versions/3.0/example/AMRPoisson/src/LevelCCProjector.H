#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELCELLCENTEREDPROJECTOR_H__
#define _LEVELCELLCENTEREDPROJECTOR_H__
#include "DisjointBoxLayout.H"
#include "Box.H"
#include "REAL.H"
#include "LevelData.H"
#include "FluxBox.H"
#include "FArrayBox.H"

#include "MACProjBC.H"
#include "LevelMACProjector.H"
#include "AMRMultiGrid.H"

///
/**
   Class to project a face-centered. velocity field on a level.
   u -= G (DG^-1)(D u)
   This is designed so you can use the same solver for all
   your level mac projectors (and your cell-centered projectors too, if you need
   them).   Boundary conditions at the domain boundary come in with the solver.
   Boundary conditions at coarse-fine boundaries are dirichlet with an input
   phi field sent in for the boundary condition.
*/
class LevelCCProjector
{
public:

  ///
  ~LevelCCProjector()
  {
  }

  ///
  /**
     a_grids:          boxes on the level
     a_dx:             grid spacing
     a_solver:         solver for Lphi = div (u) solve
     a_lbase is the level in a_solver (from the define) that phi
     lives on.
     maxAMRLevels is the number of levels in the solver.
     refToCoar only matters if lbase > 0
  */
  LevelCCProjector(const DisjointBoxLayout&                                    a_grids,
                   const ProblemDomain&                                        a_domain,
                   const RefCountedPtr< LevelMACProjector >&                   a_macProj,
                   const Real&                                                 a_dx,
                   const int&                                                  a_lbase, 
                   const int&                                                  a_numLevels,
                   int a_refToCoar = -1)
  {
    m_grids     = a_grids;
    m_domain    = a_domain;
    m_macProj   = a_macProj;
    m_dx        = a_dx;
    m_lbase     = a_lbase;
    m_numLevels = a_numLevels;
    m_refToCoar = a_refToCoar;

    m_macVel.define(m_grids, 1, IntVect::Unit);
    m_macGph.define(m_grids, 1, IntVect::Unit);
  }


  ///
  /**
     velocity--input and output as divergence-free
     assumed to have correct values at domain and coarse-fine
     boundary faces at input.
     gradient--output.  the pure gradient component of input vel.
     phiCoarse-input. boundary conditions for coarser levels, if they
     exist.
  */
  void
  project(LevelData<FArrayBox>&       a_velocity,
          LevelData<FArrayBox>&       a_gradient,
          LevelData<FArrayBox>&       a_phi,
          const LevelData<FArrayBox>* a_velCoarse = NULL,
          const LevelData<FArrayBox>* a_phiCoarse = NULL);


  ///for convergence tests
  void
  divergence(LevelData<FArrayBox>&       a_divergence,
             LevelData<FArrayBox>&       a_velocity,
             LevelData<FArrayBox>*       a_velCoarse = NULL);


protected:

  DisjointBoxLayout                     m_grids;
  ProblemDomain                         m_domain;
  RefCountedPtr< LevelMACProjector >    m_macProj;
  LevelData<FluxBox>                    m_macVel;
  LevelData<FluxBox>                    m_macGph;
  Real                                  m_dx;
  int                                   m_lbase; 
  int                                   m_numLevels;
  int                                   m_refToCoar;

private:
  
  LevelCCProjector(const LevelCCProjector& a_input)
  {
    MayDay::Error("We disallow copy construction for objects with pointered data");
  }

  void operator=(const LevelCCProjector& a_input)
  {
    MayDay::Error("We disallow assignment for objects with pointered data");
  }

  LevelCCProjector()
  {
    MayDay::Error("Weak construction is bad.");
  }
};
//reusable functions made external

///
extern void
averageCellsToFaces(LevelData<FluxBox>&           a_macVeloc,
                    const LevelData<FArrayBox>&   a_cellVeloc,
                    const DisjointBoxLayout&      a_grids,
                    const ProblemDomain&          a_domain);

///
extern void
averageCellsToFaces(FluxBox&               a_macVel,
                    const FArrayBox&       a_cellVel,
                    const Box&             a_grid,
                    const ProblemDomain&   a_domain);

///
extern void
averageFacesToCells(LevelData<FArrayBox>&         a_cellData,
                    const LevelData<FluxBox>&     a_macData,
                    const DisjointBoxLayout&      a_grids,
                    const ProblemDomain&          a_domain);

///
extern void
averageFacesToCells(FArrayBox&             a_cellData,
                    const FluxBox&         a_fluxData,
                    const Box&             a_grid,
                    const ProblemDomain&   a_domain);

///
extern void
extrapolateToDomainBoundaries(LevelData<FluxBox>&      a_macData,
                              const DisjointBoxLayout& a_grids,
                              const ProblemDomain&     a_domain);


///
extern void
extrapolateToDomainBoundaries(FluxBox&               a_faceData,
                              const Box&             a_grid,
                              const ProblemDomain&   a_domain,
                              int                    a_faceDir,
                              Side::LoHiSide         a_side);


#endif
