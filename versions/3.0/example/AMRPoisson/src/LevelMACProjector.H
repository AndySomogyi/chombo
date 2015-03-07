#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _LEVELMACPROJECTOR_H__
#define _LEVELMACPROJECTOR_H__
#include "DisjointBoxLayout.H"
#include "Box.H"
#include "REAL.H"
#include "LevelData.H"
#include "FluxBox.H"
#include "FArrayBox.H"
#include "MACProjBC.H"
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
class LevelMACProjector
{
public:

  ///
  ~LevelMACProjector()
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
  */
  LevelMACProjector(const DisjointBoxLayout&                                     a_grids,
                    const ProblemDomain&                                         a_domain,
                    const RefCountedPtr< AMRMultiGrid<LevelData<FArrayBox> > >&  a_solver,
                    const Real&                                                  a_dx,
                    const int&                                                   a_lbase, 
                    const int&                                                   a_numLevels, 
                    const VelocityBCFunction                                     a_bcVel,
                    const GradientBCFunction                                     a_bcGrad)
  {
    m_bcVel = a_bcVel;
    m_bcGrad= a_bcGrad;
    m_grids      = a_grids;
    m_domain     = a_domain;
    m_solver     = a_solver;
    m_dx         = a_dx;
    m_lbase      = a_lbase;
    m_numLevels  = a_numLevels;
    m_rhs.define(m_grids, 1, IntVect::Zero);
  }


  ///
  /**
     velocity--input and output as divergence-free
     assumed to have correct values at domain and coarse-fine
     boundary faces at input.
     gradient--output.  the pure gradient component of input vel.
     phiCoarse-input. boundary conditions for coarser levels, if they
     exist.   Leave NULL if lbase == 0 (no coarser levels).
  */
  void
  project(LevelData<FluxBox>&         a_velocity,
          LevelData<FluxBox>&         a_gradient,
          LevelData<FArrayBox>&       a_phi,
          const LevelData<FArrayBox>* a_phiCoarse = NULL);

protected:

  void solve(LevelData<FArrayBox>&            a_phi,
             const LevelData<FArrayBox>&      a_rhs,
             const LevelData<FArrayBox>*      a_phiCoarse = NULL);


  DisjointBoxLayout                                     m_grids;
  ProblemDomain                                         m_domain;
  RefCountedPtr< AMRMultiGrid<LevelData<FArrayBox> > >  m_solver;
  Real                                                  m_dx;
  int                                                   m_lbase;
  LevelData<FArrayBox>                                  m_rhs;
  int                                                   m_numLevels;
  VelocityBCFunction                                    m_bcVel;
  GradientBCFunction                                    m_bcGrad;

private:
  LevelMACProjector(const LevelMACProjector& a_input)
  {
    MayDay::Error("We disallow copy construction for objects with pointered data");
  }

  void operator=(const LevelMACProjector& a_input)
  {
    MayDay::Error("We disallow assignment for objects with pointered data");
  }

  LevelMACProjector()
  {
    MayDay::Error("weak construction is bad");
  }
};

//make reusable functions standalone

///
extern void
macGradient(LevelData<FluxBox>&          a_gradPhi,
            const LevelData<FArrayBox>&  a_phi,
            const DisjointBoxLayout&     a_dbl,
            const ProblemDomain&         a_domain,
            const Real&                  a_dx);


///
extern void
macEnforceGradientBC(LevelData<FluxBox>&          a_gradPhi,
                     const LevelData<FArrayBox>&  a_phi,
                     const DisjointBoxLayout&     a_dbl,
                     const ProblemDomain&         a_domain,
                     const Real&                  a_dx,
                     GradientBCFunction           a_bcfunc);

void
macEnforceGradientBC(FluxBox&               a_gradPhi,
                     const FArrayBox&       a_phi,
                     const Box&             a_grids,
                     const ProblemDomain&   a_domain,
                     const Real&            a_dx,
                     GradientBCFunction     a_bcFunc);

///
extern void
macGradient(FluxBox&             a_gradPhi,
            const FArrayBox&     a_phi,
            const Box&           a_box,
            const ProblemDomain& a_domain,
            const Real&          a_dx);

///
extern void
macDivergence(LevelData<FArrayBox>&              a_divVel,
              const LevelData<FluxBox>&          a_velocity,
              const DisjointBoxLayout&           a_dbl,
              const ProblemDomain&               a_domain,
              const Real&                        a_dx);

extern void
macEnforceVelocityBC(LevelData<FluxBox>&                a_velocity,
                     const DisjointBoxLayout&           a_dbl,
                     const ProblemDomain&               a_domain,
                     const Real&                        a_dx,
                     VelocityBCFunction                 a_bcfunc);


extern void
macEnforceVelocityBC(FluxBox&                a_velocity,
                     const Box&              a_dbl,
                     const ProblemDomain&    a_domain,
                     const Real&             a_dx,
                     VelocityBCFunction      a_bcfunc);


///
extern void
macDivergence(FArrayBox&             a_divVel,
              const FluxBox&         a_velocity,
              const Box&             a_box,
              const ProblemDomain&   a_domain,
              const Real&            a_dx);

#endif
