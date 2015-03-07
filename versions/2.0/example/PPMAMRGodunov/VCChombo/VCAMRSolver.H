#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#ifndef _VCAMRSolver_H_
#define _VCAMRSolver_H_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <assert.h>

#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "AMRLevelMG.H"
#include "LevelSolver.H"
#include "AMRSolver.H"

/// Class which manages grid hierarchy and composite elliptic solution
/// Yuri Omelchenko, 08.20.02
/**
   VCAMRSolver extends AMRSolver by implementing a Krylov (BiCGSTab) 
   solution to the elliptic equation with space dependent coefficients
   on a multiple level grid that satisfies certain proper nesting 
   conditions. 
   Warning: do not cast this class to AMRSolver 
   (no virtual functions defined yet).
 */
class VCAMRSolver : public AMRSolver
{
public:
  /// 
  /**
     Creates a VCAMRSolver whose internal state is  undefined.  
     Need to call define(...) function to use any of 
     the functionality of the class
  */
  VCAMRSolver();

  /// 
  void
  define(const Vector<DisjointBoxLayout>& a_gridsLevel,
         const Vector<Box>& a_domainLevel,
         const Vector<Real>& a_dxLevel,
         const Vector<int>& a_refRatio,
         int a_numLevels,
         int a_lBase,
         const Vector<LevelOp*>& a_opinv,
         int a_ncomp);	
  /// 
  /**
     Defines AMRSolvers' internal state
     \\ {\bf Inputs:} 
     \\ gridsLevel ---
     The grids at all levels.        
     Each element in the vector is a level's worth of grids.
     gridsLevel[0] are grids for level 0 and so forth.
     Vector index corresponds to level number.
     \\ domainLevel ---
     The domains at all levels.        
     Each element in the vector is a level's domain.
     domainLevel[0] is the domain for level 0 and so forth.
     Vector index corresponds to level number.
     \\ dxLevel ---
     The grid spacing at all levels.        
     Each element in the vector is a level's $\dx$.
     dxLevel[0]--- is $\dx$ for level 0 and so forth.
     Vector index corresponds to level number.
     \\ ref_ratio--- 
     The refinement ratio between all levels.        
     ref_ratio[0] is the refinement ratio between
     level 0 and level 1;
     Vector index corresponds to level number.
     \\ numlevels  
     The number of AMR levels in the calculation.  The length of the
     Vector~s has to be at least numlevels.

     \\ 
     lBase - coarsest level on which solution is to be computed. This
     needs to be set at the time of definition, in order to build the bottom
     LevelSolver.
     \\ opin_a
     The levelop  to use in the solution.
     \\ a_ncomp
     number of components
  */
  void
  define(const Vector<DisjointBoxLayout>& a_gridsLevel,
         const Vector<ProblemDomain>& a_domainLevel,
         const Vector<Real>& a_dxLevel,
         const Vector<int>& a_refRatio,
         int a_numLevels,
         int a_lBase,
         const Vector<LevelOp*>& a_opinv,
         int a_ncomp); 

  /// 
  ~VCAMRSolver();	// not virtual

  /// Solves on hierarchy to tolerance m_tolerance
  /**
     Solves the elliptic equation over the
     hierarchy of levels lBase ... lMax where $lMax = numlevels-1$.
     If  lBase > 0, then the data at level lBase - 1 is used 
     to interpolate boundary conditions at boundary cells that are not
     adjacent to the domain boundary.
     \\{\bf Inputs:} 
     \\
     phiLevel - pointers to current guess at the
     solution values for levels (lMin = max(lBase-1,0)) ... lMax.
     Vector index corresponds to level number.
     \\
     rhsLevel - pointers to right-hand side for levels lmin ... lMax.
     Vector index corresponds to level number.
     \\{\bf Outputs:}
     \\
     phiLevel -  LevelData<FArrayBox>s pointed to for levels 
     lMin,...,lMax are updated in place.
  */
  void
  solveBiCGStab(Vector<LevelData<FArrayBox> *>& a_phiLevel, 
                const Vector<LevelData<FArrayBox> *>& a_rhsLevel);

  ///
  void setNumAMRVCycle(int a_numAMRVCycle)
       {m_numAMRVCycle = a_numAMRVCycle;}

protected:

  ///
  void
  applyPreconditioner(Vector<LevelData<FArrayBox> *>& a_phiLevel, 
                const Vector<LevelData<FArrayBox> *>& a_rhsLevel,
                      Vector<LevelData<FArrayBox> *>& a_corrLevel,
                      Vector<LevelData<FArrayBox> *>& a_residLevel);
  ///
  Vector<Real>
  computeNorm(const Vector<LevelData<FArrayBox> *>& a_arrv, int a_normType);

  /// complete, 3-level operator with homogeneous BCs
  /**
     we need it because Chombo does not have it yet!
  */
  void 
  applyAMROperatorHphys(Vector<LevelData<FArrayBox>* >& a_phiLevel, 
                        LevelData<FArrayBox> & a_LofPhi,
                        int a_ilev);

  ///
  void setDefaultValues(); // not virtual

  ///
  void clear();		   // not virtual

  ///
  int m_numAMRVCycle;

private:
  
  VCAMRSolver(const VCAMRSolver&);
  VCAMRSolver& operator=(const VCAMRSolver&);


};
#endif

