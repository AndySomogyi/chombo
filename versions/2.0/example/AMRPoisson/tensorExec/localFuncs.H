#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _localFuncs_
#define _localFuncs_
#include "LevelData.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "Vector.H"
#include "SPMD.H"
#include "ProblemDomain.H"
#include "GhostBC.H"
#include "AMRSolver.H"
#include "UsingNamespace.H"

/*
  set domain boundary conditions from input file
 */
int setDomainBC(DomainGhostBC& domghostbc, bool verbose);

/* 
   set domain boundary conditions for the tangential gradients
*/
int setTanGradBC(DomainGhostBC& domghostbc, bool verbose);

/*
  Set grid hierarchy, either from input file, or based on RHS
 */
int setGrids(Vector<DisjointBoxLayout>& vectGrids,
             Vector<ProblemDomain>&     vectDomain, 
             Vector<Real>&              vectDx, 
             Vector<int>&               vectRefRatio, 
             int& numlevels, bool verbose,
             Real domainMult = 1.0);

/*
  Set RHS on hierarchy from input file
 */
int setRHS(Vector<LevelData<FArrayBox>* >& vectRhs,
           Vector<Real>&                   vectDx, 
           Vector<ProblemDomain>&          vectDomain, 
           int numlevels,
           bool verbose);

/* set exact solution
 */
void setExact(LevelData<FArrayBox>& a_phi,
              const ProblemDomain& a_domain,
              const Real a_dx);

void
tagCells(Vector<LevelData<FArrayBox>*>& vectRHS,
	 Vector<IntVectSet>& tagVect,
	 Vector<Real>& vectDx,
	 Vector<ProblemDomain>& vectDomain,
	 const Real refine_thresh,
	 const int tags_grow,
	 const int baselevel, 
	 int numlevels);

/*
   Output phi and rhs to file -- need amrSolver to compute residuals.
   VectPhi is not const because of BC issues.
 */
int outputData(Vector<LevelData<FArrayBox>* >& vectPhi,
               const Vector<LevelData<FArrayBox>* >& vectRhs,
               const Vector<DisjointBoxLayout>& vectGrids,
               const Vector<ProblemDomain>& domain,
               const Vector<int>& vectRatio,
               const Vector<Real>& vectDx,
               AMRSolver& amr_solver,
               int numlevels,
               int lBase,
               Vector<Real>& L1Old,
               Vector<Real>& L2Old,
               Vector<Real>& MaxOld,
               bool verbose);

void 
computeError(LevelData<FArrayBox>& a_error, 
             const LevelData<FArrayBox>& a_phi,
             const DisjointBoxLayout* a_fineGridsPtr,
             const ProblemDomain& a_domain,
             const Real a_dx,
             const int a_nRefFine);

#endif

