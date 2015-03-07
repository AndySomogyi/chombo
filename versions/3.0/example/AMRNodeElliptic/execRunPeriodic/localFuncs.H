#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef LOCALFUNCS_H
#define LOCALFUNCS_H
#include  <iostream>
#include  <cstdio>
#include "LevelData.H"
#include "NodeFArrayBox.H"
#include "NodePoissonBC.H"
#include "Vector.H"

/*
  set domain boundary conditions from input file
 */
int setDomainBC(DomainNodeBC& a_dombc,
                const bool a_verbose);

/*
  set homogeneous Dirichlet boundary conditions
 */
int setDirichletBC(DomainNodeBC& a_dombc,
                   const bool a_verbose);

/*
  Set RHS on hierarchy from input file
 */
int setRHS(Vector<LevelData<NodeFArrayBox>* >& a_vectRhs,
           const Vector<Real>&                 a_vectDx,
           const Vector<ProblemDomain>&        a_vectDomain, // CELL-centered
           const int a_numlevels,
           const bool a_verbose);

/*
  Set phi on hierarchy from input file
*/
int setPhiExact(Vector<LevelData<NodeFArrayBox>* >& a_phi,
                const Vector<Real>&                 a_vectDx,
                const Vector<ProblemDomain>&        a_vectDomain,
                const int a_numlevels,
                const bool a_verbose);

#endif
