#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// generalFuncs.H
// petermc, 24 March 2003

#ifndef GENERALFUNCS_H
#define GENERALFUNCS_H
#include  <iostream>
#include  <cstdio>
#include "LevelData.H"
#include "NodeFArrayBox.H"
#include "Vector.H"
#include "RealVect.H"

#include "UsingNamespace.H"

Box boxFromVector(const Vector<int>&  a_ints,
                  int                 a_start = 0);

int readGrids(Vector<DisjointBoxLayout>& a_vectGrids,
              Vector<ProblemDomain>&     a_vectDomain,
              Vector<Real>&              a_vectDx,
              Vector<int>&               a_vectRefRatio,
              RealVect&                  a_probLo,
              RealVect&                  a_probHi,
              int& a_numlevels,
              const bool a_verbose,
              int a_refined = 1);

int readGrids(Vector<DisjointBoxLayout>& a_vectGrids,
              Vector<Box>&               a_vectDomain,
              Vector<Real>&              a_vectDx,
              Vector<int>&               a_vectRefRatio,
              RealVect&                  a_probLo,
              RealVect&                  a_probHi,
              int& a_numlevels,
              const bool a_verbose,
              int a_refined = 1);

int
readDomain(Vector<ProblemDomain>&     a_vectDomain,
           Vector<Real>&              a_vectDx,
           Vector<int>&               a_vectRefRatio,
           RealVect&                  a_probLo,
           RealVect&                  a_probHi,
           int& a_numlevels,
           const bool a_verbose,
           int a_refined = 1);

int setGridsTwoLevel(const IntVect&             a_lengths,
                     const IntVect&             a_offsets,
                     const IntVect&             a_sublengths,
                     const int                  a_refRatio,
                     const Real                 a_dx,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<ProblemDomain>&     a_vectDomain,
                     Vector<Real>&              a_vectDx,
                     Vector<int>&               a_vectRefRatio,
                     const bool a_verbose);

int getGridsLevel(Vector<Box>&          a_grids,
                  const ProblemDomain&  a_domain,
                  const Vector<int>&    a_subdomainIndices,
                  const int             a_boxCount,
                  const Vector<int>&    a_boxIndices);

int assignAllGrids(Vector<DisjointBoxLayout>&    a_vectGrids,
                   const Vector<ProblemDomain>&  a_vectDomain,
                   const Vector<Vector<Box> >&   a_boxes,
                   int a_maxboxsize = 0,
                   int a_blockFactor = 1);

int assignLevelGrids(Vector<DisjointBoxLayout>&  a_vectGrids,
                     const ProblemDomain&        a_domain,
                     const Vector<Box>&          a_boxes,
                     int a_level,
                     int a_maxboxsize = 0,
                     int a_blockFactor = 1);

int assignLevelGrids(DisjointBoxLayout&          a_vectGrids,
                     const ProblemDomain&        a_domain,
                     const Vector<Box>&          a_boxes,
                     int a_maxboxsize = 0,
                     int a_blockFactor = 1);

int assignDx(Vector<Real>&              a_vectDx,
             const Vector<int>&         a_vectRefRatio,
             const Real                 a_dx0,
             const int a_numlevels);

int assignDomains(Vector<ProblemDomain>&     a_vectDomain,
                  const Vector<int>&         a_vectRefRatio,
                  const ProblemDomain&       a_baseDomain,
                  const int a_numlevels);

int getDiff(Vector<LevelData<NodeFArrayBox>* >& a_vectDiff,
            const Vector<LevelData<NodeFArrayBox>* >& a_vectPhi1,
            const Vector<LevelData<NodeFArrayBox>* >& a_vectPhi2);

int getMagnitude(Vector<LevelData<NodeFArrayBox>* >& a_vectMag,
                 const Vector<LevelData<NodeFArrayBox>* >& a_vectField);

int getTruncError(Vector<LevelData<NodeFArrayBox>* >& a_vectErr,
                  const Vector<LevelData<NodeFArrayBox>* >& a_vectLap,
                  const Vector<LevelData<NodeFArrayBox>* >& a_vectRhs,
                  const Vector<DisjointBoxLayout>& a_vectGrids,
                  const Vector<ProblemDomain>& a_vectDomain,
                  const Vector<int>& a_vectRatio,
                  const Vector<Real>& a_vectDx,
                  const int a_numlevels,
                  const bool a_verbose);

int getTruncError(Vector<LevelData<NodeFArrayBox>* >& a_vectErr,
                  const Vector<LevelData<NodeFArrayBox>* >& a_vectLap,
                  const Vector<LevelData<NodeFArrayBox>* >& a_vectRhs,
                  const Vector<DisjointBoxLayout>& a_vectGrids,
                  const Vector<Box>& a_vectDomain,
                  const Vector<int>& a_vectRatio,
                  const Vector<Real>& a_vectDx,
                  const int a_numlevels,
                  const bool a_verbose);

int project2(Vector<LevelData<NodeFArrayBox>* >& a_vectProj,
             const Vector<LevelData<NodeFArrayBox>* >& a_vectPhi,
             const Vector<DisjointBoxLayout>& a_vectGrids);


int printErrorNorms(const string& a_prefix,
                    Real a_dxFine,
                    int a_deginterp,
                    Real a_normMaxerrFine,
                    Real a_norm1errFine,
                    Real a_norm2errFine,
                    Real a_normMaxerrCoarse,
                    Real a_norm1errCoarse,
                    Real a_norm2errCoarse);

int printDiffNorms(const string& a_prefix,
                   Real a_dxFine,
                   int a_deginterp,
                   Real a_normMaxdiff,
                   Real a_norm1diff,
                   Real a_norm2diff);

#endif
