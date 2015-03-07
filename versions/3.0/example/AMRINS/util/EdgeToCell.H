#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EDGETOCELL_H_
#define _EDGETOCELL_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"

#include "UsingNamespace.H"

/// averages cell-centered cellData to edge-centered edgeData
void EdgeToCell(const LevelData<FluxBox>& a_edgeData,
                LevelData<FArrayBox>& a_cellData);

/// averages cellData->edgeData (all directions) for one FAB
/** averages cell-centered cellData to SpaceDim-dimensioned
    EdgeData -- assumes all boundary conditions already set.
    Does this for all grids. */
void EdgeToCell(const FluxBox& a_edgeData,
                FArrayBox& a_cellData);

/// averages cellData->edgeData in 1 direction, for one component in a FAB
void EdgeToCell(const FluxBox& a_edgeData, const int a_edgeComp,
                FArrayBox& a_cellData, const int a_cellComp,
                const int a_dir);

/// averages cellData->edgeData in 1 direction, for one component in a FAB, on a given cell box
void EdgeToCell(const FluxBox& a_edgeData, const int a_edgeComp,
                FArrayBox& a_cellData, const int a_cellComp,
                const Box& a_cellBox, const int a_dir);

#endif
