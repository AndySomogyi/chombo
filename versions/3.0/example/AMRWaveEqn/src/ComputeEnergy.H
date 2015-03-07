#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _COMPUTEENERGY_H_
#define _COMPUTEENERGY_H_

#include "REAL.H"
#include "Vector.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "UsingNamespace.H"

/// Returns the energy integrated over the whole domain.
/** Uses the finest box available at each cell.  See computeEnergy.ChF.
 */
Real computeEnergy(const Vector<LevelData<FArrayBox>* >& a_phi,
                   const Vector<LevelData<FArrayBox>* >& a_pi,
                   const Vector<int>&                    a_nRefFine,
                   Real                                  a_dxCrse,
                   Real                                  a_lambda);

/// Returns the energy integrated over a single box
Real computeEnergy(const LevelData<FArrayBox>& a_phi,
                   const LevelData<FArrayBox>& a_pi,
                   const DisjointBoxLayout*    a_finerGrids,
                   int                         a_nRefFine,
                   Real                        a_dx,
                   Real                        a_lambda);
#endif
