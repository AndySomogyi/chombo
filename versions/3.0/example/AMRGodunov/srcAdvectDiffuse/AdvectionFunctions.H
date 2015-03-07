#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _AdvectionFUNC_H_
#define _AdvectionFUNC_H_

#include "IntVect.H"
#include "RealVect.H"

#include "UsingNamespace.H"

///
/**
   Velocity function interface.
   point is the location of the function evaluation
   velcomp is the velocity component.
*/
///
typedef Real (*AdvectionVelocityFunction)(const RealVect& a_point,
                                          const int&      a_velComp);


/// handy example--constant flow
Real constantAdvection(const RealVect& a_point,
                       const int&      a_velComp);


/// handy example--rotating flow
Real rotatingAdvection(const RealVect& a_point,
                       const int&      a_velComp);


#endif
