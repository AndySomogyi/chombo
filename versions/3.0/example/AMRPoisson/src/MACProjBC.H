#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _MACProjFUNC_H_
#define _MACProjFUNC_H_

#include "IntVect.H"
#include "RealVect.H"
#include "LoHiSide.H"

///
/**
   Velocity function interface.
   point is the location of the function evaluation
*/
///
typedef Real (*VelocityBCFunction)(const IntVect&        a_point,
                                   const Real&           a_dx,
                                   const FluxBox  &      a_vel,
                                   const int&            a_faceDir, //= velcomp in mac land
                                   const Side::LoHiSide& a_domainSide);

////
typedef Real (*GradientBCFunction)(const IntVect&        a_point,
                                   const Real&           a_dx,
                                   const FArrayBox&      a_phi,
                                   const int&            a_faceDir, //= velcomp in mac land
                                   const Side::LoHiSide& a_domainSide);


///handy example
extern Real noFlowVelocityBC(const IntVect&        a_point,     
                             const Real&           a_dx,
                             const FluxBox  &      a_vel,
                             const int&            a_faceDir,   
                             const Side::LoHiSide& a_domainSide);

///handy example
extern Real noFlowGradientBC(const IntVect&        a_point,     
                             const Real&           a_dx,
                             const FArrayBox&      a_phi,
                             const int&            a_faceDir,   
                             const Side::LoHiSide& a_domainSide);


///handy example
extern Real inflowOutflowVelocityBC(const IntVect&        a_point,     
                                    const Real&           a_dx,
                                    const FluxBox  &      a_vel,
                                    const int&            a_faceDir,   
                                    const Side::LoHiSide& a_domainSide);

///handy example
extern Real inflowOutflowGradientBC(const IntVect&        a_point,     
                                    const Real&           a_dx,
                                    const FArrayBox&      a_phi,
                                    const int&            a_faceDir,   
                                    const Side::LoHiSide& a_domainSide);



#endif 
