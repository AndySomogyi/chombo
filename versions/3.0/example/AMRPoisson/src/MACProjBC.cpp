#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "IntVect.H"
#include "LoHiSide.H"
#include "ParmParse.H"
#include "REAL.H"
#include "FluxBox.H"
#include "FArrayBox.H"

///
Real 
noFlowVelocityBC(const IntVect&        a_point,     
                 const Real&           a_dx,
                 const FluxBox  &      a_vel,
                 const int&            a_faceDir,   
                 const Side::LoHiSide& a_domainSide)
{
  return 0.;
}

///
Real 
noFlowGradientBC(const IntVect&        a_point,     
                 const Real&           a_dx,
                 const FArrayBox&      a_phi,
                 const int&            a_faceDir,   
                 const Side::LoHiSide& a_domainSide)

{
  return 0.;
}

///
Real 
inflowOutflowVelocityBC(const IntVect&        a_point,     
                        const Real&           a_dx,
                        const FluxBox  &      a_vel,
                        const int&            a_faceDir,   
                        const Side::LoHiSide& a_domainSide)
{
  ParmParse pp;
  int  inflowDir;
  Real inflowVel;

  pp.get("inflow_vel", inflowVel);
  pp.get("inflow_dir", inflowDir);
  Real retval;
  if(      (a_faceDir == inflowDir) && (a_domainSide == Side::Lo)) //inflow side
    {
      retval = inflowVel;
    }
  else if ((a_faceDir == inflowDir) && (a_domainSide == Side::Hi)) //outflow side
    {
      //extrapolate for outflow bc      
      //a_point is the iv for the outflow face (because it is on the high side)
      //a_point -  BASISV(a_faceDir) is the face one back 
      //a_point -2*BASISV(a_faceDir) is the face two back 
      //the above -'s are - because this is the high side of the domain
      Real velClose = a_vel[a_faceDir](a_point -  BASISV(a_faceDir), 0);
      Real velFar   = a_vel[a_faceDir](a_point -2*BASISV(a_faceDir), 0);
      Real extrapVal = 2*velClose - velFar;
      retval = extrapVal;
    }
  else
    {
      retval = 0;
    }
  return retval;
}

///
Real 
inflowOutflowGradientBC(const IntVect&        a_point,     
                        const Real&           a_dx,
                        const FArrayBox&      a_phi,
                        const int&            a_faceDir,   
                        const Side::LoHiSide& a_domainSide)
{
  ParmParse pp;
  int  inflowDir;

  pp.get("inflow_dir", inflowDir);
  bool isOutflow= (a_domainSide==Side::Hi) && (a_faceDir==inflowDir);

  Real retval;
  if(isOutflow)
    {
      Real phiVal = a_phi(a_point, 0);
      //this is -(phi-phi0)/(0.5*dx) estimate of gradient at outflow
      //negative because it is at the high side of the domain 
      //so phi0 is higher in x than phi
      //oh, and phi0=0 in case that was not obvious
      retval = -2.*phiVal/a_dx;
    }
  else
    {
      retval = 0;
    }
  return retval;
}



