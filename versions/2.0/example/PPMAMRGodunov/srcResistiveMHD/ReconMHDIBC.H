#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef  RECONMHDIBC_H_
#define  RECONMHDIBC_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "PhysIBC.H"

class ReconMHDIBC:public PhysIBC
{
public:
  /// Null constructor
  /**
   */
  ReconMHDIBC();
  
  /// Constructor which defines parameters used by Fortran routines
  /**
   */
  ReconMHDIBC(const Real& a_gamma,
	      const Real& a_mu,
	      const Real& a_eta,
	      const Real& a_kappa);

  /// Destructor
  /**
   */
  ~ReconMHDIBC() {};

  /// Sets parameters in a common block used by Fortran routines
  /**
   */
  void setParameters(const Real& a_gamma,
	      const Real& a_mu,
	      const Real& a_eta,
	      const Real& a_kappa);

  /**
     Set the flag m_isFortranCommonSet to true so that new IBCs made with
     new_physIBC() will have this flag set without calling setFortranCommon()
     (this is a clumsy design and should be improved).
   */
  void setParameterSet();

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used) and m_isFortranCommonSet
     set to value of m_isFortranCommonset in the current (factory) object.
   */
  PhysIBC *new_physIBC();

  /// Set boundary fluxes
  /**
   */
  void fluxBC(FArrayBox&            a_F,
              const FArrayBox&      a_W,
              const FArrayBox&      a_Wextrap,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const Real&           a_time);
  
  virtual void diffusiveEnergyFluxBC(FArrayBox&            a_F,
				     FArrayBox&      a_W,
				     const int&            a_dir,
				     const Side::LoHiSide& a_side,
				     const Real&           a_time);

  virtual void faceMagBC(FArrayBox&            a_F,
			 const int&            a_dir,
			 const Side::LoHiSide& a_side,
			 const Real&           a_time);

  virtual void parabolicFluxBC(FArrayBox&            a_F,
			       FArrayBox&      a_W,
			       const int&            a_dir,
			       const Side::LoHiSide& a_side,
			       const Real&           a_time);

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one sided difference
     approximations.  If this function doesn't change them they will be
     used for the slopes at the boundaries.
   */
  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
                     const Real&      a_time);
  
  /// Set up initial conditions
  /**
   */
  void initialize(LevelData<FArrayBox>& a_U);

  void initializePhi(LevelData<FArrayBox>& a_U);


  void phiBC(FArrayBox& a_U,
	     const ProblemDomain& a_domain,
	     const Real& a_dxLevel);
// Set boundary fluxes
  void primBC(FArrayBox&            a_WGdnv,
	      const FArrayBox&      a_Wextrap,
	      const FArrayBox&      a_W,
	      const int&            a_dir,
	      const Side::LoHiSide& a_side,
	      const Real&           a_time);

// Adjust boundary fluxes to account for artificial viscosity
// Do nothing for this problem
virtual void artViscBC(FArrayBox&       a_F,
		       const FArrayBox& a_U,
		       const FArrayBox& a_divVel,
		       const int&       a_dir,
		       const Real&      a_time);
  
protected:
  // True if the Fortran common block has been set
  bool m_isParameterSet;
};

#endif
