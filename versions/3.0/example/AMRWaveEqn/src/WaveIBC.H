#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef  WAVEIBC_H_
#define  WAVEIBC_H_

#include  <iostream>

#include "RealVect.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "AMRIO.H"
#include "PhysIBC.H"
#include "UsingNamespace.H"


/// Initial conditions for gaussian.
/**
   Parameters used in the initial conditions
     m_rr0            - initial radius
 */
class WaveIBC:public PhysIBC
{
public:
  /// Null constructor
  /**
   */
  WaveIBC();

  /// Destructor
  /**
   */
  ~WaveIBC() {};

  /// Sets parameters
  /**
   */
  void setParams(Real a_r0, RealVect a_x0=RealVect::Zero);

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used) and m_params_are_set
     set to value of m_params_are_set in the current (factory) object.
   */
  PhysIBC *new_physIBC();

  /// Set up initial conditions
  /**
   */
  void initialize(LevelData<FArrayBox>& a_phi, LevelData<FArrayBox>& a_pi,
                  Real a_dx);

protected:
  RealVect m_x0;
  Real m_rr0;
  bool m_params_are_set;
};

#endif
