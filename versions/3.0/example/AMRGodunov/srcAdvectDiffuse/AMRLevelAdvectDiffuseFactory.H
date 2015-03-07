#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _AMRLEVELADVECTDIFFUSEFACTORY_H_
#define _AMRLEVELADVECTDIFFUSEFACTORY_H_

#include "AMRLevelFactory.H"

#include "AdvectionFunctions.H"
#include "AMRLevelAdvectDiffuse.H"

#include "UsingNamespace.H"

///
/**
 */
class AMRLevelAdvectDiffuseFactory: public AMRLevelFactory
{
public:

  ///
  /**
   */
  AMRLevelAdvectDiffuseFactory(const AdvectPhysics&        a_gphys,
                               AdvectionVelocityFunction   a_advFunc,
                               BCHolder                    a_bcFunc, //bc for diffusion solve
                               const Real&                 a_cfl,
                               const Real&                 a_domainLength,
                               const Real&                 a_refineThresh,
                               const int&                  a_tagBufferSize,
                               const Real&                 a_initialDtMultiplier,
                               const bool&                 a_useLimiting,
                               const Real&                 a_nu);


  ///
  virtual ~AMRLevelAdvectDiffuseFactory()
  {
  }


  ///
  virtual AMRLevel* new_amrlevel() const;

protected:

  Real                  m_cfl;
  Real                  m_domainLength;
  Real                  m_refineThresh;
  int                   m_tagBufferSize;
  Real                  m_initialDtMultiplier;
  bool                  m_useLimiting;
  Real                  m_nu;
  RefCountedPtr<AdvectPhysics>     m_advPhys;
  AdvectionVelocityFunction        m_advFunc;
  BCHolder                    m_bcFunc; //bc for diffusion solve

private:
  //weak construction is evil
  AMRLevelAdvectDiffuseFactory()
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  void operator=(const AMRLevelAdvectDiffuseFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelAdvectDiffuseFactory(const AMRLevelAdvectDiffuseFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
