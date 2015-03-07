#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "PhysIBC.H"
#include "NamespaceHeader.H"

// Indicate that define() hasn't been called
PhysIBC::PhysIBC()
{
  m_isDefined = false;
}

PhysIBC::~PhysIBC()
{
}

// Define the object
void PhysIBC::define(const ProblemDomain& a_domain,
                     const Real&          a_dx)
{
  m_domain = a_domain;
  m_dx     = a_dx;

  m_isDefined = true;
}
#include "NamespaceFooter.H"
