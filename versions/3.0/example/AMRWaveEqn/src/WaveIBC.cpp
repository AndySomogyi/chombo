#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <cmath>

#include "LoHiSide.H"
#include "BoxIterator.H"
#include "WaveIBC.H"

// Null constructor
WaveIBC::WaveIBC()
{
    m_params_are_set = false;
}
// Sets parameters used by initial conditions
//     a_r0            - Full width of Gaussian.
void WaveIBC::setParams(Real a_r0, RealVect a_x0)
{
    CH_assert(m_params_are_set == false);

    m_rr0 = a_r0;
    m_x0 = a_x0;

    m_params_are_set = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysIBC* WaveIBC::new_physIBC()
{
  WaveIBC* retval = new WaveIBC();
  retval -> setParams(m_rr0,m_x0);
  //[NOTE: do this even though setParams() will set it true
  //       because it might not have been set in the source
  //       object (this).  Of course, that would be a bad idea
  //       because then any object created by the factor wont
  //       be usable.  Caveat usor.  -dbs]
  retval->m_params_are_set = m_params_are_set ;

  return static_cast<PhysIBC*>(retval);
}

// Set up initial conditions
void WaveIBC::initialize(LevelData<FArrayBox>& a_phi,
                         LevelData<FArrayBox>& a_pi, Real dx)
{
  CH_assert(m_params_are_set == true);

  DataIterator dit = a_phi.boxLayout().dataIterator();
  int ncomp = a_phi.nComp();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& phi = a_phi[dit()];
    a_pi[dit()].setVal(0.);

    // Box of current grid
    Box phiBox = phi.box();
//    phiBox &= m_domain;   //?? Shouldn't this be an error check -dbs??
                        // No - ghost cells are allocated independent of
                        // whether they are in the domain. - PC.

    // compute initial values

    BoxIterator bit(phiBox);
    for (bit.begin();bit.ok();++bit)
    {
      IntVect iv = bit();
      Real dist = 0.;
      for (int idir = 0;idir < SpaceDim; idir++)
      {
        dist = dist + (m_x0[idir] + (0.5+iv[idir])*dx) * (m_x0[idir] + (0.5+iv[idir])*dx);
      }
      for (int icomp = 0;icomp < ncomp;icomp++)
      {
      phi(iv,icomp) = exp(-dist/m_rr0/m_rr0)/pow(m_rr0,SpaceDim);
      }
    }
  }
}
