#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "BoxIterator.H"

#include "EBPlanarShockSolverBC.H"
#include "NeumannViscousTensorDomainBC.H"
#include "DirichletViscousTensorDomainBC.H"
#include "EBArith.H"
#include "Stencils.H"
#include "DirichletViscousTensorEBBC.H"
#include "VoFIterator.H"

#include "NamespaceHeader.H"

bool
EBPlanarShockSolverBC::
getNeumann(const int&            a_idir, 
           const Side::LoHiSide& a_side)
{
  //velocity is neumann only at inflow or at slipwalls
  bool doingNeumann;
  if(m_slipWall)
    {
      doingNeumann = true;
    }
  else if(a_idir != m_shockNorm)
    {
      doingNeumann = false;
    }
  else
    {
      CH_assert(a_idir == m_shockNorm);
      //inflow = neumann.  outflow = dirichlet.   a_dir = shockNorm established
      bool atInflow =  ((!m_shockBackward && a_side == Side::Lo) ||
                        ( m_shockBackward && a_side == Side::Hi));

      doingNeumann = atInflow;
    }
  return doingNeumann;
}

////
void 
EBPlanarShockSolverBC::
getFaceFlux(BaseFab<Real>&        a_faceFlux,
            const BaseFab<Real>&  a_phi,
            const RealVect&       a_probLo,
            const RealVect&       a_dx,
            const int&            a_idir,
            const Side::LoHiSide& a_side,
            const DataIndex&      a_dit,
            const Real&           a_time,
            const bool&           a_useHomogeneous)
{
  //velocity is neumann only at inflow or at slipwalls
  bool doingNeumann = getNeumann(a_idir, a_side);
  if(doingNeumann)
    {
      if(m_slipWall)
        {
          MayDay::Error("need to write no flow EBVT bc ");
        }
      NeumannViscousTensorDomainBC neumannBC;
      neumannBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else
    {
      DirichletViscousTensorDomainBC diriBC;
      diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      diriBC.setValue(0.0);
      diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}


void 
EBPlanarShockSolverBC::
getFaceFlux(Real&                 a_faceFlux,
            const VolIndex&       a_vof,
            const int&            a_comp,
            const EBCellFAB&      a_phi,
            const RealVect&       a_probLo,
            const RealVect&       a_dx,
            const int&            a_idir,
            const Side::LoHiSide& a_side,
            const DataIndex&      a_dit,
            const Real&           a_time,
            const bool&           a_useHomogeneous)
{
  bool doingNeumann = getNeumann(a_idir, a_side);
  if(doingNeumann)
    {
      MayDay::Error("need to write no flow EBVT bc ");
      NeumannViscousTensorDomainBC neumannBC;
      neumannBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else
    {
      DirichletViscousTensorDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}


#include "NamespaceFooter.H"
