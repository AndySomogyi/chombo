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

#include "InflowOutflowPoissonDomainBC.H"
#include "NeumannPoissonDomainBC.H"
#include "DirichletPoissonDomainBC.H"
#include "EBArith.H"
#include "Stencils.H"
#include "DirichletPoissonEBBC.H"
#include "VoFIterator.H"

#include "NamespaceHeader.H"

void
InflowOutflowPoissonDomainBC::
getFaceVel(Real&                 a_faceFlux,
           const FaceIndex&      a_face,
           const EBFluxFAB&      a_vel,
           const RealVect&       a_probLo,
           const RealVect&       a_dx,
           const int&            a_idir,
           const int&            a_icomp,
           const Real&           a_time,
           const Side::LoHiSide& a_side,
           const bool&           a_doDivFreeOutflow)
{

  int velcomp =  DirichletPoissonEBBC::s_velComp;
  CH_assert(a_vel.nComp() == 1);
  bool isInflow = (a_side==Side::Lo) && (a_idir==m_flowDir);
  bool isOutflow= (a_side==Side::Hi) && (a_idir==m_flowDir);

  RealVect point =  EBArith::getFaceLocation(a_face, a_dx, a_probLo);
  RealVect normal = EBArith::getDomainNormal(a_idir, a_side);

  const Side::LoHiSide inSide = flip(a_side);
  const EBISBox& ebisBox= a_vel[0].getEBISBox();
  const VolIndex& vof = a_face.getVoF(inSide);

  if(isOutflow)
    {//outflow
      if(a_doDivFreeOutflow)
        {//do div free outflow bc
          a_faceFlux = EBArith::getFaceVelForDivFreeCell(a_face,vof,a_idir,a_side,a_vel,a_dx,ebisBox);
        }
      else
        {//extrapolate for outflow bc
          a_faceFlux = EBArith::extrapToDomainFace(a_face,a_side,a_idir,ebisBox.getEBGraph(),a_vel[a_idir],0);
        }
    }
  else if(isInflow && (velcomp==m_flowDir))
    {
      //input component is always zero.
      //value of face direction corresponds to velocity direciton here
      a_faceFlux = m_inflowVel;
    }
  else
    {
      //must be a solid wall
      a_faceFlux = 0;
    }
}

////
void InflowOutflowPoissonDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
                                               const BaseFab<Real>&  a_phi,
                                               const RealVect&       a_probLo,
                                               const RealVect&       a_dx,
                                               const int&            a_idir,
                                               const Side::LoHiSide& a_side,
                                               const DataIndex&      a_dit,
                                               const Real&           a_time,
                                               const bool&           a_useHomogeneous)
{
  //for phi: outflow  dirichlet. inflow neumann
  bool isOutflow= (a_side==Side::Hi) && (a_idir==m_flowDir);
  if(!isOutflow)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}


void InflowOutflowPoissonDomainBC::getFaceFlux(Real&                 a_faceFlux,
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
  bool isOutflow= (a_side==Side::Hi) && (a_idir==m_flowDir);
  if(!isOutflow)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}


void InflowOutflowPoissonDomainBC::getFaceGradPhi(Real&                 a_faceFlux,
                                                  const FaceIndex&      a_face,
                                                  const int&            a_comp,
                                                  const EBCellFAB&      a_phi,
                                                  const RealVect&       a_probLo,
                                                  const RealVect&       a_dx,
                                                  const int&            a_idir,
                                                  const Side::LoHiSide& a_side,
                                                  const DataIndex&      a_dit,
                                                  const Real&           a_time,
                                                  const bool&           a_useAreaFrac,
                                                  const RealVect&       a_centroid,
                                                  const bool&           a_useHomogeneous)
{
  bool isOutflow= (a_side==Side::Hi) && (a_idir==m_flowDir);

  if(!isOutflow)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getFaceGradPhi(a_faceFlux, a_face, a_comp,a_phi, a_probLo,
                               a_dx, a_idir, a_side, a_dit,a_time,a_useAreaFrac,a_centroid,a_useHomogeneous);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.getFaceGradPhi(a_faceFlux, a_face, a_comp, a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit,a_time,a_useAreaFrac,a_centroid,a_useHomogeneous);
    }
}


void
InflowOutflowHelmholtzDomainBC::
getFaceVel(Real&                 a_faceFlux,
           const FaceIndex&      a_face,
           const EBFluxFAB&      a_vel,
           const RealVect&       a_probLo,
           const RealVect&       a_dx,
           const int&            a_idir,
           const int&            a_icomp,
           const Real&           a_time,
           const Side::LoHiSide& a_side,
           const bool&           a_doDivFreeOutflow)
{
  MayDay::Error("not needed");
}

void InflowOutflowHelmholtzDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
                                                 const BaseFab<Real>&  a_phi,
                                                 const RealVect&       a_probLo,
                                                 const RealVect&       a_dx,
                                                 const int&            a_idir,
                                                 const Side::LoHiSide& a_side,
                                                 const DataIndex&      a_dit,
                                                 const Real&           a_time,
                                                 const bool&           a_useHomogeneous)
{
  //vel: outflow is neumann. all others dirichlet.  inflow uses inflow vel as the value
  int velcomp =  DirichletPoissonEBBC::s_velComp;
  bool isOutflow =  ((a_side==Side::Hi) && (a_idir==m_flowDir));
  bool isInflow =  ((a_side==Side::Lo) && (a_idir==m_flowDir));
  bool isSlipWall = ((a_idir!=m_flowDir) && (a_idir != velcomp) && (m_doSlipWalls));
  bool isVelNeum =  (isOutflow || isSlipWall);

  if(isVelNeum)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.);
      neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else if(isInflow)
    {
      DirichletPoissonDomainBC diriBC;
      if(velcomp==m_flowDir)
        {
          diriBC.setValue(m_inflowVel);
        }
      else
        {
          diriBC.setValue(0.0);
        }
      diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else
    {
      //wall bc no slip
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}

void InflowOutflowHelmholtzDomainBC::getFaceFlux(Real&                 a_faceFlux,
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
  //vel: outflow is Neumann. all others Dirichlet
  int velcomp =  DirichletPoissonEBBC::s_velComp;
  bool isOutflow =  ((a_side==Side::Hi) && (a_idir==m_flowDir));
  bool isInflow =  ((a_side==Side::Lo) && (a_idir==m_flowDir));
  bool isSlipWall = ((a_idir!=m_flowDir) && (a_idir != velcomp) && (m_doSlipWalls));
  bool isVelNeum =  (isOutflow || isSlipWall);

  if(isVelNeum)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.);
      neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else if(isInflow)
    {
      DirichletPoissonDomainBC diriBC;
      if(velcomp==m_flowDir)
        {
          diriBC.setValue(m_inflowVel);
        }
      else
        {
          diriBC.setValue(0.0);
        }
      diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else
    {
      //wall bc no slip
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}

void InflowOutflowHelmholtzDomainBC::getFaceGradPhi(Real&                 a_faceFlux,
                                                    const FaceIndex&       a_face,
                                                    const int&            a_comp,
                                                    const EBCellFAB&      a_phi,
                                                    const RealVect&       a_probLo,
                                                    const RealVect&       a_dx,
                                                    const int&            a_idir,
                                                    const Side::LoHiSide& a_side,
                                                    const DataIndex&      a_dit,
                                                    const Real&           a_time,
                                                    const bool&           a_useAreaFrac,
                                                    const RealVect&       a_centroid,
                                                    const bool&           a_useHomogeneous)
{
  MayDay::Error("InOuAnyHelmholtz::getFaceGradPhi: not needed");
}

#include "NamespaceFooter.H"
