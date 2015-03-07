#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "DirichletConductivityDomainBC.H"
#include "EBArithF_F.H"
#include "NamespaceHeader.H"
/*****/
DirichletConductivityDomainBC::
DirichletConductivityDomainBC()
{
}
/*****/
DirichletConductivityDomainBC::
~DirichletConductivityDomainBC()
{
}

/*****/
void
DirichletConductivityDomainBC::
setValue(Real a_value)
{
  m_bc.setValue(a_value);
}

/*****/
int
DirichletConductivityDomainBC::
whichBC(int                  a_idir,
        Side::LoHiSide       a_side)
{
  return 0;
}

/*****/
void
DirichletConductivityDomainBC::
setFunction(RefCountedPtr<BaseBCValue> a_flux)
{
  m_bc.setFunction(a_flux);
}

/*****/
void
DirichletConductivityDomainBC::
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
  m_bc.getFaceFlux(a_faceFlux,
                   a_phi,
                   a_probLo,
                   a_dx,
                   a_idir,
                   a_side,
                   a_dit,
                   a_time,
                   a_useHomogeneous);
  //again, following the odd convention of EBAMRPoissonOp
  //(because I am reusing its BC classes),
  //the input flux here is CELL centered and the input box
  //is the box adjacent to the domain boundary on the valid side.
  //because I am not insane (yet) I will just shift the flux's box
  //over and multiply by the appropriate coefficient
  a_faceFlux.shiftHalf(a_idir, -sign(a_side));
  const Box& faceBox = a_faceFlux.box();
  const BaseFab<Real>&   regCoef = (*m_bcoef)[a_dit][a_idir].getSingleValuedFAB();
  int  isrc = 0;
  int  idst = 0;
  int  inum = 1;
  FORT_MULTIPLYTWOFAB(CHF_FRA(a_faceFlux),
                      CHF_CONST_FRA(regCoef),
                      CHF_BOX(faceBox),
                      CHF_INT(isrc),CHF_INT(idst),CHF_INT(inum));

  //shift flux back to cell centered land
  a_faceFlux.shiftHalf(a_idir,  sign(a_side));
}

/*****/
void

DirichletConductivityDomainBC::
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
  m_bc.getFaceFlux(a_faceFlux,
                   a_vof,
                   a_comp,
                   a_phi,
                   a_probLo,
                   a_dx,
                   a_idir,
                   a_side,
                   a_dit,
                   a_time,
                   a_useHomogeneous);
  const EBISBox& ebisBox = a_phi.getEBISBox();
  Real bcoave = 0;
  Vector<FaceIndex> faces = ebisBox.getFaces(a_vof, a_idir, a_side);
  Real areaTot = 0;
  for(int iface = 0; iface < faces.size(); iface++)
    {
      Real areaFrac = ebisBox.areaFrac(faces[iface]);
      Real bcoFace  = (*m_bcoef)[a_dit][a_idir](faces[iface], 0);
      bcoave += areaFrac*bcoFace;
      areaTot += areaFrac;
    }
  if(areaTot > 1.0e-8)
    {
      bcoave /= areaTot;
    }
  a_faceFlux *= bcoave;
}

/*****/
void
DirichletConductivityDomainBC::
getFaceGradPhi(Real&                 a_faceFlux,
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
  m_bc.getFaceGradPhi(a_faceFlux,
                      a_face,
                      a_comp,
                      a_phi,
                      a_probLo,
                      a_dx,
                      a_idir,
                      a_side,
                      a_dit,
                      a_time,
                      a_useAreaFrac,
                      a_centroid,
                      a_useHomogeneous);
}

/*****/
void
DirichletConductivityDomainBC::
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
  m_bc.getFaceVel(a_faceFlux, a_face, a_vel, a_probLo, a_dx, a_idir,
                  a_icomp, a_time, a_side, a_doDivFreeOutflow);
}
/******/
DirichletConductivityDomainBCFactory::
DirichletConductivityDomainBCFactory()
{
  m_value = 12345.6789;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = true;
  m_isFunction = false;
}

/******/
DirichletConductivityDomainBCFactory::
~DirichletConductivityDomainBCFactory()
{
}
/******/
void
DirichletConductivityDomainBCFactory::
setValue(Real a_value)
{
  m_value = a_value;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = false;
  m_isFunction = false;
}
/******/
void
DirichletConductivityDomainBCFactory::
setFunction(RefCountedPtr<BaseBCValue> a_flux)
{
  m_value = 12345.6789;
  m_flux = a_flux;

  m_onlyHomogeneous = false;
  m_isFunction = true;
}
/******/
DirichletConductivityDomainBC*
DirichletConductivityDomainBCFactory::
create(const ProblemDomain& a_domain,
       const EBISLayout&    a_layout,
       const RealVect&      a_dx)
{
  DirichletConductivityDomainBC* newBC = new DirichletConductivityDomainBC();
  if(!m_onlyHomogeneous)
    {
      if(m_isFunction)
        {
          newBC->setFunction(m_flux);
        }
      else
        {
          newBC->setValue(m_value);
        }
    }
  return newBC;
}
#include "NamespaceFooter.H"
