#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "NeumannConductivityDomainBC.H"
#include "EBArithF_F.H"
#include "NamespaceHeader.H"
/*****/
NeumannConductivityDomainBC::
NeumannConductivityDomainBC()
{
}
/*****/
NeumannConductivityDomainBC::
~NeumannConductivityDomainBC()
{
}

/*****/
void
NeumannConductivityDomainBC::
setValue(Real a_value)
{
  m_bc.setValue(a_value);
}

/*****/
int
NeumannConductivityDomainBC::
whichBC(int                  a_idir,
        Side::LoHiSide       a_side)
{
  return 0;
}

/*****/
void
NeumannConductivityDomainBC::
setFunction(RefCountedPtr<BaseBCValue> a_flux)
{
  m_bc.setFunction(a_flux);
}

/*****/
void
NeumannConductivityDomainBC::
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
NeumannConductivityDomainBC::
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
  //already has area multiplied in
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
  Real areaTot = 0.0;
  for(int iface = 0; iface < faces.size(); iface++)
    {
      Real areaFrac = ebisBox.areaFrac(faces[iface]);
      areaTot += areaFrac;
      Real bcoFace  = (*m_bcoef)[a_dit][a_idir](faces[iface], 0);
      bcoave += areaFrac*bcoFace;
    }
  if(areaTot > 1.0e-8)
    {
      bcoave /= areaTot;
    }
  a_faceFlux *= bcoave;
}

/*****/
void
NeumannConductivityDomainBC::
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
NeumannConductivityDomainBC::
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
NeumannConductivityDomainBCFactory::
NeumannConductivityDomainBCFactory()
{
  m_value = 12345.6789;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = true;
  m_isFunction = false;
}

/******/
NeumannConductivityDomainBCFactory::
~NeumannConductivityDomainBCFactory()
{
}
/******/
void
NeumannConductivityDomainBCFactory::
setValue(Real a_value)
{
  m_value = a_value;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = false;
  m_isFunction = false;
}
/******/
void
NeumannConductivityDomainBCFactory::
setFunction(RefCountedPtr<BaseBCValue> a_flux)
{
  m_value = 12345.6789;
  m_flux = a_flux;

  m_onlyHomogeneous = false;
  m_isFunction = true;
}
/******/
NeumannConductivityDomainBC*
NeumannConductivityDomainBCFactory::
create(const ProblemDomain& a_domain,
       const EBISLayout&    a_layout,
       const RealVect&      a_dx)
{
  NeumannConductivityDomainBC* newBC = new NeumannConductivityDomainBC();
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
