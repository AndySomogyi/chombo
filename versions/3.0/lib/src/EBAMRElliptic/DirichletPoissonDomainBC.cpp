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
#include "EBArith.H"
#include "PolyGeom.H"

#include "DirichletPoissonDomainBC.H"
#include "DirichletPoissonDomainBCF_F.H"
#include "NamespaceHeader.H"


void
DirichletPoissonDomainBC::
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
  CH_assert(a_idir == a_face.direction());
  Real value;
  if (m_isFunctional)
    {
      RealVect pt;
      IntVect iv = a_face.gridIndex(Side::Hi);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(idir != a_face.direction())
            {
              Real ptval = a_dx[a_idir]*(Real(iv[idir]) + 0.5);
              pt[idir] = ptval;
            }
          else
            {
              pt[idir] = a_dx[a_idir]*(Real(iv[idir]));
            }
        }
      RealVect normal = EBArith::getDomainNormal(a_idir, a_side);

      value = m_func->value(pt, normal, a_time,a_icomp);

    }
  else
    {
      value = m_value;
    }
  a_faceFlux = value;
}

DirichletPoissonDomainBC::DirichletPoissonDomainBC()
 : m_onlyHomogeneous(true),
   m_isFunctional(false),
   m_value(12345.6789),
   m_func(RefCountedPtr<BaseBCValue>()),
   m_ebOrder(2)
{
}

DirichletPoissonDomainBC::~DirichletPoissonDomainBC()
{
}

void DirichletPoissonDomainBC::setValue(Real a_value)
{
   m_onlyHomogeneous = false;
   m_isFunctional = false;
   m_value = a_value;
   m_func = RefCountedPtr<BaseBCValue>();
}

void DirichletPoissonDomainBC::setFunction(RefCountedPtr<BaseBCValue> a_func)
{
  m_value = 12345.6789;
  m_func = a_func;

  m_onlyHomogeneous = false;
  m_isFunctional = true;
}

void DirichletPoissonDomainBC::setEBOrder(int a_ebOrder)
{
  CH_assert(a_ebOrder >= 1);
  CH_assert(a_ebOrder <= 2);

  m_ebOrder = a_ebOrder;
}

void DirichletPoissonDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
                                           const BaseFab<Real>&  a_phi,
                                           const RealVect&       a_probLo,
                                           const RealVect&       a_dx,
                                           const int&            a_idir,
                                           const Side::LoHiSide& a_side,
                                           const DataIndex&      a_dit,
                                           const Real&           a_time,
                                           const bool&           a_useHomogeneous)
{
  //CH_assert(a_phi.nComp() == 1);
  for(int comp=0; comp<a_phi.nComp(); comp++) {

    const Box& box = a_faceFlux.box();

    int iside;

    if (a_side == Side::Lo)
      {
        iside = 1;
      }
    else
      {
        iside = -1;
      }

    if (a_useHomogeneous)
      {
        Real value = 0.0;

//      if(m_currentBox.size(a_idir) >= BOXSIZE_HO)
//        {
//          FORT_SETHODIRICHLETFACEFLUX(CHF_FRA1(a_faceFlux,comp),
//                                      CHF_CONST_FRA1(a_phi,comp),
//                                      CHF_CONST_REAL(value),
//                                      CHF_CONST_REALVECT(a_dx),
//                                      CHF_CONST_INT(a_idir),
//                                      CHF_CONST_INT(iside),
//                                      CHF_BOX(box));
//        }
//      else
//        {
            FORT_SETDIRICHLETFACEFLUX(CHF_FRA1(a_faceFlux,comp),
                                      CHF_CONST_FRA1(a_phi,comp),
                                      CHF_CONST_REAL(value),
                                      CHF_CONST_REALVECT(a_dx),
                                      CHF_CONST_INT(a_idir),
                                      CHF_CONST_INT(iside),
                                      CHF_BOX(box));
//        }
      }
    else
      {
        if (m_isFunctional)
          {
            Real ihdx;

            ihdx = 2.0 / a_dx[a_idir];

            BoxIterator bit(box);

//          if(m_currentBox.size(a_idir) >= BOXSIZE_HO)
//            {
//              for (bit.begin(); bit.ok(); ++bit)
//                {
//                  IntVect iv = bit();
//                  IntVect ivNeigh = iv;
//                  ivNeigh[a_idir] += sign(a_side);
//                  const VolIndex vof      = VolIndex(iv,     0);
//                  const VolIndex vofNeigh = VolIndex(ivNeigh,0);
//                  const FaceIndex face = FaceIndex(vof,vofNeigh,a_idir);
//                  const RealVect  point = EBArith::getFaceLocation(face,a_dx,a_probLo);
//                  const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
//                  Real value = m_func->value(face,a_side,a_dit,point,normal,a_time,comp);
//                  Real phi0 = a_phi(iv,comp);
//                  iv[a_idir] += iside;
//                  Real phi1 = a_phi(iv,comp);
//                  iv[a_idir] -= iside;
//                  a_faceFlux(iv,comp) = iside*(8.0*value + phi1 - 9.0*phi0)/(3.0*a_dx[a_idir]);
//                }
//            }
//          else
//            {
                for (bit.begin(); bit.ok(); ++bit)
                  {
                    IntVect iv = bit();
                    IntVect ivNeigh = iv;
                    ivNeigh[a_idir] += sign(a_side);
                    const VolIndex vof      = VolIndex(iv,     0);
                    const VolIndex vofNeigh = VolIndex(ivNeigh,0);
                    const FaceIndex face = FaceIndex(vof,vofNeigh,a_idir);
                    const RealVect  point = EBArith::getFaceLocation(face,a_dx,a_probLo);
                    const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
                    Real value = m_func->value(face,a_side,a_dit,point,normal,a_time,comp);
                    Real phiVal = a_phi(iv,comp);
                    a_faceFlux(iv,comp) = iside * ihdx * (phiVal - value);
                  }
//            }
          }
        else
          {
            if (m_onlyHomogeneous)
              {
                MayDay::Error("DirichletPoissonDomainBC::getFaceFlux called with undefined inhomogeneous BC");
              }

            Real value = m_value;

//          if(m_currentBox.size(a_idir) >= BOXSIZE_HO)
//            {
//              FORT_SETHODIRICHLETFACEFLUX(CHF_FRA1(a_faceFlux,comp),
//                                          CHF_CONST_FRA1(a_phi,comp),
//                                          CHF_CONST_REAL(value),
//                                          CHF_CONST_REALVECT(a_dx),
//                                          CHF_CONST_INT(a_idir),
//                                          CHF_CONST_INT(iside),
//                                          CHF_BOX(box));
//            }
//          else
//            {
                FORT_SETDIRICHLETFACEFLUX(CHF_FRA1(a_faceFlux,comp),
                                          CHF_CONST_FRA1(a_phi,comp),
                                          CHF_CONST_REAL(value),
                                          CHF_CONST_REALVECT(a_dx),
                                          CHF_CONST_INT(a_idir),
                                          CHF_CONST_INT(iside),
                                          CHF_BOX(box));
//            }
          }
      }
  }
}

void DirichletPoissonDomainBC::getFaceFlux(Real&                 a_faceFlux,
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
  const EBISBox& ebisBox = a_phi.getEBISBox();

  a_faceFlux = 0.0;
  Vector<FaceIndex> faces = ebisBox.getFaces(a_vof,a_idir,a_side);
  for (int i = 0; i < faces.size(); i++)
    {
      const RealVect centroid = ebisBox.centroid(faces[i]);
      Real thisFaceFlux;
      getFaceFluxGradPhi(thisFaceFlux,faces[i],a_comp,a_phi,a_probLo,
                     a_dx,a_idir,a_side,a_dit,a_time,true,centroid,a_useHomogeneous);
      a_faceFlux += thisFaceFlux;
    }
}

void DirichletPoissonDomainBC::getFaceFluxGradPhi(Real&                 a_faceFlux,
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
  int iside = -sign(a_side);
  const Real ihdx = 2.0 / a_dx[a_idir];
  const EBISBox& ebisBox = a_phi.getEBISBox();

  Real value = -1.e99;
  if (a_useHomogeneous)
    {
      value = 0.0;
    }
  else if (m_isFunctional)
    {
      const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
      RealVect point  = a_centroid;
      point *= a_dx;
      point[a_face.direction()] = 0.0;//make this not depend on whatever ebisbox is returning for centroid in the face direction.
      point += EBArith::getFaceLocation(a_face,a_dx,a_probLo);
      value = m_func->value(a_face,a_side,a_dit,point,normal,a_time,a_comp);
    }
  else
    {
      if (m_onlyHomogeneous)
        {
          MayDay::Error("DirichletPoissonDomainBC::getFaceFlux called with undefined inhomogeneous BC");
        }
      value = m_value;
    }

  const VolIndex& vof = a_face.getVoF(flip(a_side));
//   if(m_currentBox.size(a_idir) >= BOXSIZE_HO)
//     {
//   Vector<FaceIndex> facesInsideDomain = ebisBox.getFaces(vof,a_idir,flip(a_side));
//   if(facesInsideDomain.size() == 1)
//     {
//       const VolIndex& vofNextInsideDomain = facesInsideDomain[0].getVoF(flip(a_side));
//       a_faceFlux = iside*(9.0*a_phi(vof,a_comp) - a_phi(vofNextInsideDomain,a_comp) - 8.0*value)/(3.0*a_dx[a_idir]);
//     }
//   else
//     {
//       a_faceFlux = iside * ihdx * (a_phi(vof,a_comp) - value);
//     }
//     }
//   else
//     {
  a_faceFlux = iside * ihdx * (a_phi(vof,a_comp) - value);
//     }

  if(a_useAreaFrac)
    {
      a_faceFlux *= ebisBox.areaFrac(a_face);
    }
}

void DirichletPoissonDomainBC::getFaceGradPhi(Real&                 a_faceFlux,
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
  int iside = -sign(a_side);
  const EBISBox& ebisBox = a_phi.getEBISBox();

  Real value = -1.e99;
  if (a_useHomogeneous)
    {
      value = 0.0;
    }
  else if (m_isFunctional)
    {
      const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
      RealVect point  = a_centroid;
      point *= a_dx;
      point[a_face.direction()] = 0.0;//make this not depend on whatever ebisbox is returning for centroid in the face direction.
      point += EBArith::getFaceLocation(a_face,a_dx,a_probLo);
      value = m_func->value(a_face,a_side,a_dit,point,normal,a_time,a_comp);
    }
  else
    {
      if (m_onlyHomogeneous)
        {
          MayDay::Error("DirichletPoissonDomainBC::getFaceFlux called with undefined inhomogeneous BC");
        }
      value = m_value;
    }

  const VolIndex& vof1 = a_face.getVoF(flip(a_side));
  Vector<FaceIndex> faces1 = ebisBox.getFaces(vof1,a_idir,flip(a_side));
  bool notOneSided = ((faces1.size()==1) && (!faces1[0].isBoundary()));
  if(notOneSided)
    {
      const VolIndex& vof2 = faces1[0].getVoF(flip(a_side));
      Vector<FaceIndex> faces2 = ebisBox.getFaces(vof2,a_idir,flip(a_side));
      bool linearExtrap = (faces2.size()==1 && (!faces2[0].isBoundary()));
      if(linearExtrap)
        {
          const VolIndex& vof3 = faces2[0].getVoF(flip(a_side));
          Vector<FaceIndex> faces3 = ebisBox.getFaces(vof3,a_idir,flip(a_side));
          bool quadExtrap = (faces3.size()==1 && (!faces3[0].isBoundary()));
          if(quadExtrap)
            {
              const VolIndex& vof4 = faces3[0].getVoF(flip(a_side));
              a_faceFlux = iside*(-4.*a_phi(vof1,a_comp) + 7.*a_phi(vof2,a_comp) 
                                  -4.*a_phi(vof3,a_comp) + a_phi(vof4,a_comp))/(2.*a_dx[a_idir]);
            }
          else
            {
              a_faceFlux = iside*(-3.*a_phi(vof1,a_comp) + 4.*a_phi(vof2,a_comp) 
                                  -a_phi(vof3,a_comp))/(2.*a_dx[a_idir]);
            }
        }
      else
        {
          a_faceFlux = iside*(a_phi(vof1,a_comp) - a_phi(vof2,a_comp))/a_dx[a_idir];
        }
    }
  else
    {
      Tuple<int,SpaceDim-1> tanDirs = PolyGeom::computeTanDirs(a_idir);
      IntVect iv = vof1.gridIndex();
      Real numGradFound = 0.;
      for(int itan = 0; itan < SpaceDim-1; itan++)
        {
          int tanDir = tanDirs[itan];
          for(SideIterator sit; sit.ok(); ++sit)
            {
              Vector<FaceIndex> facesTan = ebisBox.getFaces(vof1,tanDir,sit());
              if(facesTan.size()!=0 && !facesTan[0].isBoundary())
                {
                  IntVect ivTan = iv + sign(sit())*BASISV(tanDir);
                  VolIndex vofTan = ebisBox.getVoFs(ivTan)[0];
                  a_faceFlux += -sign(sit())*(a_phi(vof1,a_comp) - a_phi(vofTan,a_comp))/a_dx[a_idir];
                  numGradFound = numGradFound + 1.;
                }
            }
        }
      if(numGradFound > 0.)
        {
          a_faceFlux /= numGradFound;
        }
      else
        {
          a_faceFlux = 0.;
        }
    }
  
}
DirichletPoissonDomainBCFactory::DirichletPoissonDomainBCFactory()
: m_onlyHomogeneous(true),
  m_isFunctional(false),
  m_value(12345.6789),
  m_func(RefCountedPtr<BaseBCValue>()),
  m_ebOrder(2)
{
}

DirichletPoissonDomainBCFactory::~DirichletPoissonDomainBCFactory()
{
}

void DirichletPoissonDomainBCFactory::setValue(Real a_value)
{
  m_value = a_value;
  m_func = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = false;
  m_isFunctional = false;
}

void DirichletPoissonDomainBCFactory::setFunction(RefCountedPtr<BaseBCValue> a_func)
{
  m_value = 12345.6789;
  m_func = a_func;

  m_onlyHomogeneous = false;
  m_isFunctional = true;
}

void DirichletPoissonDomainBCFactory::setEBOrder(int a_ebOrder)
{
  CH_assert(a_ebOrder >= 1);
  CH_assert(a_ebOrder <= 2);

  m_ebOrder = a_ebOrder;
}

DirichletPoissonDomainBC*
DirichletPoissonDomainBCFactory::
create(const ProblemDomain& a_domain,
       const EBISLayout&    a_layout,
       const RealVect&      a_dx)
{
  DirichletPoissonDomainBC* newBC = new DirichletPoissonDomainBC();
  newBC->setEBOrder(m_ebOrder);
  if(!m_onlyHomogeneous)
    {
      if(m_isFunctional)
        {
          newBC->setFunction(m_func);
        }
      else
        {
          newBC->setValue(m_value);
        }
    }
  return newBC;
}
#include "NamespaceFooter.H"
