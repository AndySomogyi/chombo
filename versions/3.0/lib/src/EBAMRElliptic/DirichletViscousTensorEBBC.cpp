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

#include "BoxIterator.H"
#include "VoFIterator.H"

#include "EBArith.H"
#include "PolyGeom.H"

#include "DirichletViscousTensorEBBC.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

DirichletViscousTensorEBBC::
DirichletViscousTensorEBBC()
{
}


DirichletViscousTensorEBBC::
DirichletViscousTensorEBBC(
                           const ProblemDomain& a_domain,
                           const EBISLayout&    a_layout,
                           const RealVect&      a_dx,
                           const IntVect*       a_ghostCellsPhi,
                           const IntVect*       a_ghostCellsRhs)
  : m_ghostCellsPhi( *a_ghostCellsPhi ),
    m_ghostCellsRHS( *a_ghostCellsRhs )
{
  CH_assert( bool(a_ghostCellsPhi) == bool(a_ghostCellsRhs) ); // !xor

  m_domain = a_domain;
  m_ebisl = a_layout;

  m_dx = a_dx;

  m_isDefined = false;
}

DirichletViscousTensorEBBC::
~DirichletViscousTensorEBBC()
{
}

void
DirichletViscousTensorEBBC::
define(const LayoutData<IntVectSet>& a_cfivs,
       const Real&                   a_factor)
{
  //factor is 1 here.  ignore.
  CH_assert(m_coefSet);
  const DisjointBoxLayout& dbl = m_ebisl.getDisjointLayout();

  LayoutData<IntVectSet>   irregIVS;
  LayoutData<VoFIterator > vofItIrreg;
  vofItIrreg.define(dbl); // vofiterator cache
  irregIVS.define(dbl);
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      irregIVS[dit()] = m_ebisl[dit()].getIrregIVS(dbl.get(dit()));
      vofItIrreg[dit()].define(irregIVS[dit()],m_ebisl[dit()].getEBGraph());
    }

  //make the Dirichlet stencils
  for(int ivar = 0; ivar < SpaceDim; ivar++)
    {
      m_fluxStencil[ivar].define(dbl);
      m_fluxWeight [ivar].define(dbl);

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {

          m_fluxStencil[ivar][dit()].define(irregIVS[dit()],m_ebisl[dit()].getEBGraph(), 1);
          m_fluxWeight [ivar][dit()].define(irregIVS[dit()],m_ebisl[dit()].getEBGraph(), 1);


          for (vofItIrreg[dit()].reset(); vofItIrreg[dit()].ok(); ++vofItIrreg[dit()])
            {
              Real        fluxWeight[SpaceDim];
              VoFStencil fluxStencil[SpaceDim];
              getFluxStencil(fluxStencil, fluxWeight, dit(),
                             vofItIrreg[dit()](),m_ebisl[dit()],
                             m_dx,a_cfivs[dit()]);

              m_fluxStencil[ivar][dit()](vofItIrreg[dit()](), 0)  = fluxStencil[ivar];
              m_fluxWeight [ivar][dit()](vofItIrreg[dit()](), 0)  =  fluxWeight[ivar];
            }
        }
    }

  m_isDefined = true;
}
//full normal grad--for debugging when ebstencil stuff is turned off
void
DirichletViscousTensorEBBC::
getNormalGradSlow(Real              a_normalGrad[CH_SPACEDIM],
                  const DataIndex&  a_dit,
                  const VolIndex &  a_vof,
                  const RealVect&   a_probLo,
                  const EBISBox  &  a_ebisBox,
                  const EBCellFAB&  a_phi)
{
  for(int icomp = 0; icomp < SpaceDim; icomp++)
    {
      Real value;
      if (m_isFunction)
        {
          // Compute the bndryCentroid location in physical coordinates
          RealVect probLo = a_ebisBox.bndryCentroid(a_vof);
          probLo *= m_dx[0];
          probLo += a_probLo;
          RealVect point = EBArith::getVofLocation(a_vof, m_dx, probLo);
          value = m_func->value(point,  icomp);
        }
      else
        {
          value = m_value;
        }
      VoFStencil normalStencil;
      Real weight;
      IntVectSet cfivs;//this will break in eb-cf land
      getNormalStencil(normalStencil, weight,
                       a_vof, a_ebisBox, m_dx, cfivs);

      normalStencil.setAllVariables(icomp);

      Real   homogComp = applyVoFStencil(normalStencil, a_phi, icomp);
      Real inhomogComp = value*weight;
      a_normalGrad[icomp] = homogComp + inhomogComp;
    }
}
//inhomogeneous contribution only
void
DirichletViscousTensorEBBC::
getNormalGradFast(Real              a_normalGrad[CH_SPACEDIM],
                  const DataIndex&  a_dit,
                  const VolIndex &  a_vof,
                  const RealVect&   a_probLo,
                  const EBISBox  &  a_ebisBox)
{
  for(int icomp = 0; icomp < SpaceDim; icomp++)
    {
      Real value;
      if (m_isFunction)
        {
          // Compute the bndryCentroid location in physical coordinates
          RealVect probLo = a_ebisBox.bndryCentroid(a_vof);
          probLo *= m_dx[0];
          probLo += a_probLo;
          RealVect point = EBArith::getVofLocation(a_vof, m_dx, probLo);
          value = m_func->value(point,  icomp);
        }
      else
        {
          value = m_value;
        }
      a_normalGrad[icomp] = value*m_fluxWeight[icomp][a_dit](a_vof,0);
    }
}
//go from normal-tangential coordinates to xy coords
void
transformToNTFromXY(Real                  a_gradNT[CH_SPACEDIM][CH_SPACEDIM],
                    const Real            a_gradXY[CH_SPACEDIM][CH_SPACEDIM],
                    const RealVect&       a_normal,
                    const RealVect        a_tangen[CH_SPACEDIM-1])
{
  Real matrix[CH_SPACEDIM][CH_SPACEDIM];
  for(int irow = 0; irow < SpaceDim; irow++)
    {
      RealVect tvect;
      if(irow == 0)
        {
          tvect = a_normal;
        }
      else
        {
          tvect = a_tangen[irow-1];
        }
      for(int icol = 0; icol < SpaceDim; icol++)
        {
          matrix[irow][icol] = tvect[icol];
        }
    }

  for(int icomp = 0; icomp < SpaceDim; icomp++)
    {
      for(int irow = 0; irow < SpaceDim; irow++)
        {
          a_gradNT[icomp][irow] = 0;
          for(int icol = 0; icol < SpaceDim; icol++)
            {
              a_gradNT[icomp][irow] += matrix[irow][icol]*a_gradXY[icomp][icol];
            }
        }
    }
}
//go from xy coords to normal-tangential coordinates
void
transformToXYFromNT(Real                  a_gradXY[CH_SPACEDIM][CH_SPACEDIM],
                    const Real            a_gradNT[CH_SPACEDIM][CH_SPACEDIM],
                    const RealVect&       a_normal,
                    const RealVect        a_tangen[CH_SPACEDIM-1])
{
  //set up transformation matrix
  Real matrix[CH_SPACEDIM][CH_SPACEDIM];
  for(int irow = 0; irow < SpaceDim; irow++)
    {
      RealVect tvect;
      if(irow == 0)
        {
          tvect = a_normal;
        }
      else
        {
          tvect = a_tangen[irow-1];
        }
      for(int icol = 0; icol < SpaceDim; icol++)
        {
          //since the vectors are orthogonal and unit vectors
          //transpose = inverse
          matrix[icol][irow] = tvect[icol];
        }
    }

  for(int icomp = 0; icomp < SpaceDim; icomp++)
    {
      for(int irow = 0; irow < SpaceDim; irow++)
        {
          a_gradXY[icomp][irow] = 0;
          for(int icol = 0; icol < SpaceDim; icol++)
            {
              a_gradXY[icomp][irow] += matrix[irow][icol]*a_gradNT[icomp][icol];
            }
        }
    }
}
//given a normal vector, get spacedim-1 tangential vectors that are
//unit vectors and
void
getTangentVectors( RealVect              a_tangen[CH_SPACEDIM-1],
                   const RealVect&       a_normal)
{
#if CH_SPACEDIM==2
  a_tangen[0][1] =  a_normal[0];
  a_tangen[0][0] = -a_normal[1];
#else
  //get some vector that is not parallel to a_normal
  RealVect notNormal = a_normal + RealVect::Unit;
  Real dotProd = PolyGeom::dot(notNormal, a_normal);
  //there is one bizzare case where I have just scaled up a_normal
  //in this case, add something diff to each component.
  if(Abs(dotProd) < 1.0e-6)
    {
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          notNormal[idir] += Real(idir+1);
        }
      dotProd = PolyGeom::dot(notNormal, a_normal);
      CH_assert(Abs(dotProd) > 1.0e-6);
    }
  Real sumSquare;
  PolyGeom::unifyVector(notNormal, sumSquare);
  //make the first tangential vector = cross(notnormal, normal);
  a_tangen[0] = PolyGeom::cross(notNormal,   a_normal);
  //make the second one the cross between that and normal
  a_tangen[1] = PolyGeom::cross(a_tangen[0], a_normal);
#endif

  //make them all unit vectors
  for(int itan = 0; itan < CH_SPACEDIM-1; itan++)
    {
      Real sumSquare;
      PolyGeom::unifyVector(a_tangen[itan], sumSquare);
    }
}
void
DirichletViscousTensorEBBC::
applyEBFlux(EBCellFAB&                    a_lphi,
            const EBCellFAB&              a_phi,
            VoFIterator&                  a_vofit,
            const LayoutData<IntVectSet>& a_cfivs,
            const DataIndex&              a_dit,
            const RealVect&               a_probLo,
            const RealVect&               a_dx,
            const Real&                   a_factor,
            const bool&                   a_useHomogeneous,
            const Real&                   a_time)
{
  CH_assert(m_coefSet);
  const EBISBox&   ebisBox = a_phi.getEBISBox();
  for(a_vofit.reset(); a_vofit.ok(); ++a_vofit)
    {
      const VolIndex& vof = a_vofit();

      Real gradXY[SpaceDim][SpaceDim];
      //get boundary gradient in x-y coordinates
      getBoundaryGrad(gradXY, vof, a_dx, a_probLo, ebisBox);

      Real normalGrad[SpaceDim];
      if(m_nullReturned)
        {
          //slow way-- not using EBStencil stuff--- for debugging
          //get full normal gradient (not just inhomogeneous contribution)
          getNormalGradSlow(normalGrad, a_dit, vof, a_probLo, ebisBox, a_phi);
        }
      else
        {
          //fast way--  for when EBStencil are being used
          //this is just the contribution of the inhomogeneous term
          getNormalGradFast(normalGrad, a_dit, vof, a_probLo, ebisBox);
        }

      //get normal and tangential vectors
      RealVect       tangen[CH_SPACEDIM-1];
      RealVect       normal= ebisBox.normal(vof);
      getTangentVectors(tangen, normal);

      //transform the gradient into normal-tangential space.
      Real gradNT[SpaceDim][SpaceDim];

      //grad is of the form grad[comp][dir]
      transformToNTFromXY(gradNT, gradXY, normal, tangen);

      //replace the normal griadients with normalgrad
      //first index is normal, tan1, tan2
      for(int icomp = 0; icomp < SpaceDim; icomp++)
        {
          gradNT[icomp][0] = normalGrad[icomp];
        }
      //transform back to XY space and we can get the flux from that
      transformToXYFromNT(gradXY, gradNT, normal, tangen);

      Real flux[SpaceDim][SpaceDim];

      getFluxFromGrad(flux, gradXY, vof, a_dit);

      Real deltaLph[SpaceDim];
      getChangeInSolution(deltaLph, flux, a_dx, vof, a_dit, ebisBox);

      for(int comp = 0; comp < a_lphi.nComp(); comp++)
        {
          a_lphi(vof, comp) += deltaLph[comp];
        }

    }
}


void
DirichletViscousTensorEBBC::
getFluxStencil(VoFStencil*          a_stencils,
               Real*                a_weights,
               const DataIndex&     a_dit,
               const VolIndex&      a_vof,
               const EBISBox&       a_ebisBox,
               const RealVect&      a_dx,
               const IntVectSet&    a_cfivs)
{
  //So the boundary condition for this is that gphi = gphiNormal
  //assuming the value at the boundary is a given.   The tangential
  //derivatives are just externally applied so they have no stencil.
  //we will use the johansen stencil where to get the normal gradient and
  //the least squares thing where we cannot.
  VoFStencil normalStencil;
  Real weight;
  getNormalStencil(normalStencil, weight,
                   a_vof, a_ebisBox, a_dx, a_cfivs);

  //i do not know how to do this cleverly so I will just go step by step
  //compute the stencils for all gradients given that we have the normal
  //gradient and we assume tangential grads = 0 (fixed in apply EBFlux)
  VoFStencil gradStencil[SpaceDim][SpaceDim];
  RealVect normal = a_ebisBox.normal(a_vof);
  for(int icomp = 0; icomp < SpaceDim; icomp++)
    {
      for(int ideriv = 0; ideriv < SpaceDim; ideriv++)
        {
          gradStencil[icomp][ideriv] = normalStencil;
          gradStencil[icomp][ideriv].setAllVariables(icomp);
          gradStencil[icomp][ideriv] *= -normal[ideriv];
        }
    }
  //now get the stencil for the divergence
  VoFStencil divergenceStencil;
  for(int ideriv = 0; ideriv < SpaceDim; ideriv++)
    {
      divergenceStencil += gradStencil[ideriv][ideriv];
    }

  //now get the flux stencils
  //now for all the flux stencils Fij = eta*(gradij +gradji) + lambda*deltaij*divergence
  VoFStencil fluxStencil[SpaceDim][SpaceDim];
  Real lambdaPt = (*m_lambda)[a_dit](a_vof, 0);
  Real    etaPt =    (*m_eta)[a_dit](a_vof, 0);
  for(int icomp = 0; icomp < SpaceDim; icomp++)
    {
      for(int ideriv = 0; ideriv < SpaceDim; ideriv++)
        {
          if(icomp == ideriv)
            {
              fluxStencil[icomp][ideriv] += divergenceStencil;
              fluxStencil[icomp][ideriv] *= lambdaPt;
            }
          VoFStencil gradContrib;
          gradContrib += gradStencil[icomp ][ideriv];
          gradContrib += gradStencil[ideriv][icomp ];
          gradContrib *= etaPt;
          fluxStencil[icomp][ideriv] += gradContrib;
        }
    }
  //flux through face = sum flux[icomp][ideriv]*normal[deriv]
  for(int icomp = 0; icomp < SpaceDim; icomp++)
    {
      a_stencils[icomp].clear();
      for(int ideriv = 0; ideriv < SpaceDim; ideriv++)
        {
          VoFStencil derivContrib = fluxStencil[icomp][ideriv];
          derivContrib *= normal[ideriv];
          a_stencils[icomp] += derivContrib;
        }
      a_weights[icomp] = weight;
    }
}

void
DirichletViscousTensorEBBC::
getNormalStencil(VoFStencil&          a_stencil,
                 Real&                a_weight,
                 const VolIndex&      a_vof,
                 const EBISBox&       a_ebisBox,
                 const RealVect&      a_dx,
                 const IntVectSet&    a_cfivs)
{
  Vector<VoFStencil>  pointStencils;
  Vector<Real>        distanceAlongLine;
  bool needToDropOrder = getSecondOrderStencil(a_stencil, a_weight,
                                               pointStencils, distanceAlongLine,
                                               a_vof, a_ebisBox, a_dx, a_cfivs);
  if(needToDropOrder)
    {
      getFirstOrderStencil(a_stencil, a_weight, a_vof, a_ebisBox, a_dx, a_cfivs);
    }
}
bool
DirichletViscousTensorEBBC::
getSecondOrderStencil(VoFStencil&          a_stencil,
                      Real&                a_weight,
                      Vector<VoFStencil>&  a_pointStencils,
                      Vector<Real>&        a_distanceAlongLine,
                      const VolIndex&      a_vof,
                      const EBISBox&       a_ebisBox,
                      const RealVect&      a_dx,
                      const IntVectSet&    a_cfivs)
{

  a_stencil.clear();
  bool dropOrder = false;
  EBArith::johanStencil(dropOrder, a_pointStencils, a_distanceAlongLine,
                        a_vof, a_ebisBox, a_dx, a_cfivs);
  if(dropOrder)
    {
      return true;
    }

  //if we got this far, sizes should be at least 2
  CH_assert(a_distanceAlongLine.size() >= 2);
  CH_assert(a_pointStencils.size() >= 2);
  Real x1 = a_distanceAlongLine[0];
  Real x2 = a_distanceAlongLine[1];
  //fit quadratic function to line and find the gradient at the origin
  //grad = (x2*x2*(phi1-phi0)-x1*x1(phi2-phi0))/(x2*x2*x1 - x1*x1*x2);
  Real denom = x2*x2*x1 - x1*x1*x2;
  //not done by reference because i want point stencils checkable externally.
  VoFStencil phi1Sten = a_pointStencils[0];
  VoFStencil phi2Sten = a_pointStencils[1];
  phi1Sten *= x2*x2/denom;
  phi2Sten *=-x1*x1/denom;
  //weight is the multiplier of the inhomogeneous value (phi0)
  a_weight = x1*x1/denom - x2*x2/denom;
  a_stencil += phi1Sten;
  a_stencil += phi2Sten;
  //if we got this far, we have a second order stencil;
  return false;
}

void
DirichletViscousTensorEBBC::
getFirstOrderStencil(VoFStencil&       a_stencil,
                     Real&             a_weight,
                     const VolIndex&   a_vof,
                     const EBISBox&    a_ebisBox,
                     const RealVect&   a_dx,
                     const IntVectSet& a_cfivs)
{
  EBArith::getLeastSquaresGradSten(a_stencil, a_weight, a_vof, a_ebisBox, a_dx, m_domain);
  //we want stencil to be grad dot normal.  this returns -grad dot normal
  a_stencil *= -1.0;
  a_weight  *= -1.0;
}


DirichletViscousTensorEBBCFactory::DirichletViscousTensorEBBCFactory()
{
  m_value = 12345.6789;
  m_func = RefCountedPtr<BaseBCFuncEval>();

  m_isFunction = false;
}

DirichletViscousTensorEBBCFactory::~DirichletViscousTensorEBBCFactory()
{
}

void
DirichletViscousTensorEBBCFactory::
setValue(Real a_value)
{
  m_value = a_value;
  m_func = RefCountedPtr<BaseBCFuncEval>();

  m_isFunction = false;
}

void
DirichletViscousTensorEBBCFactory::
setFunction(RefCountedPtr<BaseBCFuncEval> a_func)
{
  m_value = 12345.6789;
  m_func = a_func;

  m_isFunction = true;
}

DirichletViscousTensorEBBC*
DirichletViscousTensorEBBCFactory::
create(
       const ProblemDomain& a_domain,
       const EBISLayout&    a_layout,
       const RealVect&      a_dx,
       const IntVect*       a_ghostCellsPhi,
       const IntVect*       a_ghostCellsRhs)
{
  CH_TIME("DirichletViscousTensorEBBC::create");
  DirichletViscousTensorEBBC* fresh;
  fresh = new DirichletViscousTensorEBBC(a_domain,a_layout,a_dx,a_ghostCellsPhi,
                                         a_ghostCellsRhs);

  if (m_isFunction)
    {
      fresh->setFunction(m_func);
    }
  else
    {
      fresh->setValue(m_value);
    }

  return fresh;
}

#include "NamespaceFooter.H"
