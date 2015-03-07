#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _DIRICHLETVISCOUSTENSOREBBC_H_
#define _DIRICHLETVISCOUSTENSOREBBC_H_

#include "RefCountedPtr.H"

#include "BaseEBBC.H"
#include "BaseBCFuncEval.H"
#include "LayoutData.H"
#include "IntVectSet.H"
#include "EBStencil.H"
#include "NamespaceHeader.H"
class VoFCache;
class VoFStencil;

///
/**
 */
class DirichletViscousTensorEBBC: public ViscousBaseEBBC
{
public:

  ///
  DirichletViscousTensorEBBC();

  ///
  /**
     If the optional a_ghostCells{Phi,Rhs} args are NULL, we won't use the
     VoF cache.
  */
  DirichletViscousTensorEBBC( const ProblemDomain& a_domain,
                              const EBISLayout&    a_layout,
                              const RealVect&      a_dx,
                              const IntVect*       a_ghostCellsPhi,
                              const IntVect*       a_ghostCellsRhs);

  ///
  /**
   */
  virtual ~DirichletViscousTensorEBBC();

  ///
  /**
   */
  virtual void define(const LayoutData<IntVectSet>& a_cfivs,
                      const Real&                   a_factor);

  virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int ivar)
  {
    CH_assert((ivar >= 0) && (ivar < SpaceDim));
    m_nullReturned = false;
    //debugging hook--turns off EBStencil stuff
    //m_nullReturned = true;
    // end debugging hook

    LayoutData<BaseIVFAB<VoFStencil> >* retval = NULL;
    if(!m_nullReturned)
      {
        //non debug mode
        return &m_fluxStencil[ivar];
      }

    return retval;
  }

  virtual void applyEBFlux(EBCellFAB&                    a_lphi,
                           const EBCellFAB&              a_phi,
                           VoFIterator&                  a_vofit,
                           const LayoutData<IntVectSet>& a_cfivs,
                           const DataIndex&              a_dit,
                           const RealVect&               a_probLo,
                           const RealVect&               a_dx,
                           const Real&                   a_factor,
                           const bool&                   a_useHomogeneous,
                           const Real&                   a_time);
  ///
  /**  deprecated interface not used by viscous operator
   */
  virtual void getEBFlux(Real&                         a_flux,
                         const VolIndex&               a_vof,
                         const LevelData<EBCellFAB>&   a_phi,
                         const LayoutData<IntVectSet>& a_cfivs,
                         const DataIndex&              a_dit,
                         const RealVect&               a_probLo,
                         const RealVect&               a_dx,
                         const bool&                   a_useHomogeneous,
                         const Real&                   a_time,
                         const pair<int,Real>*         a_cacheHint=0 )
  {
    MayDay::Error();
  }


protected:
  void  getNormalStencil(VoFStencil&          a_stencil,
                         Real&                a_weight,
                         const VolIndex&      a_vof,
                         const EBISBox&       a_ebisBox,
                         const RealVect&      a_dx,
                         const IntVectSet&    a_cfivs);

  void getFluxStencil(VoFStencil*       a_stencil,
                      Real*             a_weight,
                      const DataIndex&  a_dit,
                      const VolIndex&   a_vof,
                      const EBISBox&    a_ebisBox,
                      const RealVect&   a_dx,
                      const IntVectSet& a_cfivs);

  //return true if you need to drop order
  bool getSecondOrderStencil(VoFStencil&          a_stencil,
                             Real&                a_weight,
                             Vector<VoFStencil>&  a_pointStencil,
                             Vector<Real>&        a_distanceAlongLine,
                             const VolIndex&      a_vof,
                             const EBISBox&       a_ebisBox,
                             const RealVect&      a_dx,
                             const IntVectSet&    a_cfivs);

  void getFirstOrderStencil(VoFStencil&       a_stencil,
                            Real&             a_weight,
                            const VolIndex&   a_vof,
                            const EBISBox&    a_ebisBox,
                            const RealVect&   a_dx,
                            const IntVectSet& a_cfivs);

  //only inhomogeneous contribution
  void getNormalGradFast(Real              a_normalGrad[CH_SPACEDIM],
                         const DataIndex&  a_dit,
                         const VolIndex &  a_vof,
                         const RealVect&   a_probLo,
                         const EBISBox  &  a_ebisBox);


  //whole normal grad--- this is for debugging
  void getNormalGradSlow(Real              a_normalGrad[CH_SPACEDIM],
                         const DataIndex&  a_dit,
                         const VolIndex &  a_vof,
                         const RealVect&   a_probLo,
                         const EBISBox  &  a_ebisBox,
                         const EBCellFAB&  a_phi);


private:
  bool m_isDefined;
  //debugging hook
  bool m_nullReturned;
  //end debugging hook

  ProblemDomain m_domain;
  EBISLayout    m_ebisl;

  RealVect m_dx;


  //stencils for operator evaluation
  LayoutData<BaseIVFAB<VoFStencil> > m_fluxStencil[CH_SPACEDIM];
  LayoutData<BaseIVFAB<Real> >        m_fluxWeight[CH_SPACEDIM];
  const IntVect                   m_ghostCellsPhi;
  const IntVect                   m_ghostCellsRHS;


};

///
/**
 */
class DirichletViscousTensorEBBCFactory: public BaseEBBCFactory
{
public:
  ///
  /**
   */
  DirichletViscousTensorEBBCFactory();

  ///
  /**
   */
  virtual ~DirichletViscousTensorEBBCFactory();


  ///
  /**
   */
  virtual void setValue(Real a_value);

  ///
  /**
   */
  virtual void setFunction(RefCountedPtr<BaseBCFuncEval> a_func);

  ///
  /**
   */
  virtual DirichletViscousTensorEBBC* create(  const ProblemDomain& a_domain,
                                               const EBISLayout&    a_layout,
                                               const RealVect&      a_dx,
                                               const IntVect*       a_ghostCellsPhi,
                                               const IntVect*       a_ghostCellsRhs);

protected:

private:
  bool m_isFunction;

  Real m_value;
  RefCountedPtr<BaseBCFuncEval>   m_func;
  const IntVect                   m_ghostCellsPhi;
  const IntVect                   m_ghostCellsRHS;

};



#include "NamespaceFooter.H"
#endif
