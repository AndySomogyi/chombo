#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// VCHelmHoltzOp.cpp
// petermc, 22 Nov 2002

#include "VCHelmholtzOp.H"

// ---------------------------------------------------------
void 
VCHelmholtzOp::define(const DisjointBoxLayout& a_Ba,
                      const DisjointBoxLayout* a_baseBa,
                      Real  a_dxLevel,
                      int a_refratio,
                      const Box& a_domain,
                      bool a_homogeneousOnly,
                      int a_ncomp)
{
  ProblemDomain probdomain(a_domain);
  define(a_Ba, a_baseBa, a_dxLevel, a_refratio, probdomain, a_homogeneousOnly, a_ncomp);
}


// ---------------------------------------------------------
void 
VCHelmholtzOp::define(const DisjointBoxLayout& a_Ba,
                      const DisjointBoxLayout* a_baseBa,
                      Real  a_dxLevel,
                      int a_refratio,
                      const ProblemDomain& a_probdomain,
                      bool a_homogeneousOnly,
                      int a_ncomp)
{
  // Call this when defining your first VCHelmholtzOp.
  // Then call setDomainGhostBC with boundary conditions.

  if (m_vcdivop == NULL)
    {
      m_vcdivop = new VCDivOp();
      m_vcdivop->defineVC(a_Ba, a_ncomp);
    }
  else if ( m_vcdivop->isVCDefined() &&
            ! (a_Ba == m_vcdivop->m_diagn.getBoxes()) )
    {
      // New set of grids.
      delete m_vcdivop;
      m_vcdivop = new VCDivOp();
      m_vcdivop->defineVC(a_Ba, a_ncomp);
    }
  m_vcdivop->define(a_Ba, a_baseBa, a_dxLevel, a_refratio, a_probdomain,
                    a_homogeneousOnly, a_ncomp);
  // If m_vcdivop is non-null, then it has likely been filled,
  // so you don't call defineVC because that would erase the coefficients.
  if (! m_vcdivop->isVCDefined() )
    m_vcdivop->defineVC(a_Ba, a_ncomp);
}


// ---------------------------------------------------------
void 
VCHelmholtzOp::define(const LevelOp* a_opfine, int a_reftoFine)
{
  const VCHelmholtzOp* opfineptr = dynamic_cast<const VCHelmholtzOp*>(a_opfine);
  if (opfineptr == NULL)
    MayDay::Error("VCHelmholtzOp::define: casting failed");
  const VCHelmholtzOp& opfine = *opfineptr;

  const VCDivOp* opfineDivOp = opfineptr->m_vcdivop;
  // The coeffs in m_vcdivop are set by coarsening those in opfineDivOp.
  m_vcdivop->define(opfineDivOp, a_reftoFine);
  CH_assert(opfine.isDefined());
}


// ---------------------------------------------------------
void
VCHelmholtzOp::smooth(LevelData<FArrayBox>& a_phi,
                      const LevelData<FArrayBox>& a_rhs)
{
  CH_assert(isDefined());
  m_vcdivop->smooth(a_phi, a_rhs);
}


// ---------------------------------------------------------
void
VCHelmholtzOp::applyOpI(LevelData<FArrayBox>& a_phi, 
                        const LevelData<FArrayBox>* a_phicPtr, 
                        LevelData<FArrayBox>& a_lofPhi)
{
  CH_assert(isDefined());
  m_vcdivop->applyOpI(a_phi, a_phicPtr, a_lofPhi);
}


// ---------------------------------------------------------
void
VCHelmholtzOp::applyOpIcfHphys(LevelData<FArrayBox>& a_phi, 
                               const LevelData<FArrayBox>* a_phicPtr, 
                               LevelData<FArrayBox>& a_lofPhi)
{
  CH_assert(isDefined());
  m_vcdivop->applyOpIcfHphys(a_phi, a_phicPtr, a_lofPhi);
}


// ---------------------------------------------------------
void
VCHelmholtzOp::applyOpH(LevelData<FArrayBox>& a_phi, 
                        LevelData<FArrayBox>& a_lofPhi)
{
  CH_assert(isDefined());
  m_vcdivop->applyOpH(a_phi, a_lofPhi);
}


// ---------------------------------------------------------
void
VCHelmholtzOp::applyOpHcfIphys(LevelData<FArrayBox>& a_phi, 
                               LevelData<FArrayBox>& a_lofPhi)
{
  CH_assert(isDefined());
  m_vcdivop->applyOpHcfIphys(a_phi, a_lofPhi);
}


// ---------------------------------------------------------
void
VCHelmholtzOp::getFlux(
                       FArrayBox& a_fineFlux,
                       const FArrayBox& a_data,
                       const DataIndex& a_dit, 
                       int a_dir)
{
  CH_assert(isDefined());
  m_vcdivop->getFlux(a_fineFlux, a_data, a_dit, a_dir);
}


// ---------------------------------------------------------
void
VCHelmholtzOp::levelPreconditioner(LevelData<FArrayBox> & a_phihat,
                                   const LevelData<FArrayBox>& a_rhshat)
{
  CH_assert(isDefined());
  m_vcdivop->levelPreconditioner(a_phihat, a_rhshat);
}


// ---------------------------------------------------------
void
VCHelmholtzOp::setAlpha(const LevelData<FArrayBox>& a_diag)
{
  CH_assert(isDefined());
  m_vcdivop->setDiagCoeff(a_diag);
}


// ---------------------------------------------------------
void
VCHelmholtzOp::setBeta(const LevelData<FluxBox>& a_coeff,
                       Real a_dxLevel)
{
  CH_assert(isDefined());
  m_vcdivop->setDivCoeff(a_coeff, a_dxLevel);
}


// ---------------------------------------------------------
void 
VCHelmholtzOp::scaleAlpha(Real a_scale)
{
  CH_assert(isDefined());
  CH_assert(m_vcdivop->isVCDefined());
  LevelData<FArrayBox>& diagn = m_vcdivop->m_diagn;
  for (DataIterator dit = diagn.dataIterator(); dit.ok(); ++dit)
    diagn[dit()].mult(a_scale);

  // since you've changed alpha, you'll need to recalculate lambda.
  m_vcdivop->needNewLambda();

  // compute m_lambda
  // m_vcdivop->computeInverseDiag();
  // m_vcdivop->computeVC(*m_vcdivop);
}


// ---------------------------------------------------------
void
VCHelmholtzOp::scaleBeta(Real a_scale)
{
  CH_assert(isDefined());
  CH_assert(m_vcdivop->isVCDefined());
  DataIterator dit = (m_vcdivop->m_diagn).dataIterator();
  for (dit.begin(); dit.ok(); ++dit) 
    for (int dir=0; dir<SpaceDim; dir++) 
      {
        FArrayBox& lowerCoeff = m_vcdivop->m_lower[dir][dit()];
        FArrayBox& upperCoeff = m_vcdivop->m_upper[dir][dit()];
        FArrayBox& centrCoeff = m_vcdivop->m_centr[dir][dit()];
        lowerCoeff.mult(a_scale);
        upperCoeff.mult(a_scale);
        centrCoeff.mult(a_scale);
      }
  // since you've changed beta, you'll need to recalculate lambda.
  m_vcdivop->needNewLambda();

  // compute m_lambda
  // m_vcdivop->computeInverseDiag();
  // m_vcdivop->computeVC(*m_vcdivop);
}


// ---------------------------------------------------------
LevelOp*
VCHelmholtzOp::new_levelop() const
{
  CH_assert(isDefined());
  VCHelmholtzOp* osh_ptr = new VCHelmholtzOp();
  // if coeffs of m_vcdivop are defined, then they are copied
  // to osh_ptr->m_vcdivop.
  osh_ptr->m_vcdivop = (VCDivOp*) m_vcdivop->new_levelop();
  if (osh_ptr == NULL)
    {
      MayDay::Error("Out of Memory in VCHelmholtzOp::new_levelop");
    }
  // CH_assert(osh_ptr->isDefined());
  // m_domainbc is defined in m_vcdivop
  return static_cast<LevelOp*>(osh_ptr);
}


// ---------------------------------------------------------
void
VCHelmholtzOp::setDomainGhostBC(const DomainGhostBC& a_dombc)
{
  CH_assert(m_vcdivop != NULL);
  m_vcdivop->setDomainGhostBC(a_dombc);
}


// ---------------------------------------------------------
void
VCHelmholtzOp::CFInterp(
                        LevelData<FArrayBox>& phi,
                        const LevelData<FArrayBox>& phiCoarse
                        )
{
  CH_assert(isDefined());
  m_vcdivop->CFInterp(phi, phiCoarse);
}


// ---------------------------------------------------------
void
VCHelmholtzOp::homogeneousCFInterp(LevelData<FArrayBox>& phi)
{
  CH_assert(isDefined());
  m_vcdivop->homogeneousCFInterp(phi);
}


// ---------------------------------------------------------
bool
VCHelmholtzOp::isDefined() const
{
  if (m_vcdivop == NULL) return false;
  bool isDivOpDefined = m_vcdivop->isDefined();
  return (isDivOpDefined);
}


// ---------------------------------------------------------
void
VCHelmholtzOp::bottomSmoother(LevelData<FArrayBox>& phi,
                              const LevelData<FArrayBox>& rhs)
{
  CH_assert(isDefined());
  m_vcdivop->bottomSmoother(phi, rhs);
}



void
VCHelmholtzOp::setConvergenceMetric(Real a_metric, int a_comp)
{
  m_vcdivop->setConvergenceMetric(a_metric, a_comp);
}
