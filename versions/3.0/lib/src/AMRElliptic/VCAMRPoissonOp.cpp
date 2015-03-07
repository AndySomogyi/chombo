#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "FORT_PROTO.H"
#include "BoxIterator.H"
#include "AverageF_F.H"
#include "InterpF_F.H"
#include "LayoutIterator.H"
#include "FineInterp.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "AMRIO.H"

#include "VCAMRPoissonOp.H"

#include "VCAMRPoissonOpF_F.H"
#include "AMRPoissonOpF_F.H"

#include "NamespaceHeader.H"

void
VCAMRPoissonOp::
setAlphaAndBeta(const Real& a_alpha,
                const Real& a_beta)
{
  m_alpha = a_alpha;
  m_beta = a_beta;
  resetLambda();
}
void
VCAMRPoissonOp::
setCoefs(const RefCountedPtr<LevelData<FArrayBox> >& a_aCoef,
         const RefCountedPtr<LevelData<FluxBox  > >& a_bCoef,
         const Real& a_alpha,
         const Real& a_beta)
{
  m_alpha = a_alpha;
  m_beta  = a_beta;
  m_aCoef = a_aCoef;
  m_bCoef = a_bCoef;
  resetLambda();
}
void
VCAMRPoissonOp::
resetLambda()
{
  Real scale = 1.0 / (m_dx*m_dx);

  // Compute it box by box, point by point
  for (DataIterator dit = m_lambda.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox&       lambdaFab = m_lambda[dit];
      const FArrayBox& aCoefFab  = (*m_aCoef)[dit];
      const FluxBox&   bCoefFab  = (*m_bCoef)[dit];

      const Box& curBox = lambdaFab.box();

      // Compute the diagonal term
      lambdaFab.copy(aCoefFab);
      lambdaFab.mult(m_alpha);

      for (int dir = 0; dir < SpaceDim; dir++)
        {
          FORT_SUMFACES(CHF_FRA(lambdaFab),
                        CHF_CONST_REAL(m_beta),
                        CHF_CONST_FRA(bCoefFab[dir]),
                        CHF_BOX(curBox),
                        CHF_CONST_INT(dir),
                        CHF_CONST_REAL(scale));
        }

      // Take its reciprocal
      lambdaFab.invert(1.0);
    }
}

void VCAMRPoissonOp::define(const DisjointBoxLayout& a_grids,
                            const Real&              a_dx,
                            const ProblemDomain&     a_domain,
                            BCHolder                   a_bc)
{
  CH_TIME("VCAMRPoissonOp::define1");

  m_bc     = a_bc;
  m_domain = a_domain;
  m_dx     = a_dx;
  m_dxCrse = 2*a_dx;

  // redefined in AMRLevelOp<LevelData<FArrayBox> >::define virtual function.
  m_refToCoarser = 2;
  m_refToFiner   = 2;

  m_exchangeCopier.exchangeDefine(a_grids, IntVect::Unit);
  m_exchangeCopier.trimEdges(a_grids, IntVect::Unit);

  for (int i=0; i<CH_SPACEDIM; ++i)
    {
      LayoutData<CFIVS>& lo =  m_loCFIVS[i];
      LayoutData<CFIVS>& hi =  m_hiCFIVS[i];
      lo.define(a_grids);
      hi.define(a_grids);
      for (DataIterator dit(a_grids); dit.ok(); ++dit)
        {
          lo[dit].define(a_domain, a_grids.get(dit),a_grids, i, Side::Lo);
          hi[dit].define(a_domain, a_grids.get(dit),a_grids, i, Side::Hi);
        }
    }
}

void VCAMRPoissonOp::define(const DisjointBoxLayout& a_grids,
                            const DisjointBoxLayout& a_coarse,
                            const Real&              a_dxLevel,
                            int                      a_refRatio,
                            const ProblemDomain&     a_domain,
                            BCHolder                   a_bc)
{
  CH_TIME("VCAMRPoissonOp::define2");

  m_refToCoarser = a_refRatio;
  m_dxCrse = a_refRatio*a_dxLevel;
  m_refToFiner = 1;

  m_bc     = a_bc;
  m_domain = a_domain;
  m_dx     = a_dxLevel;
  m_dxCrse = a_refRatio*a_dxLevel;

  // these get set again after define is called

  m_exchangeCopier.exchangeDefine(a_grids, IntVect::Unit);
  m_exchangeCopier.trimEdges(a_grids, IntVect::Unit);

  for (int i=0; i<CH_SPACEDIM; ++i)
    {
      LayoutData<CFIVS>& lo =  m_loCFIVS[i];
      LayoutData<CFIVS>& hi =  m_hiCFIVS[i];
      lo.define(a_grids);
      hi.define(a_grids);
      for (DataIterator dit(a_grids); dit.ok(); ++dit)
        {
          lo[dit].define(m_domain, a_grids.get(dit), a_grids, i, Side::Lo);
          hi[dit].define(m_domain, a_grids.get(dit), a_grids, i, Side::Hi);
        }
    }
  m_interpWithCoarser.define(a_grids, &a_coarse, a_dxLevel,
                             m_refToCoarser, 1, m_domain);
}

/** full define function for AMRLevelOp with both coarser and finer levels */
void VCAMRPoissonOp::define(const DisjointBoxLayout& a_grids,
                            const DisjointBoxLayout& a_gridsFiner,
                            const DisjointBoxLayout& a_gridsCoarser,
                            const Real&              a_dxLevel,
                            int                      a_refRatio,
                            int                      a_refRatioFiner,
                            const ProblemDomain&     a_domain,
                            BCHolder                   a_bc)
{
  CH_TIME("VCAMRPoissonOp::define3");

  this->define(a_grids, a_gridsCoarser, a_dxLevel, a_refRatio, a_domain, a_bc);
  m_refToFiner = a_refRatioFiner;
}

/** full define function for AMRLevelOp<LevelData<FArrayBox> > with finer
 *  levels, but no coarser */
void VCAMRPoissonOp::define(const DisjointBoxLayout& a_grids,
                            const DisjointBoxLayout& a_gridsFiner,
                            const Real&              a_dxLevel,
                            int                      a_refRatio, // dummy arg
                            int                      a_refRatioFiner,
                            const ProblemDomain&     a_domain,
                            BCHolder                   a_bc)
{
  CH_TIME("VCAMRPoissonOp::define4");

  CH_assert(a_refRatio == 1);

  // calls the MG version of define
  this->define(a_grids, a_dxLevel, a_domain, a_bc);

  m_refToFiner = a_refRatioFiner;
}

// Compute the reciprocal of the diagonal entry of the operator matrix
void VCAMRPoissonOp::computeLambda()
{
  CH_TIME("VCAMRPoissonOp::computeLambda");

  CH_assert(!m_lambda.isDefined());

  // Define lambda
  m_lambda.define(m_aCoef->disjointBoxLayout(),m_aCoef->nComp());
  resetLambda();
}

void VCAMRPoissonOp::residual(LevelData<FArrayBox>&       a_lhs,
                              const LevelData<FArrayBox>& a_phi,
                              const LevelData<FArrayBox>& a_rhs,
                              bool                        a_homogeneous)
{
  CH_TIME("VCAMRPoissonOp::residual");

  applyOp(a_lhs, a_phi, a_homogeneous);
  incr(a_lhs, a_rhs, -1);
  scale(a_lhs, -1.0);
}

/**************************/
// this preconditioner first initializes phihat to (IA)phihat = rhshat
// (diagonization of L -- A is the matrix version of L)
// then smooths with a couple of passes of levelGSRB
void VCAMRPoissonOp::preCond(LevelData<FArrayBox>&       a_phi,
                             const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("VCAMRPoissonOp::preCond");

  // diagonal term of this operator in:
  //
  //       alpha * a(i)
  //     + beta  * sum_over_dir (b(i-1/2*e_dir) + b(i+1/2*e_dir)) / (dx*dx)
  //
  // The inverse of this is our initial multiplier.

  int ncomp = a_phi.nComp();

  CH_assert(m_lambda.isDefined());
  CH_assert(a_rhs.nComp()    == ncomp);
  CH_assert(m_bCoef->nComp() == ncomp);

  // don't need to use a Copier -- plain copy will do
  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      // also need to average and sum face-centered bCoefs to cell-centers
      Box gridBox = a_rhs[dit].box();

      // approximate inverse
      a_phi[dit].copy(a_rhs[dit]);
      a_phi[dit].mult(m_lambda[dit], gridBox, 0, 0, ncomp);
    }

  relax(a_phi, a_rhs, 2);
}

void VCAMRPoissonOp::applyOp(LevelData<FArrayBox>&       a_lhs,
                             const LevelData<FArrayBox>& a_phi,
                             bool                        a_homogeneous )
{
  CH_TIME("VCAMRPoissonOp::applyOp");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  Real dx = m_dx;
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_bc(phi[dit], dbl[dit()],m_domain, dx, a_homogeneous);
    }

  phi.exchange(phi.interval(), m_exchangeCopier);

  computeOperatorNoBCs(a_lhs, phi);
}

void VCAMRPoissonOp::applyOpNoBoundary(LevelData<FArrayBox>&       a_lhs,
                             const LevelData<FArrayBox>& a_phi)
{
  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  phi.exchange(phi.interval(), m_exchangeCopier);
  computeOperatorNoBCs(a_lhs, a_phi);
}

void VCAMRPoissonOp::create(LevelData<FArrayBox>&       a_lhs,
                            const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("VCAMRPoissonOp::create");

  m_levelOps.create(a_lhs, a_rhs);
}

void VCAMRPoissonOp::createCoarsened(LevelData<FArrayBox>&       a_lhs,
                                     const LevelData<FArrayBox>& a_rhs,
                                     const int&                  a_refRat)
{
  CH_TIME("VCAMRPoissonOp::createCoarsened");

  int ncomp = a_rhs.nComp();
  IntVect ghostVect = a_rhs.ghostVect();

  DisjointBoxLayout dbl = a_rhs.disjointBoxLayout();
  CH_assert(dbl.coarsenable(a_refRat));

  // fill ebislayout
  DisjointBoxLayout dblCoarsenedFine;
  coarsen(dblCoarsenedFine, dbl, a_refRat);

  a_lhs.define(dblCoarsenedFine, ncomp, a_rhs.ghostVect());
}

void VCAMRPoissonOp::assign(LevelData<FArrayBox>&       a_lhs,
                            const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("VCAMRPoissonOp::assign");

  m_levelOps.assign(a_lhs, a_rhs);
}

Real VCAMRPoissonOp::dotProduct(const LevelData<FArrayBox>& a_1,
                                const LevelData<FArrayBox>& a_2)
{
  CH_TIME("VCAMRPoissonOp::dotProduct");

  return m_levelOps.dotProduct(a_1, a_2);
}

void VCAMRPoissonOp::incr(LevelData<FArrayBox>&       a_lhs,
                          const LevelData<FArrayBox>& a_x,
                          Real                        a_scale)
{
  CH_TIME("VCAMRPoissonOp::incr");

  m_levelOps.incr(a_lhs, a_x, a_scale);
}

void VCAMRPoissonOp::axby(LevelData<FArrayBox>&       a_lhs,
                          const LevelData<FArrayBox>& a_x,
                          const LevelData<FArrayBox>& a_y,
                          Real                        a_a,
                          Real                        a_b)
{
  CH_TIME("VCAMRPoissonOp::axby");

  m_levelOps.axby(a_lhs, a_x, a_y, a_a, a_b);
}

void VCAMRPoissonOp::scale(LevelData<FArrayBox>& a_lhs,
                           const Real&           a_scale)
{
  CH_TIME("VCAMRPoissonOp::scale");

  m_levelOps.scale(a_lhs, a_scale);
}

Real VCAMRPoissonOp::norm(const LevelData<FArrayBox>& a_x,
                          int                         a_ord)
{
  CH_TIME("VCAMRPoissonOp::norm");

  return CH_XD::norm(a_x, a_x.interval(), a_ord);
}

void VCAMRPoissonOp::setToZero(LevelData<FArrayBox>& a_lhs)
{
  CH_TIME("VCAMRPoissonOp::setToZero");

  m_levelOps.setToZero(a_lhs);
}

void VCAMRPoissonOp::relax(LevelData<FArrayBox>&       a_e,
                           const LevelData<FArrayBox>& a_residual,
                           int                         a_iterations)
{
  CH_TIME("VCAMRPoissonOp::relax");

  for (int i=0; i<a_iterations; i++)
    {
      levelGSRB(a_e, a_residual);
      // levelJacobi(a_e, a_residual);
    }
}

void VCAMRPoissonOp::createCoarser(LevelData<FArrayBox>&       a_coarse,
                                   const LevelData<FArrayBox>& a_fine,
                                   bool                        a_ghosted)
{
  CH_TIME("VCAMRPoissonOp::createCoarser");

  // CH_assert(!a_ghosted);
  IntVect ghost = a_fine.ghostVect();
  DisjointBoxLayout dbl;
  CH_assert(dbl.coarsenable(2));
  coarsen(dbl, a_fine.disjointBoxLayout(), 2); // multigrid, so coarsen by 2
  a_coarse.define(dbl, a_fine.nComp(), ghost);
}

void VCAMRPoissonOp::restrictResidual(LevelData<FArrayBox>&       a_resCoarse,
                                      LevelData<FArrayBox>&       a_phiFine,
                                      const LevelData<FArrayBox>& a_rhsFine)
{
  CH_TIME("VCAMRPoissonOp::restrictResidual");

  homogeneousCFInterp(a_phiFine);
  const DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& phi = a_phiFine[dit];
      m_bc(phi, dblFine[dit()], m_domain, m_dx, true);
    }

  a_phiFine.exchange(a_phiFine.interval(), m_exchangeCopier);

  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox&       phi = a_phiFine[dit];
      const FArrayBox& rhs = a_rhsFine[dit];
      FArrayBox&       res = a_resCoarse[dit];

      const FArrayBox& thisACoef = (*m_aCoef)[dit];
      const FluxBox&   thisBCoef = (*m_bCoef)[dit];

      Box region = dblFine.get(dit());
      const IntVect& iv = region.smallEnd();
      IntVect civ = coarsen(iv, 2);

      res.setVal(0.0);

#if CH_SPACEDIM == 1
      FORT_RESTRICTRESVC1D(
#elif CH_SPACEDIM == 2
                           FORT_RESTRICTRESVC2D(
#elif CH_SPACEDIM == 3
                                                FORT_RESTRICTRESVC3D(
#else
                                                                     This_will_not_compile))
#endif
                           CHF_FRA_SHIFT(res, civ),
                           CHF_CONST_FRA_SHIFT(phi, iv),
                           CHF_CONST_FRA_SHIFT(rhs, iv),
                           CHF_CONST_REAL(m_alpha),
                           CHF_CONST_FRA_SHIFT(thisACoef, iv),
                           CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
                           CHF_CONST_FRA_SHIFT(thisBCoef[0], iv),
#endif
#if CH_SPACEDIM >= 2
                           CHF_CONST_FRA_SHIFT(thisBCoef[1], iv),
#endif
#if CH_SPACEDIM >= 3
                           CHF_CONST_FRA_SHIFT(thisBCoef[2], iv),
#endif
#if CH_SPACEDIM >= 4
                           This_will_not_compile!
#endif
                           CHF_BOX_SHIFT(region, iv),
                           CHF_CONST_REAL(m_dx)
                           );
    }
}

void VCAMRPoissonOp::prolongIncrement(LevelData<FArrayBox>&       a_phiThisLevel,
                                      const LevelData<FArrayBox>& a_correctCoarse)
{
  CH_TIME("VCAMRPoissonOp::prolongIncrement");

  DisjointBoxLayout dbl = a_phiThisLevel.disjointBoxLayout();
  int mgref = 2; //this is a multigrid func

  for (DataIterator dit = a_phiThisLevel.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& phi =  a_phiThisLevel[dit];
      const FArrayBox& coarse = a_correctCoarse[dit];
      Box region = dbl.get(dit());
      const IntVect& iv = region.smallEnd();
      IntVect civ=coarsen(iv, 2);

      FORT_PROLONG(CHF_FRA_SHIFT(phi, iv),
                   CHF_CONST_FRA_SHIFT(coarse, civ),
                   CHF_BOX_SHIFT(region, iv),
                   CHF_CONST_INT(mgref));

    }
}

void VCAMRPoissonOp::AMRResidual(LevelData<FArrayBox>&              a_residual,
                                 const LevelData<FArrayBox>&        a_phiFine,
                                 const LevelData<FArrayBox>&        a_phi,
                                 const LevelData<FArrayBox>&        a_phiCoarse,
                                 const LevelData<FArrayBox>&        a_rhs,
                                 bool                               a_homogeneousPhysBC,
                                 AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIME("VCAMRPoissonOp::AMRResidual");

  AMROperator(a_residual, a_phiFine, a_phi, a_phiCoarse,
              a_homogeneousPhysBC, a_finerOp);

  incr(a_residual, a_rhs, -1.0);
  scale(a_residual, -1.0);
}

void VCAMRPoissonOp::AMRResidualNC(LevelData<FArrayBox>&              a_residual,
                                   const LevelData<FArrayBox>&        a_phiFine,
                                   const LevelData<FArrayBox>&        a_phi,
                                   const LevelData<FArrayBox>&        a_rhs,
                                   bool                               a_homogeneousPhysBC,
                                   AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIME("VCAMRPoissonOp::AMRResidualNC");

  AMROperatorNC(a_residual, a_phiFine, a_phi,
                a_homogeneousPhysBC, a_finerOp);
  axby(a_residual, a_residual, a_rhs, -1.0, 1.0);
}

void VCAMRPoissonOp::AMRResidualNF(LevelData<FArrayBox>&       a_residual,
                                   const LevelData<FArrayBox>& a_phi,
                                   const LevelData<FArrayBox>& a_phiCoarse,
                                   const LevelData<FArrayBox>& a_rhs,
                                   bool                        a_homogeneousPhysBC)
{
  CH_TIME("VCAMRPoissonOp::AMRResidualNF");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  if (a_phiCoarse.isDefined())
    {
      m_interpWithCoarser.coarseFineInterp(phi, a_phiCoarse);
    }

  // apply boundary conditions
  this->residual(a_residual, a_phi, a_rhs, a_homogeneousPhysBC );
}

void VCAMRPoissonOp::AMROperator(LevelData<FArrayBox>&              a_LofPhi,
                                 const LevelData<FArrayBox>&        a_phiFine,
                                 const LevelData<FArrayBox>&        a_phi,
                                 const LevelData<FArrayBox>&        a_phiCoarse,
                                 bool                               a_homogeneousPhysBC,
                                 AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIME("VCAMRPoissonOp::AMROperator");

  CH_assert(a_phiFine  .isDefined());
  CH_assert(a_phi      .isDefined());
  CH_assert(a_phiCoarse.isDefined());
  CH_assert(a_finerOp != NULL);

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  m_interpWithCoarser.coarseFineInterp(phi, a_phiCoarse);
  applyOp(a_LofPhi, a_phi, a_homogeneousPhysBC);
  reflux(a_phiFine, a_phi,  a_LofPhi, a_finerOp);
}

void VCAMRPoissonOp::AMROperatorNC(LevelData<FArrayBox>&              a_LofPhi,
                                   const LevelData<FArrayBox>&        a_phiFine,
                                   const LevelData<FArrayBox>&        a_phi,
                                   bool                               a_homogeneousPhysBC,
                                   AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIME("VCAMRPoissonOp::AMROperatorNC");

  CH_assert(a_phiFine  .isDefined());
  CH_assert(a_phi      .isDefined());
  CH_assert(a_finerOp != NULL);

  // no coarse-fine interpolation here
  applyOp(a_LofPhi, a_phi, a_homogeneousPhysBC);
  reflux(a_phiFine, a_phi,  a_LofPhi, a_finerOp);
}

void VCAMRPoissonOp::AMROperatorNF(LevelData<FArrayBox>&       a_LofPhi,
                                   const LevelData<FArrayBox>& a_phi,
                                   const LevelData<FArrayBox>& a_phiCoarse,
                                   bool                        a_homogeneousPhysBC)
{
  CH_TIME("VCAMRPoissonOp::AMROperatorNF");

  CH_assert(a_phi      .isDefined());
  CH_assert(a_phiCoarse.isDefined());

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  if (a_phiCoarse.isDefined())
    {
      m_interpWithCoarser.coarseFineInterp(phi, a_phiCoarse);
    }
  // apply boundary conditions in applyOp
  this->applyOp(a_LofPhi, a_phi, a_homogeneousPhysBC );
}

void VCAMRPoissonOp::AMRRestrict(LevelData<FArrayBox>&       a_resCoarse,
                                 const LevelData<FArrayBox>& a_residual,
                                 const LevelData<FArrayBox>& a_correction,
                                 const LevelData<FArrayBox>& a_coarseCorrection)
{
  CH_TIME("VCAMRPoissonOp::AMRRestrict");

  LevelData<FArrayBox> r;
  create(r, a_residual);
  AMRResidualNF(r, a_correction, a_coarseCorrection, a_residual, true);
  DisjointBoxLayout dblCoar = a_resCoarse.disjointBoxLayout();
  DataIterator dit = a_residual.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& coarse = a_resCoarse[dit];
      const FArrayBox& fine = r[dit];
      const Box& b = dblCoar.get(dit());
      Box refbox(IntVect::Zero,
                 (m_refToCoarser-1)*IntVect::Unit);
      FORT_AVERAGE(CHF_FRA(coarse),
                   CHF_CONST_FRA(fine),
                   CHF_BOX(b),
                   CHF_CONST_INT(m_refToCoarser),
                   CHF_BOX(refbox));
    }
}

/** a_correction += I[2h->h](a_coarseCorrection) */
void VCAMRPoissonOp::AMRProlong(LevelData<FArrayBox>&       a_correction,
                                const LevelData<FArrayBox>& a_coarseCorrection)
{
  CH_TIME("VCAMRPoissonOp::AMRProlong");

  DisjointBoxLayout c;
  coarsen(c,  a_correction.disjointBoxLayout(), m_refToCoarser);
  LevelData<FArrayBox> eCoar(c, a_correction.nComp(),a_coarseCorrection.ghostVect());
  a_coarseCorrection.copyTo(eCoar.interval(), eCoar, eCoar.interval());

  DisjointBoxLayout dbl = a_correction.disjointBoxLayout();
  for (DataIterator dit = a_correction.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& phi =  a_correction[dit];
      const FArrayBox& coarse = eCoar[dit];

      Box region = dbl.get(dit());
      const IntVect& iv = region.smallEnd();
      IntVect civ = coarsen(iv, m_refToCoarser);

      FORT_PROLONG(CHF_FRA_SHIFT(phi, iv),
                   CHF_CONST_FRA_SHIFT(coarse, civ),
                   CHF_BOX_SHIFT(region, iv),
                   CHF_CONST_INT(m_refToCoarser));
    }
}

void VCAMRPoissonOp::AMRUpdateResidual(LevelData<FArrayBox>&       a_residual,
                                       const LevelData<FArrayBox>& a_correction,
                                       const LevelData<FArrayBox>& a_coarseCorrection)
{
  CH_TIME("VCAMRPoissonOp::AMRUpdateResidual");

  LevelData<FArrayBox> r;
  this->create(r, a_residual);
  this->AMRResidualNF(r, a_correction, a_coarseCorrection, a_residual, true);
  this->assign(a_residual, r);
}

// compute norm over all cells on coarse not covered by finer
Real VCAMRPoissonOp::AMRNorm(const LevelData<FArrayBox>& a_coarResid,
                             const LevelData<FArrayBox>& a_fineResid,
                             const int&                  a_refRat,
                             const int&                  a_ord)
{
  CH_TIME("VCAMRPoissonOp::AMRNorm");

  const DisjointBoxLayout& coarGrids = a_coarResid.disjointBoxLayout();

  // create temp and zero out under finer grids
  LevelData<FArrayBox> coarTemp;
  m_levelOps.create(coarTemp, a_coarResid);
  m_levelOps.assign(coarTemp, a_coarResid);

  if (a_fineResid.isDefined())
    {
      const DisjointBoxLayout& fineGrids = a_fineResid.disjointBoxLayout();
      int ncomp = coarTemp.nComp();
      for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
        {
          FArrayBox& coarTempFAB = coarTemp[dit()];
          LayoutIterator litFine = fineGrids.layoutIterator();
          for (litFine.reset(); litFine.ok(); ++litFine)
            {
              Box overlayBox = coarTempFAB.box();
              Box coarsenedGrid = coarsen(fineGrids[litFine()], a_refRat);

              overlayBox &= coarsenedGrid;

              if (!overlayBox.isEmpty())
                {
                  coarTempFAB.setVal(0.0,overlayBox,0, ncomp);
                }
            }
        }
    } // end if there is a finer level

  // return norm of temp
  return norm(coarTemp, a_ord);
}

void VCAMRPoissonOp::levelGSRB(LevelData<FArrayBox>&       a_phi,
                               const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("VCAMRPoissonOp::levelGSRB");

  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());
  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();

  DataIterator dit = a_phi.dataIterator();

  // do first red, then black passes
  for (int whichPass = 0; whichPass <= 1; whichPass++)
    {
      CH_TIMERS("VCAMRPoissonOp::levelGSRB::Compute");

      CH_TIMER("VCAMRPoissonOp::levelGSRB::homogeneousCFInterp",timeCFInterp);
      CH_TIMER("VCAMRPoissonOp::levelGSRB::exchange",           timeExchange);
      CH_TIMER("VCAMRPoissonOp::levelGSRB::BCs",                timeBCs);

      CH_START(timeCFInterp);
      homogeneousCFInterp(a_phi);
      CH_STOP(timeCFInterp);

      // fill in intersection of ghostcells and a_phi's boxes
      CH_START(timeExchange);
      a_phi.exchange(a_phi.interval(), m_exchangeCopier);
      CH_STOP(timeExchange);

      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& phiFab = a_phi[dit];
          const Box& region = dbl[dit()];

          const FArrayBox& thisACoef  = (*m_aCoef)[dit];
          const FluxBox&   thisBCoef  = (*m_bCoef)[dit];

          const FArrayBox& thisLambda = m_lambda[dit];

          // invoke physical BC's where necessary
          CH_START(timeBCs);
          m_bc(phiFab, region, m_domain, m_dx, true);
          CH_STOP(timeBCs);

#if CH_SPACEDIM == 1
       FORT_GSRBHELMHOLTZVC1D(
#elif CH_SPACEDIM == 2
       FORT_GSRBHELMHOLTZVC2D(
#elif CH_SPACEDIM == 3
       FORT_GSRBHELMHOLTZVC3D(
#else
                              (This_will_not_compile!)))
#endif
       CHF_FRA(phiFab),
       CHF_CONST_FRA(a_rhs[dit]),
       CHF_BOX(region),
       CHF_CONST_REAL(m_dx),
       CHF_CONST_REAL(m_alpha),
       CHF_CONST_FRA(thisACoef),
       CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
       CHF_CONST_FRA(thisBCoef[0]),
#endif
#if CH_SPACEDIM >= 2
       CHF_CONST_FRA(thisBCoef[1]),
#endif
#if CH_SPACEDIM >= 3
       CHF_CONST_FRA(thisBCoef[2]),
#endif
#if CH_SPACEDIM >= 4
       This_will_not_compile!
#endif
       CHF_CONST_FRA(thisLambda),
       CHF_CONST_INT(whichPass)
       );
        }
    } // end loop through red-black
}

void VCAMRPoissonOp::levelJacobi(LevelData<FArrayBox>&       a_phi,
                                 const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("VCAMRPoissonOp::levelJacobi");

  LevelData<FArrayBox> resid;
  create(resid, a_rhs);
  // Get the residual
  residual(resid,a_phi,a_rhs,true);

  DataIterator dit = m_lambda.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      resid[dit].mult(m_lambda[dit]);
    } // end loop over grids

  // Do the Jacobi relaxation
  incr(a_phi, resid, 0.5);
}

void VCAMRPoissonOp::homogeneousCFInterp(LevelData<FArrayBox>& a_phif)
{
  CH_TIME("VCAMRPoissonOp::homogeneousCFInterp1");

  CH_assert( a_phif.ghostVect() >= IntVect::Unit);

  DataIterator dit = a_phif.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const DataIndex& datInd = dit();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              homogeneousCFInterp(a_phif,datInd,idir,sit());
            }
        }
    }
}

void VCAMRPoissonOp::homogeneousCFInterp(LevelData<FArrayBox>& a_phif,
                                         const DataIndex&      a_datInd,
                                         int                   a_idir,
                                         Side::LoHiSide        a_hiorlo)
{
  CH_TIME("VCAMRPoissonOp::homogeneousCFInterp2");

  CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
  CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
  CH_assert( a_phif.ghostVect() >= IntVect::Unit);
  //  CH_assert (m_ncomp == a_phif.nComp());

  const CFIVS* cfivs_ptr = NULL;
  if (a_hiorlo == Side::Lo)
    cfivs_ptr = &m_loCFIVS[a_idir][a_datInd];
  else
    cfivs_ptr = &m_hiCFIVS[a_idir][a_datInd];

  if (cfivs_ptr->isPacked())
    {
      int ihiorlo = sign(a_hiorlo);
      FORT_INTERPHOMO(CHF_FRA(a_phif[a_datInd]),
                      CHF_BOX(cfivs_ptr->packedBox()),
                      CHF_CONST_REAL(m_dx),
                      CHF_CONST_REAL(m_dxCrse),
                      CHF_CONST_INT(a_idir),
                      CHF_CONST_INT(ihiorlo));
    }
  else
    {
      const IntVectSet& interp_ivs = cfivs_ptr->getFineIVS();
      if (!interp_ivs.isEmpty())
        {
          // Assuming homogenous, interpolate on fine ivs
          interpOnIVSHomo(a_phif, a_datInd, a_idir,
                          a_hiorlo, interp_ivs);
        }
    }
}

void VCAMRPoissonOp::interpOnIVSHomo(LevelData<FArrayBox>& a_phif,
                                     const DataIndex&      a_datInd,
                                     const int             a_idir,
                                     const Side::LoHiSide  a_hiorlo,
                                     const IntVectSet&     a_interpIVS)
{
  CH_TIME("VCAMRPoissonOp::interpOnIVSHomo");

  CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
  CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
  CH_assert( a_phif.ghostVect() >= IntVect::Unit);
  IVSIterator fine_ivsit(a_interpIVS);
  FArrayBox& a_phi = a_phif[a_datInd];

  // much of these scalar values can be precomputed and stored if
  // we ever need to speed-up this function (ndk)
  Real x1 = m_dx;
  Real x2 = 0.5*(3. * x1 + m_dxCrse);
  Real denom = 1.0-((x1+x2)/x1);
  Real idenom = 1/(denom); // divide is more expensive usually
  Real x = 2.*x1;
  Real xsquared = x*x;

  Real m1 = 1/(x1*x1);
  Real m2 = 1/(x1*(x1-x2));

  Real q1 = 1/(x1-x2);
  Real q2 = x1+x2;

  int ihilo = sign(a_hiorlo);
  Real pa,pb,a,b;
  IntVect ivf;
  for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit)
    {
      ivf = fine_ivsit();
      // quadratic interpolation
      for (int ivar = 0; ivar < a_phif.nComp(); ivar++)
        {
          ivf[a_idir]-=2*ihilo;
          pa = a_phi(ivf, ivar);
          ivf[a_idir]+=ihilo;
          pb = a_phi(ivf, ivar);

          a = ((pb-pa)*m1 - (pb)*m2)*idenom;
          b = (pb)*q1 - a*q2;

          ivf[a_idir]+=ihilo;
          a_phi(fine_ivsit(), ivar) = a*xsquared + b*x + pa;
        } // end loop over components
    } // end loop over fine intvects
}

void VCAMRPoissonOp::getFlux(FArrayBox&       a_flux,
                             const FArrayBox& a_data,
                             const FluxBox&   a_bCoef,
                             const Box&       a_facebox,
                             int              a_dir,
                             int              a_ref) const
{
  CH_TIME("VCAMRPoissonOp::getFlux");

  CH_assert(a_dir >= 0);
  CH_assert(a_dir <  SpaceDim);
  CH_assert(!a_data.box().isEmpty());
  CH_assert(!a_facebox.isEmpty());

  // probably the simplest way to test centering
  // a_box needs to be face-centered in the a_dir
  Box faceTestBox(IntVect::Zero, IntVect::Unit);
  faceTestBox.surroundingNodes(a_dir);
  CH_assert(a_facebox.type() == faceTestBox.type());

  const FArrayBox& bCoefDir = a_bCoef[a_dir];

  // reality check for bCoef
  CH_assert(bCoefDir.box().contains(a_facebox));

  a_flux.resize(a_facebox, a_data.nComp());
  BoxIterator bit(a_facebox);

  Real scale = m_beta * a_ref / m_dx;

  for ( bit.begin(); bit.ok(); bit.next())
    {
      IntVect iv = bit();
      IntVect shiftiv = BASISV(a_dir);
      IntVect ivlo = iv - shiftiv;
      IntVect ivhi = iv;

      CH_assert(a_data.box().contains(ivlo));
      CH_assert(a_data.box().contains(ivhi));

      for (int ivar = 0; ivar < a_data.nComp(); ivar++)
        {
          Real phihi = a_data(ivhi,ivar);
          Real philo = a_data(ivlo,ivar);
          Real gradphi = (phihi - philo ) * scale;

          a_flux(iv,ivar) = -bCoefDir(iv, ivar) * gradphi;
        }
    }
}

void VCAMRPoissonOp::reflux(const LevelData<FArrayBox>&        a_phiFine,
                            const LevelData<FArrayBox>&        a_phi,
                            LevelData<FArrayBox>&              a_residual,
                            AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIME("VCAMRPoissonOp::reflux");

  int ncomp = 1;
  ProblemDomain fineDomain = refine(m_domain, m_refToFiner);
  LevelFluxRegister levfluxreg(a_phiFine.disjointBoxLayout(),
                               a_phi.disjointBoxLayout(),
                               fineDomain,
                               m_refToFiner,
                               ncomp);

  levfluxreg.setToZero();
  Interval interv(0,a_phi.nComp()-1);

  DataIterator dit = a_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      const FArrayBox& coarfab = a_phi[dit];
      const FluxBox& coarBCoef  = (*m_bCoef)[dit];
      const Box& gridBox = a_phi.getBoxes()[dit];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          FArrayBox coarflux;
          Box faceBox = surroundingNodes(gridBox, idir);
          getFlux(coarflux, coarfab, coarBCoef , faceBox, idir);

          Real scale = 1.0;
          levfluxreg.incrementCoarse(coarflux, scale,dit(),
                                     interv,interv,idir);
        }
    }
  LevelData<FArrayBox>& p = ( LevelData<FArrayBox>&)a_phiFine;

  // has to be its own object because the finer operator
  // owns an interpolator and we have no way of getting to it
  VCAMRPoissonOp* finerAMRPOp = (VCAMRPoissonOp*) a_finerOp;
  QuadCFInterp& quadCFI = finerAMRPOp->m_interpWithCoarser;

  quadCFI.coarseFineInterp(p, a_phi);
  // p.exchange(a_phiFine.interval()); // BVS is pretty sure this is not necesary.
  IntVect phiGhost = p.ghostVect();

  DataIterator ditf = a_phiFine.dataIterator();
  const  DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (ditf.reset(); ditf.ok(); ++ditf)
    {
      const FArrayBox& phifFab = a_phiFine[ditf];
      const FluxBox& fineBCoef  = (*(finerAMRPOp->m_bCoef))[ditf];
      const Box& gridbox = dblFine.get(ditf());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          int normalGhost = phiGhost[idir];
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              Side::LoHiSide hiorlo = sit();
              Box fabbox;
              Box facebox;

              // assumption here that the stencil required
              // to compute the flux in the normal direction
              // is 2* the number of ghost cells for phi
              // (which is a reasonable assumption, and probably
              // better than just assuming you need one cell on
              // either side of the interface
              // (dfm 8-4-06)
              if (sit() == Side::Lo)
                {
                  fabbox = adjCellLo(gridbox,idir, 2*normalGhost);
                  fabbox.shift(idir, 1);
                  facebox = bdryLo(gridbox, idir,1);
                }
              else
                {
                  fabbox = adjCellHi(gridbox,idir, 2*normalGhost);
                  fabbox.shift(idir, -1);
                  facebox = bdryHi(gridbox, idir, 1);
                }

              // just in case we need ghost cells in the transverse direction
              // (dfm 8-4-06)
              for (int otherDir=0; otherDir<SpaceDim; ++otherDir)
                {
                  if (otherDir != idir)
                    {
                      fabbox.grow(otherDir, phiGhost[otherDir]);
                    }
                }
              CH_assert(!fabbox.isEmpty());

              FArrayBox phifab(fabbox, a_phi.nComp());
              phifab.copy(phifFab);

              FArrayBox fineflux;
              getFlux(fineflux, phifab, fineBCoef, facebox, idir,
                      m_refToFiner);

              Real scale = 1.0;
              levfluxreg.incrementFine(fineflux, scale, ditf(),
                                       interv, interv, idir, hiorlo);
            }
        }
    }

  Real scale =  1.0/m_dx;
  levfluxreg.reflux(a_residual, scale);
}

// computes VC operator after BC's have been set
void VCAMRPoissonOp::computeOperatorNoBCs(LevelData<FArrayBox>&       a_lhs,
                                          const LevelData<FArrayBox>& a_phi)
{
  CH_TIME("VCAMRPoissonOp::computeOperatorNoBCs");

  const DisjointBoxLayout dbl= a_phi.getBoxes();
  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisLHS = a_lhs[dit];
      const Box& region = dbl[dit()];

      const FArrayBox& thisACoef = (*m_aCoef)[dit];
      const FluxBox&   thisBCoef = (*m_bCoef)[dit];

#if CH_SPACEDIM == 1
      FORT_VCCOMPUTEOP1D(
#elif CH_SPACEDIM == 2
      FORT_VCCOMPUTEOP2D(
#elif CH_SPACEDIM == 3
      FORT_VCCOMPUTEOP3D(
#else
                         (This_will_not_compile!)))
#endif
      CHF_FRA(thisLHS),
      CHF_CONST_FRA(a_phi[dit]),
      CHF_CONST_REAL(m_alpha),
      CHF_CONST_FRA(thisACoef),
      CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
      CHF_CONST_FRA(thisBCoef[0]),
#endif
#if CH_SPACEDIM >= 2
      CHF_CONST_FRA(thisBCoef[1]),
#endif
#if CH_SPACEDIM >= 3
      CHF_CONST_FRA(thisBCoef[2]),
#endif
#if CH_SPACEDIM >= 4
      This_will_not_compile!
#endif
      CHF_BOX(region),
      CHF_CONST_REAL(m_dx)
      );
    } // end loop over boxes
}

void VCAMRPoissonOp::write(const LevelData<FArrayBox>* a_data,
                           const char*                 a_filename)
{
#ifdef CH_USE_HDF5
  writeLevelname(a_data, a_filename);
#else
  MayDay::Warning("VCAMRPoissonOp::write unimplemented since CH_USE_HDF5 undefined");
#endif
}

// Factory
VCAMRPoissonOpFactory::VCAMRPoissonOpFactory()
{
  setDefaultValues();
}

//  AMR Factory define function
void VCAMRPoissonOpFactory::define(
                                   const ProblemDomain&                           a_coarseDomain,
                                   const Vector<DisjointBoxLayout>&               a_grids,
                                   const Vector<int>&                             a_refRatios,
                                   const Real&                                    a_coarsedx,
                                   BCHolder                                         a_bc,
                                   const Real&                                    a_alpha,
                                   Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                                   const Real&                                    a_beta,
                                   Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_bCoef)
{
  setDefaultValues();

  m_domains.resize(a_grids.size());
  m_boxes=a_grids;
  m_refRatios=a_refRatios;
  m_dx.resize(a_grids.size());
  m_bc = a_bc;
  m_domains[0] = a_coarseDomain;
  m_dx[0] = a_coarsedx;
  for (int i=1; i<a_grids.size(); i++)
    {
      m_dx[i] = m_dx[i-1]/m_refRatios[i-1];
      m_domains[i] = m_domains[i-1];
      m_domains[i].refine(m_refRatios[i-1]);
    }

  m_alpha = a_alpha;
  m_aCoef = a_aCoef;

  m_beta  = a_beta;
  m_bCoef = a_bCoef;
}

VCAMRPoissonOp* VCAMRPoissonOpFactory::MGnewOp(const ProblemDomain& a_indexSpace,
                                               int                  a_depth,
                                               bool                 a_homoOnly)
{
  // CH_assert(m_boxes.size()>a_depth);

  int lev;
  for (lev = 0; lev < m_domains.size(); lev++)
    {
      if (a_indexSpace.domainBox() == m_domains[lev].domainBox()) break;
    }
  CH_assert(lev !=  m_domains.size()); // didn't find domain

  DisjointBoxLayout layout(m_boxes[lev]);
  ProblemDomain domain(m_domains[lev]);
  Real dx = m_dx[lev];
  int coarsenRatio = 1;

  for (int i = 0; i < a_depth; i++)
    {
      if (!layout.coarsenable(4))
        {
          return NULL;
        }
      DisjointBoxLayout dbl;
      coarsen_dbl(dbl, layout, 2);
      layout = dbl;
      dx*=2;
      coarsenRatio *= 2;
      domain.coarsen(2);
    }

  VCAMRPoissonOp* newOp = new VCAMRPoissonOp;
  newOp->define(layout, dx, domain, m_bc);

  if (a_depth == 0)
    {
      newOp->m_alpha = m_alpha;
      newOp->m_beta  = m_beta;

      // don't need to coarsen anything for this
      newOp->m_aCoef = m_aCoef[lev];
      newOp->m_bCoef = m_bCoef[lev];
    }
  else
    {
      // need to coarsen coefficients
      RefCountedPtr<LevelData<FArrayBox> > aCoef( new LevelData<FArrayBox> );
      RefCountedPtr<LevelData<FluxBox> > bCoef( new LevelData<FluxBox> );
      aCoef->define(layout, m_aCoef[lev]->nComp(), m_aCoef[lev]->ghostVect());
      bCoef->define(layout, m_bCoef[lev]->nComp(), m_bCoef[lev]->ghostVect());

      // average coefficients to coarser level
      // for now, do this with a CoarseAverage --
      // may want to switch to harmonic averaging at some point
      CoarseAverage averager(m_aCoef[lev]->getBoxes(),
                             layout, aCoef->nComp(), coarsenRatio);

      CoarseAverageFace faceAverager(m_bCoef[lev]->getBoxes(),
                                     bCoef->nComp(), coarsenRatio);

      if (m_coefficient_average_type == CoarseAverage::arithmetic)
        {
          averager.averageToCoarse(*aCoef, *(m_aCoef[lev]));
          faceAverager.averageToCoarse(*bCoef, *(m_bCoef[lev]));
        }
      else if (m_coefficient_average_type == CoarseAverage::harmonic)
        {
          averager.averageToCoarseHarmonic(*aCoef, *(m_aCoef[lev]));
          faceAverager.averageToCoarseHarmonic(*bCoef, *(m_bCoef[lev]));
        }
      else
        {
          MayDay::Error("VCAMRPoissonOpFactory::MGNewOp -- bad averagetype");
        }

      newOp->m_alpha = m_alpha;
      newOp->m_beta  = m_beta;

      newOp->m_aCoef = aCoef;
      newOp->m_bCoef = bCoef;
    }

  newOp->computeLambda();

  Real dxCrse = -1;
  if (lev > 0)
    {
      dxCrse = m_dx[lev-1];
    }
  newOp->m_dxCrse  = dxCrse;

  return newOp;
}

VCAMRPoissonOp* VCAMRPoissonOpFactory::AMRnewOp(const ProblemDomain& a_indexSpace)
{
  VCAMRPoissonOp* newOp = new VCAMRPoissonOp;
  int ref = 0;
  for (;ref< m_domains.size(); ref++)
    {
      if (a_indexSpace.domainBox() == m_domains[ref].domainBox()) break;
    }

  if (ref == 0)
    {
      // coarsest AMR level
      if (m_domains.size() == 1)
        {
          // no finer level
          newOp->define(m_boxes[0],  m_dx[0],
                        a_indexSpace, m_bc);
        }
      else
        {
          // finer level exists but no coarser
          int dummyRat = 1;  // argument so compiler can find right function
          int refToFiner = m_refRatios[0]; // actual refinement ratio
          newOp->define(m_boxes[0],  m_boxes[1], m_dx[0],
                        dummyRat, refToFiner,
                        a_indexSpace, m_bc);
        }
    }
  else if (ref ==  m_domains.size()-1)
    {
      // finest AMR level
      newOp->define(m_boxes[ref], m_boxes[ref-1], m_dx[ref],
                    m_refRatios[ref-1],
                    a_indexSpace, m_bc);
    }
  else if ( ref == m_domains.size())
    {
      MayDay::Error("Did not find a domain to match AMRnewOp(const ProblemDomain& a_indexSpace)");
    }
  else
    {
      // intermediate AMR level, full define
      newOp->define(m_boxes[ref], m_boxes[ref+1], m_boxes[ref-1], m_dx[ref],
                    m_refRatios[ref-1], m_refRatios[ref], a_indexSpace, m_bc);
    }

  int lev = ref;

  newOp->m_alpha = m_alpha;
  newOp->m_beta  = m_beta;

  newOp->m_aCoef = m_aCoef[lev];
  newOp->m_bCoef = m_bCoef[lev];

  newOp->computeLambda();

  Real dxCrse = -1;
  if (ref > 0)
    {
      dxCrse = m_dx[ref-1];
    }
  newOp->m_dxCrse = dxCrse;

  return newOp;
}

int VCAMRPoissonOpFactory::refToFiner(const ProblemDomain& a_domain) const
{
  int retval = -1;
  bool found = false;
  for (int ilev = 0; ilev < m_domains.size(); ilev++)
    {
      if (m_domains[ilev].domainBox() == a_domain.domainBox())
        {
          retval = m_refRatios[ilev];
          found = true;
        }
    }
  if (!found)
    {
      MayDay::Error("Domain not found in AMR hierarchy");
    }
  return retval;
}

void VCAMRPoissonOpFactory::setDefaultValues()
{
  // Default to Laplacian operator
  m_alpha = 0.0;
  m_beta = -1.0;

  m_coefficient_average_type = CoarseAverage::arithmetic;
}

#include "NamespaceFooter.H"
