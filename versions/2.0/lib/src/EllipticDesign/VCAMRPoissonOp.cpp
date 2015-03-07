#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "VCAMRPoissonOp.H"
#include "FORT_PROTO.H"
#include "VCAMRPoissonOpF_F.H"
#include "BoxIterator.H"
#include "AverageF_F.H"
#include "InterpF_F.H"
#include "LayoutIterator.H"
#include "FineInterp.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "NamespaceHeader.H"

void VCAMRPoissonOp::createCoarsened(LevelData<FArrayBox>&       a_lhs,
                                     const LevelData<FArrayBox>& a_rhs,
                                     const int &                 a_refRat)
{
  int ncomp = a_rhs.nComp();
  IntVect ghostVect = a_rhs.ghostVect();

  DisjointBoxLayout dbl = a_rhs.disjointBoxLayout();
  CH_assert(dbl.coarsenable(a_refRat));

  //fill ebislayout
  DisjointBoxLayout dblCoarsenedFine;
  coarsen(dblCoarsenedFine, dbl, a_refRat);

  a_lhs.define(dblCoarsenedFine, ncomp, a_rhs.ghostVect());
  
}
/** full define function for AMRLevelOp with both coarser and finer levels */
void VCAMRPoissonOp::define(const DisjointBoxLayout& a_grids,
                            const DisjointBoxLayout& a_gridsFiner,
                            const DisjointBoxLayout& a_gridsCoarser,
                            const Real&              a_dxLevel,
                            int                      a_refRatio,
                            int                      a_refRatioFiner,
                            const ProblemDomain&     a_domain,
                            BCFunc                   a_bc)
{
  this->define(a_grids, a_gridsCoarser, a_dxLevel, a_refRatio, a_domain, a_bc);
  m_refToFiner = a_refRatioFiner;
}

  /** full define function for AMRLevelOp<LevelData<FArrayBox> > with finer levels, but no coarser */
void VCAMRPoissonOp::define(const DisjointBoxLayout& a_grids,
                            const DisjointBoxLayout& a_gridsFiner,
                            const Real&              a_dxLevel,
                            int                      a_refRatio, // dummy arg
                            int                      a_refRatioFiner,
                            const ProblemDomain&     a_domain,
                            BCFunc                   a_bc)
{
  CH_assert(a_refRatio == 1);
  this->define(a_grids, a_dxLevel, a_domain, a_bc); //calls the MG version of define
  m_refToFiner = a_refRatioFiner;
}

void VCAMRPoissonOp::define(const DisjointBoxLayout& a_grids,
                            const DisjointBoxLayout& a_coarse,
                            const Real&              a_dxLevel,
                            int                      a_refRatio,
                            const ProblemDomain&     a_domain,
                            BCFunc                   a_bc)
{
  m_refToCoarser = a_refRatio;
  m_dxCrse = a_refRatio*a_dxLevel;
  m_refToFiner = 1;

  m_bc     = a_bc;
  m_domain = a_domain;
  m_dx     = a_dxLevel;
  m_dxCrse = a_refRatio*a_dxLevel;

  //these get set again after define is called
    
  m_exchangeCopier.define(a_grids, a_grids, IntVect::Unit, true);
  for(int i=0; i<CH_SPACEDIM; ++i)
    {
      LayoutData<CFIVS>& lo =  m_loCFIVS[i];
      LayoutData<CFIVS>& hi =  m_hiCFIVS[i];
      lo.define(a_grids);
      hi.define(a_grids);
      for(DataIterator dit(a_grids); dit.ok(); ++dit)
        {
          lo[dit].define(m_domain, a_grids.get(dit), a_grids, i, Side::Lo);
          hi[dit].define(m_domain, a_grids.get(dit), a_grids, i, Side::Hi);
        }
    }
  m_interpWithCoarser.define(a_grids, &a_coarse, a_dxLevel,
                             m_refToCoarser, 1, m_domain);  
}

void VCAMRPoissonOp::define(const DisjointBoxLayout& a_grids,
                            const Real&              a_dx,
                            const ProblemDomain&     a_domain, 
                            BCFunc                   a_bc)
{
  m_bc     = a_bc;
  m_domain = a_domain;
  m_dx     = a_dx;
  m_dxCrse = 2*a_dx;
  m_refToCoarser = 2; // redefined in AMRLevelOp<LevelData<FArrayBox> >::define virtual function.
  m_refToFiner   = 2;

    
  m_exchangeCopier.define(a_grids, a_grids, IntVect::Unit, true);
  for(int i=0; i<CH_SPACEDIM; ++i)
    {
      LayoutData<CFIVS>& lo =  m_loCFIVS[i];
      LayoutData<CFIVS>& hi =  m_hiCFIVS[i];
      lo.define(a_grids);
      hi.define(a_grids);
      for(DataIterator dit(a_grids); dit.ok(); ++dit)
        {
          lo[dit].define(a_domain, a_grids.get(dit),a_grids, i, Side::Lo);
          hi[dit].define(a_domain, a_grids.get(dit),a_grids, i, Side::Hi);
        }
    }
}

void VCAMRPoissonOp::residual(LevelData<FArrayBox>& a_lhs, 
                              const LevelData<FArrayBox>& a_phi,
                              const LevelData<FArrayBox>& a_rhs, 
                              bool a_homogeneous)
{
  applyOp(a_lhs, a_phi, a_homogeneous);
  incr(a_lhs, a_rhs, -1);
  scale(a_lhs, -1.0);
}

/**************************/
// this preconditioner first initializes phihat to (IA)phihat = rhshat
// (diagonization of L -- A is the matrix version of L)
// then smooths with a couple of passes of levelGSRB
void VCAMRPoissonOp::preCond(LevelData<FArrayBox>& a_phi, 
                             const LevelData<FArrayBox>& a_rhs)
{
  // diagonal term of this operator is 4/h/h in 2D, 6/h/h in 3D,
  // so inverse of this is our initial multiplier
  Real mult = -m_dx*m_dx;
  //Interval comps = a_phi.interval();
  int ncomp = a_phi.nComp();
  CH_assert(a_phi.nComp() == a_rhs.nComp());
  CH_assert(m_beta->nComp() == ncomp);

  // don't need to use a Copier -- plain copy will do
  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      // also need to average and sum face-centered betas to cell-centers
      Box gridBox = a_rhs[dit].box();
      FArrayBox multFAB(gridBox, ncomp);
      multFAB.copy((*m_alpha)[dit]);
      const FluxBox& thisBeta = (*m_beta)[dit];      
      Real scale = 1.0;
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FORT_SUMFACES(CHF_FRA(multFAB),
                        CHF_CONST_FRA(thisBeta[dir]),
                        CHF_BOX(gridBox),
                        CHF_INT(dir),
                        CHF_REAL(scale));
        }
      multFAB.invert(mult);

      a_phi[dit].copy(a_rhs[dit]);

      a_phi[dit].mult(multFAB, gridBox, 0, 0, ncomp);
    }

  relax(a_phi, a_rhs, 2);
}

void VCAMRPoissonOp::applyOp(LevelData<FArrayBox>& a_lhs, 
                             const LevelData<FArrayBox>& a_phi,
                             bool a_homogeneous )
{

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


// computes VC operator after BC's have been set
void
VCAMRPoissonOp::computeOperatorNoBCs(LevelData<FArrayBox>& a_lhs,
                                     const LevelData<FArrayBox>& a_phi)
{

  Real dx = m_dx;
  const DisjointBoxLayout dbl= a_phi.getBoxes();
  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisLHS = a_lhs[dit];
      const Box& region = dbl[dit()];
      // first set with alpha*I
      if (m_alpha != NULL)
        {          
          thisLHS.copy(a_phi[dit], region);
          thisLHS.mult((*m_alpha)[dit], region,0,0,a_lhs.nComp());
        }
      
      // need to loop over directions to do the Div(Beta*Grad) part,
      // at least for now 
      const FluxBox& thisBeta = (*m_beta)[dit];

      for (int idir=0; idir<SpaceDim; idir++)
        {
          
          FORT_OPERATORLAPVC(CHF_FRA(thisLHS),
                             CHF_CONST_FRA(a_phi[dit]),
                             CHF_CONST_FRA(thisBeta[idir]),
                             CHF_BOX(region),
                             CHF_CONST_REAL(dx),
                             CHF_CONST_INT(idir));
        }
    } // end loop over boxes
}

///compute norm over all cells on coarse not covered by finer
Real
VCAMRPoissonOp::AMRNorm(const LevelData<FArrayBox>& a_coarResid,
                        const LevelData<FArrayBox>& a_fineResid,
                        const int& a_refRat,
                        const int& a_ord)

{
  const DisjointBoxLayout& coarGrids = a_coarResid.disjointBoxLayout();

  //create temp and zero out under finer grids
  LevelData<FArrayBox> coarTemp;
  m_levelOps.create(coarTemp, a_coarResid);
  m_levelOps.assign(coarTemp, a_coarResid);

  if (a_fineResid.isDefined())
    {
      const DisjointBoxLayout& fineGrids = a_fineResid.disjointBoxLayout();
      int ncomp = coarTemp.nComp();
      for(DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
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
void VCAMRPoissonOp::create(LevelData<FArrayBox>& a_lhs, 
                            const LevelData<FArrayBox>& a_rhs)
{
  m_levelOps.create(a_lhs, a_rhs);
}

void VCAMRPoissonOp::assign(LevelData<FArrayBox>& a_lhs, 
                            const LevelData<FArrayBox>& a_rhs)
{
  m_levelOps.assign(a_lhs, a_rhs);
}

Real VCAMRPoissonOp::dotProduct(const LevelData<FArrayBox>& a_1, 
                                const LevelData<FArrayBox>& a_2)
{
  return m_levelOps.dotProduct(a_1, a_2);
}
void VCAMRPoissonOp::incr( LevelData<FArrayBox>& a_lhs, 
                           const LevelData<FArrayBox>& a_x, Real a_scale)
{
  m_levelOps.incr(a_lhs, a_x, a_scale);
}
void VCAMRPoissonOp::axby( LevelData<FArrayBox>& a_lhs, 
                           const LevelData<FArrayBox>& a_x,
                           const LevelData<FArrayBox>& a_y, 
                           Real a_a, Real a_b)
{
  m_levelOps.axby(a_lhs, a_x, a_y, a_a, a_b);
}

void VCAMRPoissonOp::scale(LevelData<FArrayBox>& a_lhs, 
                           const Real& a_scale)
{
  m_levelOps.scale(a_lhs, a_scale);
}
void VCAMRPoissonOp::setToZero(LevelData<FArrayBox>& a_lhs)
{
  m_levelOps.setToZero(a_lhs);
}

void VCAMRPoissonOp::relax(LevelData<FArrayBox>& a_e,
                           const LevelData<FArrayBox>& a_residual,
                           int a_iterations)
{
  for(int i=0; i<a_iterations; i++)
    {      
      levelGSRB(a_e, a_residual);
      //levelJacobi(a_e, a_residual);
    }
}

void VCAMRPoissonOp::createCoarser(LevelData<FArrayBox>& a_coarse,
                                   const LevelData<FArrayBox>& a_fine,
                                   bool ghosted)
{
  // CH_assert(!ghosted);
  IntVect ghost = a_fine.ghostVect();
  DisjointBoxLayout dbl;
  CH_assert(dbl.coarsenable(2));
  coarsen(dbl, a_fine.disjointBoxLayout(), 2); //multigrid, so coarsen by 2
  a_coarse.define(dbl, a_fine.nComp(), ghost);
}

void VCAMRPoissonOp::restrictResidual(LevelData<FArrayBox>& a_resCoarse,
                                      LevelData<FArrayBox>& a_phiFine,
                                      const LevelData<FArrayBox>& a_rhsFine)
{
  homogeneousCFInterp(a_phiFine);
  const DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for(DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& phi = a_phiFine[dit];
      m_bc(phi,  dblFine[dit()], m_domain,  m_dx, true);
    }

  a_phiFine.exchange(a_phiFine.interval());

  // temp storage
  LevelData<FArrayBox> LofPhi;
  create(LofPhi, a_rhsFine);
  computeOperatorNoBCs(LofPhi, a_phiFine);

  for(DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit)
    {
      const FArrayBox& rhs = a_rhsFine[dit];
      FArrayBox& res = a_resCoarse[dit];

      Box region = dblFine.get(dit());

      res.setVal(0.0);
      FORT_RESTRICTRESVC(CHF_CONST_FRA(LofPhi[dit]),
                         CHF_CONST_FRA(rhs),
                         CHF_FRA(res),
                         CHF_BOX(region),
                         CHF_CONST_REAL(m_dx));

    }
}

Real VCAMRPoissonOp::norm(const LevelData<FArrayBox>& a_x, int a_ord)
{
  return CH_XD::norm(a_x, a_x.interval(), a_ord);
}
/***/

void VCAMRPoissonOp::prolongIncrement(LevelData<FArrayBox>& a_phiThisLevel,
                                      const LevelData<FArrayBox>& a_correctCoarse)
{
  DisjointBoxLayout dbl = a_phiThisLevel.disjointBoxLayout();
  int mgref = 2; //this is a multigrid func
  for(DataIterator dit = a_phiThisLevel.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& phi =  a_phiThisLevel[dit];
      const FArrayBox& coarse = a_correctCoarse[dit];
      Box region = dbl.get(dit());
      Box cBox = coarsen(region, m_refToCoarser);

      FORT_PROLONGVC(CHF_FRA(phi),
                     CHF_CONST_FRA(coarse),
                     CHF_BOX(region),
                     CHF_CONST_INT(mgref));

    }
}

/***/
void VCAMRPoissonOp::
levelGSRB(LevelData<FArrayBox>&       a_phi,
          const LevelData<FArrayBox>& a_rhs)
{

  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());
  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();

  // simplest thing to do is to allocate a temporary and apply 
  // the operator all at once.  We can revisit this later if 
  // necessary (DFM 12/6/2006)
  LevelData<FArrayBox> LofPhi;
  create(LofPhi, a_rhs);

  // will also need the cell-averaged beta to compute the
  // proper GSRB coefficient
  LevelData<FArrayBox> sumBeta;
  create(sumBeta, a_rhs);
  DataIterator dit = a_phi.dataIterator();
  Real scale = 1.0;
  for (dit.begin(); dit.ok(); ++dit)
    {
      sumBeta[dit].setVal(0.0);
      FluxBox& thisBeta = (*m_beta)[dit];
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FORT_SUMFACES(CHF_FRA(sumBeta[dit]),
                        CHF_FRA(thisBeta[dir]),
                        CHF_BOX(sumBeta[dit].box()),
                        CHF_INT(dir),
                        CHF_REAL(scale));
        }
    }
                            

  // do first red, then black passes
  for (int whichPass =0; whichPass <= 1; whichPass++)
    {

      homogeneousCFInterp(a_phi);

      //fill in intersection of ghostcells and a_phi's boxes
      a_phi.exchange(a_phi.interval(), m_exchangeCopier);
      
      for (dit.begin(); dit.ok(); ++dit)
        {
          // invoke physical BC's where necessary
          m_bc(a_phi[dit],  dbl[dit()], m_domain,  m_dx, true);
        }

      // note that this isn't ideal, since we're computing 
      // over 2x as many points as necessary.  However, it's 
      // simpler (and safer) than re-writing the operator at this point 
      // we will probably revisit this once things are stable.
      computeOperatorNoBCs(LofPhi, a_phi);

      // now step through grids...
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& region = dbl.get(dit());

          FORT_GSRBLEVELLAPVC(CHF_FRA(a_phi[dit]),
                              CHF_CONST_FRA(LofPhi[dit]),
                              CHF_CONST_FRA(a_rhs[dit]),
                              CHF_CONST_FRA((*m_alpha)[dit]),
                              CHF_CONST_FRA(sumBeta[dit]),
                              CHF_BOX(region),
                              CHF_CONST_REAL(m_dx),
                              CHF_CONST_INT(whichPass));

        } // end loop through grids

    } // end loop through red-black
}

void VCAMRPoissonOp::levelJacobi(LevelData<FArrayBox>& a_phi,
                                 const LevelData<FArrayBox>& a_rhs)
{
  LevelData<FArrayBox> resid;
  create(resid, a_rhs);
  // Get the residual
  residual(resid,a_phi,a_rhs,true);

  LevelData<FArrayBox> weights;
  create(weights, a_rhs);

  // Create the weight
  m_alpha->copyTo(m_alpha->interval(), weights, weights.interval());

  DataIterator dit = weights.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisBeta = (*m_beta)[dit];
      FArrayBox& thisWeight = weights[dit];
      const Box thisBox = thisWeight.box();
      Real scale = -2.0/(m_dx*m_dx);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          //weight += -2.0 * m_beta/ (m_dx*m_dx);
          FORT_SUMFACES(CHF_FRA(thisWeight),
                        CHF_FRA(thisBeta[idir]),
                        CHF_BOX(thisBox),
                        CHF_INT(idir),
                        CHF_REAL(scale));
        } // end loop over directions
      // scale residual by 1/weight
      Real invertScale = 1.0;
      thisWeight.invert(invertScale);
      resid[dit].mult(thisWeight);
    } // end loop over grids
  // Do the Jacobi relaxation
  
  incr(a_phi, resid, 0.5);
}

void
VCAMRPoissonOp::homogeneousCFInterp(LevelData<FArrayBox>& a_phif)
{

  CH_assert( a_phif.ghostVect() >= IntVect::Unit);

  DataIterator dit = a_phif.dataIterator();
  for(dit.begin(); dit.ok(); ++dit)
    {
      const DataIndex& datInd = dit();
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          SideIterator sit;
          for(sit.begin(); sit.ok(); sit.next())
            {
              homogeneousCFInterp(a_phif,datInd,idir,sit());
            }
        }
    }
}

void
VCAMRPoissonOp::homogeneousCFInterp(LevelData<FArrayBox>& a_phif,
                                    const DataIndex& a_datInd,
                                    int a_idir,
                                    Side::LoHiSide a_hiorlo)
{

  CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
  CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
  CH_assert( a_phif.ghostVect() >= IntVect::Unit);
  //  CH_assert (m_ncomp == a_phif.nComp());

  const CFIVS* cfivs_ptr = NULL;
  if(a_hiorlo == Side::Lo)
    cfivs_ptr = &m_loCFIVS[a_idir][a_datInd];
  else
    cfivs_ptr = &m_hiCFIVS[a_idir][a_datInd];

  if(cfivs_ptr->isPacked())
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
      if(!interp_ivs.isEmpty())
        {
          // Assuming homogenous, interpolate on fine ivs
          interpOnIVSHomo(a_phif, a_datInd, a_idir,
                          a_hiorlo, interp_ivs);
        }
    }
}

void
VCAMRPoissonOp::interpOnIVSHomo(LevelData<FArrayBox>& a_phif,
                                const DataIndex& a_datInd,
                                const int a_idir,
                                const Side::LoHiSide a_hiorlo,
                                const IntVectSet& a_interpIVS)
{

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
  for(fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit)
    {
      ivf = fine_ivsit();
      // quadratic interpolation
      for(int ivar = 0; ivar < a_phif.nComp(); ivar++)
        {
          ivf[a_idir]-=2*ihilo;
          pa = a_phi(ivf, ivar);
          ivf[a_idir]+=ihilo;
          pb = a_phi(ivf, ivar);

          a = ((pb-pa)*m1 - (pb)*m2)*idenom;
          b = (pb)*q1 - a*q2;

          ivf[a_idir]+=ihilo;
          a_phi(fine_ivsit(), ivar) = a*xsquared + b*x + pa;

        } //end loop over components
    } //end loop over fine intvects
}

void VCAMRPoissonOp::AMRResidualNF(LevelData<FArrayBox>& a_residual,
                                   const LevelData<FArrayBox>& a_phi,
                                   const LevelData<FArrayBox>& a_phiCoarse,
                                   const LevelData<FArrayBox>& a_rhs,
                                   bool a_homogeneousPhysBC)
{
  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  if (a_phiCoarse.isDefined())
    {
      m_interpWithCoarser.coarseFineInterp(phi, a_phiCoarse);
    }
  this->residual(a_residual, a_phi, a_rhs, a_homogeneousPhysBC ); //apply boundary conditions
}


void VCAMRPoissonOp::AMROperatorNF(LevelData<FArrayBox>& a_LofPhi,
                                   const LevelData<FArrayBox>& a_phi,
                                   const LevelData<FArrayBox>& a_phiCoarse,
                                   bool a_homogeneousPhysBC)
{
  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  if (a_phiCoarse.isDefined())
    {
      m_interpWithCoarser.coarseFineInterp(phi, a_phiCoarse);
    }
  //apply boundary conditions in applyOp
  this->applyOp(a_LofPhi, a_phi, a_homogeneousPhysBC );
}

void VCAMRPoissonOp::getFlux(FArrayBox&       a_flux,
                             const FArrayBox& a_data,
                             const FluxBox&   a_beta,
                             const Box& a_facebox,
                             int              a_dir,
                             int ref) const
{

  CH_assert(a_dir >= 0);
  CH_assert(a_dir <  SpaceDim);
  CH_assert(!a_data.box().isEmpty());
  CH_assert(!a_facebox.isEmpty());

  // probably the simplest way to test centering
  // a_box needs to be face-centered in the a_dir
  Box faceTestBox(IntVect::Zero, IntVect::Unit);
  faceTestBox.surroundingNodes(a_dir);
  CH_assert(a_facebox.type() == faceTestBox.type());

  const FArrayBox& betaDir = a_beta[a_dir];


  // reality check for beta
  CH_assert(betaDir.box().contains(a_facebox));
  
  a_flux.resize(a_facebox, a_data.nComp());
  BoxIterator bit(a_facebox);

  Real dx(m_dx);
  dx/=ref;
  for( bit.begin(); bit.ok(); bit.next())
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
          Real gradphi =(phihi-philo)/dx;

          a_flux(iv,ivar) = betaDir(iv, ivar)*gradphi;
        }
    }
}

void VCAMRPoissonOp::reflux(const LevelData<FArrayBox>& a_phiFine,
                            const LevelData<FArrayBox>& a_phi,
                            LevelData<FArrayBox>& residual,
                            AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
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
      const FluxBox& coarBeta = (*m_beta)[dit];
      const Box& gridBox = a_phi.getBoxes()[dit];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          FArrayBox coarflux;
          Box faceBox = surroundingNodes(gridBox, idir);
          getFlux(coarflux, coarfab, coarBeta, faceBox, idir);

          Real scale = 1.0;
          levfluxreg.incrementCoarse(coarflux, scale,dit(), 
                                     interv,interv,idir);
        }
    }
  LevelData<FArrayBox>& p = ( LevelData<FArrayBox>&)a_phiFine;

  //has to be its own object because the finer operator 
  //owns an interpolator and we have no way of getting to it
  VCAMRPoissonOp* finerAMRPOp = (VCAMRPoissonOp*) a_finerOp;
  QuadCFInterp& quadCFI = finerAMRPOp->m_interpWithCoarser;
  
  quadCFI.coarseFineInterp(p, a_phi);
  p.exchange(a_phiFine.interval());
  IntVect phiGhost = p.ghostVect();

  DataIterator ditf = a_phiFine.dataIterator();
  const  DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (ditf.reset(); ditf.ok(); ++ditf)
    {
      const FArrayBox& phifFab = a_phiFine[ditf];
      const FluxBox& fineBeta = (*(finerAMRPOp->m_beta))[ditf];
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
              getFlux(fineflux, phifab, fineBeta, facebox, idir,
                      m_refToFiner);
              
              Real scale = 1.0;
              levfluxreg.incrementFine(fineflux, scale, ditf(),
                                       interv, interv, idir, hiorlo);
            }
        }
    }

  Real scale =  1.0/m_dx;
  levfluxreg.reflux(residual, scale);
}

void VCAMRPoissonOp::AMRResidual(LevelData<FArrayBox>& a_residual,
                                 const LevelData<FArrayBox>& a_phiFine,
                                 const LevelData<FArrayBox>& a_phi,
                                 const LevelData<FArrayBox>& a_phiCoarse,
                                 const LevelData<FArrayBox>& a_rhs,
                                 bool a_homogeneousPhysBC,
                                 AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{

  AMROperator(a_residual, a_phiFine, a_phi, a_phiCoarse,
              a_homogeneousPhysBC, a_finerOp);
  incr(a_residual, a_rhs, -1.0);

  scale(a_residual, -1.0);
    
}

void VCAMRPoissonOp::AMRResidualNC(LevelData<FArrayBox>& a_residual,
                                   const LevelData<FArrayBox>& a_phiFine,
                                   const LevelData<FArrayBox>& a_phi,
                                   const LevelData<FArrayBox>& a_rhs,
                                   bool a_homogeneousPhysBC,
                                   AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  
  AMROperatorNC(a_residual, a_phiFine, a_phi, 
                a_homogeneousPhysBC, a_finerOp);
  axby(a_residual, a_residual, a_rhs, -1.0, 1.0);
}


void VCAMRPoissonOp::AMROperator(LevelData<FArrayBox>& a_LofPhi,
                                 const LevelData<FArrayBox>& a_phiFine,
                                 const LevelData<FArrayBox>& a_phi,
                                 const LevelData<FArrayBox>& a_phiCoarse,
                                 bool a_homogeneousPhysBC,
                                 AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  if (a_phiCoarse.isDefined())
    {
      m_interpWithCoarser.coarseFineInterp(phi, a_phiCoarse);
    }
  applyOp(a_LofPhi, a_phi, a_homogeneousPhysBC);
  if (a_finerOp != NULL)
    {
      reflux(a_phiFine, a_phi,  a_LofPhi, a_finerOp);
    }  
}

void VCAMRPoissonOp::AMROperatorNC(LevelData<FArrayBox>& a_LofPhi,
                                   const LevelData<FArrayBox>& a_phiFine,
                                   const LevelData<FArrayBox>& a_phi,
                                   bool a_homogeneousPhysBC,
                                   AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  //no coarse-fine interpolation here
  applyOp(a_LofPhi, a_phi, a_homogeneousPhysBC);
  reflux(a_phiFine, a_phi,  a_LofPhi, a_finerOp);
}

void VCAMRPoissonOp::AMRRestrict(LevelData<FArrayBox>& a_resCoarse,
                                 const LevelData<FArrayBox>& a_residual,
                                 const LevelData<FArrayBox>& a_correction,
                                 const LevelData<FArrayBox>& a_coarseCorrection)
{
  
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
      FORT_AVERAGE( CHF_FRA(coarse),
                    CHF_CONST_FRA(fine),
                    CHF_BOX(b),
                    CHF_CONST_INT(m_refToCoarser),
                    CHF_BOX(refbox)
                    );
    }
}

/** a_correction += I[2h->h](a_coarseCorrection) */
void VCAMRPoissonOp::AMRProlong(LevelData<FArrayBox>& a_correction,
                                const LevelData<FArrayBox>& a_coarseCorrection)
{
  DisjointBoxLayout c;
  coarsen(c,  a_correction.disjointBoxLayout(), m_refToCoarser);
  LevelData<FArrayBox> eCoar(c, a_correction.nComp(),a_coarseCorrection.ghostVect());
  a_coarseCorrection.copyTo(eCoar.interval(), eCoar, eCoar.interval());

  DisjointBoxLayout dbl = a_correction.disjointBoxLayout();
  for(DataIterator dit = a_correction.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& phi =  a_correction[dit];
      const FArrayBox& coarse = eCoar[dit];
      Box region = dbl.get(dit());

      FORT_PROLONGVC(CHF_FRA(phi),
                     CHF_CONST_FRA(coarse),
                     CHF_BOX(region),
                     CHF_CONST_INT(m_refToCoarser));

    }
}

void VCAMRPoissonOp::AMRUpdateResidual(LevelData<FArrayBox>& a_residual,
                                       const LevelData<FArrayBox>& a_correction,
                                       const LevelData<FArrayBox>& a_coarseCorrection)
{
  LevelData<FArrayBox> r;
  this->create(r, a_residual);
  this->AMRResidualNF(r, a_correction, a_coarseCorrection, a_residual, true);
  this->assign(a_residual, r);
}

//==========================================================
//
//==========================================================

// MultiGrid define function

VCAMRPoissonOpFactory::VCAMRPoissonOpFactory()
{
  setDefaultValues();
}

void VCAMRPoissonOpFactory::define(const ProblemDomain& a_domain,
                                   const DisjointBoxLayout& a_grid,
                                   const Real& a_dx,
                                   BCFunc a_bc,
                                   int maxDepth)
{
  setDefaultValues();

  Vector<DisjointBoxLayout> grids(1, a_grid);
  Vector<int> refRatio(1, 2);
  define(a_domain, grids, refRatio, a_dx, a_bc);
}

VCAMRPoissonOp*
VCAMRPoissonOpFactory::MGnewOp(const ProblemDomain& a_indexSpace,
                               int depth,
                               bool homoOnly)
{
  // CH_assert(m_boxes.size()>depth);

  int lev=0;
  for(;lev< m_domains.size(); lev++)
    {
      if(a_indexSpace.domainBox() == m_domains[lev].domainBox()) break;
    }
  CH_assert(lev !=  m_domains.size()); // didn't find domain

  DisjointBoxLayout layout(m_boxes[lev]);
  ProblemDomain domain(m_domains[lev]);
  Real dx = m_dx[lev];
  int coarsenRatio = 1;

  for(int i=0; i<depth; i++)
    {
      if(!layout.coarsenable(4)) return NULL;
      DisjointBoxLayout dbl;
      coarsen_dbl(dbl, layout, 2);
      layout = dbl;
      dx*=2;
      coarsenRatio *= 2;
      domain.coarsen(2);
    }

  VCAMRPoissonOp* newOp = new VCAMRPoissonOp;
  newOp->define(layout, dx, domain, m_bc);
  if (depth == 0)
    {
      // don't need to coarsen anything for this
      newOp->m_alpha = m_alpha[lev];
      newOp->m_beta  = m_beta[lev];
    } 
  else
    {
      // need to coarsen coefficients
      RefCountedPtr<LevelData<FArrayBox> > alpha( new LevelData<FArrayBox> );
      RefCountedPtr<LevelData<FluxBox> > beta( new LevelData<FluxBox> );
      alpha->define(layout, m_alpha[lev]->nComp(), m_alpha[lev]->ghostVect());
      beta->define(layout, m_beta[lev]->nComp(), m_beta[lev]->ghostVect());

      // average coefficients to coarser level
      // for now, do this with a CoarseAverage -- 
      // may want to switch to harmonic averaging at some point
      CoarseAverage averager(m_alpha[lev]->getBoxes(),
                             layout, alpha->nComp(), coarsenRatio);

      CoarseAverageFace faceAverager(m_beta[lev]->getBoxes(),
                                     beta->nComp(), coarsenRatio);

      if (m_coefficient_average_type == CoarseAverage::arithmetic) 
        {
          averager.averageToCoarse(*alpha, *(m_alpha[lev]));
          faceAverager.averageToCoarse(*beta, *(m_beta[lev]));
        }
      else if (m_coefficient_average_type == CoarseAverage::harmonic) 
        {
          averager.averageToCoarseHarmonic(*alpha, *(m_alpha[lev]));
          faceAverager.averageToCoarseHarmonic(*beta, *(m_beta[lev]));
        }
      else 
        { 
          MayDay::Error("VCAMRPoissonOpFactory::MGNewOp -- bad averagetype");
        }
      
      newOp->m_alpha = alpha;
      newOp->m_beta = beta;
    }

  Real dxCrse = -1;
  if(lev > 0)
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
  for(;ref< m_domains.size(); ref++)
    {
      if(a_indexSpace.domainBox() == m_domains[ref].domainBox()) break;
    }

  if(ref == 0) 
    {
      // coarsest AMR level
      if(m_domains.size() == 1)
        {
          //no finer level
          newOp->define(m_boxes[0],  m_dx[0],
                        a_indexSpace, m_bc);
        }
      else
        {
          //finer level exists but no coarser 
          int dummyRat = 1;  //argument so compiler can find right function
          int refToFiner = m_refRatios[0]; //actual refinement ratio
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
  newOp->m_alpha = m_alpha[lev];
  newOp->m_beta  = m_beta[lev];
  Real dxCrse = -1;
  if(ref > 0)
    {
      dxCrse = m_dx[ref-1];
    }
  newOp->m_dxCrse = dxCrse;
  return newOp;
}

//  AMR Factory define function
void VCAMRPoissonOpFactory::define(const ProblemDomain& a_coarseDomain,
                                   const Vector<DisjointBoxLayout>& a_grids,
                                   const Vector<int>& a_refRatios,
                                   const Real& a_coarsedx,
                                   BCFunc a_bc)
{
  setDefaultValues();
  
  m_domains.resize(a_grids.size());
  m_boxes=a_grids;
  m_refRatios=a_refRatios;
  m_dx.resize(a_grids.size());
  m_bc = a_bc;
  m_domains[0] = a_coarseDomain;
  m_dx[0] = a_coarsedx;
  for(int i=1; i<a_grids.size(); i++)
    {
      m_dx[i] = m_dx[i-1]/m_refRatios[i] ;
      m_domains[i] = m_domains[i-1];
      m_domains[i].refine(m_refRatios[i-1]);
    }

}

/// set diagonal coefficient
void
VCAMRPoissonOpFactory::setAlpha(Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_alpha)
{
  m_alpha = a_alpha;
}
  
/// only sets one level
void
VCAMRPoissonOpFactory::setAlpha(RefCountedPtr<LevelData<FArrayBox> >& a_alpha, 
                                int a_level)
{
  if (m_alpha.size() < a_level+1)
    {
      m_alpha.resize(a_level+1);
    }
  m_alpha[a_level] = a_alpha;
}

///
void
VCAMRPoissonOpFactory::setBeta(Vector<RefCountedPtr<LevelData<FluxBox> > >& a_beta)
{
  m_beta = a_beta;
}

///
void
VCAMRPoissonOpFactory::setBeta(RefCountedPtr<LevelData<FluxBox> >& a_beta,
                               int a_level)
{
  if (m_beta.size() < a_level+1)
    {
      m_beta.resize(a_level+1);
    }
  
  m_beta[a_level] = a_beta;
}




int VCAMRPoissonOpFactory::refToFiner(const ProblemDomain& a_domain) const
  
{

  int retval = -1;
  bool found = false;
  for(int ilev = 0; ilev < m_domains.size(); ilev++)
    {
      if(m_domains[ilev].domainBox() == a_domain.domainBox())
        {
          retval = m_refRatios[ilev];
          found = true;
        }
    }
  if(!found)
    {
      MayDay::Error("Domain not found in AMR hierarchy");
    }
  return retval;
}

void
VCAMRPoissonOpFactory::setDefaultValues()
{
  m_coefficient_average_type = CoarseAverage::arithmetic;
}

#include "NamespaceFooter.H"
