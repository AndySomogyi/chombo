#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
// DTGraves, Weds, July 21, 1999

/// This is purely a workbench class encapsulating the most common
/// stuff used by derived operators (Yuri Omelchenko)
#include "ParentLevelOp.H"
#include "DotProduct.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
#include "BiCGStabSmoother.H"

///
void 
ParentLevelOp::setDomainGhostBC(const DomainGhostBC& a_dombcIn)
{
  m_domghostbc = a_dombcIn;
  m_isBCDefined = true;
}

///
void
ParentLevelOp::setBottomSmoother(const BaseBottomSmoother& a_bottomSmoother)
{
  if (m_bottom_smoother_ptr != NULL)
    {
      delete m_bottom_smoother_ptr;
    }
  m_bottom_smoother_ptr = a_bottomSmoother.new_bottomSmoother();
}

///
bool 
ParentLevelOp::isDefined() const
{
  return (m_isDefined && m_isBCDefined);
}

// set default values
void
ParentLevelOp::setDefaultValues()
{
  m_dxLevel = -1.0;
  m_refRatio = -1;
  m_isDefined = false;
  m_ihcfiEnabled = false;
  m_isBCDefined = false;
  // use BiCGSTab as bottom smoother by default
  m_bottom_smoother_ptr = new BiCGStabSmoother();
}

// return to undefined state
void
ParentLevelOp::clearMemory()
{
  m_dxLevel = -1.0;
  m_refRatio = -1;
  m_isDefined = false;
  m_ihcfiEnabled = false;
  // don't touch bottom smoother here, since we want it to be able
  // to pass through the define statement unaltered

}

// default constructor
ParentLevelOp::ParentLevelOp(): m_bottom_smoother_ptr(NULL)
{
  ParentLevelOp::setDefaultValues();
}

///
ParentLevelOp::~ParentLevelOp() 
{
  ParentLevelOp::clearMemory();
  if(m_bottom_smoother_ptr != NULL)
    {
      delete m_bottom_smoother_ptr;
      m_bottom_smoother_ptr = NULL;
    }
}

///
void
ParentLevelOp::define(
                  const DisjointBoxLayout& a_grids,
                  const DisjointBoxLayout* a_baseBAPtr,
                  Real  a_dxLevel, 
                  int a_refRatio,
                  const Box& a_domain,
                  bool a_homogeneousOnly,
                  int a_ncomp)
{
  ProblemDomain probdomain(a_domain);
  define(a_grids, a_baseBAPtr, a_dxLevel, a_refRatio, probdomain,
         a_homogeneousOnly, a_ncomp);
}

///
void
ParentLevelOp::define(
                  const DisjointBoxLayout& a_grids,
                  const DisjointBoxLayout* a_baseBAPtr,
                  Real  a_dxLevel, 
                  int a_refRatio,
                  const ProblemDomain& a_domain,
                  bool a_homogeneousOnly,
                  int a_ncomp)
{
  ParentLevelOp::clearMemory();
  m_isDefined = true;
  CH_assert(!a_domain.isEmpty());
  CH_assert(a_grids.checkPeriodic(a_domain));

  m_dxLevel = a_dxLevel;
  m_domain = a_domain;
  m_grids = a_grids;
  m_exchangeCopier.define(a_grids, a_grids, IntVect::Unit);
  m_refRatio = a_refRatio;
  m_dxCrse = m_refRatio*m_dxLevel;
  m_ihcfiEnabled = (!a_homogeneousOnly);
  m_ncomp = a_ncomp;
  if(m_ihcfiEnabled)
    m_quadCFI.define(a_grids,a_baseBAPtr, a_dxLevel, 
                     a_refRatio, m_ncomp, m_domain);
  if(a_baseBAPtr != NULL)
    {
      m_baseBA = *a_baseBAPtr;
    }
  DataIterator lit = m_grids.dataIterator();
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      m_loCFIVS[idir].define(m_grids);
      m_hiCFIVS[idir].define(m_grids);
      for(lit.begin(); lit.ok(); ++lit)
        {
         m_loCFIVS[idir][lit()].define(m_domain, m_grids.get(lit()), 
                                       m_grids, idir,Side::Lo);
         m_hiCFIVS[idir][lit()].define(m_domain, m_grids.get(lit()), 
                                       m_grids, idir,Side::Hi);
        }
    }
}

//
//reftofine_a is the refinement between this operator and opfine_a
//

void 
ParentLevelOp::define(const LevelOp* a_opfine,
                      int a_reftoFine)
{
  ParentLevelOp::clearMemory();
  const ParentLevelOp* opfineptr = 
    dynamic_cast<const ParentLevelOp*>(a_opfine);
  if(opfineptr == NULL)
    MayDay::Error("ParentLevelOp::define(): casting failed");
  const ParentLevelOp& opfine = *opfineptr; 
  CH_assert(opfine.isDefined());
  CH_assert(a_reftoFine > 0);
  m_isDefined = true;
  m_ihcfiEnabled = false;
  setDomainGhostBC(opfine.m_domghostbc);

  m_dxLevel = (opfine.m_dxLevel)*a_reftoFine;
  m_domain = coarsen(opfine.m_domain,a_reftoFine);
  coarsen(m_grids, opfine.m_grids, a_reftoFine);
  m_exchangeCopier.define(m_grids, m_grids, IntVect::Unit);
  m_dxCrse = opfine.m_dxCrse;
  m_ncomp = opfine.m_ncomp;

  //this refratio is the ratio between this level and the
  //next level in the amr hierarchy
  m_refRatio = -1;
  //don't need to define the cfinterp since inhomogeneous cfinterpolaton
  //is never done here.  We shall leave m_quadCFI and base_ba undefined and
  //and just leave the refratio -1.  We DO need fineinterp_ivs for doing
  //homogeneous cfinterpolation
  DataIterator lit = m_grids.dataIterator();
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      m_loCFIVS[idir].define(m_grids);
      m_hiCFIVS[idir].define(m_grids);
      for(lit.begin(); lit.ok(); ++lit)
        {
         m_loCFIVS[idir][lit()].define(m_domain, m_grids.get(lit()), 
                                       m_grids, idir,Side::Lo);
         m_hiCFIVS[idir][lit()].define(m_domain, m_grids.get(lit()), 
                                       m_grids, idir,Side::Hi);
        }
    }
  setBottomSmoother(*(opfine.m_bottom_smoother_ptr));
}

/***********************/
//  apply coarse-fine boundary conditions -- assume that phi grids
//  are grown by one
/***********************/
void
ParentLevelOp::CFInterp(LevelData<FArrayBox>& a_phif,
                    const LevelData<FArrayBox>& a_phic)
{
  CH_assert(isDefined());
  CH_assert(m_ihcfiEnabled);
  CH_assert( a_phif.ghostVect() >= IntVect::Unit);
  CH_assert(m_quadCFI.isDefined());
  m_quadCFI.coarseFineInterp(a_phif, a_phic);
}


/***********************/
// does homogeneous coarse/fine interpolation
/***********************/
void
ParentLevelOp::homogeneousCFInterp(LevelData<FArrayBox>& a_phif)
{
  CH_assert(isDefined());
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

/***********************/
// does homogeneous coarse/fine interpolation
/***********************/
void
ParentLevelOp::homogeneousCFInterp(LevelData<FArrayBox>& a_phif, 
                               const DataIndex& a_datInd, 
                               int a_idir, 
                               Side::LoHiSide a_hiorlo)
{
  CH_assert(isDefined());
  CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
  CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
  CH_assert( a_phif.ghostVect() >= IntVect::Unit);
  //  CH_assert (m_ncomp == a_phif.nComp());

  const CFIVS* cfivs_ptr = NULL;
  if(a_hiorlo == Side::Lo)
    cfivs_ptr = &m_loCFIVS[a_idir][a_datInd];
  else
    cfivs_ptr = &m_hiCFIVS[a_idir][a_datInd];

  const IntVectSet& interp_ivs = cfivs_ptr->getFineIVS();
  if(!interp_ivs.isEmpty())
    {
      int ihilo = sign(a_hiorlo);
      Box phistarbox = interp_ivs.minBox();
      phistarbox.shift(a_idir, ihilo);
      FArrayBox phistar(phistarbox, m_ncomp);
      //hence the homogeneous...
      phistar.setVal(0.);

      //given phistar, interpolate on fine ivs
      interpOnIVS(a_phif, phistar, 
                  a_datInd, a_idir, a_hiorlo,
                  interp_ivs);
    }
}

void
ParentLevelOp::interpOnIVS(LevelData<FArrayBox>& a_phif, 
                       const FArrayBox& a_phistar, 
                       const DataIndex& a_datInd,
                       const int a_idir,
                       const Side::LoHiSide a_hiorlo,
                       const IntVectSet& a_interpIVS)
{
  CH_assert(isDefined());
  CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
  CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
  CH_assert( a_phif.ghostVect() >= IntVect::Unit);
  IVSIterator fine_ivsit(a_interpIVS);
  FArrayBox& a_phi = a_phif[a_datInd];
  CH_assert (m_ncomp = a_phi.nComp());

  Real x1 = m_dxLevel;
  Real x2 = 0.5*(3.*m_dxLevel+m_dxCrse);
  Real denom = 1.0-((x1+x2)/x1);
  Real idenom = 1/(denom); // divide is more expensive usually
  Real x = 2.*m_dxLevel;
  Real xsquared = x*x;

  Real m1 = 1/(x1*x1);
  Real m2 = 1/(x1*(x1-x2));

  Real q1 = 1/(x1-x2);
  Real q2 = x1+x2;

  int ihilo = sign(a_hiorlo);
  IntVect ai = -2*ihilo*BASISV(a_idir);
  IntVect bi = -  ihilo*BASISV(a_idir);
  IntVect ci =    ihilo*BASISV(a_idir);


  for(fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit)
    {
      IntVect ivf = fine_ivsit();
      // quadradic interpolation 
      for(int ivar = 0; ivar < m_ncomp; ivar++)
        {
          Real pa =      a_phi(ivf + ai, ivar);
          Real pb =      a_phi(ivf + bi, ivar);
          Real pc =  a_phistar(ivf + ci, ivar);
          //phi = ax**2 + bx + c, x = 0 at pa
          Real a = (pb-pa)*m1 - (pb-pc)*m2;
          a *= idenom;
          Real b = (pb-pc)*q1 - a*q2;
          Real c = pa;
          a_phi(ivf,ivar) = a*xsquared + b*x + c;
        } //end loop over components
    } //end loop over fine intvects
}

/***********************/
// This does a GSRB Pre/Conditioned BiCGStab on a level
// for the bottom solver.  
/***********************/


void 
ParentLevelOp::bottomSmoother(LevelData<FArrayBox>& a_phi, 
                          const LevelData<FArrayBox>& a_rhs)
{
  CH_assert(isDefined());
  CH_assert(m_bottom_smoother_ptr != NULL);
  m_bottom_smoother_ptr->doBottomSmooth(a_phi, a_rhs, this);
}


void 
ParentLevelOp::setConvergenceMetric(Real a_metric, int a_comp)
{
  m_bottom_smoother_ptr->setConvergenceMetric(a_metric, a_comp);
}
