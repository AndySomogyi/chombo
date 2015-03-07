#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#include "VCDivOp.H"	
#include "LayoutIterator.H"
#include "DotProduct.H"
#include "BoxIterator.H"
#include "VCDivOpF_F.H"   

//////////////////////////////////////////////////////
// Yuri Omelchenko, 08.20.02
//////////////////////////////////////////////////////
/// 
VCDivOp::VCDivOp()
{
  setDefaultValues();
}

///
VCDivOp::~VCDivOp() 
{
}

///
void
VCDivOp::setDefaultValues()
{
  ParentLevelOp::setDefaultValues();
  m_isVCDefined = false;
  m_isDiagSet = false;
  m_isDivSet = false;
  m_isLambdaSet = false;
}

///
void
VCDivOp::clearMemory()
{
  ParentLevelOp::clearMemory();
  m_isVCDefined = false;
  m_isDiagSet = false;
  m_isDivSet = false;
  m_isLambdaSet = false;
}

///
/// virtual constructor workaround :
/// Here's a problem. VCs are defined and initialized without
/// a new operator being defined in LevelOp's context ! 
///

LevelOp* 
VCDivOp::new_levelop() const
{
  CH_assert(m_isBCDefined);
  VCDivOp *newop = new VCDivOp();
  if(newop == NULL)
    {
      MayDay::Error("Out of Memory in VCDivOp::new_levelop");
    }
  newop->setDomainGhostBC(m_domghostbc);

  // need to define and copy newop's coefficients  
  if (m_isDiagSet && m_isDivSet) // if(m_isVCDefined) 
  {
  // do not use m_grids and m_ncomp since they may not be defined yet!
  const DisjointBoxLayout &grids = m_diagn.getBoxes();
  int ncomp = m_diagn.nComp();
  newop->defineVC(grids, ncomp);

  DataIterator dit= m_diagn.dataIterator();
  Interval intComp= m_diagn.interval();
  for(int idir=0; idir < SpaceDim; idir++)
     {
     m_lower[idir].copyTo(intComp,(newop->m_lower)[idir],intComp);
     m_centr[idir].copyTo(intComp,(newop->m_centr)[idir],intComp);
     m_upper[idir].copyTo(intComp,(newop->m_upper)[idir],intComp);
     }
  newop->m_isDivSet = true;

  m_diagn.copyTo(intComp,newop->m_diagn,intComp);
  newop->m_isDiagSet = true;

  newop->m_isLambdaSet = m_isLambdaSet;
  if (m_isLambdaSet)
    m_lambda.copyTo(intComp,newop->m_lambda,intComp);
  }

  if(m_bottom_smoother_ptr != NULL)
    {
      newop->setBottomSmoother(*m_bottom_smoother_ptr);
    }
  return static_cast<LevelOp*>(newop);
}

//////////////////////////////////////////////////////
//// other stuff
//////////////////////////////////////////////////////

bool
VCDivOp::isVCDefined() const
{
  return m_isVCDefined;
}

bool
VCDivOp::isDiagSet() const
{
  return m_isDiagSet;
}

bool
VCDivOp::isDivSet() const
{
  return m_isDivSet;
}

void
VCDivOp::needNewLambda()
{
  m_isLambdaSet = false;
}

/// view VC operator (coefficients)
void 
VCDivOp::viewOp(int idir, const Box& viewbox) const
{
  CH_assert(m_isVCDefined);
  CH_assert(m_isDiagSet);
  CH_assert(m_isDivSet);
  // can't do this because viewOp is const.
  // if (! m_isLambdaSet) computeInverseDiag();
  cout << "View VCs =" << " in dir=" << idir << endl;

  DataIterator dit= m_diagn.dataIterator();
     for(dit.begin(); dit.ok(); ++dit)
     {
      const FArrayBox& diagn = m_diagn[dit()];
      const FArrayBox& lower = m_lower[idir][dit()];
      const FArrayBox& centr = m_centr[idir][dit()];
      const FArrayBox& upper = m_upper[idir][dit()];
      const FArrayBox& lambda = m_lambda[dit()];

      Box vbox= viewbox & diagn.box();
      BoxIterator bit(vbox);
      for(bit.begin(); bit.ok(); ++bit)
          cout << bit() << " diagn=" << diagn(bit())
                        << " lower=" << lower(bit())
                        << " centr=" << centr(bit())
                        << " upper=" << upper(bit()) 
                        << " lambda="<< lambda(bit()) << endl;
     }
     cout << flush;
}

/// compute default (Poisson) coefficients
void 
VCDivOp::computeDefaultVC(Real a_dxLevel)
{
  CH_assert(isVCDefined());
  const Real dxinv2 = 1./(a_dxLevel*a_dxLevel);
  DataIterator dit= m_diagn.dataIterator();
  for(int idir=0; idir < SpaceDim; idir++)
     for(dit.begin(); dit.ok(); ++dit)
     {
     Real xl = 1.; // +0.*(drand48()-0.5);
     Real yl = 1.1; // +0.*(drand48()-0.5);
     m_lower[idir][dit()].setVal(xl*dxinv2);
     m_upper[idir][dit()].setVal(yl*dxinv2);
     m_centr[idir][dit()].setVal(0.);
     m_centr[idir][dit()].plus(m_lower[idir][dit()],-1.);
     m_centr[idir][dit()].plus(m_upper[idir][dit()],-1.);
     if(idir == 0) m_diagn[dit()].setVal(0.);
     }
  // this->computeInverseDiag(); // => m_lambda for GSRB
  m_isDiagSet = true;
  m_isDivSet = true;
  m_isLambdaSet = false;
}

/// compute arbitrary coefficients
void 
VCDivOp::computeVC(void (*a_levelVC)(VCDivOp& a_vcdivop))
{
  CH_assert(isVCDefined());
  a_levelVC(*this);
  // this->computeInverseDiag(); // => m_lambda for GSRB
  m_isDiagSet = true;
  m_isDivSet = true;
  m_isLambdaSet = false;
}

// set diagonal coefficients
void
VCDivOp::setDiagCoeff(const LevelData<FArrayBox>& a_diag)
{
  CH_assert(isVCDefined());
  CH_assert(a_diag.getBoxes() == m_diagn.getBoxes());
  a_diag.copyTo(a_diag.interval(), m_diagn, m_diagn.interval());

  // computeInverseDiag();  // (compute m_lambda)
  m_isDiagSet = true;
  m_isLambdaSet = false; // lambda must be recalculated
}

// set div coefficients (face centered beta for div(beta grad)phi
void
VCDivOp::setDivCoeff(const LevelData<FluxBox>& a_coeff,
                     Real a_dxLevel)
{
  CH_assert(isVCDefined());
  CH_assert(a_coeff.getBoxes() == m_centr[0].getBoxes());

  // this is a bit odd, but bear with me. Face-centered coefficients
  // are both lower and upper coefficients, but shifted properly
  
  // first loop over grids
  DataIterator dit = a_coeff.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) 
    {
      const FluxBox& thisCoeff = a_coeff[dit()];
      for (int dir=0; dir<SpaceDim; dir++) 
        {
          const FArrayBox& thisCoeffDir = thisCoeff[dir];
          FArrayBox& lowerCoeff = m_lower[dir][dit()];
          FArrayBox& upperCoeff = m_upper[dir][dit()];
          // need to shift lower and upper so that they match
          // locations of face-centered coeff
          lowerCoeff.shiftHalf(dir,-1);
          upperCoeff.shiftHalf(dir,1);
          // now copy coefficients
          lowerCoeff.copy(thisCoeffDir);
          upperCoeff.copy(thisCoeffDir);
          // now shift back
          lowerCoeff.shiftHalf(dir,1);
          upperCoeff.shiftHalf(dir,-1);

          // need to divide by dx^2
          Real dxInvSqr = 1./(a_dxLevel*a_dxLevel);
          lowerCoeff *= dxInvSqr;
          upperCoeff *= dxInvSqr;

          // now define center
          m_centr[dir][dit()].setVal(0.);
          m_centr[dir][dit()].plus(m_lower[dir][dit()],-1.);
          m_centr[dir][dit()].plus(m_upper[dir][dit()],-1.);
          
        }
    }
  // computeInverseDiag();  // (compute m_lambda)
  m_isDivSet = true;
  m_isLambdaSet = false; // lambda must be recalculated
}

  

///
void
VCDivOp::defineVC(const DisjointBoxLayout& a_grids, int a_ncomp)
{
  m_diagn.define(a_grids,a_ncomp);
  for(int idir=0; idir < SpaceDim; idir++)
    {
    m_lower[idir].define(a_grids,a_ncomp);
    m_upper[idir].define(a_grids,a_ncomp);
    m_centr[idir].define(a_grids,a_ncomp);
    }
  // private (scratch storage) 
  m_lofphi.define(a_grids,a_ncomp);
  m_lambda.define(a_grids,a_ncomp);

  m_isVCDefined = true; 
}

///
void
VCDivOp::define(
                  const DisjointBoxLayout& a_grids,
                  const DisjointBoxLayout* a_baseBAPtr,
                  Real  a_dxLevel,
                  int a_refRatio,
                  const Box& a_domain,
                  bool a_homogeneousOnly,
                  int a_ncomp)
{
  ProblemDomain probdomain(a_domain);
  VCDivOp::
  define(a_grids, a_baseBAPtr, a_dxLevel, a_refRatio, probdomain,
         a_homogeneousOnly, a_ncomp);
}

///
void
VCDivOp::define(
                  const DisjointBoxLayout& a_grids,
                  const DisjointBoxLayout* a_baseBAPtr,
                  Real  a_dxLevel,
                  int a_refRatio,
                  const ProblemDomain& a_probdomain,
                  bool a_homogeneousOnly,
                  int a_ncomp)
{

  ParentLevelOp::
  define(a_grids, a_baseBAPtr, a_dxLevel, a_refRatio, a_probdomain, 
         a_homogeneousOnly, a_ncomp);
}
//
void 
VCDivOp::define(const LevelOp *a_opfine, int a_reftoFine)
{
  CH_assert (this != a_opfine); 
//a_reftofine is the refinement between this operator and a_opfine
  const VCDivOp *opfinePtr = 
       dynamic_cast<const VCDivOp*>(a_opfine);
  if(opfinePtr == NULL)
       MayDay::Error("VCDivOp::define(a_opfine): casting failed");
  CH_assert(opfinePtr->isDefined() && opfinePtr->isVCDefined());

  // first define this as a coarse ParentLevelOp 
  ParentLevelOp::define(opfinePtr,a_reftoFine);
  // now define coarse coefficients (m_grids and m_ncomp have just been defined)
  this->defineVC(m_grids,m_ncomp);

/**
  now get coarse grid coefficients by averaging the fine values
  assuming piece-wise constant MG restriction/prolongation:
  du_coarse = const for all the overlying fine cells,
  {L,rhs_coarse} = (1/nref^SpaceDim) sum {L,rhs_fine} 
**/

  Box refbox(IntVect::Zero,
             (a_reftoFine-1)*IntVect::Unit);

  // coarsen diagonal coefficient (diagn*phi) 
  DataIterator dit= m_diagn.dataIterator();
  for(dit.begin(); dit.ok(); ++dit)
     {
     FORT_VCDIVOP2DEPOSITD(
                   CHF_FRA(m_diagn[dit()]),
                   CHF_CONST_FRA(opfinePtr->m_diagn[dit()]),
                   CHF_BOX(m_grids.get(dit())),
                   CHF_CONST_INT(a_reftoFine),
                   CHF_BOX(refbox));
     }
  m_isDiagSet = true;

  // coarsen other coefficients 
  for(int idir=0; idir < SpaceDim; idir++)
     for(dit.begin(); dit.ok(); ++dit)
     {
     FORT_VCDIVOP2DEPOSITLCU(
                   CHF_FRA(m_lower[idir][dit()]),
                   CHF_FRA(m_upper[idir][dit()]),
                   CHF_FRA(m_centr[idir][dit()]),
                   CHF_CONST_FRA(opfinePtr->m_lower[idir][dit()]),
                   CHF_CONST_FRA(opfinePtr->m_upper[idir][dit()]),
                   CHF_CONST_FRA(opfinePtr->m_centr[idir][dit()]),
                   CHF_BOX(m_grids.get(dit())),
                   CHF_CONST_INT(a_reftoFine),
                   CHF_BOX(refbox),
                   CHF_CONST_INT(idir));
    }
  // this->computeInverseDiag(); // => m_lambda for GSRB
  m_isDivSet = true;

  m_isLambdaSet = false; // lambda must be recalculated
}

///
/// compute approximate inverse (diagonal) operator
///
void 
VCDivOp::computeInverseDiag() 
{
  CH_assert(isVCDefined());
  CH_assert(isDiagSet() && isDivSet());
  DataIterator dit= m_diagn.dataIterator();
  for(int idir=0; idir < SpaceDim; idir++)
     for(dit.begin(); dit.ok(); ++dit)
     {
     FORT_VCDIVOP2INVERSEADD(CHF_FRA(m_lambda[dit()]),
                   CHF_CONST_FRA(m_diagn[dit()]),
                   CHF_CONST_FRA(m_centr[idir][dit()]),
                   CHF_BOX(m_diagn[dit()].box()),
                   CHF_CONST_INT(idir));
     }
  m_isLambdaSet = true;
}

/**
     Smoother. 
     Assumes that problem has already been put in 
     residual correction form, 
     so that C/F boundary conditions are homogeneous.
*/
void
VCDivOp::smooth(LevelData<FArrayBox>& a_phi,
                  const LevelData<FArrayBox>& a_rhs)
{
  CH_assert(isDefined() && isVCDefined());
  CH_assert(isDiagSet() && isDivSet());
  levelGSRB(a_phi,a_rhs);
}

/**
// this preconditioner first initializes phihat to 
// phihat = (IA)^-1 * rhshat and then smoothes it with a couple
// of passes of levelGSRB. It assumes homogeneous BC and CFI.
*/

void
VCDivOp::levelPreconditioner(LevelData<FArrayBox> & a_phihat,
                              const LevelData<FArrayBox>& a_rhshat)                    
{
// diagonal form of this operator is D == (IA)
// phi = {D^-1}*rhs

  CH_assert(isDefined() && isVCDefined());
  CH_assert(isDiagSet() && isDivSet());
  CH_assert(m_ncomp == a_rhshat.nComp());
  CH_assert(m_ncomp == a_phihat.nComp());

  if (! m_isLambdaSet) computeInverseDiag(); // get m_lambda

  DataIterator dit = a_phihat.dataIterator();
  for(dit.begin(); dit.ok(); ++dit)
      {
      a_phihat[dit()].copy(a_rhshat[dit()]);
      a_phihat[dit()].mult(m_lambda[dit()]);
      }
  levelGSRB(a_phihat, a_rhshat);
  levelGSRB(a_phihat, a_rhshat);
}

/***********************/
//  levelGSRB does GSRB on a level.  it has no knowledge of overlying
//  finer grids, and assumes coarse/fine BC's are homogeneous
/***********************/
void
VCDivOp::levelGSRB(LevelData<FArrayBox>& a_phi, 
                     const LevelData<FArrayBox>& a_rhs)
{
  CH_assert(isDefined() && isVCDefined());
  CH_assert(isDiagSet() && isDivSet());
  CH_assert(m_ncomp == a_phi.nComp());
  CH_assert(m_ncomp == a_rhs.nComp());

  if (! m_isLambdaSet) computeInverseDiag(); // get m_lambda

  // do first red, then black passes
  DataIterator dit = a_phi.dataIterator();
  for (int whichPass =0; whichPass <= 1; whichPass++) 
    {
      // do coarse/fine and copy bc's
      //should be done patch by patch.
      homogeneousCFInterp(a_phi);
      //fill in intersection of ghostcells and a_phi's boxes
      a_phi.exchange(a_phi.interval(), m_exchangeCopier);

      // now step through grids...
      for(dit.begin(); dit.ok(); ++dit)
        {
#ifndef NDEBUG
          FArrayBox& thisPhi = a_phi[dit()];
          const FArrayBox& thisRhs = a_rhs[dit()];
          const Box& thisBox = m_grids.get(dit());
          CH_assert(thisRhs.box().contains(thisBox));
          CH_assert(thisPhi.box().contains(grow(thisBox,1)));
#endif
          // invoke physical BC's where necessary
          // can do this here (w/o a following exchange)
          // because we use a 5-pt stencil
          m_domghostbc.applyHomogeneousBCs(a_phi[dit()],
                                           m_domain,
                                           m_dxLevel);
          // evaluate dimensionally split operator 
          for(int idir=0; idir < SpaceDim; idir++)
          {
          FORT_VCDIVOP2LOPHIADD(CHF_FRA(m_lofphi[dit()]),
                        CHF_CONST_FRA(a_phi[dit()]),
                        CHF_CONST_FRA(m_diagn[dit()]),
                        CHF_CONST_FRA(m_lower[idir][dit()]),
                        CHF_CONST_FRA(m_centr[idir][dit()]),
                        CHF_CONST_FRA(m_upper[idir][dit()]),
                        CHF_BOX(m_grids.get(dit())),
                        CHF_CONST_INT(idir),
                        CHF_CONST_INT(whichPass));
          }

          // do GSRB relaxation
          FORT_VCDIVOP2GSRB(CHF_FRA(a_phi[dit()]),
                        CHF_CONST_FRA(a_rhs[dit()]),
                        CHF_CONST_FRA(m_lofphi[dit()]),
                        CHF_CONST_FRA(m_lambda[dit()]),
                        CHF_BOX(m_grids.get(dit())),
                        CHF_CONST_INT(whichPass));
        } // end loop through grids
    } // end loop through red-black

/*
    // do final coarse/fine interpolation and copy bc's
    homogeneousCFInterp(a_phi);
    //fill in intersection of ghostcells and a_phi's boxes
    a_phi.exchange(a_phi.interval(), m_exchangeCopier);
    for(dit.begin(); dit.ok(); ++dit)
        {
        m_domghostbc.applyHomogeneousBCs(a_phi[dit()],
                                           m_domain,
                                           m_dxLevel);
        }
*/
}

/***********************/
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
/***********************/
void
VCDivOp::applyOpI(LevelData<FArrayBox>& a_phi, 
                      const LevelData<FArrayBox>* a_phicPtr, 
                      LevelData<FArrayBox>& a_lofPhi) 
{
  CH_assert(isDefined() && isVCDefined());
  CH_assert(isDiagSet() && isDivSet());
  CH_assert(m_ihcfiEnabled);
  CH_assert (a_phi.nComp() == a_lofPhi.nComp());
  CH_assert (a_phi.nComp() == m_ncomp);

  if(a_phicPtr != NULL)
    {
      // apply C/F boundary conditions...
      CFInterp(a_phi,*a_phicPtr);
    }
  
  //fill in intersection of ghostcells and phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
/*
#ifndef NDEBUG
      FArrayBox& thisPhi = a_phi[dit()];
      FArrayBox& thisLPhi = a_lofPhi[dit()];
      const Box& thisBox = m_grids.get(dit());
      CH_assert(thisPhi.box().contains(grow(thisBox,1)));
      CH_assert(thisLPhi.box().contains(thisBox));
#endif
*/

      // can do this here (w/o a following exchange)
      // because we use a 5-pt stencil
      m_domghostbc.applyInhomogeneousBCs(a_phi[dit()], 
                                         m_domain, 
                                         m_dxLevel);

      // evaluate dimensionally split operator 
      const int fullop = -1;	// evaluate full operator
      for(int idir=0; idir < SpaceDim; idir++)
      {
      FORT_VCDIVOP2LOPHIADD(CHF_FRA(a_lofPhi[dit()]),
                    CHF_CONST_FRA(a_phi[dit()]),
                    CHF_CONST_FRA(m_diagn[dit()]),
                    CHF_CONST_FRA(m_lower[idir][dit()]),
                    CHF_CONST_FRA(m_centr[idir][dit()]),
                    CHF_CONST_FRA(m_upper[idir][dit()]),
                    CHF_BOX(m_grids.get(dit())),
                    CHF_CONST_INT(idir),
                    CHF_CONST_INT(fullop));
      }
    }
}
/***********************/
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
// homogeneous physical bcs
/***********************/
void
VCDivOp::applyOpIcfHphys(LevelData<FArrayBox>& a_phi, 
                             const LevelData<FArrayBox>* a_phicPtr, 
                             LevelData<FArrayBox>& a_lofPhi) 
{
  CH_assert(isDefined() && isVCDefined());
  CH_assert(isDiagSet() && isDivSet());
  CH_assert(m_ihcfiEnabled);
  CH_assert (a_phi.nComp() == a_lofPhi.nComp());
  CH_assert (a_phi.nComp() == m_ncomp);

  if(a_phicPtr != NULL)
    {
      // apply C/F boundary conditions...
      CFInterp(a_phi,*a_phicPtr);
    }

  //fill in intersection of ghostcells and phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);
  
  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
#ifndef NDEBUG
      FArrayBox& thisPhi = a_phi[dit()];
      FArrayBox& thisLPhi = a_lofPhi[dit()];
      const Box& thisBox = m_grids.get(dit());
      CH_assert(thisPhi.box().contains(grow(thisBox,1)));
      CH_assert(thisLPhi.box().contains(thisBox));
#endif

      // can do this here (w/o a following exchange)
      // because we use a 5-pt stencil
      m_domghostbc.applyHomogeneousBCs(a_phi[dit()], 
                                       m_domain, 
                                       m_dxLevel);

      // evaluate dimensionally split Helmholtz operator 
      const int fullop = -1;	// evaluate full operator
      for(int idir=0; idir < SpaceDim; idir++)
      {
      FORT_VCDIVOP2LOPHIADD(CHF_FRA(a_lofPhi[dit()]),
                    CHF_CONST_FRA(a_phi[dit()]),
                    CHF_CONST_FRA(m_diagn[dit()]),
                    CHF_CONST_FRA(m_lower[idir][dit()]),
                    CHF_CONST_FRA(m_centr[idir][dit()]),
                    CHF_CONST_FRA(m_upper[idir][dit()]),
                    CHF_BOX(m_grids.get(dit())),
                    CHF_CONST_INT(idir),
                    CHF_CONST_INT(fullop));
      }
    }
}
/***********************/
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
// homogeneous physical bcs
/***********************/
void
VCDivOp::applyOpH(LevelData<FArrayBox>& a_phi, 
                      LevelData<FArrayBox>& a_lofPhi)
{
  CH_assert(isDefined() && isVCDefined());
  CH_assert(isDiagSet() && isDivSet());
  CH_assert (a_phi.nComp() == a_lofPhi.nComp());
  CH_assert (a_phi.nComp() == m_ncomp);

  homogeneousCFInterp(a_phi);
  //fill in intersection of ghostcells and phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
#ifndef NDEBUG
      FArrayBox& thisPhi = a_phi[dit()];
      FArrayBox& thisLPhi = a_lofPhi[dit()];
      const Box& thisBox = m_grids.get(dit());
      CH_assert(thisPhi.box().contains(grow(thisBox,1)));
      CH_assert(thisLPhi.box().contains(thisBox));
#endif

      // can do this here (w/o a following exchange)
      // because we use a 5-pt stencil
      m_domghostbc.applyHomogeneousBCs(a_phi[dit()], 
                                       m_domain, 
                                       m_dxLevel);

      // evaluate dimensionally split operator 
      const int fullop = -1;	// evaluate full operator
      for(int idir=0; idir < SpaceDim; idir++)
      {
      FORT_VCDIVOP2LOPHIADD(CHF_FRA(a_lofPhi[dit()]),
                    CHF_CONST_FRA(a_phi[dit()]),
                    CHF_CONST_FRA(m_diagn[dit()]),
                    CHF_CONST_FRA(m_lower[idir][dit()]),
                    CHF_CONST_FRA(m_centr[idir][dit()]),
                    CHF_CONST_FRA(m_upper[idir][dit()]),
                    CHF_BOX(m_grids.get(dit())),
                    CHF_CONST_INT(idir),
                    CHF_CONST_INT(fullop));
      }
    }
}
/***********************/
// evaluate operator, homogeneous C/F boundary conditions
// note -- grids for phi need not correspond to bx
// inhomogeneous physical boundary conditions
/***********************/
void
VCDivOp::applyOpHcfIphys(LevelData<FArrayBox>& a_phi, 
                             LevelData<FArrayBox>& a_lofPhi)
{
  CH_assert(isDefined() && isVCDefined());
  CH_assert(isDiagSet() && isDivSet());
  CH_assert (a_phi.nComp() == a_lofPhi.nComp());
  CH_assert (a_phi.nComp() == m_ncomp);

  homogeneousCFInterp(a_phi);
  //fill in intersection of ghostcells and phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
#ifndef NDEBUG
      FArrayBox& thisPhi = a_phi[dit()];
      FArrayBox& thisLPhi = a_lofPhi[dit()];
      const Box& thisBox = m_grids.get(dit());
      CH_assert(thisPhi.box().contains(grow(thisBox,1)));
      CH_assert(thisLPhi.box().contains(thisBox));
#endif

      // can do this here (w/o a following exchange)
      // because we use a 5-pt stencil
      m_domghostbc.applyInhomogeneousBCs(a_phi[dit()], 
                                         m_domain, 
                                         m_dxLevel);

      // evaluate dimensionally split operator 
      const int fullop = -1;	// evaluate full operator
      for(int idir=0; idir < SpaceDim; idir++)
      {
      FORT_VCDIVOP2LOPHIADD(CHF_FRA(a_lofPhi[dit()]),
                    CHF_CONST_FRA(a_phi[dit()]),
                    CHF_CONST_FRA(m_diagn[dit()]),
                    CHF_CONST_FRA(m_lower[idir][dit()]),
                    CHF_CONST_FRA(m_centr[idir][dit()]),
                    CHF_CONST_FRA(m_upper[idir][dit()]),
                    CHF_BOX(m_grids.get(dit())),
                    CHF_CONST_INT(idir),
                    CHF_CONST_INT(fullop));
      }
    }
}

/** 
    get flux( == flux at THIS level) 
    The fluxes live on the cell faces with direction dir.
    Fluxes are computed for all interior edges of data.
    The flux fab is resized inside the routine.
    
    We assume: 

    L[i] = d/dx Flux = (1/dx){Flux[i+1/2]-Flux[i-1/2]}  
    Flux[i+1/2] = (1/dx){alpha[i+1/2]*f[i]+beta[i+1/2]*f[i+1]} 
    for each component.

    Therefore:

    alpha[i+1/2] = -dx^2*lower[i+1]
    beta[i+1/2]  =  dx^2*upper[i]
    centr[i+1/2] = (alpha[i+1/2]-beta[i-1/2])/dx^2

*/
void 
VCDivOp::getFlux(
                  FArrayBox& a_fineFlux,
                  const FArrayBox& a_data,
	          const DataIndex& a_dit, 
                  int a_dir)
{
  CH_assert(isDefined() && isVCDefined());
  CH_assert(isDiagSet() && isDivSet());
  CH_assert(a_dir >= 0 && a_dir  < SpaceDim);
  CH_assert(!a_data.box().isEmpty());

// old stuff
// Box edgebox = surroundingNodes(a_data.box(),a_dir);
// edgebox.grow(a_dir, -1);

// new algorithm for flux calculation (Yuri Omelchenko/Dan Martin, 1.10.02)
// make a global flux box (needed to process domain edges correctly)
   Box edgebox= surroundingNodes(m_grids.get(a_dit),a_dir);

// clip data to prevent accessing undefined coefficients in the buffer cells
   Box clipbox= m_grids.get(a_dit);
   clipbox.grow(a_dir,1);
   clipbox &= a_data.box();

// form a flux calculation box
   Box fluxbox= surroundingNodes(clipbox,a_dir);
   fluxbox.grow(a_dir, -1);
   
//  cout << "a_dir="<<a_dir<<endl;
//  cout << "gridbox="<<m_grids.get(a_dit)<<endl;
//  cout << "databox="<<a_data.box()<<endl;
//  cout << "edgebox="<<edgebox<<endl;

// if this fails, the data box was too small (one cell wide, in fact)
  CH_assert(!fluxbox.isEmpty());
  a_fineFlux.resize(fluxbox, a_data.nComp());
  FORT_VCDIVOP2GETFLUX(CHF_FRA(a_fineFlux),
                CHF_CONST_FRA(a_data),
                CHF_CONST_FRA(m_lower[a_dir][a_dit]),
                CHF_CONST_FRA(m_centr[a_dir][a_dit]),
                CHF_CONST_FRA(m_upper[a_dir][a_dit]),
                CHF_BOX(fluxbox),
                CHF_BOX(edgebox),
                CHF_CONST_REAL(m_dxLevel),
                CHF_CONST_INT(a_dir));
}


