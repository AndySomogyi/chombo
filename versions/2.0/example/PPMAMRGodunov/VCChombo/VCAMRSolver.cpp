#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#include "VCAMRSolver.H"
#include "AMRLevelTools.H"
#include "LayoutIterator.H"
#include "parstream.H"
using std::cout;
using std::endl;

/// Implementation of class VCAMRSolver (Yuri Omelchenko, 08.20.02)
/*****************/
/*****************/
void
VCAMRSolver::setDefaultValues()
{
m_verbose = 0;
m_numAMRVCycle = 2;
}
/*****************/
/*****************/
void 
VCAMRSolver::clear()
{
}
/*****************/
/*****************/
VCAMRSolver::VCAMRSolver() 
{
  setDefaultValues();
}
/*****************/
/*****************/
VCAMRSolver::~VCAMRSolver()
{
  clear();
}
/*****************/
/*****************/
void
VCAMRSolver::define(const Vector<DisjointBoxLayout>& a_gridsLevel,
                    const Vector<Box>& a_domainLevel,
                    const Vector<Real>& a_dxLevel,
                    const Vector<int>& a_refRatio,
                    int a_numLevels,
                    int a_lBase,
                    const Vector<LevelOp*>& a_opinv,
                    int a_ncomp) 
{
  Vector<ProblemDomain> physdomains(a_domainLevel.size());

  for(int lev=0; lev<physdomains.size(); lev++)
     physdomains[lev] = ProblemDomain(a_domainLevel[lev]);

  define(a_gridsLevel, physdomains, a_dxLevel, 
         a_refRatio, a_numLevels, a_lBase, a_opinv, a_ncomp);
}

/*****************/
/*****************/
void
VCAMRSolver::define(const Vector<DisjointBoxLayout>& a_gridsLevel,
                    const Vector<ProblemDomain>& a_domainLevel,
                    const Vector<Real>& a_dxLevel,
                    const Vector<int>& a_refRatio,
                    int a_numLevels,
                    int a_lBase,
                    const Vector<LevelOp*>& a_opinv, 
                    int a_ncomp)
{
  clear();
  setDefaultValues();
  m_isDefined = true;

  m_gridsLevel = a_gridsLevel;
  m_domainLevel = a_domainLevel;
  m_numLevels = a_numLevels;
  m_finestLevel = a_numLevels - 1;
  m_refRatio = a_refRatio;
  m_dxLevel = a_dxLevel;
  m_lBase = a_lBase;
  m_ncomp = a_ncomp;
  ///define amrlevelmgs
  m_amrmgLevel.resize(m_numLevels, NULL);
  // (dfm -- 12/10/01) if lBase > 1, then we only need to define levels 
  // from lBase-1 (needed for c-f boundary conditions)
  int startLev= a_lBase;
  //if (startLev>0) startLev--;

  for(int ilev = startLev; ilev < m_numLevels; ilev++)
     m_amrmgLevel[ilev]=new AMRLevelMG(this,ilev,a_opinv[ilev],m_ncomp);

  CH_assert(m_lBase >= 0);

  //define levelsolver stuff
  const DisjointBoxLayout& levsolvGrids= m_gridsLevel[m_lBase];
  const ProblemDomain&     levsolvDom  = m_domainLevel[m_lBase];
  const Real&     levsolvDx   = m_dxLevel[m_lBase];
  const DisjointBoxLayout* levsolvBase = NULL;
  int levsolvRef = 2;
  if(m_lBase > 0) 
    {
      levsolvBase = &m_gridsLevel[m_lBase-1];
      levsolvRef = m_refRatio[m_lBase-1];
    }
  m_levelSolver.define(levsolvGrids, levsolvBase, 
                       levsolvDom,   levsolvDx, 
                       levsolvRef,   a_opinv[m_lBase], m_ncomp);
  m_levelSolver.setMaxIter(m_numVCyclesBottom);
}

/*****************/
/*****************/
///
/// complete, 3-level operator with homogeneous BCs
void 
VCAMRSolver::applyAMROperatorHphys(Vector<LevelData<FArrayBox>* >& a_phiLevel, 
                                   LevelData<FArrayBox> & a_LofPhi,
                                   int a_ilev)
{
  CH_assert(isDefined());
  CH_assert(a_ilev <= m_finestLevel);
  CH_assert(a_ilev >= 0);
  CH_assert(a_phiLevel.size() > m_finestLevel);
  m_amrmgLevel[a_ilev]->applyAMROperatorHphys(a_phiLevel,a_LofPhi);

}

/*****************/
/*****************/
///
///

void 
VCAMRSolver::applyPreconditioner(Vector<LevelData<FArrayBox>* >& a_phiv, 
                           const Vector<LevelData<FArrayBox>* >& a_rhsv,
                                 Vector<LevelData<FArrayBox>* >& a_corrv,
                                 Vector<LevelData<FArrayBox>* >& a_residv)
{
  for(int ilev=m_lBase; ilev <= m_finestLevel; ilev++)
    {
    DataIterator dit = (*a_phiv[ilev]).dataIterator();
    for(dit.begin(); dit.ok(); ++dit)
       {
       (*a_phiv[ilev])[dit()].setVal(0.);
       (*a_corrv[ilev])[dit()].setVal(0.);
       (*a_residv[ilev])[dit()].setVal(0.);
       }
    }
   for(int pass=0; pass < m_numAMRVCycle; pass++)
   { 
   computeResidualNorm(a_residv,a_phiv,a_rhsv,0); 

   // compute correction
   AMRVCycleMG(a_corrv,a_residv); 

   // correct solution, zero out correction
   for (int ilev = m_lBase; ilev <= m_finestLevel; ilev++)
        {
         DataIterator dit = a_phiv[ilev]->dataIterator();
         for (dit.begin(); dit.ok(); ++dit)
            {
            (*a_phiv[ilev])[dit()] += (*a_corrv[ilev])[dit()];
            (*a_corrv[ilev])[dit()].setVal(0.);
            }
        }
    } // end pass
}

/**
    Calculate norm of multilevel data
    Does not include data covered by finer levels.
    */
Vector<Real>
VCAMRSolver::computeNorm(const Vector<LevelData<FArrayBox> *>& a_arrv,
                         int a_normType)
{
  CH_assert(isDefined());
  CH_assert(a_arrv.size() > m_finestLevel);
  CH_assert((a_normType >= 0) && (a_normType <= 2));

  int ncomp = a_arrv[m_lBase]->nComp();
  if (ncomp > m_ncomp) 
    {
    MayDay::Warning("AMRSolver::computeNorm -- arr.nComp > solver.ncomp");
    ncomp = m_ncomp;
    }

  Vector<Real> normTot(ncomp,0.);
  Vector<Real> normLevel(ncomp,0.);

  for (int ilev = m_finestLevel; ilev >= m_lBase; ilev--) 
    {
      normLevel = m_amrmgLevel[ilev]->computeNorm(*a_arrv[ilev],a_normType);
      for (int comp=0; comp<ncomp; comp++)
        {
          if (a_normType == 0) 
            {
              normTot[comp] = Max(normTot[comp], normLevel[comp]);
            } 
          else 
            {
              normTot[comp] += normLevel[comp];
            }
        }
    } // end loop over levels

  if (a_normType == 2) 
    for (int comp=0; comp<ncomp; comp++)
         normTot[comp] = sqrt(normTot[comp]);

  return normTot;
}
/*****************/
/*****************/
///
/// This does a AMRVCycleMG Pre/Conditioned BiCGStab on all levels
///

void 
VCAMRSolver::solveBiCGStab(Vector<LevelData<FArrayBox>* >& a_phiv, 
                     const Vector<LevelData<FArrayBox>* >& a_rhsv)
{
  int maxNumberPass = 2;
  Real small= 1e-60, 
       small2= 1e-8;
  Real fault= 1.0-m_operatorTolerance;

  CH_assert(isDefined());
  CH_assert(a_phiv.size() > m_finestLevel);
  CH_assert(a_rhsv.size() > m_finestLevel);

  int ncomp = a_phiv[m_lBase]->nComp();
  if (ncomp > m_ncomp) 
    {
    MayDay::Warning("AMRSolver::solveBiCGStab -- phi.nComp > solver.ncomp");
    ncomp = m_ncomp;
    }

  //initial guess at the solution is zero
  for(int lev = m_lBase; lev <= m_finestLevel; lev++)
    {
      LevelData<FArrayBox>& a_phi = *a_phiv[lev]; 
      DataIterator dit = a_phi.dataIterator();
      for(dit.reset(); dit.ok(); ++dit)
         a_phi[dit()].setVal(0.);
    }

  // allocate BiCGStab work arrays
  Vector<LevelData<FArrayBox>* > corremfv(m_numLevels);
  Vector<LevelData<FArrayBox>* > residmfv(m_numLevels);
  Vector<LevelData<FArrayBox>* > rtwidmfv(m_numLevels);
  Vector<LevelData<FArrayBox>* > shatmfv(m_numLevels);
  Vector<LevelData<FArrayBox>* > phatmfv(m_numLevels);
  Vector<LevelData<FArrayBox>* > pmfv(m_numLevels);
  Vector<LevelData<FArrayBox>* > vmfv(m_numLevels);
  Vector<LevelData<FArrayBox>* > smfv(m_numLevels);
  Vector<LevelData<FArrayBox>* > tmfv(m_numLevels);

  // these are needed in the preconditioner
  Vector<LevelData<FArrayBox>* > t_corrv(m_numLevels);
  Vector<LevelData<FArrayBox>* > t_residv(m_numLevels);

  // coarsened finer grids for AMR Dot
  Vector<DisjointBoxLayout*> gridsCFPtrv(m_numLevels);

  for(int lev=m_lBase; lev <= m_finestLevel; lev++)
    {
    const DisjointBoxLayout& grids= m_gridsLevel[lev]; 
    const IntVect unit= IntVect::Unit;
    const IntVect zero= IntVect::Zero;
    
    corremfv[lev]= new LevelData<FArrayBox>(grids,ncomp,unit);
    if(corremfv[lev] == NULL) MayDay::Error("Out of Memory in solveBiCGStab");

    residmfv[lev]= new LevelData<FArrayBox>(grids,ncomp,zero);
    if(residmfv[lev] == NULL) MayDay::Error("Out of Memory in solveBiCGStab");

    rtwidmfv[lev]= new LevelData<FArrayBox>(grids,ncomp,zero);
    if(rtwidmfv[lev] == NULL) MayDay::Error("Out of Memory in solveBiCGStab");

    shatmfv[lev]= new LevelData<FArrayBox>(grids,ncomp,unit);
    if(shatmfv[lev] == NULL) MayDay::Error("Out of Memory in solveBiCGStab");

    phatmfv[lev]= new LevelData<FArrayBox>(grids,ncomp,unit);
    if(phatmfv[lev] == NULL) MayDay::Error("Out of Memory in solveBiCGStab");

    pmfv[lev]= new LevelData<FArrayBox>(grids,ncomp,zero);
    if(pmfv[lev] == NULL) MayDay::Error("Out of Memory in solveBiCGStab");

    vmfv[lev]= new LevelData<FArrayBox>(grids,ncomp,zero);
    if(vmfv[lev] == NULL) MayDay::Error("Out of Memory in solveBiCGStab");

    smfv[lev]= new LevelData<FArrayBox>(grids,ncomp,zero);
    if(smfv[lev] == NULL) MayDay::Error("Out of Memory in solveBiCGStab");

    tmfv[lev]= new LevelData<FArrayBox>(grids,ncomp,zero);
    if(tmfv[lev] == NULL) MayDay::Error("Out of Memory in solveBiCGStab");

    t_corrv[lev]= new LevelData<FArrayBox>(grids,ncomp,unit);
    if(t_corrv[lev] == NULL) MayDay::Error("Out of Memory in solveBiCGStab");

    t_residv[lev]= new LevelData<FArrayBox>(grids,ncomp,zero);
    if(t_residv[lev] == NULL) MayDay::Error("Out of Memory in solveBiCGStab");

    if(lev == m_finestLevel)
      {
      gridsCFPtrv[lev]= NULL;   // no finer grids to coarsen
      }
    else
      {
      gridsCFPtrv[lev]= new DisjointBoxLayout();
      if(gridsCFPtrv[lev]==NULL) MayDay::Error("Out of Memory in solveBiCGStab");
      coarsen(*gridsCFPtrv[lev], m_gridsLevel[lev+1], m_refRatio[lev]);
      }
    } // endfor lev

  // initialize all arrays to zero
  for(int lev=m_finestLevel; lev >= m_lBase; lev--)
    {
    const LevelData<FArrayBox>& a_rhs= *a_rhsv[lev];
    LevelData<FArrayBox>& corremf= *corremfv[lev];
    LevelData<FArrayBox>& residmf= *residmfv[lev];
    LevelData<FArrayBox>& rtwidmf= *rtwidmfv[lev];
    LevelData<FArrayBox>& shatmf= *shatmfv[lev];
    LevelData<FArrayBox>& phatmf= *phatmfv[lev];
    LevelData<FArrayBox>& pmf= *pmfv[lev];
    LevelData<FArrayBox>& vmf= *vmfv[lev];
    LevelData<FArrayBox>& smf= *smfv[lev];
    LevelData<FArrayBox>& tmf= *tmfv[lev];

    DataIterator dit = residmf.dataIterator();
    for(dit.begin(); dit.ok(); ++dit)
      {
      corremf[dit()].setVal(0.0);
      residmf[dit()].setVal(0.0);
      rtwidmf[dit()].setVal(0.0);
      shatmf[dit()].setVal(0.0);
      phatmf[dit()].setVal(0.0);
      pmf[dit()].setVal(0.0);
      vmf[dit()].setVal(0.0);
      smf[dit()].setVal(0.0);
      tmf[dit()].setVal(0.0);

      //(*t_corrv[lev])[dit()].setVal(0.0);
      //(*t_residv[lev])[dit()].setVal(0.0);
      }
    } // enfor lev

  Vector<Real> bnorm(ncomp,0.);
  Vector<Real> rnorm(ncomp,0.), oldnorm(ncomp,0.); 

  // calculate norm(rhs)
  bnorm = this->computeNorm(a_rhsv,0);

 // not sure if this is the right thing
 //  for(int comp=0; comp<ncomp; comp++)
 //    if(bnorm[comp] == 0.) bnorm[comp] = 1.;

  int iter = 0; // iteration counter (per pass)
  int nBiCGStabPass = 0; // BiCGStab pass counter
  do
  {
  // compute residual residmf = rhs-LofPhi [lev]
  for(int lev=m_finestLevel; lev >= m_lBase; lev--)
     this->computeAMRResidual(a_phiv, a_rhsv, *residmfv[lev], lev); 

  // compute norm of residual 
  rnorm = this->computeNorm(residmfv,0);

  // finish if norm(residual)/norm(rhs) < m_tolerance
  // for all components
  bool finished = true;
  for(int comp=0; comp<ncomp; comp++)
    {
    if(nBiCGStabPass == 0 && bnorm[comp] == 0.)
       bnorm[comp] = rnorm[comp];
    if(rnorm[comp] > m_tolerance*bnorm[comp]) 
       finished = false;
    if(nBiCGStabPass == 0 && m_verbose)
       pout() << "BiCGStab solver: initNorm[" << comp
              << "] = " << rnorm[comp]
              << " rhsNorm[" << comp
              << "] = " << bnorm[comp] << endl;
    }

  if(finished) break; // out of BiCGSTab pass loop => P.S.

  // BiCGStab coefficients for all components
  Vector<Real> rho1(ncomp,0.), rho2(ncomp,0.);
  Vector<Real> alpha(ncomp,0.), beta(ncomp,0.), omega(ncomp,0.);

  // not sure if this is the right thing to do
  // option 1: if any rhoDots are bad, add in current correction and exit
  // option 2: keep going with other components? -- not implemented

  iter = 1; // iteration counter (per pass)
  for(; (iter <= m_maxIter) && !finished; iter++)
    {
    oldnorm = rnorm;

    if(iter == 1)  
      {
       // r~ = r(0)
       for(int lev=m_lBase; lev <= m_finestLevel; lev++)
          {
           Interval intComp = residmfv[lev]->interval();
           residmfv[lev]->copyTo(intComp, *rtwidmfv[lev], intComp);
           }
      }
 
    bool badDotProd = false;
    for(int comp=0; comp<ncomp; comp++) 
       {
       rho1[comp] = 0.;
       Interval compInt(comp,comp);
       // rho(i-1) = Dot(r,r~)
       for(int lev=m_lBase; lev <= m_finestLevel; lev++)
          {
          rho1[comp] += AMRDotProduct(*residmfv[lev], *rtwidmfv[lev], 
                                      m_gridsLevel[lev], gridsCFPtrv[lev],
                                      compInt);
          }
       if(Abs(rho1[comp]) < small) badDotProd = true;
       } // endfor comp

    if(badDotProd)
       {
       if(m_verbose)
          pout() << "VCAMRSolver::solveBiCGStab failed after nPass="
                 <<  nBiCGStabPass << " and iter= " << iter << endl;
       break; // out of iteration loop
       }

    // if we're restarting the computation, act as if it was
    // the first iteration
    if(iter == 1)
      {
        // p(i) = r(i-1)
        for(int lev=m_lBase; lev <= m_finestLevel; lev++)
          {
           Interval intComp = residmfv[lev]->interval();
           residmfv[lev]->copyTo(intComp,*pmfv[lev],intComp);
          }
      }
    else
      {
        // beta(i-1)=[rho(i-1)/rho(i-2)]*[alpha(i-1)/omega(i-1)]
        for(int comp=0; comp<ncomp; comp++)
          beta[comp]=(rho1[comp]/rho2[comp])*(alpha[comp]/omega[comp]);

        //p(i) = r(i-1) + beta(i-1)*(p(i-1) -omega(i-1)*v(i-1))
        for(int lev=m_lBase; lev <= m_finestLevel; lev++)
          {
              DataIterator dit = pmfv[lev]->dataIterator();
              for(dit.begin(); dit.ok(); ++dit)
                  {
                  FArrayBox& pfab = (*pmfv[lev])[dit()];
                  FArrayBox& vfab = (*vmfv[lev])[dit()];
                  FArrayBox& rfab = (*residmfv[lev])[dit()];
                  Box fabbox = m_gridsLevel[lev].get(dit());
                  FArrayBox temp(fabbox, 1);
                  for(int comp=0; comp<ncomp; comp++)
                    {
                    temp.copy(pfab,comp,0,1);
                    temp.plus(vfab,-omega[comp],comp,0,1);
                    pfab.copy(rfab,comp,comp,1);
                    pfab.plus(temp,beta[comp],0,comp,1);
                    } // endfor comp
                  } // endfor dit
          } // endfor lev
      }// endif iter ==  1 

      ////////////////////////////////////////////////////
      //solve M(phat) = p(i) => phat
      ////////////////////////////////////////////////////
      this->applyPreconditioner(phatmfv,pmfv,
                                t_corrv,t_residv);

      //v(i) = A(phat)
      for(int lev=m_finestLevel; lev >= m_lBase; lev--)
         this->applyAMROperatorHphys(phatmfv, *vmfv[lev], lev);

      //alpha(i) = rho(i-1)/(rtwid^T v(i))
      bool bad_correction = false;
      for(int comp=0; comp<ncomp; comp++) 
         {
         Interval onecomp(comp,comp);
         Real denom = 0.0;
         for(int lev=m_lBase; lev <= m_finestLevel; lev++)
            {
            denom += AMRDotProduct(*vmfv[lev], *rtwidmfv[lev],
                                   m_gridsLevel[lev], gridsCFPtrv[lev],
                                   onecomp);
            } // endfor lev

         alpha[comp] = (Abs(denom) > small2*Abs(rho1[comp])) ?
                        rho1[comp]/denom : 0.;
         rho2[comp] = rho1[comp];
         if(alpha[comp] == 0.) bad_correction = true; 
         } // endfor comp

      if(bad_correction) break; // out of iteration loop

      //s = r(i-1) - alpha(i)*v(i)
      for(int lev=m_lBase; lev <= m_finestLevel; lev++)
          {
            DataIterator dit = vmfv[lev]->dataIterator();
            for(dit.begin(); dit.ok(); ++dit)
               {
               FArrayBox& rfab = (*residmfv[lev])[dit()];
               FArrayBox& vfab = (*vmfv[lev])[dit()];
               FArrayBox& sfab = (*smfv[lev])[dit()];
               for(int comp=0; comp<ncomp; comp++)
                 {
                 sfab.copy(rfab,comp,comp,1);
                 sfab.plus(vfab,-alpha[comp],comp,comp,1);
                 } 
               } // endfor dit 
          } // endfor lev 

      rnorm = this->computeNorm(smfv,0);
      finished = true;
      for(int comp=0; comp<ncomp; comp++) 
         if(rnorm[comp] > m_tolerance*bnorm[comp]) finished = false; 
        
      if(finished)
      {
      // increment correction by alpha(i)*phat and exit iteration
      for(int lev=m_lBase; lev <= m_finestLevel; lev++)
          {
            LevelData<FArrayBox>& corremf = *corremfv[lev];
            LevelData<FArrayBox>& phatmf = *phatmfv[lev];
            DataIterator dit = corremf.dataIterator();
            for(dit.begin(); dit.ok(); ++dit)
               for(int comp=0; comp<ncomp; comp++)
                 corremf[dit()].plus(phatmf[dit()],alpha[comp],comp,comp,1);
          } // endfor lev
       break; // out of iteration loop
      }
      ////////////////////////////////////////////////////
      //solve Mshat = s => shatmf 
      ////////////////////////////////////////////////////
      this->applyPreconditioner(shatmfv,smfv,
                                t_corrv,t_residv);
      //t(i) = A(shat)
      for(int lev=m_finestLevel; lev >= m_lBase; lev--)
         this->applyAMROperatorHphys(shatmfv, *tmfv[lev], lev);
      
      //w(i) = tTresid/tTt
      for(int comp=0; comp<ncomp; comp++) 
        {
        Real denomw = 0.0, numerw = 0.0;
        for(int lev=m_lBase; lev <= m_finestLevel; lev++)
           {
           Interval onecomp(comp,comp);
           numerw += AMRDotProduct(*tmfv[lev], *smfv[lev],
                                  m_gridsLevel[lev], gridsCFPtrv[lev],
                                  onecomp);
           denomw += AMRDotProduct(*tmfv[lev], *tmfv[lev],
                                  m_gridsLevel[lev], gridsCFPtrv[lev],
                                  onecomp);
           } // endfor lev
        omega[comp] = (Abs(denomw)>small) ? numerw/denomw : 0.0;
        } // endfor comp

      //corr(i) = corr(i-1) + alpha(i)*phat + omega(i)*shat
      //resid(i) = s - omega(i)*t
      for(int lev=m_lBase; lev <= m_finestLevel; lev++)
         {
         DataIterator dit = phatmfv[lev]->dataIterator();
         for(dit.begin(); dit.ok(); ++dit)
            {
            FArrayBox& corrfab = (*corremfv[lev])[dit()];
            FArrayBox& shatfab = (*shatmfv[lev])[dit()];
            FArrayBox& phatfab = (*phatmfv[lev])[dit()];
            FArrayBox& rfab = (*residmfv[lev])[dit()];
            FArrayBox& sfab = (*smfv[lev])[dit()];
            FArrayBox& tfab = (*tmfv[lev])[dit()];
            for(int comp=0; comp<ncomp; comp++) 
              {
              corrfab.plus(phatfab,alpha[comp],comp,comp,1);
              corrfab.plus(shatfab,omega[comp],comp,comp,1);
              rfab.copy(sfab,comp,comp,1);
              rfab.plus(tfab,-omega[comp],comp,comp,1);
              } // endfor comp
            } // endfor dit 
         } // endfor lev

       // now check residual.
       rnorm = this->computeNorm(residmfv,0);
       finished = true;
       for(int comp=0; comp<ncomp; comp++) 
          {
          if( (iter > m_minIter) && (rnorm[comp]>fault*oldnorm[comp]) )  
              {finished = true; break;}
          if( (rnorm[comp] > m_tolerance*bnorm[comp])
              && (Abs(omega[comp]) > m_tolerance) ) finished = false; 
          }
       if(m_verbose)
       for(int comp=0; comp<ncomp; comp++) 
          pout() << "BiCGStab pass= " << nBiCGStabPass << " iter= " << iter 
                 << " norm[" << comp << "] = " << rnorm[comp] << endl;
       if(finished) break; // not to increment iter 
    } // endfor main loop over BiCGStab iterations
      
    // correct solution and reset correction to zero
    for(int lev=m_lBase; lev <= m_finestLevel; lev++)
       {
          LevelData<FArrayBox>& corremf = *corremfv[lev];
          LevelData<FArrayBox>& a_phi = *a_phiv[lev];
          DataIterator dit = corremf.dataIterator();
          for(dit.begin(); dit.ok(); ++dit)
             {
             a_phi[dit()].plus(corremf[dit()],0,0,ncomp);
             corremf[dit()].setVal(0.);
             }
       } // endfor lev
  nBiCGStabPass++;
  }
  while(nBiCGStabPass <= maxNumberPass);

  if(m_verbose)  // check convergence criteria
    for (int comp=0; comp < ncomp; comp++)
      {
      if(rnorm[comp] <= m_tolerance*bnorm[comp]) 
          {
          pout() << "VCAMRSolver::solveBiCGStab OK:\n"; 
          }
      else if(rnorm[comp]>fault*oldnorm[comp])
          {
          pout() << "VCAMRSolver::solveBiCGStab REACHED HANG POINT:\n"; 
          }
      else
          {
          pout() << "VCAMRSolver::solveBiCGStab DIVERGED:\n"; 
          }
      pout() << "nPass= " << nBiCGStabPass << " iterations= " << iter 
             << " final rnorm[" << comp << "] = " << rnorm[comp]
             << endl;
      }

  // clear memory
  for(int lev=m_lBase; lev <= m_finestLevel; lev++)
    {
    delete corremfv[lev];
    delete residmfv[lev];
    delete rtwidmfv[lev];
    delete shatmfv[lev];
    delete phatmfv[lev];
    delete pmfv[lev];
    delete vmfv[lev];
    delete smfv[lev];
    delete tmfv[lev];

    delete t_corrv[lev];
    delete t_residv[lev];
    if(lev != m_finestLevel) delete gridsCFPtrv[lev];
    }
}

