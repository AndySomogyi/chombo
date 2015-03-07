#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <iostream>
#include <fstream>
#include <string>
using std::fstream;
using std::string;
#include <cstdio>

#include "localFuncs.H"
#include "parstream.H"
#include "BoxIterator.H"
#include "TensorProbF_F.H"
#include "LoadBalance.H"
#include "AMRIO.H"
#include "PoissonBC.H"
#include "BRMeshRefine.H"
#include "computeNorm.H"
#include "LayoutIterator.H"
#include "UsingNamespace.H"

int outputData(Vector<LevelData<FArrayBox>* >& vectPhi,
               const Vector<LevelData<FArrayBox>* >& vectRhs,
               const Vector<DisjointBoxLayout>& vectGrids,
               const Vector<ProblemDomain>& vectDomain,
               const Vector<int>& vectRatio,
               const Vector<Real>& vectDx,
               AMRSolver& amr_solver,
               int numlevels,
               int lBase,
               Vector<Real>& L1Old,
               Vector<Real>& L2Old,
               Vector<Real>& MaxOld,
               bool verbose)
{
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  int ncomp = vectPhi[0]->nComp();

  int numCompPerVariable = 4;

  string phiNameBase("phi");
  string rhsNameBase("rhs");
  string residNameBase("resid");
  Vector<string> vectName(numCompPerVariable*ncomp);
  for (int comp=0; comp<ncomp; comp++) 
    {
      char labelString[80];
      sprintf(labelString, "%s%d", phiNameBase.c_str(), comp);
      string phiName(labelString);
      vectName[comp] = phiName;
      sprintf(labelString, "%s%d", rhsNameBase.c_str(), comp);
      string rhsName(labelString);
      vectName[ncomp+comp] = rhsName;
      sprintf(labelString, "%s%d", "error", comp);
      string errName(labelString);
      vectName[2*ncomp+comp] = errName;
      sprintf(labelString, "%s%d", residNameBase.c_str(), comp);
      string residName(labelString);
      vectName[3*ncomp+comp] = residName;
    }
  Box domain = vectDomain[0].domainBox();
  Real dt = 1.;
  Real time = 1.;
  Vector<LevelData<FArrayBox>* > vectPhiAndRHS(numlevels, NULL);
  Vector<LevelData<FArrayBox>* > vectError(numlevels, NULL);
  Vector<LevelData<FArrayBox>* > vectResid(numlevels, NULL);
  for(int ilev = 0; ilev < numlevels; ilev++)
    {
      vectPhiAndRHS[ilev] = new LevelData<FArrayBox>(vectGrids[ilev], 
                                                     numCompPerVariable*ncomp);


      Interval phiInterval(0,ncomp-1);
      Interval rhsInterval(ncomp, 2*ncomp-1);
      vectPhi[ilev]->copyTo(vectPhi[ilev]->interval(), 
                            *vectPhiAndRHS[ilev],
                            phiInterval);
      vectRhs[ilev]->copyTo(vectRhs[ilev]->interval(), 
                            *vectPhiAndRHS[ilev],
                            rhsInterval);

      
      vectError[ilev] = new LevelData<FArrayBox>(vectGrids[ilev], ncomp);
      LevelData<FArrayBox>& levelError = *vectError[ilev];
      const DisjointBoxLayout* fineGridsPtr = NULL;
      if (ilev < (numlevels-1)) fineGridsPtr = &vectGrids[ilev+1];

      computeError(levelError, *vectPhi[ilev], fineGridsPtr,
                   vectDomain[0], vectDx[ilev], vectRatio[ilev]);
      Interval errorInterval(2*ncomp, 3*ncomp-1);
      levelError.copyTo(levelError.interval(), *vectPhiAndRHS[ilev],
                        errorInterval);

      // now compute residual
      vectResid[ilev] = new LevelData<FArrayBox>(vectGrids[ilev], ncomp);
      LevelData<FArrayBox>& levelResid = *vectResid[ilev];

      LevelData<FArrayBox>& levelRHS = *vectRhs[ilev];
      if (ilev >= lBase)
        {
          amr_solver.applyAMROperator(vectPhi, levelResid, ilev);
          DataIterator levelDit = levelResid.dataIterator();
          for (levelDit.begin(); levelDit.ok(); ++levelDit)
            {
              levelResid[levelDit()] -= levelRHS[levelDit()];
              levelResid[levelDit()].negate();
            }
        }
      else 
        {
          // if we're not computing on this level, set resid to 0
          DataIterator levelDit = levelResid.dataIterator();
          for (levelDit.begin(); levelDit.ok(); ++levelDit)
            {
              levelResid[levelDit()].setVal(0.0);
            }
        }
      Interval residInterval(3*ncomp, 4*ncomp-1);
      levelResid.copyTo(levelResid.interval(),*vectPhiAndRHS[ilev],
                        residInterval);

    }

  int gridIndex = vectDomain[0].domainBox().size(0);
  cout << "Size = " << gridIndex << ":" << endl;
  // also compute error here
  for (int comp=0; comp<ncomp; comp++) 
    {
      Real MaxErr=0, L1Err=0, L2Err=0;
      Real MaxRate=0, L1Rate=0, L2Rate=0;
      Interval comps(comp,comp);
      
      MaxErr = computeNorm(vectError, vectRatio, vectDx[0], comps, 0);
      L1Err = computeNorm(vectError, vectRatio, vectDx[0], comps, 1);
      L2Err = computeNorm(vectError, vectRatio, vectDx[0], comps, 2);

      Real log_2 = log(2.0);
      if (MaxOld[comp] != 0.0)
        MaxRate = log(MaxOld[comp]/MaxErr)/log_2;
      if (L1Old[comp] != 0.0)
        L1Rate = log(L1Old[comp]/L1Err)/log_2;
      if (L2Old[comp] != 0.0)
        L2Rate = log(L2Old[comp]/L2Err)/log_2;

      MaxOld[comp] = MaxErr;
      L1Old[comp] = L1Err;
      L2Old[comp] = L2Err;
      
      cout << "comp[" << comp << "]: Max(Err) = " << MaxErr 
           << " L1(Err) = " << L1Err 
           << " L2(Err) = " << L2Err << endl;

      cout << "convergence:       Max = " << MaxRate
           << " L1 = " << L1Rate 
           << " L2 = " << L2Rate << endl;


    }



#ifdef CH_USE_HDF5
  string filename("solverOut.hdf5");
  WriteAMRHierarchyHDF5(filename,                             
                        vectGrids, 
                        vectPhiAndRHS, 
                        vectName,                 
                        domain,                               
                        vectDx[0], dt, time,                      
                        vectRatio,                           
                        numlevels);
#endif

  for(int ilev = 0; ilev < numlevels; ilev++)
    {
      delete vectPhiAndRHS[ilev];
      delete vectError[ilev];
      delete vectResid[ilev];
    }

  return 0;
}

/*
  Set RHS on hierarchy from input file
 */
int setRHS(Vector<LevelData<FArrayBox>* >& vectRhs,
           Vector<Real>&                   vectDx, 
           Vector<ProblemDomain>&          vectDomain, 
           int numlevels, bool verbose)
{

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  Real rhono, rno;
  int iprob;

  // problem 1 is better for dirichlet problem, while 
  // problem 2 is for the periodic case
  if (vectDomain[0].isPeriodic())
    {
      iprob = 2;
    } else {
      iprob = 1;
      //iprob = 2;
    }

  rhono = 0.75;
  rno = 0.5;
  if(verbose)
    pout() 
      << " rhono  = " << rhono
      << " rno  = " << rno
      << " iprob  = " << iprob
      << endl;
  for(int ilev = 0; ilev < numlevels; ilev++)
    {
      if(verbose)
        pout() 
          << " ilev = " << ilev 
          << " dom = " << vectDomain[ilev]
          << " dx  = " << vectDx[ilev]
        << endl;
      Real dxlev = vectDx[ilev];
      Box domlev =vectDomain[ilev].domainBox();
      LevelData<FArrayBox>& rhsLD = *vectRhs[ilev];
      DataIterator dit =  rhsLD.dataIterator();
      for(dit.reset(); dit.ok(); ++dit)
        {
          FArrayBox& rhsFab = rhsLD[dit()];
          //kluge to get things the same in 0.2
          rhsFab.setVal(7.0);
          /**/
          FORT_GETRHSTENSOR(CHF_FRA(rhsFab),
                            CHF_BOX(rhsFab.box()),
                            CHF_BOX(domlev),
                            CHF_CONST_REAL(dxlev),
                            CHF_CONST_REAL(rhono),
                            CHF_CONST_REAL(rno),
                            CHF_CONST_INT(iprob));
        }
#ifdef CH_MPI
      MPI_Barrier(Chombo_MPI::comm);
#endif
    }
  return 0;
}

/*
  Set grid hierarchy from input file
 */
int setGrids(Vector<DisjointBoxLayout>& vectGrids,
             Vector<ProblemDomain>&     vectDomain, 
             Vector<Real>&              vectDx, 
             Vector<int>&               vectRefRatio, 
             int& numlevels, bool verbose,
             Real domainMult)
{
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  int max_level;
  Real fillRat = 0.77;
  std::vector<int> ancells(SpaceDim), refs;
  std::vector<Real>  prob_loa(SpaceDim);
  std::vector<Real>  prob_hia(SpaceDim);
  // set default to non-periodic
  std::vector<int> is_periodica(SpaceDim,0);
  bool is_periodic[SpaceDim];
  int temp;
  bool predefinedGrids = false;
  int maxboxsize = 32;
  Real refine_thresh = 0.5;

  ParmParse pp("main");
  pp.get("max_level", max_level);
  numlevels = max_level + 1;
  //kluge to get things the same in 0.2
  //max_level = 1;
  pp.query("fill_ratio", fillRat);
  pp.getarr("n_cell", ancells, 0, SpaceDim);
  pp.getarr("ref_ratio", refs, 0, numlevels);
  vectRefRatio = refs;
  pp.getarr("prob_lo",prob_loa,0,SpaceDim);
  pp.getarr("prob_hi",prob_hia,0,SpaceDim);
  pp.queryarr("is_periodic", is_periodica,0, SpaceDim);
  pp.get("maxboxsize",maxboxsize);

  int ncomp = 2;
  pp.query("ncomp", ncomp);

  numlevels = max_level + 1;

  temp = predefinedGrids;
  pp.query("presetGrids", temp);
  predefinedGrids = (temp == 1);

  pp.query("refineoment_threshold", refine_thresh);

  // grid generation parameters
  int maxrat = 2;
  for(int irat = 0; irat < vectRefRatio.size(); irat++)
    maxrat = Max(vectRefRatio[irat], maxrat);
  int blockFactor = maxrat*2;
  int nesting_radius = 1;


  IntVect ivlo = IntVect::Zero;
  IntVect ivhi;
  for(int idir = 0; idir < SpaceDim; idir++)
    ivhi[idir] = ancells[idir] - 1;

  Box basedom = Box(ivlo, ivhi);

  basedom.refine(domainMult);

  vectGrids.resize(numlevels);
  vectDomain.resize(numlevels);
  vectDx.resize(numlevels);


  vectDx.resize(numlevels,0.0);
  vectDx[0] = (prob_hia[0]-prob_loa[0])/ancells[0];
  vectDx[0] = vectDx[0]/domainMult;

  for(int ilev = 1; ilev < numlevels; ilev++)
    {
      CH_assert(vectRefRatio[ilev-1] > 0);
      vectDx[ilev] = vectDx[ilev-1]/vectRefRatio[ilev-1];
    }


  // convert periodic from int->bool
  for (int dim=0; dim<SpaceDim; dim++) 
    {
      is_periodic[dim] = (is_periodica[dim] == 1);
      if (is_periodic[dim] && verbose && procID() == 0)
        pout() << "Using Periodic BCs in direction: " << dim << endl;
    }

  vectDomain[0] = ProblemDomain(basedom,is_periodic);
  for(int ilev = 1;ilev < numlevels; ilev++)
    {
      vectDomain[ilev] = refine(vectDomain[ilev-1],vectRefRatio[ilev-1]);
    }

  int maxLevel = numlevels-1;
  Vector<Vector<Box> > newBoxes(numlevels);
  Vector<Vector<Box> > oldBoxes(numlevels);

  // use predefined grid configuration designed to test quadCFInterp configs
  if (predefinedGrids) 
    {
      int nc = ancells[0]*((int)domainMult);
      int nmi = nc/2;//16
      int nqu = nc/4;//8
      int ntf = (nc*3)/4;  //24
      int nte = (nc*3)/8; //12
      int nfe = (nc*5)/8; //20
#if(CH_SPACEDIM ==2)
      Box centeredBox(IntVect(nqu,nqu), IntVect(ntf-1, ntf-1));
      Box boxf1(IntVect(0, nqu), IntVect(nmi-1,ntf-1));
      Box boxf2(IntVect(nmi,nte), IntVect(ntf-1,nfe-1));
      Box boxf3(IntVect(nqu,0  ), IntVect(nfe-1,nqu-1));
      Box boxf4(IntVect(nfe,nqu), IntVect(nc -1,nte-1));
#else
      Box centeredBox(IntVect(nqu,nqu,nqu), IntVect(ntf-1,ntf-1,ntf-1));
      Box boxf1(IntVect(0, nqu,nqu), IntVect(nmi-1,ntf-1,ntf-1));
      Box boxf2(IntVect(nmi,nte,nte), IntVect(ntf-1,nfe-1,nfe-1));
      Box boxf3(IntVect(nqu,0,0  ), IntVect(nfe-1,nqu-1,nqu-1));
      Box boxf4(IntVect(nfe,nqu,nqu), IntVect(nc -1,nte-1,nte-1));
#endif
      //comment out for kluge 
      /**/
      IntVectSet tags;
      tags |= centeredBox;
      //tags |= boxf1;
      //tags |= boxf2;
      //tags |= boxf3;
      //tags |= boxf4;
      
      for(int ilev = 0; ilev <numlevels; ilev++)
        {
          oldBoxes[ilev].push_back(vectDomain[ilev].domainBox());
        }
      int baseLevel = 0;
      int topLevel  = numlevels - 2;
      int eekflag = 0;
      if(topLevel >= 0) 
        {
          int newFinestLev;
          BRMeshRefine meshrefine(vectDomain[0],vectRefRatio, fillRat,
                                  blockFactor, nesting_radius, maxboxsize);
          newFinestLev = meshrefine.regrid(newBoxes, tags, baseLevel, 
                                           topLevel, oldBoxes);
        } else {
          newBoxes = oldBoxes;
        }
      
      Vector< Vector<int> > procAssign;
      Real effRatio = 0.75;
      Vector< Vector<long> > loads(numlevels);
      for(int ilev = 0; ilev <numlevels; ilev++)
        {
          loads[ilev].resize(newBoxes[ilev].size());
          for(int ibox = 0; ibox < newBoxes[ilev].size() ; ibox++)
            {
              loads[ilev][ibox] = newBoxes[ilev][ibox].numPts();
            }
        }
      LoadBalance(procAssign, effRatio, newBoxes, loads, vectRefRatio);
      
      if(eekflag != 0)
        {
          cerr << "setGrids: loadBalance returned error code " 
               << eekflag << endl;
          return(1);
        }
      for(int ilev = 0; ilev <numlevels; ilev++)
        {
          vectGrids[ilev].define(newBoxes[ilev], procAssign[ilev],
                                 vectDomain[ilev]);
          vectGrids[ilev].close();
        }
      
    } else {
      // determine grids dynamically, based on grad(RHS)
      // will need temp storage for RHS
      Vector<LevelData<FArrayBox>* > vectRHS(maxLevel+1,NULL);

      // define base level first
      Vector< Vector<int> > procAssign(maxLevel+1);
      domainSplit(vectDomain[0], oldBoxes[0], maxboxsize, blockFactor);
      procAssign[0].resize(oldBoxes[0].size());
      LoadBalance(procAssign[0],oldBoxes[0]);
      
      vectGrids[0].define(oldBoxes[0],procAssign[0],vectDomain[0]);
      
      vectRHS[0] = new LevelData<FArrayBox>(vectGrids[0], ncomp, 
                                            IntVect::Zero);

      int topLevel = 0;

      bool moreLevels = (maxLevel > 0);
      
      // create grid generation object
      BRMeshRefine meshrefine(vectDomain[0], vectRefRatio, fillRat,
                              blockFactor, nesting_radius, maxboxsize);

      while (moreLevels) 
        {
          // default is moreLevels = false
          // (only repeat loop in the case where a new level 
          // is generated which is still less than maxLevel)
          moreLevels = false;

          int baseLevel = 0;
          int oldTopLevel = topLevel;

          // now initialize RHS for this existing hierarchy
          setRHS(vectRHS, vectDx, vectDomain, topLevel+1, verbose);

          Vector<IntVectSet> tagVect(topLevel+1);
          int tags_grow = 1;
          tagCells(vectRHS, tagVect, vectDx, vectDomain, refine_thresh,
                   tags_grow, baseLevel, topLevel+1);

          if (procID() == uniqueProc(SerialTask::compute))
            {
              if (topLevel >= 0 && ! tagVect[topLevel].isEmpty())
                {
                  
                  meshrefine.regrid(newBoxes, tagVect, baseLevel, 
                                    topLevel, oldBoxes);

                  topLevel++;

                } // end if there was a new level generated
            } // end if proc is serial node
          broadcast(topLevel, uniqueProc(SerialTask::compute) );
          broadcast(newBoxes, uniqueProc(SerialTask::compute) );

          oldBoxes = newBoxes;

          //  no need to do this for the base level (already done)
          for (int lev=1; lev<= topLevel; lev++) 
            {
              // do load balancing
              procAssign[lev].resize(newBoxes[lev].size());
              LoadBalance(procAssign[lev], newBoxes[lev]);
              const DisjointBoxLayout newDBL(newBoxes[lev], procAssign[lev],
                                             vectDomain[lev]);
              vectGrids[lev] = newDBL;
              if (vectRHS[lev] != NULL)
                {
                  delete vectRHS[lev];
                  vectRHS[lev] = NULL;
                }
              vectRHS[lev] = new LevelData<FArrayBox>(vectGrids[lev], ncomp,
                                                      IntVect::Zero);
            } // end loop over levels for initialization
          
          // figure out whether we need another pass through grid generation
          if ((topLevel<maxLevel) && (topLevel > oldTopLevel))
            moreLevels = true;

        } // end while moreLevels loop

    } // end if grids generated dynamically
  
  return 0;
}

/*
  tag cells for refinement based on magnitude(RHS)
*/
void
tagCells(Vector<LevelData<FArrayBox>* >& vectRHS,
         Vector<IntVectSet>& tagVect,
         Vector<Real>& vectDx,
         Vector<ProblemDomain>& vectDomain,
         const Real refine_thresh, 
         const int tags_grow,
         const int baseLevel, 
         int numLevels) 
{
  for (int lev=baseLevel; lev!= numLevels; lev++) 
    {
      IntVectSet local_tags;
      LevelData<FArrayBox> & levelRhs = *vectRHS[lev];
      DisjointBoxLayout level_domain = levelRhs.getBoxes();
      DataIterator dit = levelRhs.dataIterator();

      Real maxRHS = 0;
      Real thisMax = 0;
      for (dit.reset(); dit.ok(); ++dit) 
        {
          thisMax = levelRhs[dit()].norm(0,0,1);
          if (thisMax > maxRHS) maxRHS = thisMax;
        }

      Real tagVal = maxRHS * refine_thresh;
      
      // now loop through grids and tag cells where RHS > tagVal
      for (dit.reset(); dit.ok(); ++dit) 
        {
          const Box thisBox = level_domain.get(dit());
          const FArrayBox& thisRhs = levelRhs[dit()];
          BoxIterator bit(thisBox);
          for (bit.begin(); bit.ok(); ++bit) 
            {
              const IntVect& iv = bit();
              if (abs(thisRhs(iv)) >= tagVal)
                local_tags |= iv;
            }
        } // end loop over grids on this level
      
      local_tags.grow(tags_grow);
      const Box& domainBox = vectDomain[lev].domainBox();
      local_tags &= domainBox;
      int numTags = local_tags.numPts();
      Vector<IntVectSet> all_tags;
      const int dest_proc = uniqueProc(SerialTask::compute);
      gather (all_tags, local_tags, dest_proc);
      if (procID() == uniqueProc (SerialTask::compute)) {
        for (int i=0; i< all_tags.size(); ++i) 
          {
            tagVect[lev] |= all_tags[i];
          }
      }
      numTags = tagVect[lev].numPts();
    } // end loop over levels
}  


/*
  set domain boundary conditions from input file
  note that periodic BC's will override whichever physical
  boundary conditions are set here.
 */
int
setDomainBC(DomainGhostBC& domghostbc, bool verbose)
{
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  std::vector<int> ibclo;
  std::vector<int> ibchi;
  ParmParse pp("main");
  pp.getarr("bc_lo", ibclo, 0, SpaceDim);
  pp.getarr("bc_hi", ibchi, 0, SpaceDim);
  int ncomp;
  pp.get("ncomp", ncomp);
  Interval comps(0,ncomp-1);

  for(int idir = 0; idir < SpaceDim; idir++)
    {
      //lo bc
      {
        NeumannBC neumbc(idir, Side::Lo,comps);
        DirichletBC dircbc(idir, Side::Lo,comps);
        if(ibclo[idir] == 0) 
          {
            if(verbose&& procID()==0)
              pout() << "Dirichlet bcs in low direction for side " << idir << endl;
            domghostbc.setBoxGhostBC(dircbc);
          }
        else if(ibclo[idir] == 1)
          {
            if(verbose&& procID()==0)
              pout() << "neumann bcs in low direction for side " << idir << endl;
            domghostbc.setBoxGhostBC(neumbc);
          }
        else
          {
            if(verbose)
              cerr << "setDomainBC:: bogus input bc_lo flag" << endl;
            return(1);
          }
      }
      //hi bc
      {
        NeumannBC neumbc(idir, Side::Hi,comps);
        DirichletBC dircbc(idir, Side::Hi,comps);
        if(ibchi[idir] == 0) 
          {
            if(verbose&& procID()==0)
              pout() << "Dirichlet bcs in high direction for side " << idir << endl;
            domghostbc.setBoxGhostBC(dircbc);
          }
        else if(ibchi[idir] == 1)
          {
            if(verbose&& procID()==0)
              pout() << "neumann bcs in high direction for side " << idir << endl;
            domghostbc.setBoxGhostBC(neumbc);
          }
        else
          {
            if(verbose)
              cerr << "setDomainBC:: bogus input bc_hi flag" << endl;
            return(2);
          }
      }
    }
  return(0);
}


int
setTanGradBC(DomainGhostBC& a_domghostBC, bool verbose)
{

  Interval BCcomps(0,(SpaceDim*SpaceDim-1));
  // for now, assume Dirichlet BC's
  for (int idir=0; idir<SpaceDim; idir++) 
    {
      // lo bc
      {
        DirichletBC thisbc(idir, Side::Lo, BCcomps);
        //NoOpBC thisbc(idir, Side::Lo, BCcomps);
        a_domghostBC.setBoxGhostBC(thisbc);
      }
      
      // hi bc
      {
        DirichletBC thisbc(idir, Side::Hi, BCcomps);
        //scalarDirichletBC thisbc(idir, Side::Hi, BCcomps);
        //thisbc.setVal(1.0);
        //NoOpBC thisbc(idir, Side::Hi, BCcomps);
        a_domghostBC.setBoxGhostBC(thisbc);
      }
    }
  return(0);
}


void 
computeError(LevelData<FArrayBox>& a_error, 
             const LevelData<FArrayBox>& a_phi,
             const DisjointBoxLayout* a_fineGridsPtr,
             const ProblemDomain& a_domain,
             const Real a_dx,
             const int a_nRefFine)
{


  setExact(a_error, a_domain, a_dx);
  
  
  DataIterator dit = a_error.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisErr = a_error[dit()];

      // now subtract computed solution
      thisErr -= a_phi[dit()];
      
    }
  
  // now white out finer grids
  if (a_fineGridsPtr != NULL)
    {
      const DisjointBoxLayout& crseGrids = a_error.getBoxes();
      LayoutIterator fineLit = a_fineGridsPtr->layoutIterator();
      for (fineLit.begin(); fineLit.ok(); ++fineLit)
        {
          Box coarsenedFineBox(a_fineGridsPtr->get(fineLit()));
          coarsenedFineBox.coarsen(a_nRefFine);
          for (dit.begin(); dit.ok(); ++dit) 
            {
              Box intersectBox(crseGrids.get(dit()));
              intersectBox &= coarsenedFineBox;
              if (!intersectBox.isEmpty())
                {
                  a_error[dit()].setVal(0.0, intersectBox, 0, 
                                        a_error.nComp());
                } // end if intersection with fine grid
            } // end loop over this level's grids
        } // end loop over fine grid
    } // end if there are fine grids

}

void
setExact(LevelData<FArrayBox>& a_phi, 
         const ProblemDomain& a_domain,
         const Real a_dx)
{
  int iprob;
  // problem 1 is better for dirichlet problem, while 
  // problem 2 is for the periodic case
  if (a_domain.isPeriodic())
    {
      iprob = 2;
    } else {
      iprob = 1;
    }
  
  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisPhi = a_phi[dit()];
      const Box& thisBox = thisPhi.box();
      FORT_EXACT(CHF_FRA(thisPhi),
                 CHF_BOX(thisBox),
                 CHF_REAL(a_dx),
                 CHF_INT(iprob));
      
    }
  
}
