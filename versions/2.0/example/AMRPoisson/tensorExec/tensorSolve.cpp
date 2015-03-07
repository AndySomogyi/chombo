#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include  <iostream>
#include "LevelOp.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "PoissonBC.H"
#include "ParmParse.H"
#include "AMRSolver.H"
#include "TensorOp.H"
#include "Vector.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "localFuncs.H"
#include "computeSum.H"
#include "TensorProbF_F.H"
#include "ProblemDomain.H"
#include "UsingNamespace.H"

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  //scoping trick
  {
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  if (argc < 2) 
    {
      cerr<< " usage " << argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
  char* in_file = argv[1];
  ParmParse pp(argc-2,argv+2,NULL,in_file);

  int lbase;

  int iverbose;
  ParmParse pp2("main");
  pp2.get("lbase", lbase);
  pp2.get("verbose", iverbose);
  bool verbose = (iverbose == 1);
  int ncomp = 2;
  pp2.query("ncomp", ncomp);

  int numPasses = 1;
  pp2.query("num_passes", numPasses);
  Real domainMult = 1.0;
      
  DomainGhostBC domghostbc;
  int eekflag = setDomainBC(domghostbc, verbose);
  if(eekflag != 0)
    {
      cerr << "error: setDomainBC returned error code  " 
           << eekflag << endl;
      exit(1);
    }
  
  DomainGhostBC tanGradBC;
  eekflag = setTanGradBC(tanGradBC, verbose);
  if(eekflag != 0)
    {
      cerr << "error: setTanGradBC returned error code  " 
           << eekflag << endl;
      exit(1);
    }  

  Vector<Real> L1Old(ncomp,0);
  Vector<Real> L2Old(ncomp,0);
  Vector<Real> MaxOld(ncomp,0);

  for (int passNo = 0; passNo<numPasses; passNo++) 
    {
      
      TensorOp levelop;
      levelop.setDomainGhostBC(domghostbc);
      levelop.setTanGradBC(tanGradBC);
      
      int numlevels;
      Vector<DisjointBoxLayout> vectGrids;
      Vector<ProblemDomain> vectDomain;
      Vector<int> vectRefRatio;
      Vector<Real> vectDx;
      eekflag = setGrids(vectGrids, vectDomain, vectDx, 
                         vectRefRatio, numlevels, verbose,
                         domainMult);
      if(eekflag != 0)
        {
          cerr << "error: setGrids returned error code " 
               << eekflag << endl;
          exit(2);
        }
      
      Vector<LevelData<FArrayBox>* > phi(numlevels, NULL);
      Vector<LevelData<FArrayBox>* > rhs(numlevels, NULL);
      for(int ilev = 0; ilev < numlevels; ilev++)
        {
          phi[ilev] = new LevelData<FArrayBox>(vectGrids[ilev], ncomp, IntVect::Unit);
          rhs[ilev] = new LevelData<FArrayBox>(vectGrids[ilev], ncomp, IntVect::Zero);
          LevelData<FArrayBox>& phifabs= *phi[ilev];
#if 0
          Real rval = Real(3*ilev);
          for(DataIterator dit = phifabs.dataIterator(); dit.ok(); ++dit)
            {
              rval += 1.0;
              phifabs[dit()].setVal(rval);
            }
#endif
          setExact(phifabs, vectDomain[ilev], vectDx[ilev]);
          
        }
      eekflag = setRHS(rhs, vectDx, vectDomain, numlevels, verbose);
      if(eekflag != 0)
        {
          cerr << "error: setRHS returned error code " 
               << eekflag << endl;
          exit(3);
        }
      
      // now compute sum(RHS) as a diagnostic
      // (should be 0 for all-Neumann or periodic BC's for solvability).
      Interval comps(0,ncomp-1);
      Real sumRHS = computeSum(rhs,vectRefRatio,
                               vectDx[0], comps, 0);
      
      if (verbose && procID()==0)
        pout() << "Sum(RHS) = " << sumRHS << endl;
     
#if 0 
      // probably not the best way to go about this
      Real sumTol = 1.0e-12;
      if (abs(sumRHS) > sumTol)
        {
          pout () << "Rescaling..." << endl;
          for (int lev=0; lev<rhs.size(); lev++)
            {
              Real scal = vectDx[lev]*vectDx[lev];
              if (SpaceDim == 3) scal *= vectDx[lev];
              scal *= sumRHS;
              pout() << "level " << lev << ", scal = " << scal << endl;

              LevelData<FArrayBox>& levelRHS = *rhs[lev];
              DataIterator rhsDit = levelRHS.dataIterator();
              for (rhsDit.begin(); rhsDit.ok(); ++rhsDit)
                {
                  levelRHS[rhsDit] -= scal;
                }
              
            }
          
          sumRHS = computeSum(rhs,vectRefRatio,
                              vectDx[0], comps, 0);
          
          if (verbose && procID()==0)
            pout() << "Sum(RHS) = " << sumRHS << endl;
          
        } // end if rescaling for solvability
#endif

      AMRSolver amrSolver;

      int numSmoothUp = 4;
      int numSmoothDown = 4;
      pp2.query("num_smooth_up", numSmoothUp);
      pp2.query("num_smooth_down", numSmoothDown);
      amrSolver.setNumSmoothUp(numSmoothUp);
      amrSolver.setNumSmoothDown(numSmoothDown);

      amrSolver.define(vectGrids,  vectDomain,
                       vectDx,     vectRefRatio, 
                       numlevels, lbase, &levelop,
                       ncomp);


      amrSolver.setVerbose(verbose);
      int maxiter;
      pp2.get("max_iterations", maxiter);
      amrSolver.setMaxIter(maxiter);
      int numvbot;
      pp2.get("num_v_cycles_bottom", numvbot);
      amrSolver.setNumVCyclesBottom(numvbot);

      // Add the ability to loop over the solve (the same solve) for benchmarking 
      //  purposes only.  (ndk)  note that it will default to one time thru...
      int solve_iterations = 1;
      pp2.get("iterations", solve_iterations);
      for(int i=0; i<solve_iterations; i++) {
        amrSolver.solveAMR(phi, rhs);
      }
      eekflag = outputData(phi, rhs, vectGrids, vectDomain,
                           vectRefRatio, vectDx, amrSolver, 
                           numlevels, lbase,
                           L1Old, L2Old, MaxOld, verbose);
      if(eekflag != 0)
        {
          cerr << "error: outputData returned error code " 
               << eekflag << endl;
          exit(4);
        }
      for(int ilev = 0; ilev < numlevels; ilev++)
        {
          delete phi[ilev];
          delete rhs[ilev];
        }
      domainMult *= 2.0;
      
    } // end loop over convergence-testing passes
  }//end scoping trick
#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
  
  return(0);
  }
