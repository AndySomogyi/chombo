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
#include "PoissonOp.H"
#include "Vector.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "localFuncs.H"
#include "computeSum.H"
#include "PoissProbF_F.H"
#include "ProblemDomain.H"
#include "UsingNamespace.H"

#ifdef CH_MPI
void dumpmemoryatexit();
#endif

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  //scoping trick
  {
  barrier();
  if (argc < 2) 
    {
      cerr<< " usage " << argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
  char* in_file = argv[1];
  ParmParse  pp(argc-2,argv+2,NULL,in_file);

  int lbase;

  int iverbose;
  ParmParse pp2("main");
  pp2.get("lbase", lbase);
  pp2.get("verbose", iverbose);
  bool verbose = (iverbose <= 1);
  int ncomp = 2;
  pp2.query("ncomp", ncomp);


  DomainGhostBC domghostbc;
  int eekflag = setDomainBC(domghostbc, verbose,ncomp);
  if(eekflag != 0)
    {
      cerr << "error: setDomainBC returned error code  " 
           << eekflag << endl;
      exit(1);
    }
  PoissonOp levelop;
  levelop.setDomainGhostBC(domghostbc);


  int numlevels;
  Vector<DisjointBoxLayout> vectGrids;
  Vector<ProblemDomain> vectDomain;
  Vector<int> vectRefRatio;
  Vector<Real> vectDx;
  eekflag = setGrids(vectGrids, vectDomain, vectDx, 
                     vectRefRatio, numlevels, verbose);
  if (iverbose >= 2)
    {
      for (int ilev = 0; ilev < vectGrids.size(); ilev++)
        {
          pout() << "grid at level " << ilev << " = " << vectGrids[ilev] << endl;
        }
    }

  if(eekflag != 0)
    {
      cerr << "error: setGrids returned error code " 
           << eekflag << endl;
      exit(2);
    }

  Vector<LevelData<FArrayBox>* > phi(numlevels, NULL);
  Vector<LevelData<FArrayBox>* > rhs(numlevels, NULL);

  bool constPhi = false;
  Real initialPhi = 0.0;
  if (pp2.contains("initialGuess"))
    {
      constPhi = true;
      pp2.get("initialGuess", initialPhi);
    }

  for(int ilev = 0; ilev < numlevels; ilev++)
    {
      phi[ilev] = new LevelData<FArrayBox>(vectGrids[ilev], ncomp, IntVect::Unit);
      rhs[ilev] = new LevelData<FArrayBox>(vectGrids[ilev], ncomp, IntVect::Zero);
      LevelData<FArrayBox>& phifabs= *phi[ilev];
      Real rval = Real(3*ilev);
      for(DataIterator dit = phifabs.dataIterator(); dit.ok(); ++dit)
      {
        rval += 1.0;
        if (constPhi)
          {
            // option to set phi to constant-value initial guess 
            phifabs[dit()].setVal(initialPhi);
          }
        else
          {          
          }
      }
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
  for (int comp=0; comp<ncomp; comp++) 
    {
      Interval comps(comp,comp);
      Real sumRHS = computeSum(rhs,vectRefRatio,
                               vectDx[0], comps, 0);
      
      if (verbose && procID()==0)
        pout() << "Component " << comp << ",  Sum(RHS) = " << sumRHS << endl;
      
    } // end loop over components

      

  AMRSolver amrSolver(vectGrids,  vectDomain,
                      vectDx,     vectRefRatio, 
                      numlevels, lbase, &levelop,
                      ncomp);

  amrSolver.setNumSmoothUp(4);
  amrSolver.setNumSmoothDown(4);

  amrSolver.setVerbose(verbose);
  int maxiter;
  pp2.get("max_iterations", maxiter);
  amrSolver.setMaxIter(maxiter);
  int numvbot;
  pp2.get("num_v_cycles_bottom", numvbot);
  amrSolver.setNumVCyclesBottom(numvbot);
  if (pp2.contains("solverTol"))
  {
    Real solverTol;
    pp2.get("solverTol", solverTol);
    amrSolver.setTolerance(solverTol);
  }

  if (pp2.contains("numSmoothUp"))
    {
      int numSmoothUp;
      pp2.get("numSmoothUp", numSmoothUp);
      amrSolver.setNumSmoothUp(numSmoothUp);
      pout() << "numSmoothUp = " << numSmoothUp << endl;
    }


  if (pp2.contains("numSmoothDown"))
    {
      int numSmoothDown;
      pp2.get("numSmoothDown", numSmoothDown);
      amrSolver.setNumSmoothDown(numSmoothDown);
      pout() << "numSmoothDown = " << numSmoothDown << endl;
    }


  if (pp2.contains("normType"))
    {
      int normType = 0;
      pp2.query("normType", normType);
      amrSolver.setNormType(normType);
    }


  // Add the ability to loop over the solve (the same solve) for benchmarking 
  //  purposes only.  (ndk)  note that it will default to one time thru...
  int solve_iterations = 1;
  pp2.get("iterations", solve_iterations);
  for(int i=0; i<solve_iterations; i++) 
    {
      for(int ilev = 0; ilev < numlevels; ilev++)
        {
          for(DataIterator dit = phi[ilev]->dataIterator(); dit.ok(); ++dit)
            {
              (*phi[ilev])[dit()].setVal(0.0);
            }
        }
      amrSolver.solveAMR(phi, rhs, false);
    }

  eekflag = outputData(phi, rhs, vectGrids, vectDomain,
                       vectRefRatio, vectDx[0], numlevels, verbose);
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
  }//end scoping trick

#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif

  return(0);
}
