#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// testNodePoisson.cpp
// petermc, 6 Aug 2003

/**
   Solves
   L(phi) = rho
   where rho is set in setRHS.
 */

#include <string>
using std::string;

#include  <iostream>
#include "LevelData.H"
#include "NodeFArrayBox.H"
#include "ParmParse.H"
#include "AMRNodeSolver.H"
#include "NodePoissonOp.H"
#include "Vector.H"
#include "NodeAMRIO.H"
#include "FABView.H"
#include "SPMD.H"
#include "Norms.H"
#include "RealVect.H"
#include "generalFuncs.H"
#include "localFuncs.H"
#include "testProb_F.H"
// for parallel debugging
#include "CH_Attach.H"
using std::cerr;
using std::cout;
using std::endl;
#include "CH_Timer.H"
#include "memusage.H"

OldTimer TM_Everything    ("Everything", 0);
OldTimer TM_Solve         ("Solve",      TM_Everything);

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  //scoping trick
  {
    int rank, number_procs;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank=0;
    number_procs=1;
#endif

    print_memory_line("before anything");
    get_memory_usage_from_OS();

#ifndef CH_NTIMER
    OldTimer::TimerInit(rank);

    TM_Everything.start();
#endif

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#ifndef NDEBUG
    // for parallel debugging
    // registerDebugger();
#endif
#endif

    if(argc < 2)
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

    int maxIter = 20;
    pp2.query("max_iterations", maxIter);

    int deginterp = 2; // default quadratic coarse/fine interpolation
    pp2.query("deg_interp", deginterp);

    bool writeOutput = false; // default is NOT to write hdf5 output
    {
       string writeOutputString;
       int gotWriteOutput = pp2.query("write_output", writeOutputString);
       if (gotWriteOutput != 0)
         { // the write_output line was present
           if (writeOutputString == "none")
             writeOutput = false;
           else if (writeOutputString == "all")
             writeOutput = true;
           else
             {
               cerr << "main.write_output must be \"none\" or \"all\"" << endl;
               MayDay::Error("returning");
             }
         }
    }

    DomainNodeBC dombc;
    // setDomainBC is in localFuncs.  It reads in the input file
    // and calls setBoxNodeBC on dombc.
    int eekflag;
    eekflag = setDomainBC(dombc, verbose);
    if (eekflag != 0)
      {
        cerr << "error: setDomainBC returned error code  "
             << eekflag << endl;
        exit(1);
      }

    NodePoissonOp levelop;
    levelop.setDomainNodeBC(dombc);
    levelop.setInterpolationDegree(deginterp);
    int numlevels;
    Vector<DisjointBoxLayout> vectGrids;
    Vector<ProblemDomain> vectDomain;
    Vector<int> vectRefRatio;
    Vector<Real> vectDx;
    RealVect probLo, probHi;

    eekflag = readGrids(vectGrids, vectDomain, vectDx,
                        vectRefRatio, probLo, probHi,
                        numlevels, verbose);

    if (eekflag != 0)
      {
        cerr << "error: first readGrids returned error code "
             << eekflag << endl;
        exit(2);
      }

    // read parameters for right-hand side, and set common block.

    RealVect center = (probLo + probHi);
    center *= 0.5;

    Vector<LevelData<NodeFArrayBox>* > phi(numlevels, NULL);
    Vector<LevelData<NodeFArrayBox>* > rhs(numlevels, NULL);
    for (int ilev = 0; ilev < numlevels; ilev++)
      {
        // DisjointBoxLayout vectGrids[ilev], 1 component, 1 ghost cell layer
        // Need a ghost layer so that you can evaluate the operator
        // on the interior boundary cells.
        phi[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids[ilev], 1, IntVect::Unit);

        // If there is a coarser level, then the depth of the ghost cell layer
        // of rhs should be half the refinement ratio to the coarser level.
        // This is so that you can do averaging on the residual at the
        // interior boundary cells; in order to get the residual, you
        // need rhs.
        IntVect ghostExtent;
        if (ilev > 0)
          ghostExtent = (vectRefRatio[ilev-1] / 2) * IntVect::Unit;
        else
          ghostExtent = IntVect::Zero;

        ghostExtent = IntVect::Zero; // added 29 Aug 2002

        rhs[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids[ilev], 1, ghostExtent);
      }
    eekflag = setRHS(rhs, vectDx, vectDomain, numlevels, verbose);
    if (eekflag != 0)
      {
        cerr << "error: setRHS returned error code "
             << eekflag << endl;
        exit(3);
      }

    // Define amrNodeSolver.

    if (verbose)
      pout() << "*** setting up first AMRNodeSolver *** " << endl;
    AMRNodeSolver amrNodeSolver;
    amrNodeSolver.setVerbose(verbose);
    amrNodeSolver.define(vectGrids,  vectDomain,
                         vectDx,     vectRefRatio,
                         numlevels, lbase, &levelop);
    amrNodeSolver.setMaxIter(maxIter);
    amrNodeSolver.setBottomMaxIter(0);
    if (verbose)
      pout() << "*** done with first AMRNodeSolver setup *** " << endl;
#ifndef CH_NTIMER
    TM_Solve.start();
#endif
    amrNodeSolver.solveAMR(phi, rhs);
#ifndef CH_NTIMER
    TM_Solve.stop();
#endif
    if (verbose)
      pout() << "*** done with first solve *** " << endl;

    Real dx0 = vectDx[0];
    Interval myInterval = phi[0]->interval();

    /*
      Write output.
    */

    if (writeOutput)
      {
#ifdef CH_USE_HDF5
#if (CH_SPACEDIM == 2)
        string filename("out2.hdf5");
#else
        string filename("out3.hdf5");
#endif
        Vector<Vector<LevelData<NodeFArrayBox>* > > outvars;
        Vector<string> outnames;
        outvars.push_back(phi); outnames.push_back(string("phi"));
        outvars.push_back(rhs); outnames.push_back(string("rhs"));
        cout << " writing amr hierarchy on " << filename << endl;
        WriteAMRHierarchyHDF5(filename,
                              vectGrids, outvars, outnames,
                              vectDomain[0].domainBox(), vectDx[0],
                              0.0, 0.0, // dt and time, not used
                              vectRefRatio, numlevels);
        cout << " wrote amr hierarchy on " << filename << endl;
#endif // CH_USE_HDF5
      }

    if (verbose)
      pout() << " about to delete stuff" << endl;

    for(int ilev = 0; ilev < numlevels; ilev++)
      {
        delete phi[ilev];
        delete rhs[ilev];
      }
    if (verbose)
      pout() << " just finished deleting stuff " << endl;

#if !defined(CH_NTIMER) && defined(CH_MPI)
    // Gather peak memory from procs and write a single line to screen.  that's all.
    Real end_memory = get_memory_usage_from_OS();
    Real avg_memory, min_memory, max_memory;
    gather_memory_from_procs(end_memory, avg_memory, min_memory, max_memory);
#endif

#ifndef CH_NTIMER
    TM_Everything.stop();

    OldTimer::TimerSummary();
#endif

  }//end scoping trick
#ifdef CH_MPI
  MPI_Finalize();
#endif

  return(0);
}
