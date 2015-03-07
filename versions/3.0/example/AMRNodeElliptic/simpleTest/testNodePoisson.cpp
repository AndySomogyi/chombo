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
#include "DebugDump.H"
#include "testProb_F.H"
// for parallel debugging
#include "CH_Attach.H"
using std::cerr;
using std::cout;
using std::endl;

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

//      int iprob; pp2.get("iprob", iprob); // 1
//      Real rhono; pp2.get("rhono", rhono); // 0.75
//      Real rno; pp2.query("rno", rno); // 0.5

//      if (verbose)
//        cout
//          << " rhono  = " << rhono
//          << " rno  = " << rno
//          << " iprob  = " << iprob
//          << endl;

//      FORT_INITPOLY(CHF_CONST_REAL(rhono),
//                    CHF_CONST_REAL(rno),
//                    CHF_CONST_INT(iprob),
//                    CHF_REALVECT(center));

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
    amrNodeSolver.solveAMR(phi, rhs);
    if (verbose)
      pout() << "*** done with first solve *** " << endl;



#if (CH_SPACEDIM == 2)
    string filename("out.2d.hdf5");
#else
    string filename("out.3d.hdf5");
#endif
    Vector<Vector<LevelData<NodeFArrayBox>* > > outvars;
    Vector<string> outnames;
    outvars.push_back(phi); outnames.push_back(string("phi"));
    outvars.push_back(rhs); outnames.push_back(string("rhs"));

    if (verbose)
      pout() << " writing amr hierarchy" << endl;
#ifdef CH_USE_HDF5
    WriteAMRHierarchyHDF5(filename,
                          vectGrids, outvars, outnames,
                          vectDomain[0].domainBox(), vectDx[0],
                          0.0, 0.0, // dt and time, not used
                          vectRefRatio, numlevels);
    if (verbose)
      pout() << " wrote amr hierarchy" << endl;
#endif // CH_USE_HDF5
    if (verbose)
      pout() << " about to delete stuff" << endl;

    for(int ilev = 0; ilev < numlevels; ilev++)
      {
        delete phi[ilev];
        delete rhs[ilev];
      }
    if (verbose)
      pout() << " just finished deleting stuff " << endl;
  }//end scoping trick
#ifdef CH_MPI
  MPI_Finalize();
#endif

  return(0);
}
