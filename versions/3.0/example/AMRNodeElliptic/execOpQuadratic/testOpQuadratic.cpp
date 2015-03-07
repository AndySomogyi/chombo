#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// testOpQuadratic.cpp
// copy from Chombo/example/AMRPoisson/exec/poissonSolve.cpp
// petermc, 1 Nov 2000
// redone by petermc, 19 Jun 2002

/**
   Computes
   L(phi)
   where phi is a quadratic function, x*x + y*y + z*z
   Writes out truncation error, which is L[exact](phi) - L[calculated](phi).
   Run with inputs1op4 in 2-D, input1op in 3-D.
   Output opOut2.hdf5 or opOut3.hdf5:
   "phi":  phi
   "lap":  L[calculated](phi)
   "err":  error in Laplacian, L[calculated](phi) - 2*SpaceDim
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
#include "SPMD.H"
#include "RealVect.H"
#include "Norms.H"
#include "generalFuncs.H"
#include "localFuncs.H"
#include "testProb_F.H"
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
        cerr << "error: setGrids returned error code "
             << eekflag << endl;
        exit(2);
      }

    Vector<LevelData<NodeFArrayBox>* > phi(numlevels, NULL);
    Vector<LevelData<NodeFArrayBox>* > lofPhi(numlevels, NULL);
    Vector<LevelData<NodeFArrayBox>* > rhsExact(numlevels, NULL);

    for (int ilev = 0; ilev < numlevels; ilev++)
      {
        // DisjointBoxLayout vectGrids[ilev], 1 component, 1 ghost cell layer
        phi[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids[ilev], 1, IntVect::Unit);
        // not sure about ghost cells on lofPhi.
        lofPhi[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids[ilev], 1, IntVect::Zero);
        rhsExact[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids[ilev], 1, IntVect::Zero);
      }

    // Set phi.

    eekflag = setPhiQuadratic(phi, vectDx, vectDomain, numlevels, verbose);
    if (verbose)
      cout << "Did setPhiQuadratic" << endl;
    if (eekflag != 0)
      {
        cerr << "error: setPhiQuadratic returned error code "
             << eekflag << endl;
        exit(4);
      }

#ifdef CH_USE_HDF5
    writeNodeLevel(phi[0]);
#endif // CH_USE_HDF5

    // Define amrNodeSolver.

    AMRNodeSolver amrNodeSolver;
    amrNodeSolver.setVerbose(verbose);
    amrNodeSolver.define(vectGrids,  vectDomain,
                         vectDx,     vectRefRatio,
                         numlevels, lbase, &levelop);

    // Apply operator.

    // applyAMROperator() overwrites phi[ilev] with phi[ilev+1].
    // If you start with the finest level, then you'll get some overwriting
    // at each level.

    for (int ilev = numlevels-1; ilev >= 0; ilev--)
      // for (int ilev = 0; ilev < numlevels; ilev++)
      {
        LevelData<NodeFArrayBox>* lofPhiLevel = lofPhi[ilev];
        // lofPhiLevel should be same shape as phi[ilev].
        amrNodeSolver.applyAMROperator(*lofPhiLevel, phi, ilev);

        LevelData<NodeFArrayBox>* rhsLevel = rhsExact[ilev];
        for (DataIterator dit(rhsLevel->dataIterator()); dit.ok(); ++dit)
          // laplacian of quadratic
          rhsLevel->operator[](dit()).getFab().setVal(Real(2 * SpaceDim));
      }

    if (verbose)
      pout() << "Getting truncation error" << endl;

    Vector<LevelData<NodeFArrayBox>* > err(numlevels, NULL);
    eekflag = getTruncError(err, lofPhi, rhsExact,
                            vectGrids, vectDomain, vectRefRatio, vectDx,
                            numlevels, verbose);

    Interval errInterval(0, 0);
    Real normError = maxnorm(err, vectDomain, vectRefRatio, errInterval, 0, true);
    cout << "Max error: " << normError << endl;

    // Write output.
#if (CH_SPACEDIM == 2)
    string filename("opOut2.hdf5");
#else
    string filename("opOut3.hdf5");
#endif
    Vector<Vector<LevelData<NodeFArrayBox>* > > outvars;
    Vector<string> outnames;
    outvars.push_back(phi); outnames.push_back(string("phi"));
    outvars.push_back(lofPhi); outnames.push_back(string("L(phi)"));
    outvars.push_back(rhsExact); outnames.push_back(string("exact_L(phi)"));
    outvars.push_back(err); outnames.push_back(string("error_in_L(phi)"));

#ifdef CH_USE_HDF5
    if (verbose)
      pout() << "Writing output" << endl;

    WriteAMRHierarchyHDF5(filename,
                          vectGrids, outvars, outnames,
                          vectDomain[0].domainBox(), vectDx[0],
                          0.0, 0.0, // dt and time, not used
                          vectRefRatio, numlevels);
#endif // CH_USE_HDF5

    for(int ilev = 0; ilev < numlevels; ilev++)
      {
        delete phi[ilev];
        delete lofPhi[ilev];
        delete err[ilev];
        delete rhsExact[ilev];
      }
  }//end scoping trick
#ifdef CH_MPI
  MPI_Finalize();
#endif

  return(0);
}
