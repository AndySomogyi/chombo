#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// testOpExp.cpp
// copy from Chombo/example/AMRPoisson/exec/poissonSolve.cpp
// petermc, 1 Nov 2000
// redone by petermc, 19 Jun 2002

/**
   Computes
   L(phi)
   where phi is an exponential function, exp(-(x*x + y*y + z*z)/(4*s*s))
   with s = 0.125.
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
#include "FABView.H"
#include "SPMD.H"
#include "RealVect.H"
#include "generalFuncs.H"
#include "localFuncs.H"
#include "Norms.H"
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
    int deginterp = 2;
    ParmParse pp2("main");
    pp2.get("lbase", lbase);
    pp2.get("verbose", iverbose);
    pp2.query("deg_interp", deginterp);
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

    Vector<DisjointBoxLayout> vectGrids2;
    Vector<ProblemDomain> vectDomain2;
    Vector<Real> vectDx2;
    // vectRefRatio, numlevels, probLo, probHi are the same.

    eekflag = readGrids(vectGrids2, vectDomain2, vectDx2,
                        vectRefRatio, probLo, probHi,
                        numlevels, verbose, 2);

    if (eekflag != 0)
      {
        cerr << "error: second readGrids returned error code "
             << eekflag << endl;
        exit(2);
      }

    Vector<LevelData<NodeFArrayBox>* > phi(numlevels, NULL);
    Vector<LevelData<NodeFArrayBox>* > phi2(numlevels, NULL);
    Vector<LevelData<NodeFArrayBox>* > lofPhi(numlevels, NULL);
    Vector<LevelData<NodeFArrayBox>* > lofPhi2(numlevels, NULL);
    Vector<LevelData<NodeFArrayBox>* > rhs(numlevels, NULL);
    Vector<LevelData<NodeFArrayBox>* > rhs2(numlevels, NULL);

    for (int ilev = 0; ilev < numlevels; ilev++)
      {
        // DisjointBoxLayout vectGrids[ilev], 1 component, 1 ghost cell layer
        phi[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids[ilev], 1, IntVect::Unit);
        phi2[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids2[ilev], 1, IntVect::Unit);

        lofPhi[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids[ilev], 1, IntVect::Zero);
        lofPhi2[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids2[ilev], 1, IntVect::Zero);

        rhs[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids[ilev], 1, IntVect::Zero);
        rhs2[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids2[ilev], 1, IntVect::Zero);
      }

    // Set parameters in operator.

    RealVect center = (probLo + probHi);
    center *= 0.5;
    Real sigma = 0.125;
    FORT_INITEXP(CHF_CONST_REAL(sigma),
                 CHF_REALVECT(center));

    Real dx0 = vectDx[0];
    Interval myInterval = phi[0]->interval();

    // Define amrNodeSolver.

    if (verbose)
      cout << "*** setting up first AMRNodeSolver *** " << endl;
    AMRNodeSolver amrNodeSolver;
    amrNodeSolver.setVerbose(verbose);
    amrNodeSolver.define(vectGrids,  vectDomain,
                         vectDx,     vectRefRatio,
                         numlevels, lbase, &levelop);
    if (verbose)
      cout << "*** done with AMRNodeSolver setup *** " << endl;

    // Set phi.

    eekflag = setPhiExp(phi, vectDx, vectDomain, numlevels, verbose);
    if (verbose)
      cout << "Did setPhiExp on phi" << endl;
    if (eekflag != 0)
      {
        cerr << "error: setPhiExp returned error code "
             << eekflag << endl;
        exit(4);
      }

    eekflag = setLaplacianExact(rhs, vectDx, vectDomain, numlevels, verbose);
    if (verbose)
      cout << "Did setLaplacianExact on rhs" << endl;
    if (eekflag != 0)
      {
        cerr << "error: setLaplacianExact returned error code "
             << eekflag << endl;
        exit(4);
      }

    if (verbose)
      cout << "*** computing operator for first AMRNodeSolver ***" << endl;
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
      }

    // Get difference between computed and exact answers.
    Vector<LevelData<NodeFArrayBox>* > err(numlevels, NULL);
    getTruncError(err, lofPhi, rhs,
                  vectGrids, vectDomain, vectRefRatio, vectDx,
                  numlevels, verbose);
    // getTruncError(err, phi, lofPhi,
    // vectGrids, vectDomain, vectDx, vectRefRatio, numlevels, verbose);
    // getTruncError does NOT zero out the error on covered nodes.
    // Rather, norm() ignores the covered nodes.
    Real norm0err = norm(err, vectDomain, vectRefRatio,
                         dx0, myInterval, 0, 0, verbose);
    Real norm1err = norm(err, vectDomain, vectRefRatio,
                         dx0, myInterval, 1, 0, verbose);
    Real norm2err = norm(err, vectDomain, vectRefRatio,
                         dx0, myInterval, 2, 0, verbose);
    cout << "NORM OF ERROR IN OPERATOR (base h=" << dx0
         << "):  ";
    printf("Linf=%.5E, L1=%.5E, L2=%.5E\n", norm0err, norm1err, norm2err);

    if (verbose)
      cout << "*** setting up second AMRNodeSolver *** " << endl;

    NodePoissonOp levelop2;
    levelop2.setDomainNodeBC(dombc);
    levelop2.setInterpolationDegree(deginterp);
    AMRNodeSolver amrNodeSolver2;
    amrNodeSolver2.setVerbose(verbose);
    amrNodeSolver2.define(vectGrids2,  vectDomain2,
                          vectDx2,     vectRefRatio,
                          numlevels, lbase, &levelop2);
    if (verbose)
      cout << "*** done with AMRNodeSolver setup *** " << endl;

    eekflag = setPhiExp(phi2, vectDx2, vectDomain2, numlevels, verbose);
    if (verbose)
      cout << "Did setPhiExp on phi2" << endl;
    if (eekflag != 0)
      {
        cerr << "error: setPhiExp returned error code "
             << eekflag << endl;
        exit(4);
      }

    eekflag = setLaplacianExact(rhs2, vectDx2, vectDomain2, numlevels, verbose);
    if (verbose)
      cout << "Did setLaplacianExact on rhs2" << endl;
    if (eekflag != 0)
      {
        cerr << "error: setLaplacianExact returned error code "
             << eekflag << endl;
        exit(4);
      }

    if (verbose)
      cout << "*** computing operator for second AMRNodeSolver ***" << endl;
    for (int ilev = numlevels-1; ilev >= 0; ilev--)
      // for (int ilev = 0; ilev < numlevels; ilev++)
      {
        LevelData<NodeFArrayBox>* lofPhiLevel2 = lofPhi2[ilev];
        // lofPhiLevel2 should be same shape as phi2[ilev].
        amrNodeSolver2.applyAMROperator(*lofPhiLevel2, phi2, ilev);
      }

    // Get difference between computed and exact answers.
    Vector<LevelData<NodeFArrayBox>* > err2(numlevels, NULL);
    Real dx20 = vectDx2[0];
    getTruncError(err2, lofPhi2, rhs2,
                  vectGrids2, vectDomain2, vectRefRatio, vectDx2,
                  numlevels, verbose);
    // getError does NOT zero out the error on covered nodes.
    // Rather, norm() ignores the covered nodes.
    Real norm0err2 = norm(err2, vectDomain2, vectRefRatio,
                          dx20, myInterval, 0, 0, verbose);
    Real norm1err2 = norm(err2, vectDomain2, vectRefRatio,
                          dx20, myInterval, 1, 0, verbose);
    Real norm2err2 = norm(err2, vectDomain2, vectRefRatio,
                          dx20, myInterval, 2, 0, verbose);
    cout << "NORM OF ERROR IN OPERATOR (base h=" << dx20
         << "):  ";
    printf("Linf=%.5E, L1=%.5E, L2=%.5E\n", norm0err2, norm1err2, norm2err2);

    /*
      Get difference between computed answers.
    */

    Vector<LevelData<NodeFArrayBox>* > projLofPhi2(numlevels, NULL);
    eekflag = project2(projLofPhi2, lofPhi2, vectGrids);
    if (verbose)
      cout << "Did project2 on lofPhi2" << endl;
    if (eekflag != 0)
      {
        cerr << "error: project2 returned error code "
             << eekflag << endl;
        exit(5);
      }

    Vector<LevelData<NodeFArrayBox>* > diff(numlevels, NULL);
    eekflag = getDiff(diff, lofPhi, projLofPhi2);
    if (verbose)
      cout << "Did getDiff" << endl;
    if (eekflag != 0)
      {
        cerr << "error: getDiff returned error code "
             << eekflag << endl;
        exit(6);
      }

    Real norm0diff = norm(diff, vectDomain, vectRefRatio,
                          dx0, myInterval, 0, 0, verbose);
    Real norm1diff = norm(diff, vectDomain, vectRefRatio,
                          dx0, myInterval, 1, 0, verbose);
    Real norm2diff = norm(diff, vectDomain, vectRefRatio,
                          dx0, myInterval, 2, 0, verbose);

    cout << "NORM OF DIFF (base h " << dx20 << " to " << dx0
         << "):  ";
    printf("Linf=%.5E, L1=%.5E, L2=%.5E\n", norm0diff, norm1diff, norm2diff);

    printErrorNorms("NORME", dx20, deginterp,
                    norm0err2, norm1err2, norm2err2,
                    norm0err, norm1err, norm2err);

    printDiffNorms("NORMD", dx20, deginterp,
                   norm0diff, norm1diff, norm2diff);

    // Write output.
#if (CH_SPACEDIM == 2)
    string filename("opOut2.hdf5");
    string filename2("opOut2fine.hdf5");
#else
    string filename("opOut3.hdf5");
    string filename2("opOut3fine.hdf5");
#endif
    Vector<Vector<LevelData<NodeFArrayBox>* > > outvars;
    Vector<string> outnames;
    outvars.push_back(phi); outnames.push_back(string("phi"));
    outvars.push_back(lofPhi); outnames.push_back(string("L(phi)"));
    outvars.push_back(rhs); outnames.push_back(string("exact_L(phi)"));
    outvars.push_back(err); outnames.push_back(string("error_in_L(phi)"));
    outvars.push_back(diff); outnames.push_back(string("L(phi)_-_L(refined_phi)"));

#ifdef CH_USE_HDF5
    WriteAMRHierarchyHDF5(filename,
                          vectGrids, outvars, outnames,
                          vectDomain[0].domainBox(), vectDx[0],
                          0.0, 0.0, // dt and time, not used
                          vectRefRatio, numlevels);
#endif // CH_USE_HDF5
    Vector<Vector<LevelData<NodeFArrayBox>* > > outvars2;
    Vector<string> outnames2;
    outvars2.push_back(phi2); outnames2.push_back(string("phi"));
    outvars2.push_back(lofPhi2); outnames2.push_back(string("L(phi)"));
    outvars2.push_back(rhs2); outnames2.push_back(string("exact_L(phi)"));
    outvars2.push_back(err2); outnames2.push_back(string("error_in_L(phi)"));

#ifdef CH_USE_HDF5
    WriteAMRHierarchyHDF5(filename2,
                          vectGrids2, outvars2, outnames2,
                          vectDomain2[0].domainBox(), vectDx2[0],
                          0.0, 0.0, // dt and time, not used
                          vectRefRatio, numlevels);
#endif // CH_USE_HDF5

    for(int ilev = 0; ilev < numlevels; ilev++)
      {
        delete phi[ilev];
        delete phi2[ilev];
        delete lofPhi[ilev];
        delete lofPhi2[ilev];
        delete projLofPhi2[ilev];
        delete diff[ilev];
        delete err[ilev];
        delete err2[ilev];
        delete rhs[ilev];
        delete rhs2[ilev];
      }
  }//end scoping trick
#ifdef CH_MPI
  MPI_Finalize();
#endif

  return(0);
}
