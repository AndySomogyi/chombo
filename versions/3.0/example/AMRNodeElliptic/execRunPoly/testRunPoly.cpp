#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// testRunPoly.cpp
// petermc, 17 July 2002

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
    int deginterp = 2; // default quadratic coarse/fine interpolation
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

    // read parameters for right-hand side, and set common block.

    RealVect center = (probLo + probHi);
    center *= 0.5;

    int iprob; pp2.get("iprob", iprob); // 1
    Real rhono; pp2.get("rhono", rhono); // 0.75
    Real rno; pp2.query("rno", rno); // 0.5

    if (verbose)
      cout
        << " rhono  = " << rhono
        << " rno  = " << rno
        << " iprob  = " << iprob
        << endl;

    FORT_INITPOLY(CHF_CONST_REAL(rhono),
                  CHF_CONST_REAL(rno),
                  CHF_CONST_INT(iprob),
                  CHF_REALVECT(center));

    Vector<LevelData<NodeFArrayBox>* > phi(numlevels, NULL);
    Vector<LevelData<NodeFArrayBox>* > phiExact(numlevels, NULL);
    Vector<LevelData<NodeFArrayBox>* > rhs(numlevels, NULL);
    Vector<LevelData<NodeFArrayBox>* > gradExact(numlevels, NULL);
    for (int ilev = 0; ilev < numlevels; ilev++)
      {
        // DisjointBoxLayout vectGrids[ilev], 1 component, 1 ghost cell layer
        // Need a ghost layer so that you can evaluate the operator
        // on the interior boundary cells.
        phi[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids[ilev], 1, IntVect::Unit);

        phiExact[ilev] = new
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

        gradExact[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids[ilev], SpaceDim, IntVect::Zero);
      }
    eekflag = setRHS(rhs, vectDx, vectDomain, numlevels, verbose);
    if (eekflag != 0)
      {
        cerr << "error: setRHS returned error code "
             << eekflag << endl;
        exit(3);
      }
    // Set exact phi.
    eekflag = setPhiExact(phiExact, vectDx, vectDomain, numlevels, verbose);
    if (verbose)
      cout << "Did setPhiExact on phiExact" << endl;
    if (eekflag != 0)
      {
        cerr << "error: setPhiExact returned error code "
             << eekflag << endl;
        exit(4);
      }
    // Set exact gradient.
    eekflag = setGradExact(gradExact, vectDx, vectDomain, numlevels, verbose);
    if (verbose)
      cout << "Did setGradExact on gradExact" << endl;
    if (eekflag != 0)
      {
        cerr << "error: setGradExact returned error code "
             << eekflag << endl;
        exit(5);
      }

    // Define amrNodeSolver.

    if (verbose)
      cout << "*** setting up first AMRNodeSolver *** " << endl;
    AMRNodeSolver amrNodeSolver;
    amrNodeSolver.setVerbose(verbose);
    amrNodeSolver.define(vectGrids,  vectDomain,
                         vectDx,     vectRefRatio,
                         numlevels, lbase, &levelop);
    amrNodeSolver.setMaxIter(20);
    amrNodeSolver.setBottomMaxIter(0);
    if (verbose)
      cout << "*** done with first AMRNodeSolver setup *** " << endl;
    amrNodeSolver.solveAMR(phi, rhs);
    if (verbose)
      cout << "*** done with first solve *** " << endl;

    // Get difference between computed and exact solutions.
    Vector<LevelData<NodeFArrayBox>* > err(numlevels, NULL);
    eekflag = getDiff(err, phi, phiExact);
    if (verbose)
      cout << "Did getDiff" << endl;
    if (eekflag != 0)
      {
        cerr << "error: getDiff returned error code "
             << eekflag << endl;
        exit(6);
      }

    Real dx0 = vectDx[0];
    Interval myInterval = phi[0]->interval();
    Real norm0err = norm(err, vectDomain, vectRefRatio,
                         dx0, myInterval, 0, 0, verbose);
    Real norm1err = norm(err, vectDomain, vectRefRatio,
                         dx0, myInterval, 1, 0, verbose);
    Real norm2err = norm(err, vectDomain, vectRefRatio,
                         dx0, myInterval, 2, 0, verbose);
    cout << "NORM OF ERROR IN SOLUTION (base h=" << dx0
         << "):  ";
    printf("Linf=%.5E, L1=%.5E, L2=%.5E\n", norm0err, norm1err, norm2err);

    /*
      Calculate error in operator.
      This part below added by petermc, 22 July 2002
    */

    if (verbose)
      cout << "*** computing operator for first AMRNodeSolver ***" << endl;
    // Apply operator.
    // applyAMROperator() overwrites phiExact[ilev] with phiExact[ilev+1].
    // If you start with the finest level, then you'll get some overwriting
    // at each level.
    Vector<LevelData<NodeFArrayBox>* > lofPhi(numlevels, NULL);

    for (int ilev = numlevels-1; ilev >= 0; ilev--)
      // for (int ilev = 0; ilev < numlevels; ilev++)
      {
        lofPhi[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids[ilev], 1, IntVect::Zero);
        LevelData<NodeFArrayBox>* lofPhiLevel = lofPhi[ilev];
        // lofPhiLevel should be same shape as phi[ilev].
        amrNodeSolver.applyAMROperator(*lofPhiLevel, phiExact, ilev);
      }

    // Get difference between computed and exact answers.
    Vector<LevelData<NodeFArrayBox>* > errOp(numlevels, NULL);
    getTruncError(errOp, lofPhi, rhs,
                  vectGrids, vectDomain, vectRefRatio, vectDx,
                  numlevels, verbose);
    // getTruncError does NOT zero out the error on covered nodes.
    // Rather, norm() ignores the covered nodes.
    Real norm0errOp = norm(errOp, vectDomain, vectRefRatio,
                           dx0, myInterval, 0, 0, verbose);
    Real norm1errOp = norm(errOp, vectDomain, vectRefRatio,
                           dx0, myInterval, 1, 0, verbose);
    Real norm2errOp = norm(errOp, vectDomain, vectRefRatio,
                           dx0, myInterval, 2, 0, verbose);
    cout << "NORM OF ERROR IN OPERATOR (base h=" << dx0
         << "):  ";
    printf("Linf=%.5E, L1=%.5E, L2=%.5E\n", norm0errOp, norm1errOp, norm2errOp);

    /*
      Calculate gradient.
      This part below added by petermc, 28 Apr 2003
    */

    if (verbose)
      cout << "*** computing gradient for first AMRNodeSolver ***" << endl;
    // Apply operator.
    // applyAMROperator() overwrites phiExact[ilev] with phiExact[ilev+1].
    // If you start with the finest level, then you'll get some overwriting
    // at each level.
    Vector<LevelData<NodeFArrayBox>* > gradPhi(numlevels, NULL);

    for (int ilev = numlevels-1; ilev >= 0; ilev--)
      // for (int ilev = 0; ilev < numlevels; ilev++)
      {
        gradPhi[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids[ilev], SpaceDim, IntVect::Zero);
        LevelData<NodeFArrayBox>* gradPhiLevel = gradPhi[ilev];
        // lofPhiLevel should be same shape as phi[ilev].
        amrNodeSolver.applyAMRGradient(*gradPhiLevel, phi, ilev);
      }

    // Get difference between computed and exact answers.
    Vector<LevelData<NodeFArrayBox>* > errGrad(numlevels, NULL);
    getTruncError(errGrad, gradPhi, gradExact,
                  vectGrids, vectDomain, vectRefRatio, vectDx,
                  numlevels, verbose);
    // getTruncError does NOT zero out the error on covered nodes.
    // Rather, norm() ignores the covered nodes.

    Vector<LevelData<NodeFArrayBox>* > errGradMag(numlevels, NULL);
    getMagnitude(errGradMag, errGrad);

    Real norm0errGrad = norm(errGradMag, vectDomain, vectRefRatio,
                             dx0, myInterval, 0, 0, verbose);
    Real norm1errGrad = norm(errGradMag, vectDomain, vectRefRatio,
                             dx0, myInterval, 1, 0, verbose);
    Real norm2errGrad = norm(errGradMag, vectDomain, vectRefRatio,
                             dx0, myInterval, 2, 0, verbose);

    cout << "NORM OF ERROR IN GRADIENT (base h=" << dx0
         << "):  ";
    printf("Linf=%.5E, L1=%.5E, L2=%.5E\n", norm0errGrad, norm1errGrad, norm2errGrad);

    /*
      Now grids that are factor of 2 finer.
    */

    NodePoissonOp levelop2;
    levelop2.setDomainNodeBC(dombc);
    levelop2.setInterpolationDegree(deginterp);
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

    Vector<LevelData<NodeFArrayBox>* > phi2(numlevels, NULL);
    Vector<LevelData<NodeFArrayBox>* > phiExact2(numlevels, NULL);
    Vector<LevelData<NodeFArrayBox>* > rhs2(numlevels, NULL);
    Vector<LevelData<NodeFArrayBox>* > gradExact2(numlevels, NULL);
    for (int ilev = 0; ilev < numlevels; ilev++)
      {
        // DisjointBoxLayout vectGrids[ilev], 1 component, 1 ghost cell layer
        // Need a ghost layer so that you can evaluate the operator
        // on the interior boundary cells.
        phi2[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids2[ilev], 1, IntVect::Unit);

        phiExact2[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids2[ilev], 1, IntVect::Unit);

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

        rhs2[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids2[ilev], 1, ghostExtent);

        gradExact2[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids2[ilev], SpaceDim, IntVect::Zero);
      }
    eekflag = setRHS(rhs2, vectDx2, vectDomain2, numlevels, verbose);
    if (eekflag != 0)
      {
        cerr << "error: setRHS returned error code "
             << eekflag << endl;
        exit(3);
      }
    // Set exact phi.
    eekflag = setPhiExact(phiExact2, vectDx2, vectDomain2, numlevels, verbose);
    if (verbose)
      cout << "Did setPhiExact on phiExact2" << endl;
    if (eekflag != 0)
      {
        cerr << "error: setPhiExact returned error code "
             << eekflag << endl;
        exit(4);
      }
    // Set exact gradient.
    eekflag = setGradExact(gradExact2, vectDx2, vectDomain2, numlevels, verbose);
    if (verbose)
      cout << "Did setGradExact on gradExact2" << endl;
    if (eekflag != 0)
      {
        cerr << "error: setGradExact returned error code "
             << eekflag << endl;
        exit(5);
      }

    // Define amrNodeSolver.

    if (verbose)
      cout << "*** setting up second AMRNodeSolver *** " << endl;
    AMRNodeSolver amrNodeSolver2;
    amrNodeSolver2.setVerbose(verbose);
    amrNodeSolver2.define(vectGrids2,  vectDomain2,
                          vectDx2,     vectRefRatio,
                          numlevels, lbase, &levelop2);
    amrNodeSolver2.setMaxIter(20);
    amrNodeSolver2.setBottomMaxIter(0);
    if (verbose)
      cout << "*** done with second AMRNodeSolver setup *** " << endl;
    amrNodeSolver2.solveAMR(phi2, rhs2);
    if (verbose)
      cout << "*** done with second solve *** " << endl;

    // Get difference between computed and exact solutions.
    Vector<LevelData<NodeFArrayBox>* > err2(numlevels, NULL);
    eekflag = getDiff(err2, phi2, phiExact2);
    if (verbose)
      cout << "Did getDiff" << endl;
    if (eekflag != 0)
      {
        cerr << "error: getDiff returned error code "
             << eekflag << endl;
        exit(6);
      }

    Real dx20 = vectDx2[0];
    Real norm0err2 = norm(err2, vectDomain2, vectRefRatio,
                          dx20, myInterval, 0, 0, verbose);
    Real norm1err2 = norm(err2, vectDomain2, vectRefRatio,
                          dx20, myInterval, 1, 0, verbose);
    Real norm2err2 = norm(err2, vectDomain2, vectRefRatio,
                          dx20, myInterval, 2, 0, verbose);
    cout << "NORM OF ERROR IN SOLUTION (base h=" << dx20
         << "):  ";
    printf("Linf=%.5E, L1=%.5E, L2=%.5E\n", norm0err2, norm1err2, norm2err2);

    /*
      Get difference between computed answers.
    */

    Vector<LevelData<NodeFArrayBox>* > projPhi2(numlevels, NULL);
    eekflag = project2(projPhi2, phi2, vectGrids);
    if (verbose)
      cout << "Did project2 on phi2" << endl;
    if (eekflag != 0)
      {
        cerr << "error: project2 returned error code "
             << eekflag << endl;
        exit(5);
      }

    Vector<LevelData<NodeFArrayBox>* > diff(numlevels, NULL);
    eekflag = getDiff(diff, phi, projPhi2);
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

    cout << "NORM OF DIFF (base " << dx20 << " to " << dx0
         << "):  ";
    printf("Linf=%.5E, L1=%.5E, L2=%.5E\n", norm0diff, norm1diff, norm2diff);

    /*
      Calculate error in operator.
      This part below added by petermc, 22 July 2002
    */

    if (verbose)
      cout << "*** computing operator for second AMRNodeSolver ***" << endl;
    // Apply operator.
    // applyAMROperator() overwrites phiExact2[ilev] with phiExact2[ilev+1].
    // If you start with the finest level, then you'll get some overwriting
    // at each level.
    Vector<LevelData<NodeFArrayBox>* > lofPhi2(numlevels, NULL);

    for (int ilev = numlevels-1; ilev >= 0; ilev--)
      // for (int ilev = 0; ilev < numlevels; ilev++)
      {
        lofPhi2[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids2[ilev], 1, IntVect::Zero);
        LevelData<NodeFArrayBox>* lofPhiLevel = lofPhi2[ilev];
        // lofPhiLevel should be same shape as phi[ilev].
        amrNodeSolver2.applyAMROperator(*lofPhiLevel, phiExact2, ilev);
      }

    // Get difference between computed and exact answers.
    Vector<LevelData<NodeFArrayBox>* > err2Op(numlevels, NULL);
    getTruncError(err2Op, lofPhi2, rhs2,
                  vectGrids2, vectDomain2, vectRefRatio, vectDx2,
                  numlevels, verbose);
    // getTruncError does NOT zero out the error on covered nodes.
    // Rather, norm() ignores the covered nodes.
    Real norm0err2Op = norm(err2Op, vectDomain2, vectRefRatio,
                            dx20, myInterval, 0, 0, verbose);
    Real norm1err2Op = norm(err2Op, vectDomain2, vectRefRatio,
                            dx20, myInterval, 1, 0, verbose);
    Real norm2err2Op = norm(err2Op, vectDomain2, vectRefRatio,
                            dx20, myInterval, 2, 0, verbose);
    cout << "NORM OF ERROR IN OPERATOR (base h=" << dx20
         << "):  ";
    printf("Linf=%.5E, L1=%.5E, L2=%.5E\n", norm0err2Op, norm1err2Op, norm2err2Op);

    /*
      Calculate gradient.
      This part below added by petermc, 28 Apr 2003
    */

    if (verbose)
      cout << "*** computing gradient for second AMRNodeSolver ***" << endl;
    // Apply operator.
    // applyAMROperator() overwrites phiExact2[ilev] with phiExact2[ilev+1].
    // If you start with the finest level, then you'll get some overwriting
    // at each level.
    Vector<LevelData<NodeFArrayBox>* > gradPhi2(numlevels, NULL);

    for (int ilev = numlevels-1; ilev >= 0; ilev--)
      // for (int ilev = 0; ilev < numlevels; ilev++)
      {
        gradPhi2[ilev] = new
          LevelData<NodeFArrayBox>(vectGrids2[ilev], SpaceDim, IntVect::Zero);
        LevelData<NodeFArrayBox>* gradPhiLevel = gradPhi2[ilev];
        // lofPhiLevel should be same shape as phi[ilev].
        amrNodeSolver2.applyAMRGradient(*gradPhiLevel, phi2, ilev);
      }

    // Get difference between computed and exact answers.
    Vector<LevelData<NodeFArrayBox>* > err2Grad(numlevels, NULL);
    getTruncError(err2Grad, gradPhi2, gradExact2,
                  vectGrids2, vectDomain2, vectRefRatio, vectDx2,
                  numlevels, verbose);
    // getTruncError does NOT zero out the error on covered nodes.
    // Rather, norm() ignores the covered nodes.

    Vector<LevelData<NodeFArrayBox>* > err2GradMag(numlevels, NULL);
    getMagnitude(err2GradMag, err2Grad);

    Real norm0err2Grad = norm(err2GradMag, vectDomain2, vectRefRatio,
                              dx20, myInterval, 0, 0, verbose);
    Real norm1err2Grad = norm(err2GradMag, vectDomain2, vectRefRatio,
                              dx20, myInterval, 1, 0, verbose);
    Real norm2err2Grad = norm(err2GradMag, vectDomain2, vectRefRatio,
                              dx20, myInterval, 2, 0, verbose);
    cout << "NORM OF ERROR IN GRADIENT (base h=" << dx20
         << "):  ";
    printf("Linf=%.5E, L1=%.5E, L2=%.5E\n", norm0err2Grad, norm1err2Grad, norm2err2Grad);

    // print norms of error

    cout << "NORM OF ERROR IN SOLUTION (base h=" << dx0
         << "):  ";
    printf("Linf=%.5E, L1=%.5E, L2=%.5E\n", norm0err, norm1err, norm2err);
    cout << "NORM OF ERROR IN OPERATOR (base h=" << dx0
         << "):  ";
    printf("Linf=%.5E, L1=%.5E, L2=%.5E\n", norm0errOp, norm1errOp, norm2errOp);
    cout << "NORM OF ERROR IN SOLUTION (base h=" << dx20
         << "):  ";
    printf("Linf=%.5E, L1=%.5E, L2=%.5E\n", norm0err2, norm1err2, norm2err2);
    cout << "NORM OF DIFF (base " << dx20 << " to " << dx0
         << "):  ";
    printf("Linf=%.5E, L1=%.5E, L2=%.5E\n", norm0diff, norm1diff, norm2diff);
    cout << "NORM OF ERROR IN OPERATOR (base h=" << dx20
         << "):  ";
    printf("Linf=%.5E, L1=%.5E, L2=%.5E\n", norm0err2Op, norm1err2Op, norm2err2Op);

    // norm of solution error

    printErrorNorms("NORMS", dx20, deginterp,
                    norm0err2, norm1err2, norm2err2,
                    norm0err, norm1err, norm2err);

    // norm of difference between phi and refined phi

    printDiffNorms("NORMD", dx20, deginterp,
                   norm0diff, norm1diff, norm2diff);

    // norm of truncation error

    printErrorNorms("NORMO", dx20, deginterp,
                    norm0err2Op, norm1err2Op, norm2err2Op,
                    norm0errOp, norm1errOp, norm2errOp);

    // norm of error in gradient
    printErrorNorms("NORMG", dx20, deginterp,
                    norm0err2Grad, norm1err2Grad, norm2err2Grad,
                    norm0errGrad, norm1errGrad, norm2errGrad);

    /*
      Write output.
    */

#if (CH_SPACEDIM == 2)
    string filename("out2.hdf5");
    string filename2("out2fine.hdf5");
#else
    string filename("out3.hdf5");
    string filename2("out3fine.hdf5");
#endif
    Vector<Vector<LevelData<NodeFArrayBox>* > > outvars;
    Vector<string> outnames;
    outvars.push_back(phi); outnames.push_back(string("phi"));
    outvars.push_back(rhs); outnames.push_back(string("rhs"));
    outvars.push_back(phiExact); outnames.push_back(string("exact_phi"));
    outvars.push_back(gradPhi); outnames.push_back(string("E_x")); outnames.push_back(string("E_y"));
#if (CH_SPACEDIM == 3)
    outnames.push_back(string("E_z"));
#endif
    outvars.push_back(gradExact); outnames.push_back(string("exact_E_x")); outnames.push_back(string("exact_E_y"));
#if (CH_SPACEDIM == 3)
    outnames.push_back(string("exact_E_z"));
#endif
    outvars.push_back(err); outnames.push_back(string("error_in_phi"));
    outvars.push_back(diff); outnames.push_back(string("phi_-_refined_phi"));
    outvars.push_back(errOp); outnames.push_back(string("error_in_L(phi)"));
    outvars.push_back(errGradMag); outnames.push_back(string("|error_in_grad(phi)|"));

#ifdef CH_USE_HDF5
    if (verbose)
      cout << " writing amr hierarchy" << endl;
    WriteAMRHierarchyHDF5(filename,
                          vectGrids, outvars, outnames,
                          vectDomain[0].domainBox(), vectDx[0],
                          0.0, 0.0, // dt and time, not used
                          vectRefRatio, numlevels);
    if (verbose)
      cout << " wrote amr hierarchy" << endl;
#endif // CH_USE_HDF5

    Vector<Vector<LevelData<NodeFArrayBox>* > > outvars2;
    Vector<string> outnames2;
    outvars2.push_back(phi2); outnames2.push_back(string("phi"));
    outvars2.push_back(rhs2); outnames2.push_back(string("rhs"));
    outvars2.push_back(phiExact2); outnames2.push_back(string("exact_phi"));
    outvars2.push_back(gradPhi2); outnames2.push_back(string("E_x")); outnames2.push_back(string("E_y"));
#if (CH_SPACEDIM == 3)
    outnames2.push_back(string("E_z"));
#endif
    outvars2.push_back(gradExact2); outnames2.push_back(string("exact_E_x")); outnames2.push_back(string("exact_E_y"));
#if (CH_SPACEDIM == 3)
    outnames2.push_back(string("exact_E_z"));
#endif
    outvars2.push_back(err2); outnames2.push_back(string("error_in_phi"));
    outvars2.push_back(err2Op); outnames2.push_back(string("error_in_L(phi)"));
    outvars2.push_back(err2GradMag); outnames2.push_back(string("|error_in_grad(phi)|"));

#ifdef CH_USE_HDF5
    if (verbose)
      cout << " writing second amr hierarchy" << endl;
    WriteAMRHierarchyHDF5(filename2,
                          vectGrids2, outvars2, outnames2,
                          vectDomain2[0].domainBox(), vectDx2[0],
                          0.0, 0.0, // dt and time, not used
                          vectRefRatio, numlevels);
#endif // CH_USE_HDF5

    if (verbose)
      cout << " about to delete stuff" << endl;

    for(int ilev = 0; ilev < numlevels; ilev++)
      {
        delete phi[ilev];
        delete phiExact[ilev];
        delete gradPhi[ilev];
        delete gradExact[ilev];
        delete rhs[ilev];
        delete lofPhi[ilev];
        delete err[ilev];
        delete errOp[ilev];
        delete phi2[ilev];
        delete phiExact2[ilev];
        delete gradPhi2[ilev];
        delete gradExact2[ilev];
        delete rhs2[ilev];
        delete lofPhi2[ilev];
        delete err2[ilev];
        delete diff[ilev];
        delete err2Op[ilev];
        delete projPhi2[ilev];
      }
    if (verbose)
      cout << " just finished deleting stuff " << endl;
  }//end scoping trick
#ifdef CH_MPI
  MPI_Finalize();
#endif

  return(0);
}
