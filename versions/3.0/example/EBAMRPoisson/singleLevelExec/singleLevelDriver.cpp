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
using std::cerr;

#include "ParmParse.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"
#include "PoissonUtilities.H"

#include "EBFABView.H"
#include "EBDebugDump.H"

#include "EBPoissonOp.H"
#include "EBPoissonOpFactory.H"
#include "EBLevelDataOps.H"
#include "BaseBCValue.H"
#include "BaseDomainBC.H"
#include "NeumannPoissonDomainBC.H"
#include "DirichletPoissonDomainBC.H"
#include "BaseEBBC.H"
#include "DirichletPoissonEBBC.H"
#include "NeumannPoissonEBBC.H"

#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBLevelDataOps.H"
#include "EBSimpleSolver.H"
#include "BiCGStabSolver.H"
#include "MultiGrid.H"
#include "BRMeshRefine.H"
#include "EBEllipticLoadBalance.H"
#include "EBLoadBalance.H"
#include "LoadBalance.H"
#include "EBLevelGrid.H"
#include "CH_Timer.H"

/********/
void solve(const PoissonParameters&  a_params)
{
  CH_TIMERS("uber_solve");
  CH_TIMER("grid_generation", t1);
  CH_TIMER("load_balance", t2);
  CH_TIMER("solver_define", t3);
  CH_TIMER("solver_solve", t4);
  CH_TIMER("file_output", t5);


  CH_START(t1);
  int nvar = 1;
  int nghost = 4;
  Vector<Box> boxes;
  Vector<int> procs;
  domainSplit(a_params.coarsestDomain.domainBox(), boxes, a_params.maxGridSize);
  CH_STOP(t1);

  CH_START(t2);
  ParmParse pp;
  int ibalance;
  pp.get("which_loadbalance", ibalance);
  if(ibalance == 0)
    {
      pout() << "using non-eb loadbalance" << endl;
      LoadBalance(procs, boxes);
    }
  else if(ibalance == 1)
    {
      pout() << "using EBLoadBalance" << endl;
      EBLoadBalance(procs, boxes, a_params.coarsestDomain, true);
    }
  else
    {
      pout() << "using EBEllipticLoadBalance" << endl;
      EBEllipticLoadBalance(procs, boxes, a_params.coarsestDomain, true);
    }
  DisjointBoxLayout  grids(boxes, procs);
  CH_STOP(t2);

  EBLevelGrid eblg(grids, a_params.coarsestDomain, nghost, Chombo_EBIS::instance());

  pout() << "grids= " << grids << endl;

  //define  data
  IntVect nghostPhi = 4*IntVect::Unit;
  IntVect nghostRHS = IntVect::Zero;
  EBCellFactory ebcf(eblg.getEBISL());
  LevelData<EBCellFAB> phi(grids, nvar, nghostPhi, ebcf);
  LevelData<EBCellFAB> rhs(grids, nvar, nghostRHS, ebcf);;

  //for now just set phi to zero and rhs to -1.
  EBLevelDataOps::setVal(phi, 0.0);
  EBLevelDataOps::setVal(rhs, 1.0);

  CH_START(t3);
  //create the solver
  RefCountedPtr<BaseDomainBCFactory> domainBCFactory;
  RefCountedPtr<BaseEBBCFactory>     ebbcFactory;

  getBCFactories(domainBCFactory, ebbcFactory, a_params);

  int orderEB =2;
  int relaxType = 1;//gsrb
  int numPreCond;
  Real alpha, beta;
  pp.get("alpha", alpha);
  pp.get("beta",  beta);
  pp.get("num_pre_cond",  numPreCond);

  EBPoissonOpFactory factory(eblg, a_params.coarsestDx, RealVect::Zero, orderEB, numPreCond, relaxType,
                             domainBCFactory, ebbcFactory, alpha, beta, nghostPhi, nghostRHS);

  BiCGStabSolver<LevelData<EBCellFAB> > bottomSolver;
  bottomSolver.m_verbosity = 0;
  MultiGrid<LevelData<EBCellFAB> > solver;
  pout() << "defining  solver" << endl;
  solver.define(factory, &bottomSolver, a_params.coarsestDomain);

  CH_STOP(t3);
  pout() << "solving " << endl;
  int verbosity = 4;
  int maxIter, numSmooth;
  Real tolerance;
  pp.get("mg_eps", tolerance);
  pp.get("mg_iter_max", maxIter);
  pp.get("mg_num_smooths", numSmooth);
  solver.m_pre  = numSmooth;
  solver.m_post = numSmooth;
  //solve the equation
  CH_START(t4);
  solver.solve(phi, rhs, tolerance, maxIter, verbosity);
  CH_STOP(t4);


  CH_START(t5);
#ifdef CH_USE_HDF5
  pout() << "outputting the answer to file" << endl;
  //output the answer
  char charstr[100];
  sprintf(charstr, "phi.%dd.hdf5", SpaceDim);
  writeEBLevelname(&phi, charstr);
#endif
  CH_STOP(t5);
}
/******/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // Scoping trick
  {
#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif

    if (argc < 2)
      {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }

    char* inFile = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,inFile);

    PoissonParameters params;

    //read params from file -- true means force single level
    getPoissonParameters(params, true);

    //define geometry from given params
    definePoissonGeometry(params);

    //solve the stinking problem and output everything
    solve(params);

  }
  // End scoping trick

#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}
