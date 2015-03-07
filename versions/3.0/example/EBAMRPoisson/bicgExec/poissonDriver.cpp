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

#include "EBAMRPoissonOp.H"
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

/********/
void solve(const PoissonParameters&  a_params)
{
  int nvar = 1;
  Vector<DisjointBoxLayout> grids;
  Vector<EBISLayout>        ebisl;
  getAllIrregRefinedLayouts(grids, ebisl, a_params);
  for(int ilev = 0; ilev < grids.size(); ilev++)
    {
      pout() << "grids[" << ilev << "]= " << grids[ilev] << endl;
    }
  //define  data
  Vector<LevelData<EBCellFAB>* > phi(a_params.numLevels);
  Vector<LevelData<EBCellFAB>* > rhs(a_params.numLevels);

  for(int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      EBCellFactory factory(ebisl[ilev]);
      phi[ilev] = new LevelData<EBCellFAB>(grids[ilev],nvar, a_params.ghostPhi, factory);
      rhs[ilev] = new LevelData<EBCellFAB>(grids[ilev],nvar, a_params.ghostRHS, factory);

      //for now just set phi to zero and rhs to -1.
      EBLevelDataOps::setVal(*phi[ilev], 0.0);
      EBLevelDataOps::setVal(*rhs[ilev], 1.0);
    }

  pout() << "defining  solver" << endl;

  BiCGStabSolver<Vector<LevelData<EBCellFAB>* > > solver;
  Real alpha, beta;
  ParmParse pp;
  pp.get("alpha", alpha);
  pp.get("beta",  beta);
  Real time = 0;
  MultilevelLinearOp<EBCellFAB> mlOp;  //if this goes out of scope before solve all hell breakth loos
  defineMGBCGSolver(solver, mlOp, grids, ebisl,  a_params, time, alpha, beta);
  pout() << "solving " << endl;
  //solve the equation
  solver.solve(phi, rhs);


#ifdef CH_USE_HDF5
  pout() << "outputting the answer to file" << endl;
  //output the answer
  char charstr[100];
  sprintf(charstr, "phi.%dd.hdf5", SpaceDim);
  string filename(charstr);
  Vector<string> names(1, string("phi"));
  bool replaceCovered;
  Vector<Real> coveredValues(1, -1.0);
  Real dumReal =  1.0;
  Vector<LevelData<EBCellFAB> *> phidumbp(phi.size(), NULL);
  for(int ilev = 0; ilev < phi.size(); ilev++)
    {
      phidumbp[ilev] = &(*phi[ilev]);
    }
  writeEBHDF5(filename, grids, phidumbp, names,
              a_params.coarsestDomain.domainBox(), dumReal, dumReal, dumReal,
              a_params.refRatio, a_params.numLevels,
              replaceCovered, coveredValues);
#endif
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

    //read params from file
    getPoissonParameters(params);

    //define geometry from given params
    definePoissonGeometry(params);

    //solve the stinking problem and output everything
    solve(params);

  }
  // End scoping trick

#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}
