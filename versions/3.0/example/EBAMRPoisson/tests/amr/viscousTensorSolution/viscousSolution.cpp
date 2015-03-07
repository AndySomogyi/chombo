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
#include "EBAMRPoissonOp.H"

#include "EBFABView.H"
#include "EBDebugDump.H"
#include "DebugDump.H"
#include "EBSimpleSolver.H"
#include "AMRMultiGrid.H"
#include "CH_Attach.H"
#include "FaceIterator.H"

#define CHECK_FLOATING_PT 0

#if CHECK_FLOATING_PT==1
#include <fenv.h>
#endif

/******/
void getError(Vector< LevelData<EBCellFAB>* >&     a_error,
              const Vector< DisjointBoxLayout >&   a_grids,
              const Vector< EBISLayout >&          a_ebisl,
              const PoissonParameters&             a_params)
{
  ParmParse pp;

  int nlevels = a_params.numLevels;
  a_error.resize(nlevels);
  Vector<LevelData<EBCellFAB>* > velExac(nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > klpExac(nlevels, NULL);

  RealVect dxLev = RealVect::Unit;
  dxLev *=a_params.coarsestDx;
  pout() << "coarsest dx = " << dxLev << endl;
  ProblemDomain domLev = a_params.coarsestDomain;
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      a_error[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,a_params.ghostPhi, ebcellfact);
      klpExac[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,a_params.ghostRHS, ebcellfact);
      velExac[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,a_params.ghostPhi, ebcellfact);
      EBLevelDataOps::setToZero(*(a_error[ilev]));
      EBLevelDataOps::setToZero(*(velExac[ilev]));
      EBLevelDataOps::setToZero(*(klpExac[ilev]));

      //set phi = velExact, rhs=lphiexact  This makes AMRResidual return lphiexact-Lphi
      setVelViscous(*velExac[ilev], dxLev[0], a_params);
      setKLVViscous(*klpExac[ilev], dxLev[0], a_params);

      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);
    }

  //create the solver
  EBSimpleSolver bottomSolver;
  int numSmooths;
  pp.get("num_smooth", numSmooths);
  bottomSolver.setNumSmooths(numSmooths);

  AMRMultiGrid<LevelData<EBCellFAB> > solver;
  pout() << "defining  AMRMultiGrid solver" << endl;
  defineViscousTensorSolver(solver, a_grids, a_ebisl, bottomSolver, a_params);

  pout() << "solving with AMRMultiGrid" << endl;
  //solve the equation ---- a_error now = velCalc
  solver.solve(a_error, klpExac, a_params.maxLevel, 0);

  const int exitStatus = solver.m_exitStatus;
  if(exitStatus!=1)
    {
      MayDay::Warning("non-zero exit status");
    }

  //subtract off velExact so  a_error  = velCalc-velExact
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBLevelDataOps::incr(*(a_error[ilev]),*(velExac[ilev]),-1.0);
    }

  //delete the local news.  error must survive outside this scope
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete klpExac[ilev];
      delete velExac[ilev];
    }

}
/***************/
void outputError(const Vector< LevelData<EBCellFAB>* >&   a_errorFine,
                 const Vector< LevelData<EBCellFAB>* >&   a_errorCoar,
                 const Vector< DisjointBoxLayout >&       a_gridsFine,
                 const Vector< DisjointBoxLayout >&       a_gridsCoar,
                 const Vector< EBISLayout >&              a_ebislFine,
                 const Vector< EBISLayout >&              a_ebislCoar,
                 const PoissonParameters&                 a_paramsFine,
                 const PoissonParameters&                 a_paramsCoar,
                 const string& a_fileFine,
                 const string& a_fileCoar)
{
#ifdef CH_USE_HDF5
  int nvar = SpaceDim;
  Vector<string> names(SpaceDim);
  for(int ivar = 0; ivar < SpaceDim; ivar++)
    {
      char charstr[100];
      sprintf(charstr, "var%d", ivar);
      names[ivar] = string(charstr);
    }
  bool replaceCovered = true;
  Vector<Real> coveredValues(nvar, 0.0);
  //values that don't matter in output file
  Real dxFine = 1.0;
  Real dxCoar = 1.0;
  Real time   = 1.0;
  Real dt     = 1.0;

  ParmParse pp;

  Vector<int> refRatio = a_paramsFine.refRatio;
  int numlevels = a_paramsFine.numLevels;
  ProblemDomain domainFine = a_paramsFine.coarsestDomain;
  ProblemDomain domainCoar = a_paramsCoar.coarsestDomain;

  writeEBHDF5(a_fileFine, a_gridsFine, a_errorFine, names,
              domainFine, dxFine, dt, time, refRatio, numlevels,
              replaceCovered, coveredValues);

  writeEBHDF5(a_fileCoar, a_gridsCoar, a_errorCoar, names,
              domainCoar, dxCoar, dt, time, refRatio, numlevels,
              replaceCovered, coveredValues);
#endif
}
int main(int argc, char* argv[])
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
  // registerDebugger();
#endif
  // Scoping trick
  {
#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif

#if CHECK_FLOATING_PT==1
    //    int except =  FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW |  FE_INVALID ;
    int except =  FE_DIVBYZERO | FE_OVERFLOW |  FE_INVALID ;
    feenableexcept(except);
#endif

    if (argc < 2)
      {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }

    char* inFile = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,inFile);

    PoissonParameters paramFine, paramCoar;
    Vector<DisjointBoxLayout> gridsFine, gridsCoar;
    Vector<EBISLayout>        ebislFine, ebislCoar;

    //read params from file
    getPoissonParameters(paramFine);
    paramCoar = paramFine;
    paramCoar.coarsen(2);

    int nlevels = paramCoar.numLevels;
    Vector<LevelData<EBCellFAB>* > errorFine(nlevels, NULL);
    Vector<LevelData<EBCellFAB>* > errorCoar(nlevels, NULL);

    //define geometry from given params
    pout() << "defining fine geometry" << endl;
    definePoissonGeometry(paramFine);

    getAllIrregRefinedLayouts(gridsFine, ebislFine, paramFine);

    pout() << "generating fine error" << endl;
    getError(errorFine, gridsFine, ebislFine, paramFine);

    pout() << "defining coarse geometry" << endl;
    definePoissonGeometry(paramCoar);

    getCoarseLayoutsFromFine(gridsCoar, ebislCoar, gridsFine, paramCoar);

    pout() << "generating coarse error" << endl;
    getError(errorCoar, gridsCoar, ebislCoar, paramCoar);

    compareError(errorFine,   errorCoar,
                 gridsFine,   gridsCoar,
                 ebislFine,   ebislCoar,
                 paramFine,   paramCoar);

#if CH_SPACEDIM==2
    string fileFine("pltFineError.2d.hdf5");
    string fileCoar("pltCoarError.2d.hdf5");
#else
    string fileFine("pltFineError.3d.hdf5");
    string fileCoar("pltCoarError.3d.hdf5");
#endif
    int dofileout;
    pp.get("do_error_output", dofileout);
    if(dofileout == 1)
      {
        outputError(errorFine,   errorCoar,
                    gridsFine,   gridsCoar,
                    ebislFine,   ebislCoar,
                    paramFine,   paramCoar,
                    fileFine,    fileCoar);
      }

    for(int ilev = 0; ilev < nlevels; ilev++)
      {
        delete errorFine[ilev];
        delete errorCoar[ilev];
        errorFine[ilev] = NULL;
        errorCoar[ilev] = NULL;
      }

  }// End scoping trick
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
}
