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
  int nvar = 1;
  Vector<string> names(1, string("var0"));
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
  Box domainFine = a_paramsFine.coarsestDomain.domainBox();
  Box domainCoar = a_paramsCoar.coarsestDomain.domainBox();

  writeEBHDF5(a_fileFine, a_gridsFine, a_errorFine, names,
              domainFine, dxFine, dt, time, refRatio, numlevels,
              replaceCovered, coveredValues);

  writeEBHDF5(a_fileCoar, a_gridsCoar, a_errorCoar, names,
              domainCoar, dxCoar, dt, time, refRatio, numlevels,
              replaceCovered, coveredValues);
#endif
}
/******/
void getError(Vector< LevelData<EBCellFAB>* >&     a_error,
              const Vector< DisjointBoxLayout >&   a_grids,
              const Vector< EBISLayout >&          a_ebisl,
              const PoissonParameters&             a_params,
              bool isFineSolve)
{
  ParmParse pp;

  int nlevels = a_params.numLevels;
  a_error.resize(nlevels);
  Vector<LevelData<EBCellFAB>* > phiExac(nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > klpExac(nlevels, NULL);

  //set alpha and beta for the operator (kappa*alpha*I + kappa*beta*Lap)
  Real alpha = 0;
  Real beta = 1;
  pout() << "forcing alpha = 0 and beta = 1" << endl;

  RealVect dxLev = RealVect::Unit;
  dxLev *=a_params.coarsestDx;
  pout() << "coarsest dx = " << dxLev << "; dx/dz = " << dxLev[0]/dxLev[SpaceDim-1] << std::endl;
  ProblemDomain domLev(a_params.coarsestDomain);
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      a_error[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], 1, a_params.ghostPhi, ebcellfact);
      klpExac[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], 1, a_params.ghostRHS, ebcellfact);
      phiExac[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], 1, a_params.ghostPhi, ebcellfact);

      //set phi = phiExact, rhs=lphiexact  This makes AMRResidual return lphiexact-Lphi
      setMarshaPhi(         *phiExac[ilev], dxLev, a_params);
      setMarshaKappaLOfPhi (*klpExac[ilev], dxLev, a_params);
      EBLevelDataOps::setToZero(*(a_error[ilev]));

      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);
    }

  //DEBUG START
  //EBAMRDataOps::setVal(a_error,1.0);
  //DEBUG END

  //create the solver
  EBSimpleSolver bottomSolver;
  int numSmooths;
  pp.get("mg_num_smooths", numSmooths);
  bottomSolver.setNumSmooths(numSmooths);

  AMRMultiGrid<LevelData<EBCellFAB> > solver;
  pout() << "defining  AMRMultiGrid solver" << endl;
  defineSolver(solver, a_grids, a_ebisl, bottomSolver, a_params, 0.0,alpha, beta);

  pout() << "solving with AMRMultiGrid" << endl;
  //solve the equation ---- a_error now = phiCalc
  solver.solve(a_error, klpExac, a_params.maxLevel, 0);
  //solver.relaxOnly(a_error, klpExac, a_params.maxLevel, 0);

  string fileFine, fileCoar;
  if(isFineSolve)
    {
#if CH_SPACEDIM==2
      fileFine = string("phiCalcFine.2d.hdf5");
      fileCoar = string("phiExacFine.2d.hdf5");
#else
      fileFine = string("phiCalcFine.3d.hdf5");
      fileCoar = string("phiExacFine.3d.hdf5");
#endif
    }
  else
    {
#if CH_SPACEDIM==2
      fileFine = string("phiCalcCoar.2d.hdf5");
      fileCoar = string("phiExacCoar.2d.hdf5");
#else
      fileFine = string("phiCalcCoar.3d.hdf5");
      fileCoar = string("phiExacCoar.3d.hdf5");
#endif
    }

    int dofileout;
    pp.get("do_error_output", dofileout);
    if(dofileout == 1)
      {
        outputError(a_error,     phiExac,
                    a_grids,     a_grids,
                    a_ebisl,     a_ebisl,
                    a_params,    a_params,
                    fileFine,    fileCoar);
      }

  const int exitStatus = solver.m_exitStatus;
  CH_assert(exitStatus==1);

  //subtract off phiExact so  a_error  = phi-phiExact
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBLevelDataOps::incr(*(a_error[ilev]),*(phiExac[ilev]),-1.0);
    }

  //delete the local news.  error must survive outside this scope
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete klpExac[ilev];
      delete phiExac[ilev];
    }

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
    getError(errorFine, gridsFine, ebislFine, paramFine, true);

    pout() << "defining coarse geometry" << endl;
    definePoissonGeometry(paramCoar);

    getCoarseLayoutsFromFine(gridsCoar, ebislCoar, gridsFine, paramCoar);

    pout() << "generating coarse error" << endl;
    getError(errorCoar, gridsCoar, ebislCoar, paramCoar, false);

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
