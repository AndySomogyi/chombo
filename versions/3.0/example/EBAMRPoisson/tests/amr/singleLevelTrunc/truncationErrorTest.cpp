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
#include "EBPoissonOp.H"
#include "EBPoissonOpFactory.H"

#include "EBFABView.H"
#include "EBDebugDump.H"
#include "DebugDump.H"



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
  int nvar = 1;

  DisjointBoxLayout grids = a_grids[0];
  EBISLayout        ebisl = a_ebisl[0];

  EBLevelGrid eblg(grids, ebisl, a_params.coarsestDomain);

  pout() << "grids= " << grids << endl;

  //define  data
  IntVect nghostPhi = 4*IntVect::Unit;
  IntVect nghostRHS = IntVect::Zero;
  EBCellFactory ebcf(eblg.getEBISL());
  LevelData<EBCellFAB> phiExac(grids, nvar, nghostPhi, ebcf);
  LevelData<EBCellFAB> klpExac(grids, nvar, nghostRHS, ebcf);;
  a_error.resize(1, NULL);
  a_error[0] = new LevelData<EBCellFAB>(grids, nvar, nghostRHS, ebcf);
  LevelData<EBCellFAB>& error = *a_error[0];

  Real alpha,beta;
  ParmParse pp;
  pp.get("alpha", alpha);
  pp.get("beta", beta);

  //set phi = phiExact, rhs=lphiexact  This makes AMRResidual return lphiexact-Lphi
  setTrigPhi(         phiExac, a_params.coarsestDx, a_params);
  setTrigKappaLOfPhi (klpExac, a_params.coarsestDx, a_params, alpha, beta);

  //set alpha and beta for the operator (kappa*alpha*I + kappa*beta*Lap)

  //create the solver
  RefCountedPtr<BaseDomainBCFactory> domainBCFactory;
  RefCountedPtr<BaseEBBCFactory>     ebbcFactory;
  getBCFactories(domainBCFactory, ebbcFactory, a_params);
  int orderEB =2;
  int relaxType = 1;//gsrb
  int numPreCond;

  EBPoissonOpFactory factory(eblg, a_params.coarsestDx, RealVect::Zero, orderEB, numPreCond, relaxType,
                             domainBCFactory, ebbcFactory, alpha, beta, nghostPhi, nghostRHS);
  RefCountedPtr<EBPoissonOp>  ebpo = RefCountedPtr<EBPoissonOp>(factory.MGnewOp(a_params.coarsestDomain, 0, false));

  ebpo->residual(error, phiExac, klpExac, false);
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
  int nvar = 1;
  Vector<string> names(1, string("var0"));
  bool replaceCovered = true;
  Vector<Real> coveredValues(nvar, 0.0);
  //values that don't matter in output file
  Real dxFine = 1.0;
  Real dxCoar = 1.0;
  Real time   = 0.0;
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
    getPoissonParameters(paramFine, true);

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

    definePoissonGeometry(paramCoar);

    getCoarseLayoutsFromFine(gridsCoar, ebislCoar, gridsFine, paramCoar);

    pout() << "generating coarse error" << endl;
    getError(errorCoar, gridsCoar, ebislCoar, paramCoar);

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

    compareError(errorFine,   errorCoar,
                 gridsFine,   gridsCoar,
                 ebislFine,   ebislCoar,
                 paramFine,   paramCoar);

    for(int ilev = 0; ilev < nlevels; ilev++)
      {
        delete errorFine[ilev];
        delete errorCoar[ilev];
        errorFine[ilev] = NULL;
        errorCoar[ilev] = NULL;
      }

  }// End scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif
}
