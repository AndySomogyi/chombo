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
#include "EBViscousTensorOpFactory.H"
#include "EBViscousTensorOp.H"

#include "EBFABView.H"
#include "EBDebugDump.H"
#include "DebugDump.H"
#include "EBAMRDataOps.H"


/******/
void getError(Vector< LevelData<EBCellFAB>* >&     a_error,
              const Vector< DisjointBoxLayout >&   a_grids,
              const Vector< EBISLayout >&          a_ebisl,
              const PoissonParameters&             a_params)
{
  ParmParse pp;
  //set alpha and beta for the operator (kappa*alpha*I + kappa*beta*Lap)
  Real alpha;
  pp.get("alpha", alpha);

  int nlevels = a_params.numLevels;
  a_error.resize(nlevels);
  Vector<LevelData<EBCellFAB>* > velExac(nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > klbExac(nlevels, NULL);
  Vector<EBViscousTensorOp*>     ebamrpo(nlevels, NULL);

  RealVect dxLev = RealVect::Unit;
  dxLev *=a_params.coarsestDx;
  ProblemDomain domLev(a_params.coarsestDomain);
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      a_error[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,a_params.ghostRHS, ebcellfact);
      klbExac[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,a_params.ghostRHS, ebcellfact);
      velExac[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,a_params.ghostPhi, ebcellfact);

      //set vel = velExact, rhs=lvelexact  This makes AMRResidual return lvelexact-Lvel
      setVelViscous(  *velExac[ilev], dxLev[0], a_params);
      setKLVViscous(  *klbExac[ilev], dxLev[0], a_params);
      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);
    }


  Vector<RefCountedPtr<LevelData<EBFluxFAB       > > >  eta,      lambda;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >  etaIrreg, lambdaIrreg;
  Vector<RefCountedPtr<LevelData<EBCellFAB       > > >  acoeffi(nlevels);
  RefCountedPtr<EBViscousTensorOpFactory> opFactory;
  defineViscousTensorCoef(   acoeffi, eta, lambda, etaIrreg, lambdaIrreg,  a_grids, a_ebisl, a_params);
  getEBVTOFactory(opFactory, acoeffi, eta, lambda, etaIrreg, lambdaIrreg,   a_grids, a_ebisl, a_params);

  domLev = a_params.coarsestDomain;
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      ebamrpo[ilev] = opFactory->AMRnewOp(domLev);
      domLev.refine(a_params.refRatio[ilev]);
    }

  EBAMRDataOps::exchangeAll(velExac);
  EBAMRDataOps::exchangeAll(klbExac);

  LevelData<EBCellFAB> dummy;
  //create calculated data and error
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      const LevelData<EBCellFAB>*   phiFinePtr = &dummy;
      const LevelData<EBCellFAB>*   phiCoarPtr = &dummy;
      if(ilev > 0)
        {
          phiCoarPtr =  velExac[ilev-1];
        }
      if(ilev < nlevels - 1)
        {
          phiFinePtr =  velExac[ilev+1];
        }

      const LevelData<EBCellFAB>& phiFine = *phiFinePtr;
      const LevelData<EBCellFAB>& phiCoar = *phiCoarPtr;

      //inhomogeneous bcs
      if(ilev == (nlevels -1))
        {
          ebamrpo[ilev]->AMRResidual(*a_error[ilev], phiFine, *velExac[ilev], phiCoar, *klbExac[ilev], false, NULL);
        }
      else
        {
          ebamrpo[ilev]->AMRResidual(*a_error[ilev], phiFine, *velExac[ilev], phiCoar, *klbExac[ilev], false, ebamrpo[ilev+1]);
        }

    }

  //delete the local news.  error must survive outside this scope
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete klbExac[ilev];
      delete velExac[ilev];
      delete ebamrpo[ilev];

      klbExac[ilev] = NULL;
      velExac[ilev] = NULL;
      ebamrpo[ilev] = NULL;
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
  Vector<string> names(nvar);
  for(int ivar = 0; ivar < nvar; ivar++)
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

    getCoarseLayoutsFromFine(gridsCoar, ebislCoar, gridsFine, paramCoar);

    pout() << "generating coarse error" << endl;
    getError(errorCoar, gridsCoar, ebislCoar, paramCoar);

    pout() << "generating fine error" << endl;
    getError(errorFine, gridsFine, ebislFine, paramFine);


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
