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
#include "EBAMRDataOps.H"



#define CHECK_FLOATING_PT 0

#if CHECK_FLOATING_PT==1
#include <fenv.h>
#endif

/******/
void getError(Vector< LevelData<EBCellFAB>* >&     a_error,
              const Vector< DisjointBoxLayout >&   a_grids,
              const Vector< EBISLayout >&          a_ebisl,
              const PoissonParameters&             a_params,
              bool a_isFine)
{
  ParmParse pp;
  //set alpha and beta for the operator (kappa*alpha*I + kappa*beta*Lap)
  Real alpha,beta;
  pp.get("alpha", alpha);
  pp.get("beta", beta);

  int nlevels = a_params.numLevels;
  a_error.resize(nlevels);
  Vector<LevelData<EBCellFAB>* > phiExac(nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > klpExac(nlevels, NULL);
  Vector<EBAMRPoissonOp*>        ebamrpo(nlevels, NULL);

  RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > levelOpFactory;
  //placeholder---no relaxation, no preconditioning here.
  int relaxtype = 0;
  int numPreCondIters = 1;
  getEBAMRPFactory(levelOpFactory, a_grids, a_ebisl, a_params,
                   numPreCondIters,relaxtype, 0.0,alpha, beta);

  RealVect dxLev = RealVect::Unit;
  dxLev *=a_params.coarsestDx;
  ProblemDomain domLev(a_params.coarsestDomain);
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      a_error[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], 1,a_params.ghostRHS, ebcellfact);
      klpExac[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], 1,a_params.ghostRHS, ebcellfact);
      phiExac[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], 1,a_params.ghostPhi, ebcellfact);

      ebamrpo[ilev] = (EBAMRPoissonOp*)(levelOpFactory->AMRnewOp(domLev));

      //set phi = phiExact, rhs=lphiexact  This makes AMRResidual return lphiexact-Lphi
      if ( (a_params.ebBcType ==2 || a_params.ebBcType ==3) &&
           (a_params.domBcType==2 || a_params.domBcType==3) )
        {
          setTrigPhi(         *phiExac[ilev], dxLev, a_params);
          setTrigKappaLOfPhi (*klpExac[ilev], dxLev, a_params, alpha, beta);
        }
      else if ( (a_params.ebBcType ==8 || a_params.ebBcType ==9) &&
                (a_params.domBcType==8 || a_params.domBcType==9) )
        {
          setSphericalHarmonicPhi(         *phiExac[ilev], dxLev, a_params);
          setSphericalHarmonicKappaLOfPhi (*klpExac[ilev], dxLev, a_params, alpha, beta);
        }
      else
        {
          MayDay::Error("bogus EB or Domain BC type in input file");
        }

      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);
    }

  EBAMRDataOps::exchangeAll(phiExac);
  EBAMRDataOps::exchangeAll(klpExac);

#ifdef CH_USE_HDF5
  if (a_isFine)
    {
      writeEBAMRname(&phiExac,"phiExacFine.hdf5");
      writeEBAMRname(&klpExac,"klpExacFine.hdf5");
    }
  else
    {
      writeEBAMRname(&phiExac,"phiExacCoar.hdf5");
      writeEBAMRname(&klpExac,"klpExacCoar.hdf5");
    }
#endif

  LevelData<EBCellFAB> dummy;
  //create calculated data and error
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      const LevelData<EBCellFAB>*   phiFinePtr = &dummy;
      const LevelData<EBCellFAB>*   phiCoarPtr = &dummy;
      if(ilev > 0)
        {
          phiCoarPtr =  phiExac[ilev-1];
        }
      if(ilev < nlevels - 1)
        {
          phiFinePtr =  phiExac[ilev+1];
        }

      const LevelData<EBCellFAB>& phiFine = *phiFinePtr;
      const LevelData<EBCellFAB>& phiCoar = *phiCoarPtr;

      //inhomogeneous bcs
      if(ilev == (nlevels -1))
        {
          ebamrpo[ilev]->AMRResidual(*a_error[ilev], phiFine, *phiExac[ilev], phiCoar, *klpExac[ilev], false, NULL);
        }
      else
        {
          ebamrpo[ilev]->AMRResidual(*a_error[ilev], phiFine, *phiExac[ilev], phiCoar, *klpExac[ilev], false, ebamrpo[ilev+1]);
        }

    }

  //delete the local news.  error must survive outside this scope
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete klpExac[ilev];
      delete phiExac[ilev];
      delete ebamrpo[ilev];
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
  int nvar = 1;
  Vector<string> names(1, string("var0"));
  bool replaceCovered = true;
  Vector<Real> coveredValues(nvar, 0.0);
  //values that don't matter in output file
  Real dxFine = a_paramsFine.coarsestDx[0];
  Real dxCoar = a_paramsCoar.coarsestDx[0];
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

    getCoarseLayoutsFromFine(gridsCoar, ebislCoar, gridsFine, paramCoar);

    pout() << "generating coarse error" << endl;
    getError(errorCoar, gridsCoar, ebislCoar, paramCoar, false);

    pout() << "generating fine error" << endl;
    getError(errorFine, gridsFine, ebislFine, paramFine, true);


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
