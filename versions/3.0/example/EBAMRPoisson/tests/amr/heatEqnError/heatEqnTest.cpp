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
#include "AMRTGA.H"
#include "EBAMRDataOps.H"
#include "EBLevelDataOps.H"

#include "EBFABView.H"
#include "EBDebugDump.H"
#include "DebugDump.H"
#include "CH_Attach.H"





#define CHECK_FLOATING_PT 0

#if CHECK_FLOATING_PT==1
#include <fenv.h>
#endif

void makeSource(Vector< LevelData<EBCellFAB>* >& a_src,
                const PoissonParameters&         a_params,
                const Real&                      a_time,
                const Real&                      a_diffConst)
{//source term
  RealVect dxLev = a_params.coarsestDx;
  int nlevels = a_src.size();
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      setTrigSource(*a_src[ilev],
                    dxLev,
                    a_params,
                    a_diffConst,
                    a_time);
      dxLev /= a_params.refRatio[ilev];
    }
}

void interiorValue(Vector< LevelData<EBCellFAB>* >& a_phi,
                   const PoissonParameters&         a_params,
                   const Real&                      a_time)
{//set phi:: phi_t = diffConst*lap(phi) + src
  RealVect dxLev = a_params.coarsestDx;
  int nlevels = a_phi.size();
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      setTrigPhi(*a_phi[ilev],
                 dxLev,
                 a_params,
                 a_time);//initial conditions
      dxLev /= a_params.refRatio[ilev];
    }
}

/******/
void getError(Vector< LevelData<EBCellFAB>* >&     a_error,
              const Vector< DisjointBoxLayout >&   a_grids,
              const Vector< EBISLayout >&          a_ebisl,
              const PoissonParameters&             a_params,
              const int&                           a_nsteps)
{
  ParmParse pp;

  int nlevels = a_params.numLevels;
  a_error.resize(nlevels);
  Vector<LevelData<EBCellFAB>* > phiOld(nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > phiNew(nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > source(nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > phiExact(nlevels, NULL);

  //set the diffusion constant (phi_t = diffConst*lap(phi) + src)
  Real diffConst;
  pp.get("diffusionConstant", diffConst);

  //set initial, final, and currentTime
  Real initTime = 0.0;
  Real finalTime = 1.0;
  Real currentTime = initTime;
  Real dt = (finalTime - initTime)/a_nsteps;

  RealVect dxLev = RealVect::Unit;
  dxLev *=a_params.coarsestDx;
  ProblemDomain domLev(a_params.coarsestDomain);
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      a_error[ilev]    = new LevelData<EBCellFAB>(a_grids[ilev], 1,4*IntVect::Unit, ebcellfact);
      source[ilev]     = new LevelData<EBCellFAB>(a_grids[ilev], 1,4*IntVect::Unit, ebcellfact);
      phiOld[ilev]     = new LevelData<EBCellFAB>(a_grids[ilev], 1,4*IntVect::Unit, ebcellfact);
      phiNew[ilev]     = new LevelData<EBCellFAB>(a_grids[ilev], 1,4*IntVect::Unit, ebcellfact);
      phiExact[ilev]   = new LevelData<EBCellFAB>(a_grids[ilev], 1,4*IntVect::Unit, ebcellfact);

      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);
    }

  //initialize all to zero
  EBAMRDataOps::setToZero(a_error);
  EBAMRDataOps::setToZero(phiNew);
  EBAMRDataOps::setToZero(phiOld);
  EBAMRDataOps::setToZero(phiExact);
  EBAMRDataOps::setToZero(source);

  //initialize phiNew to phi exact at initTime
  interiorValue(phiNew,a_params,initTime);

  //create the solver
  int coarsestLevel = 0;
  BiCGStabSolver<LevelData<EBCellFAB> >          bottomSolver;
  RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > > TGA = newTGASolver(a_grids,
                                                                   a_ebisl,
                                                                   a_params,
                                                                   bottomSolver,
                                                                   coarsestLevel,
                                                                   diffConst,
                                                                   4*IntVect::Unit,
                                                                   4*IntVect::Unit);

  //Main loop
  int iter = 0;
  bool stepMore = true;
  int maxLev = phiNew.size() - 1;
  while(stepMore)
    {
      pout() << "Iteration #" << iter << "  Start Time: " << currentTime << "  To Time: " << currentTime + dt << "  dt: " << dt << endl;

      //copy the second argument onto the first
      EBAMRDataOps::assign(phiOld,phiNew);

      //make source at currentTime
      {
        Real time_avg = currentTime + dt/2.0;
        makeSource(source,
                   a_params,
                   time_avg,
                   diffConst);
      }

      EBAMRDataOps::setToZero(phiNew);

      TGA->oneStep(phiNew,
                   phiOld,
                   source,
                   dt, 0, maxLev);

      EBAMRDataOps::setCoveredVal(phiNew,0.0);
      Vector<EBLevelGrid> eblg(a_grids.size());
      ProblemDomain domLev = a_params.coarsestDomain;
      for(int ilev = 0; ilev < a_grids.size(); ilev++)
        {

          eblg[ilev] = EBLevelGrid(a_grids[ilev], a_ebisl[ilev], domLev);
        }
      EBAMRDataOps::setCoveredAMRVal(phiNew,eblg,a_params.refRatio,0.0);

      Real newTime = currentTime + dt;
      interiorValue(phiExact,
                    a_params,
                    newTime);

      //compute error
      //subtract off phiExact so  a_error  = phiNew - phiExact
      EBAMRDataOps::axby(a_error,phiNew,phiExact,1.0,-1.0);
      EBAMRDataOps::setCoveredVal(a_error,0.0);
      //update time
      currentTime += dt;
      iter += 1;

      if(iter == a_nsteps)
        {
          stepMore = false;
        }
    }

  //delete the local news.  error must survive outside this scope
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete source[ilev];
      delete phiOld[ilev];
      delete phiNew[ilev];
      delete phiExact[ilev];
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

    //number of timesteps to use for conv study
    int nCoarSteps;
    pp.get("numCoarSteps", nCoarSteps);
    int nFineSteps;
    if(nCoarSteps<0)
      {
        nFineSteps = 1;
        nCoarSteps = 1;
      }
    else if (nCoarSteps==0)
      {
        MayDay::Error("bad numCoarSteps in heatEqnTest.cpp");
      }
    else
      {
        nFineSteps = nCoarSteps*2;
      }


    int nlevels = paramCoar.numLevels;
    Vector<LevelData<EBCellFAB>* > errorFine(nlevels, NULL);
    Vector<LevelData<EBCellFAB>* > errorCoar(nlevels, NULL);

    //define geometry from given params
    pout() << "defining fine geometry" << endl;
    definePoissonGeometry(paramFine);

    getAllIrregRefinedLayouts(gridsFine, ebislFine, paramFine);

    pout() << "generating fine error" << endl;
    getError(errorFine, gridsFine, ebislFine, paramFine, nFineSteps);

    pout() << "defining coarse geometry" << endl;
    definePoissonGeometry(paramCoar);

    getCoarseLayoutsFromFine(gridsCoar, ebislCoar, gridsFine, paramCoar);

    pout() << "generating coarse error" << endl;
    getError(errorCoar, gridsCoar, ebislCoar, paramCoar, nCoarSteps);

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
