#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <cmath>
#include <iostream>
using std::cerr;

#include "ParmParse.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"
#include "EBFluxFAB.H"
#include "EBFluxFactory.H"

#include "PoissonUtilities.H"
#include "EBAMRPoissonOp.H"

#include "EBFABView.H"
#include "EBDebugDump.H"
#include "DebugDump.H"
#include "EBSimpleSolver.H"
#include "EBLevelCCProjector.H"



#define CHECK_FLOATING_PT 0

#if CHECK_FLOATING_PT==1
#include <fenv.h>
#endif
Real getVelExact(const VolIndex& a_vof, const Real& a_dx, const Real& a_domLen, const Real& a_freq,
                 const Real& a_magnitude, int a_idir)
{
  Real pi = 4.*atan(1.0);
  Real exactVal;
  IntVect ivFace = a_vof.gridIndex();
  Real xval = a_dx*ivFace[a_idir];
  exactVal = sin(a_freq*pi*xval/a_domLen);

  return exactVal;
}
/******/
void setExactVeloc(LevelData<EBCellFAB>&                 a_veloc,
                   const DisjointBoxLayout&              a_grids,
                   const EBISLayout&                     a_ebisl,
                   const Real&                           a_dx,
                   const PoissonParameters&              a_params)
{
  ParmParse pp;
  Vector<Real> frequencies(SpaceDim, 1.0);
  Vector<Real> magnitudes(SpaceDim, 1.0);
  pp.getarr("velocity_frequencies", frequencies, 0, SpaceDim);
  pp.getarr("velocity_magnitudes", magnitudes, 0, SpaceDim);
  for(DataIterator dit = a_grids.dataIterator(); dit.ok();++dit)
    {
      EBCellFAB& vel = a_veloc[dit()];
      vel.setVal(0.);
      IntVectSet ivsBox = IntVectSet(vel.getRegion());
      const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
      for(VoFIterator vofit(ivsBox, ebgraph); vofit.ok(); ++vofit)
        {
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              Real domainLength = a_params.domainLength[idir];
              Real velexact = getVelExact(vofit(), a_dx, domainLength, frequencies[idir], magnitudes[idir], idir);
              vel(vofit(), idir) = velexact;
            }
        }
    }
}
/******/
void getError(Vector< LevelData<EBCellFAB>* >&     a_error,
              const Vector< DisjointBoxLayout >&   a_grids,
              const Vector< EBISLayout >&          a_ebisl,
              const PoissonParameters&             a_params)
{
  int nlevels = a_params.numLevels;
  a_error.resize(nlevels);
  Vector<LevelData<EBCellFAB>* > velo(nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > gphi(nlevels, NULL);

  RealVect dxLevCoarsest = RealVect::Unit;
  dxLevCoarsest *=a_params.coarsestDx;
  ProblemDomain domLevCoarsest(a_params.coarsestDomain);

  RealVect dxLev = dxLevCoarsest;
  ProblemDomain domLev = domLevCoarsest;
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      pout() << "creating data for level " << ilev << endl;
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      a_error[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], 1       ,  IntVect::Unit,   ebcellfact);
      velo[ilev]    = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,  2*IntVect::Unit, ebcellfact);
      gphi[ilev]    = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,  2*IntVect::Unit, ebcellfact);

      //set phi = phiExact, rhs=lphiexact  This makes AMRResidual return lphiexact-Lphi
      setExactVeloc(*velo[ilev], a_grids[ilev], a_ebisl[ilev], dxLev[0], a_params);

      EBLevelDataOps::setToZero(*(a_error[ilev]));
      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);

    }

  dxLev = dxLevCoarsest;
  domLev = domLevCoarsest;
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      pout() << "creating cell-centered projector for level " << ilev << endl;
      DisjointBoxLayout gridsCoar;
      EBISLayout        ebislCoar;
      bool hasCoarser= false;
      int refToCoar = 1;
      LevelData<EBCellFAB> dummyLevelData;
      LevelData<EBCellFAB>* velCoarPtr = &dummyLevelData;
      if(ilev > 0)
        {
          hasCoarser = true;
          gridsCoar  = a_grids[ilev-1];
          ebislCoar  = a_ebisl[ilev-1];
          refToCoar  = a_params.refRatio[ilev-1];
          velCoarPtr = velo[ilev-1];
          //add gradient back into velocity so cf boundary conditions will make sense
          EBLevelDataOps::incr(*(velo[ilev-1]), *(gphi[ilev-1]), 1.0);
        }
      LevelData<EBCellFAB>& velCoar  = *velCoarPtr;

      NeumannPoissonEBBCFactory *        derebbcPhi = new NeumannPoissonEBBCFactory();
      NeumannPoissonDomainBCFactory *   derdombcPhi = new NeumannPoissonDomainBCFactory();
      DirichletPoissonDomainBCFactory * derdombcVel = new DirichletPoissonDomainBCFactory();
      derebbcPhi->setValue(0.);
      derdombcPhi->setValue(0.);
      derdombcVel->setValue(0.);

      RefCountedPtr<BaseEBBCFactory>      ebbcPhi = RefCountedPtr<BaseEBBCFactory    >(derebbcPhi );
      RefCountedPtr<BaseDomainBCFactory> dombcPhi = RefCountedPtr<BaseDomainBCFactory>(derdombcPhi);
      RefCountedPtr<BaseDomainBCFactory> dombcVel = RefCountedPtr<BaseDomainBCFactory>(derdombcVel);
      int maxDepth = -1;
      int numSmooth = 4;
      int iterMax   = 27;
      int mgCycle   = 1;
      Real tolerance = 1.0e-10;
      BiCGStabSolver<LevelData<EBCellFAB> > bottomSolver;
      bottomSolver.m_verbosity = 0;
      ParmParse pp;
      EBLevelCCProjector projectinator(a_grids[ilev], gridsCoar,
                                       a_ebisl[ilev], ebislCoar,
                                       domLev, dxLev, RealVect::Zero,
                                       refToCoar, hasCoarser, bottomSolver,
                                       ebbcPhi,  dombcPhi, dombcVel,numSmooth, mgCycle,
                                       iterMax, tolerance, maxDepth, 1.e99,false,
                                       2*IntVect::Unit, IntVect::Zero);


      pout() << "projecting velocity for level " << ilev << endl;
      projectinator.project(*velo[ilev], *gphi[ilev], velCoar, NULL);

      pout() << "computing error = kappa*div(u) for level " << ilev << endl;
      projectinator.kappaDivergence(*a_error[ilev], *velo[ilev],  velCoar, NULL);

      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);
    }

  //delete the stuff that locally newed.  error must survive outside this scope
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete velo[ilev];
      delete gphi[ilev];
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
  Vector<string> names(1, string("divu"));
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

#ifdef CH_MPI
  MPI_Finalize();
#endif
}
