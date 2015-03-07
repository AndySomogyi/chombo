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
#include "EBLevelMACProjector.H"



#define CHECK_FLOATING_PT 0

#if CHECK_FLOATING_PT==1
#include <fenv.h>
#endif
Real getVelExact(const FaceIndex& a_face, const Real& a_dx, const Real& a_domLen, const Real& a_freq,
                 const Real& a_magnitude)
{
  Real pi = 4.*atan(1.0);
  Real exactVal;
  IntVect ivFace = a_face.gridIndex(Side::Hi);
  Real xval = a_dx*ivFace[a_face.direction()];
  exactVal = sin(a_freq*pi*xval/a_domLen);

  //debug
  if(a_face.direction() == 0)
    exactVal = xval;
  else
    exactVal = 0.0;
  //end debug

  return exactVal;
}
/******/
void setExactVeloc(LevelData<EBFluxFAB>&                 a_veloc,
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
  int dummyIter = 0;
  for(DataIterator dit = a_grids.dataIterator(); dit.ok();++dit)
    {
      EBFluxFAB& velFluxFAB = a_veloc[dit()];
      IntVectSet ivsBox = IntVectSet(velFluxFAB.getRegion());
      const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          EBFaceFAB& vel = velFluxFAB[idir];
          vel.setVal(0.); //has to be zero on domain bndries
          for(FaceIterator faceit(ivsBox, ebgraph, idir, FaceStop::SurroundingNoBoundary); faceit.ok(); ++faceit)
            {
              Real domainLength = a_params.domainLength[idir];
              Real velexact = getVelExact(faceit(), a_dx, domainLength, frequencies[idir], magnitudes[idir]);
              vel(faceit(), 0) = velexact;
            }
          dummyIter++;
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
  Vector<LevelData<EBFluxFAB>* > velo(nlevels, NULL);
  Vector<LevelData<EBFluxFAB>* > gphi(nlevels, NULL);

  RealVect dxLev = RealVect::Unit;
  dxLev *=a_params.coarsestDx;
  ProblemDomain domLev(a_params.coarsestDomain);
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      pout() << "creating data for level " << ilev << endl;
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      EBFluxFactory ebfluxfact(a_ebisl[ilev]);
      a_error[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], 1,  IntVect::Unit, ebcellfact);

      velo[ilev] = new LevelData<EBFluxFAB>(a_grids[ilev], 1,  2*IntVect::Unit, ebfluxfact);
      gphi[ilev] = new LevelData<EBFluxFAB>(a_grids[ilev], 1,  2*IntVect::Unit, ebfluxfact);

      //set phi = phiExact, rhs=lphiexact  This makes AMRResidual return lphiexact-Lphi
      setExactVeloc(*velo[ilev], a_grids[ilev], a_ebisl[ilev], dxLev[0], a_params);

      EBLevelDataOps::setToZero(*(a_error[ilev]));

      pout() << "creating MAC projector for level " << ilev << endl;

      //solid wall bcs assumed for now
      NeumannPoissonEBBCFactory *        derebbcPhi = new NeumannPoissonEBBCFactory();
      NeumannPoissonDomainBCFactory *   derdombcPhi = new NeumannPoissonDomainBCFactory();
      DirichletPoissonDomainBCFactory * derdombcVel = new DirichletPoissonDomainBCFactory();
      derdombcPhi->setValue(0.);
      derebbcPhi->setValue(0.);
      derdombcVel->setValue(0.);

      RefCountedPtr<BaseEBBCFactory>      ebbcPhi = RefCountedPtr<BaseEBBCFactory>(derebbcPhi);
      RefCountedPtr<BaseDomainBCFactory> dombcPhi = RefCountedPtr<BaseDomainBCFactory>(derdombcPhi);
      RefCountedPtr<BaseDomainBCFactory> dombcVel = RefCountedPtr<BaseDomainBCFactory>(derdombcVel);
      int maxDepth = -1;
      int numSmooth, iterMax, mgCycle;
      Real tolerance;
      ParmParse pp;
      pp.get("mac_num_smooth", numSmooth);
      pp.get("mac_max_iter", iterMax);
      pp.get("mac_mg_cycle", mgCycle);
      pp.get("mac_tolerance", tolerance);
      BiCGStabSolver<LevelData<EBCellFAB> > bottomSolver;

      EBLevelMACProjector projectinator(a_grids[ilev], a_ebisl[ilev], domLev, dxLev, RealVect::Zero,bottomSolver,
                                        ebbcPhi,  dombcPhi, dombcVel, numSmooth, mgCycle,
                                        iterMax, tolerance, maxDepth, 1.e99,2*IntVect::Unit, IntVect::Zero);


      pout() << "projecting velocity for level " << ilev << endl;
      projectinator.project(*velo[ilev], *gphi[ilev], NULL);

      pout() << "computing error = kappa*div(u) for level " << ilev << endl;
      macKappaDivergence(*a_error[ilev],  *velo[ilev],
                         a_grids[ilev], a_ebisl[ilev], domLev,
                         dxLev, NULL);

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
