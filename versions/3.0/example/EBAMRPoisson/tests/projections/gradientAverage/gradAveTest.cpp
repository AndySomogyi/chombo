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
#include <cmath>
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
#include "EBCompositeCCProjector.H"



#define CHECK_FLOATING_PT 0

#if CHECK_FLOATING_PT==1
#include <fenv.h>
#endif
/******/
Real gradExact(const RealVect& a_xloc, const int& a_idir)
{
  Real retval;
  retval = a_xloc[a_idir]*a_xloc[a_idir]*a_xloc[a_idir];
  return retval;
}
Real gradExact(const VolIndex& a_vof, const RealVect& a_dx,  int a_idir)
{
  RealVect vofLoc = EBArith::getVofLocation(a_vof, a_dx, RealVect::Zero);

  Real retval = gradExact(vofLoc, a_idir);
  return retval;
}
Real gradExact(const FaceIndex& a_face, const RealVect& a_dx)
{
  RealVect faceLoc = EBArith::getFaceLocation(a_face, a_dx, RealVect::Zero);

  Real retval = gradExact(faceLoc, a_face.direction());
  return retval;
}
/******/
void setExactGrad(LevelData<EBCellFAB>&                 a_gradCell,
                  LevelData<EBFluxFAB>&                 a_gradFlux,
                  const DisjointBoxLayout&              a_grids,
                  const EBISLayout&                     a_ebisl,
                  const RealVect&                       a_dx,
                  const PoissonParameters&              a_params)
{
  for(DataIterator dit = a_grids.dataIterator(); dit.ok();++dit)
    {
      EBCellFAB& gradCell = a_gradCell[dit()];
      EBFluxFAB& gradFlux = a_gradFlux[dit()];
      gradCell.setVal(0.);
      gradFlux.setVal(0.);

      IntVectSet ivsBox = IntVectSet(gradCell.getRegion());
      const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
      for(VoFIterator vofit(ivsBox, ebgraph); vofit.ok(); ++vofit)
        {
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              Real gradExactVoF = gradExact(vofit(), a_dx,  idir);
              gradCell(vofit(), idir) = gradExactVoF;
              for(SideIterator sit; sit.ok(); ++sit)
                {
                  Vector<FaceIndex> faces = ebgraph.getFaces(vofit(), idir, sit());
                  for(int iface = 0; iface < faces.size(); iface++)
                    {
                      Real gradExactFace = gradExact(faces[iface], a_dx);
                      gradFlux[idir](faces[iface], 0) = gradExactFace;
                    }
                }
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
  Vector<LevelData<EBCellFAB>* > gphi(nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > gphiExactCell(nlevels, NULL);
  Vector<LevelData<EBFluxFAB>* > gphiExactFlux(nlevels, NULL);

  RealVect dxLevCoarsest = RealVect::Unit;
  dxLevCoarsest *=a_params.coarsestDx;

  RealVect dxLev = dxLevCoarsest;
  ProblemDomain domLev = (a_params.coarsestDomain);
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      pout() << "creating data for level " << ilev << endl;
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      EBFluxFactory ebfluxfact(a_ebisl[ilev]);
      a_error[ilev]          = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,  IntVect::Unit,   ebcellfact);
      gphi[ilev]             = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,  IntVect::Zero, ebcellfact);
      gphiExactCell[ilev]    = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,  IntVect::Zero, ebcellfact);
      gphiExactFlux[ilev]    = new LevelData<EBFluxFAB>(a_grids[ilev],        1, 4*IntVect::Unit, ebfluxfact);

      setExactGrad(*gphiExactCell[ilev],
                   *gphiExactFlux[ilev],
                   a_grids[ilev], a_ebisl[ilev], dxLev, a_params);

      ccpAverageFaceToCells(*gphi[ilev],
                            *gphiExactFlux[ilev],
                            a_grids[ilev],
                            a_ebisl[ilev],
                            domLev,
                            dxLev);

      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);

    }



  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBLevelDataOps::setToZero(*(a_error[ilev]));

      EBLevelDataOps::incr(*(a_error[ilev]), *(gphi[ilev]), 1.0);
      EBLevelDataOps::incr(*(a_error[ilev]), *(gphiExactCell[ilev]), -1.0);
    }


  //delete the stuff that locally newed.  error must survive outside this scope
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete gphi[ilev];
      delete gphiExactCell[ilev];
      delete gphiExactFlux[ilev];
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

  Vector<string> names(SpaceDim);
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      char charstr[80];
      sprintf(charstr, "gradAverageErr%d", idir);
      names[idir]  = string(charstr);
    }
  bool replaceCovered = true;
  Vector<Real> coveredValues(SpaceDim, 0.0);
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
