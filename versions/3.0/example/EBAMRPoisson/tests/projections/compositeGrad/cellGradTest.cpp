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
Real getPhiExact(const VolIndex& a_vof, const Real& a_dx, const Real& a_domLen,
                 const Vector<Real>& a_freq)
{
  Real retval = 0;
  Real pi = 4.*atan(1.0);
  IntVect iv = a_vof.gridIndex();
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      Real xval = a_dx*(0.5 + Real(iv[idir]));
      Real factor = a_freq[idir]*pi/a_domLen;
      Real dirVal = sin(factor*xval);
      retval += dirVal;
    }

  //debug
  //  Real xval = a_dx*(Real(iv[0]) + 0.5);
  //  retval = 0.5*xval*xval;
  //end debug

  return retval;
}
/******/
Real getGradExact(const VolIndex& a_vof, const Real& a_dx, const Real& a_domLen, const Real& a_freq, int a_idir)
{
  Real pi = 4.*atan(1.0);
  Real retval;
  IntVect iv = a_vof.gridIndex();
  Real xval = a_dx*(0.5+ Real(iv[a_idir]));
  Real factor = a_freq*pi/a_domLen;
  retval = factor*cos(factor*xval);

  //debug

  //  if(a_idir == 0)
  //    retval = xval;
  //  else
  //    retval = 0.0;

  //end debug

  return retval;
}
/******/
void setExactGrad(LevelData<EBCellFAB>&                 a_grad,
                  const DisjointBoxLayout&              a_grids,
                  const EBISLayout&                     a_ebisl,
                  const Real&                           a_dx,
                  const PoissonParameters&              a_params)
{
  ParmParse pp;
  Vector<Real> frequencies(SpaceDim, 1.0);
  pp.getarr("phi_frequencies", frequencies, 0, SpaceDim);
  for(DataIterator dit = a_grids.dataIterator(); dit.ok();++dit)
    {
      EBCellFAB& grad = a_grad[dit()];
      grad.setVal(0.);
      IntVectSet ivsBox = IntVectSet(grad.getRegion());
      const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
      for(VoFIterator vofit(ivsBox, ebgraph); vofit.ok(); ++vofit)
        {
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              Real domainLength = a_params.domainLength[idir];
              Real gradexact = getGradExact(vofit(), a_dx, domainLength, frequencies[idir], idir);
              grad(vofit(), idir) = gradexact;
            }
        }
    }
}

/******/
void setExactPhi(LevelData<EBCellFAB>&                 a_phi,
                 const DisjointBoxLayout&              a_grids,
                 const EBISLayout&                     a_ebisl,
                 const Real&                           a_dx,
                 const PoissonParameters&              a_params)
{
  ParmParse pp;
  Vector<Real> frequencies(SpaceDim, 1.0);
  pp.getarr("phi_frequencies", frequencies, 0, SpaceDim);
  for(DataIterator dit = a_grids.dataIterator(); dit.ok();++dit)
    {
      EBCellFAB& phi = a_phi[dit()];
      phi.setVal(0.);
      IntVectSet ivsBox = IntVectSet(phi.getRegion());
      const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
      for(VoFIterator vofit(ivsBox, ebgraph); vofit.ok(); ++vofit)
        {
          Real domainLength = a_params.domainLength[0];
          Real phiexact = getPhiExact(vofit(), a_dx, domainLength, frequencies);
          phi(vofit(), 0) = phiexact;
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
  Vector<LevelData<EBCellFAB>* >  phi(nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > gphi(nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > gphiExact(nlevels, NULL);

  RealVect dxLevCoarsest = RealVect::Unit;
  dxLevCoarsest *=a_params.coarsestDx;
  ProblemDomain domLevCoarsest(a_params.coarsestDomain);

  RealVect dxLev = dxLevCoarsest;
  ProblemDomain domLev = domLevCoarsest;
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      pout() << "creating data for level " << ilev << endl;
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      a_error[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim   ,  IntVect::Unit,   ebcellfact);
      phi[ilev]    = new LevelData<EBCellFAB>(a_grids[ilev], 1,  2*IntVect::Unit, ebcellfact);
      gphi[ilev]    = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,  IntVect::Zero, ebcellfact);
      gphiExact[ilev]    = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,  IntVect::Zero, ebcellfact);

      //set phi = phiExact,gphi =gphiexact  This makes AMRResidual return lphiexact-Lphi
      setExactPhi(*phi[ilev], a_grids[ilev], a_ebisl[ilev], dxLev[0], a_params);
      setExactGrad(*gphiExact[ilev], a_grids[ilev], a_ebisl[ilev], dxLev[0], a_params);

      EBLevelDataOps::setToZero(*(a_error[ilev]));
      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);

    }


  Real domVal = 0.0;
  NeumannPoissonDomainBCFactory* domBCPhi = new NeumannPoissonDomainBCFactory();
  RefCountedPtr<BaseDomainBCFactory> baseDomainBCPhi = RefCountedPtr<BaseDomainBCFactory>(domBCPhi);
  domBCPhi->setValue(domVal);
  NeumannPoissonEBBCFactory* EBBCPhi = new NeumannPoissonEBBCFactory();
  RefCountedPtr<BaseEBBCFactory> baseEBBCPhi = RefCountedPtr<BaseEBBCFactory>(EBBCPhi);
  EBBCPhi->setValue(domVal);
  DirichletPoissonDomainBCFactory* domBCVel = new DirichletPoissonDomainBCFactory();
  RefCountedPtr<BaseDomainBCFactory> baseDomainBCVel = RefCountedPtr<BaseDomainBCFactory>(domBCVel);
  domBCVel->setValue(domVal);


  ParmParse pp;
  Vector<EBLevelGrid>                      eblg   (a_grids.size());
  Vector<RefCountedPtr<EBQuadCFInterp> >   quadCFI(a_grids.size());
  domLev = domLevCoarsest;
  for(int ilev = 0; ilev < a_grids.size(); ilev++)
    {
      int nvar = 1;
      int nref = a_params.refRatio[ilev];
      eblg[ilev] = EBLevelGrid(a_grids[ilev], a_ebisl[ilev], domLev);
      if(ilev > 0)
        {
          int nrefOld = a_params.refRatio[ilev-1];
          ProblemDomain domLevCoar = coarsen(domLev, nrefOld);
          quadCFI[ilev] = RefCountedPtr<EBQuadCFInterp>(new EBQuadCFInterp(a_grids[ilev  ],
                                                                           a_grids[ilev-1],
                                                                           a_ebisl[ilev  ],
                                                                           a_ebisl[ilev-1],
                                                                           domLevCoar,
                                                                           nrefOld, nvar,
                                                                           *(eblg[ilev].getCFIVS())));

        }
      domLev.refine(nref);
    }


  EBCompositeCCProjector projectinator(eblg,  a_params.refRatio, quadCFI,
                                       a_params.coarsestDx,
                                       RealVect::Zero,
                                       baseDomainBCVel,
                                       baseDomainBCPhi,
                                       baseEBBCPhi,
                                       true, -1, 3 ,40,1.e99, 1, 0);


  pout() << "taking gradient  "  << endl;
  projectinator.gradient(gphi, phi);
  pout() << "subtracting gradient  "  << endl;
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBLevelDataOps::incr(*(a_error[ilev]), *(gphi[ilev]), 1.0);
      EBLevelDataOps::incr(*(a_error[ilev]), *(gphiExact[ilev]), -1.0);
    }

  ProblemDomain domainLev = a_params.coarsestDomain;
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*a_error[ilev])[dit()].setCoveredCellVal(0.0, 0);
        }
      if(ilev < nlevels-1)
        {
          //zero out stuff covered by finer levels
          const EBIndexSpace* const  ebisPtr = Chombo_EBIS::instance();
          DisjointBoxLayout gridsCoarsenedFine;
          coarsen(gridsCoarsenedFine, a_grids[ilev+1], a_params.refRatio[ilev]);

          int numGhost = 0;
          EBISLayout ebislCoarsenedFine;
          ebisPtr->fillEBISLayout(ebislCoarsenedFine,
                                  gridsCoarsenedFine,
                                  domainLev,
                                  numGhost);

          EBCellFactory ebcellfact(ebislCoarsenedFine);
          LevelData<EBCellFAB> zeroLD(gridsCoarsenedFine, 1, IntVect::Zero, ebcellfact);
          for(DataIterator dit = a_grids[ilev+1].dataIterator(); dit.ok(); ++dit)
            {
              zeroLD[dit()].setVal(0.);
            }
          Interval interv(0,0);
          zeroLD.copyTo(interv, (*a_error[ilev]), interv);
        }
      domainLev.refine(a_params.refRatio[ilev]);
    }

  //delete the stuff that locally newed.  error must survive outside this scope
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete phi[ilev];
      delete gphi[ilev];
      delete gphiExact[ilev];
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
      sprintf(charstr, "graderr%d", idir);
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
