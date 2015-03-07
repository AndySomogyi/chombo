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
#include "BaseIVFAB.H"
#include "BaseIVFactory.H"



#define CHECK_FLOATING_PT 0

#if CHECK_FLOATING_PT==1
#include <fenv.h>
#endif
/******/
RealVect getXValVoF(const VolIndex a_vof, const Real& a_dx,
                    const RealVect& a_centroid)
{
  RealVect retval;
  IntVect iv = a_vof.gridIndex();
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      Real xval = a_dx*(0.5 + a_centroid[idir]+ Real(iv[idir]));
      retval[idir] = xval;
    }
  return retval;
}
Real getDivGradExact(const VolIndex& a_vof, const Real& a_dx, const Real& a_domLen,
                     const Vector<Real>& a_freq)
{
  Real retval = 0;
  Real pi = 4.*atan(1.0);
  RealVect xval = getXValVoF(a_vof, a_dx, RealVect::Zero);
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      Real factor = a_freq[idir]*pi/a_domLen;
      Real dirVal = factor*cos(factor*xval[idir]);
      retval += dirVal;
    }

  ////debug
  //retval = 1.0 - 2.0*xval[0];
  ////end debug

  return retval;
}
/******/
Real getGradExact(RealVect a_xval, const Real& a_domLen, const Real& a_freq, int a_idir)
{
  Real pi = 4.*atan(1.0);
  Real retval;
  Real factor = a_freq*pi/a_domLen;
  retval = sin(factor*a_xval[a_idir]);

  //debug

  //if(a_idir == 0)
  //  retval = (1.0 - a_xval[a_idir])*a_xval[a_idir];
  //else
  //  retval = 0.0;

  //end debug

  return retval;
}
/******/
void setExactGradient(LevelData<EBCellFAB>&                 a_grad,
                      LevelData<BaseIVFAB<Real> >&          a_bndryGrad,
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
          RealVect xval = getXValVoF(vofit(), a_dx, RealVect::Zero);
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              Real domainLength = a_params.domainLength[idir];
              Real gradexact = getGradExact(xval, domainLength, frequencies[idir], idir);
              grad(vofit(), idir) = gradexact;
            }
        }
      IntVectSet ivsIrreg = ebgraph.getIrregCells(a_grids.get(dit()));
      BaseIVFAB<Real>& bndryGrad = a_bndryGrad[dit()];
      bndryGrad.setVal(0.);
      for(VoFIterator vofit(ivsIrreg, ebgraph); vofit.ok(); ++vofit)
        {
          RealVect bndryCent = a_ebisl[dit()].bndryCentroid(vofit());
          RealVect xval = getXValVoF(vofit(), a_dx, bndryCent);
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              Real domainLength = a_params.domainLength[idir];
              Real gradexact = getGradExact(xval, domainLength, frequencies[idir], idir);
              bndryGrad(vofit(), idir) = gradexact;
            }
        }
    }
}

/******/
void setExactKappaDiv(LevelData<EBCellFAB>&                 a_kappaDivGrad,
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
      EBCellFAB& kdivgrad = a_kappaDivGrad[dit()];
      kdivgrad.setVal(0.);
      IntVectSet ivsBox = IntVectSet(kdivgrad.getRegion());
      const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
      for(VoFIterator vofit(ivsBox, ebgraph); vofit.ok(); ++vofit)
        {
          Real domainLength = a_params.domainLength[0];
          Real divgradexact = getDivGradExact(vofit(), a_dx, domainLength, frequencies);
          Real kappa = a_ebisl[dit()].volFrac(vofit());
          kdivgrad(vofit(), 0) = kappa*divgradexact;
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
  Vector<LevelData<EBCellFAB>* >          divGrad(nlevels, NULL);
  Vector<LevelData<EBCellFAB>* >             gphi(nlevels, NULL);
  Vector<LevelData<EBCellFAB>* >     divGradExact(nlevels, NULL);
  Vector<LevelData<BaseIVFAB<Real> >* > irregGrad(nlevels, NULL);
  Vector<LayoutData<IntVectSet>* >      irregSets(nlevels, NULL);
  RealVect dxLevCoarsest = RealVect::Unit;
  dxLevCoarsest *=a_params.coarsestDx;
  ProblemDomain domLevCoarsest(a_params.coarsestDomain);

  RealVect dxLev = dxLevCoarsest;
  ProblemDomain domLev = domLevCoarsest;
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      pout() << "creating irregular sets level " << ilev << endl;
      irregSets[ilev] = new LayoutData<IntVectSet>(a_grids[ilev]);
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*irregSets[ilev])[dit()] = a_ebisl[ilev][dit()].getIrregIVS(a_grids[ilev].get(dit()));
        }

      pout() << "creating data for level " << ilev << endl;
      BaseIVFactory<Real>  baseivfact(a_ebisl[ilev], *irregSets[ilev]);
      EBCellFactory        ebcellfact(a_ebisl[ilev]);
      irregGrad[ilev]    = new LevelData<BaseIVFAB<Real> >(a_grids[ilev], SpaceDim,   IntVect::Zero, baseivfact);
      a_error[ilev]      = new LevelData<EBCellFAB>       (a_grids[ilev], 1       ,   IntVect::Zero, ebcellfact);
      divGrad[ilev]      = new LevelData<EBCellFAB>       (a_grids[ilev], 1       ,   IntVect::Zero, ebcellfact);
      divGradExact[ilev] = new LevelData<EBCellFAB>       (a_grids[ilev], 1       ,   IntVect::Zero, ebcellfact);
      gphi[ilev]         = new LevelData<EBCellFAB>       (a_grids[ilev], SpaceDim, 2*IntVect::Unit, ebcellfact);

      //set phi = phiExact,gphi =gphiexact  This makes AMRResidual return lphiexact-Lphi
      setExactKappaDiv(*divGradExact[ilev],           a_grids[ilev], a_ebisl[ilev], dxLev[0], a_params);
      setExactGradient(*gphi[ilev], *irregGrad[ilev], a_grids[ilev], a_ebisl[ilev], dxLev[0], a_params);

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
          quadCFI[ilev] = RefCountedPtr<EBQuadCFInterp> (new EBQuadCFInterp(a_grids[ilev  ],
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

  pout() << "taking divergence  "  << endl;

  projectinator.kappaDivergence(divGrad, gphi, &irregGrad);

  pout() << "subtracting exact answer  "  << endl;

  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBLevelDataOps::incr(*(a_error[ilev]),      *(divGrad[ilev]),  1.0);
      EBLevelDataOps::incr(*(a_error[ilev]), *(divGradExact[ilev]), -1.0);
    }

  pout() << "zeroing invalid data in error data" << endl;
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

          int numGhost = 4;
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
      delete gphi[ilev];
      delete divGradExact[ilev];
      delete divGrad[ilev];
      delete irregSets[ilev];
      delete irregGrad[ilev];
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

  Vector<string> names(1, string("divuError"));

  bool replaceCovered = true;
  Vector<Real> coveredValues(1, 0.0);
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
