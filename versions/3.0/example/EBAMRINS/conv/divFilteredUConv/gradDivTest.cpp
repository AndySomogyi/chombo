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

#include "ParmParse.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"

#include "AMRINSUtils.H"
#include "EBAMRPoissonOp.H"

#include "EBFABView.H"
#include "EBDebugDump.H"
#include "DebugDump.H"
#include "EBGradDivFilter.H"
#include "EBCoarseAverage.H"
#include "NeumannPoissonDomainBC.H"
#include "DirichletPoissonDomainBC.H"
#include "EBCompositeCCProjector.H"
#include "EBAMRDataOps.H"
#include "EBArith.H"
#include "NeumannPoissonEBBC.H"
#include "UsingNamespace.H"

using std::cerr;

#define CHECK_FLOATING_PT 0

#if CHECK_FLOATING_PT==1
#include <fenv.h>
#endif

Real getVelExact(const VolIndex& a_vof, const Real& a_dx, const Real& a_freq,
                 const Real& a_magnitude, int a_idir)
{
  Real pi = 4.*atan(1.0);
  Real exactVal;
  IntVect ivVoFCC = a_vof.gridIndex();
  Real xval;
  xval = a_dx*(0.5 + Real(ivVoFCC[a_idir]));

  exactVal = sin(a_freq*pi*xval);

  return exactVal;
}
/******/
void setExactVeloc(LevelData<EBCellFAB>&                 a_veloc,
                   LevelData<EBFluxFAB>&                 a_fluxVel,
                   const DisjointBoxLayout&              a_grids,
                   const EBISLayout&                     a_ebisl,
                   const Real&                           a_dx,
                   const AMRParameters&                  a_params)
{
  ParmParse pp;
  Vector<Real> frequencies(SpaceDim, 1.0);
  Vector<Real> magnitudes(SpaceDim, 1.0);
  pp.getarr("velocity_frequencies", frequencies, 0, SpaceDim);
  pp.getarr("velocity_magnitudes", magnitudes, 0, SpaceDim);
  int checkerBoard;
  pp.get("use_checkerboard_velocity", checkerBoard);
  const bool useCheckerBoard = (checkerBoard==1);

  for(DataIterator dit = a_grids.dataIterator(); dit.ok();++dit)
    {
      EBCellFAB& vel = a_veloc[dit()];
      EBFluxFAB& fluxVel = a_fluxVel[dit()];
      //set interior to zero,
      //boundary to exact value
      fluxVel.setVal(0.);
      vel.setVal(0.);
      IntVectSet ivsBox = IntVectSet(vel.getRegion());
      const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
      for(VoFIterator vofit(ivsBox, ebgraph); vofit.ok(); ++vofit)
        {
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              Real velexact;
              if(useCheckerBoard)
                {
                  const IntVect& iv = vofit().gridIndex();
                  int isRedBlackTest = 0;
                  for (int idir = 0;idir<SpaceDim;idir++)
                    {
                      isRedBlackTest += iv[idir];
                    }
                  if(isRedBlackTest%2 == 0)
                    {//red
                      velexact = 1.0;
                    }
                  else
                    {//black
                      velexact = -1.0;
                    }
                }
              else
                {
                  velexact = getVelExact(vofit(), a_dx,  frequencies[idir], magnitudes[idir], idir);
                }
              vel(vofit(), idir) = velexact;
            }
        }
    }
}
/****/
void getNormDiv(Real                                       a_normVal[3],
                Vector<LevelData<EBCellFAB>* >&            a_error,
                Vector< LevelData<EBCellFAB>* >&           a_veloc,
                const Vector< DisjointBoxLayout >&         a_grids,
                const Vector< EBISLayout >&                a_ebisl,
                const ProblemDomain&                       a_level0Domain,
                const Real&                                a_level0Dx,
                const AMRParameters&                       a_params)
{
  EBAMRDataOps::setToZero(a_error);
  int nlevels = a_grids.size();
  Real domVal = 0.0;
  NeumannPoissonDomainBCFactory* domBCPhi = new NeumannPoissonDomainBCFactory();
  RefCountedPtr<BaseDomainBCFactory> baseDomainBCPhi = RefCountedPtr<BaseDomainBCFactory>(domBCPhi);
  domBCPhi->setValue(domVal);
  DirichletPoissonDomainBCFactory* domBCVel = new DirichletPoissonDomainBCFactory();
  RefCountedPtr<BaseDomainBCFactory> baseDomainBCVel = RefCountedPtr<BaseDomainBCFactory>(domBCVel);
  domBCVel->setValue(domVal);


  RealVect dxVect = a_level0Dx*RealVect::Unit;
  Vector<LevelData<EBCellFAB>*> rhoinv;
  NeumannPoissonEBBCFactory*      ebBCPhi = new NeumannPoissonEBBCFactory();
  ebBCPhi->setValue(domVal);
  RefCountedPtr<BaseEBBCFactory>     baseEBBCPhi     = RefCountedPtr<BaseEBBCFactory>(ebBCPhi);

  Vector<EBLevelGrid>                      eblg   (a_grids.size());
  Vector<RefCountedPtr<EBQuadCFInterp> >   quadCFI(a_grids.size());
  ProblemDomain  domLev = a_level0Domain;
  for(int ilev = 0; ilev < a_grids.size(); ilev++)
    {
      int nvar = 1;
      int nref = a_params.m_refRatio[ilev];
      eblg[ilev] = EBLevelGrid(a_grids[ilev], a_ebisl[ilev], domLev);
      if(ilev > 0)
        {
          int nrefOld = a_params.m_refRatio[ilev-1];
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


  const int bottomSolverType = 1;

  EBCompositeCCProjector projectinator(eblg,  a_params.m_refRatio, quadCFI,
                                       dxVect,RealVect::Zero,
                                       baseDomainBCVel,
                                       baseDomainBCPhi,
                                       baseEBBCPhi,
                                       true, -1, 3 ,40,1.e99, 1,
                                       bottomSolverType);

  projectinator.kappaDivergence(a_error, a_veloc);

  domLev = a_level0Domain;
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
          coarsen(gridsCoarsenedFine, a_grids[ilev+1], a_params.m_refRatio[ilev]);

          int numGhost = 0;
          EBISLayout ebislCoarsenedFine;
          ebisPtr->fillEBISLayout(ebislCoarsenedFine,
                                  gridsCoarsenedFine,
                                  domLev.domainBox(),
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
      domLev.refine(a_params.m_refRatio[ilev]);
    }
  int icomp = 0;
  for(int inorm = 0; inorm < 3; inorm++)
    {
      a_normVal[inorm] = EBArith::norm(a_error, a_grids, a_ebisl, a_params.m_refRatio, icomp, inorm);
    }
}

/****/
void getDivergence(Vector< LevelData<EBCellFAB>* >&           a_error,
                   const Vector< DisjointBoxLayout >&         a_grids,
                   const Vector< EBISLayout >&                a_ebisl,
                   const ProblemDomain&                       a_level0Domain,
                   const Real&                                a_level0Dx,
                   const AMRParameters&                       a_params)
{
  int nlevels = a_grids.size();
  Vector<LevelData<EBCellFAB>* > veloc(nlevels, NULL);
  Vector<LevelData<EBFluxFAB>* > fluxVel(nlevels, NULL);
  a_error.resize(nlevels, NULL);
  Real dxLev = a_level0Dx;
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      EBFluxFactory ebfluxfact(a_ebisl[ilev]);
      a_error[ilev]  = new LevelData<EBCellFAB>(a_grids[ilev],       1 ,  IntVect::Unit,   ebcellfact);
      veloc[  ilev]  = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,4*IntVect::Unit,   ebcellfact);
      fluxVel[ilev]  = new LevelData<EBFluxFAB>(a_grids[ilev],       1 ,  IntVect::Zero,   ebfluxfact);
      setExactVeloc(*veloc[ilev], *fluxVel[ilev], a_grids[ilev], a_ebisl[ilev], dxLev, a_params);
      dxLev /= a_params.m_refRatio[ilev];
    }


  int numFilterIterations;
  ParmParse pp;


  //run it over with the filter
  pp.get("num_filter_iterations", numFilterIterations);
  Real normVal[3];
  getNormDiv(normVal, a_error, veloc, a_grids, a_ebisl, a_level0Domain, a_level0Dx, a_params);
  pout() << "before filtering kappa divergence(vel) l_inf, l_1, l2= " ;
  pout() << normVal[0] << ", " <<  normVal[1] << ", " << normVal[2] << endl;
  for(int ifilter = 0; ifilter < numFilterIterations; ifilter++)
    {
      dxLev = a_level0Dx;
      ProblemDomain domLev = a_level0Domain;
      for(int ilev = 0; ilev < nlevels; ilev++)
        {
          const DisjointBoxLayout*           coarGridsPtr = NULL;
          const EBISLayout*                  coarEBISLPtr = NULL;
          const LevelData<EBCellFAB>*        coarVelocPtr = NULL;
          if(ilev > 0)
            {
              coarGridsPtr = &a_grids[ilev-1];
              coarEBISLPtr = &a_ebisl[ilev-1];
              coarVelocPtr =  veloc[ilev-1];
            }
          int refRat = a_params.m_refRatio[ilev];
          EBGradDivFilter gdFilt(a_grids[ilev], coarGridsPtr, a_ebisl[ilev], coarEBISLPtr, domLev, dxLev*RealVect::Unit, refRat);

          gdFilt.filter(*veloc[ilev], *fluxVel[ilev], coarVelocPtr);

          dxLev /= refRat;
          domLev.refine(refRat);
        }
      getNormDiv(normVal, a_error, veloc, a_grids, a_ebisl, a_level0Domain, a_level0Dx, a_params);
      pout() << "after iteration "  << ifilter << ", kappa divergence(vel) l_inf, l_1, l2= " ;
      pout() << normVal[0] << ", " <<  normVal[1] << ", " << normVal[2] << endl;
    }

  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete veloc[ilev];
      delete fluxVel[ilev];
    }
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

    ProblemDomain domainCoar, domainFine, domainMedi;
    AMRParameters params;
    getAMRINSParameters(params, domainCoar);
    domainMedi = refine(domainCoar, 2);
    domainFine = refine(domainMedi, 2);
    Real dxFine = 1.0/(Real(domainFine.size(0)));
    Real dxMedi = 2.*dxFine;

    Vector<DisjointBoxLayout> gridsFine, gridsCoar, gridsMedi;
    Vector<EBISLayout>        ebislFine, ebislCoar, ebislMedi;

    int nlevels = params.m_maxLevel + 1;
    Vector<LevelData<EBCellFAB>* > divuFine( nlevels, NULL);
    Vector<LevelData<EBCellFAB>* > divuMedi( nlevels, NULL);


    pout() << "generating geometry only on finest level" << endl;
    AMRINSGeometry(params, domainFine);

    Vector<Vector<Box> > fineBoxes, mediBoxes, coarBoxes;
    getFixedLayouts(gridsFine,   gridsMedi,  gridsCoar,
                    ebislFine,   ebislMedi,  ebislCoar,
                    domainFine, domainMedi, domainCoar, params);


    pout() << "generating fine divu " << endl;
    getDivergence(divuFine, gridsFine, ebislFine, domainFine, dxFine, params);

    pout() << "generating medium divu" << endl;
    getDivergence(divuMedi, gridsMedi, ebislMedi, domainMedi, dxMedi, params);

    Vector<string> names(1, string("divu"));
    string testName("Filter error");
    compareError(divuFine,      divuMedi,
                 gridsFine,    gridsMedi,
                 ebislFine,    ebislMedi,
                 domainFine,  domainMedi,
                 names, testName, params);

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

        Vector<string> names(1, string("var0"));
        outputError(divuFine,     divuMedi,
                    gridsFine,    gridsMedi,
                    domainFine,   domainMedi,
                    fileFine,     fileCoar,
                    names,
                    params);
      }

    for(int ilev = 0; ilev < nlevels; ilev++)
      {
        delete divuFine[ilev];
        delete divuMedi[ilev];
        divuFine[ilev] = NULL;
        divuMedi[ilev] = NULL;
      }

  }// End scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif
}

