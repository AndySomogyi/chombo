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

#include "AMRINSUtils.H"
#include "EBAMRPoissonOp.H"

#include "EBFABView.H"
#include "EBDebugDump.H"
#include "DebugDump.H"
#include "EBAMRNoSubcycle.H"
#include "EBCoarseAverage.H"
#include "InflowOutflowIBC.H"
#include "UsingNamespace.H"

#define CHECK_FLOATING_PT 0

#if CHECK_FLOATING_PT==1
#include <fenv.h>
#endif
const Real pi = 4.*atan(1.0);

/****/
Real getScalar(const RealVect& a_xval)
{
  Real retval = 0.0;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      Real x = a_xval[idir];
      retval += x*x*x;
    }
  //debug
  //  retval = a_xval[0];
  //end debug

  return retval;
}

/****/
Real getDsDx(const RealVect& a_xval, int a_derivDir)
{
  Real retval = 0.0;
  retval = 3.*a_xval[a_derivDir]*a_xval[a_derivDir];
  //debug
  //if(a_derivDir==0)
  //  retval = 1.0;
  //else
  //  retval = 0.0;
  //end debug
  return retval;
}
/****/
Real getAdvVel(const RealVect& a_xval, int a_velDir)
{
  int variableDir = Max(1 -a_velDir, 0);
  Real retval = sin(pi*a_xval[variableDir]);

  //debug
  //retval = 1.0;
  //end debug
  return retval;
}
/****/
Real getUDelSExact(const RealVect& a_xval)
{
  Real retval = 0.0;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      Real dsdx = getDsDx(a_xval, idir);
      Real vel =  getAdvVel(a_xval, idir);
      retval += vel*dsdx;
    }
  return retval;
}
/****/
void
setExactStuff( EBCellFAB                  &  a_udelsExact,
               EBFluxFAB                  &  a_macAdvVel,
               EBFluxFAB                  &  a_macScalar,
               Vector< BaseIVFAB<Real> * >&  a_coveredAdvVelLo,
               Vector< BaseIVFAB<Real> * >&  a_coveredAdvVelHi,
               Vector< BaseIVFAB<Real> * >&  a_coveredScalarLo,
               Vector< BaseIVFAB<Real> * >&  a_coveredScalarHi,
               Vector< Vector<VolIndex>  >&  a_coveredFaceLo,
               Vector< Vector<VolIndex>  >&  a_coveredFaceHi,
               const Box                  &  a_grid,
               const EBISBox              &  a_ebisBox,
               const ProblemDomain        &  a_domain,
               const Real                 &  a_dx)
{
  a_udelsExact.setVal(0.);
  for(int faceDir = 0; faceDir < SpaceDim; ++faceDir)
    {
      a_macAdvVel[faceDir].setVal(0.);
      a_macScalar[faceDir].setVal(0.);
    }

  //first the non-covered stuff
  IntVectSet ivsGrid(a_grid);
  for(VoFIterator vofit(ivsGrid, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect xvalVoF;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          xvalVoF[idir] = a_dx*(Real(vof.gridIndex()[idir]) + 0.5);
        }
      a_udelsExact(vof, 0) = getUDelSExact(xvalVoF);

      for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          for(SideIterator sit; sit.ok(); ++sit)
            {
              Vector<FaceIndex> faces = a_ebisBox.getFaces(vof, faceDir, sit());
              for(int iface = 0; iface < faces.size(); iface++)
                {
                  const FaceIndex& face = faces[iface];
                  RealVect xvalFace;
                  for(int idir = 0; idir < SpaceDim; idir++)
                    {
                      if(idir != faceDir)
                        {
                          xvalFace[idir] = a_dx*(Real(face.gridIndex(Side::Hi)[idir]) + 0.5);
                        }
                      else
                        {
                          xvalFace[idir] = a_dx*(Real(face.gridIndex(Side::Hi)[idir]));
                        }
                    }
                  a_macAdvVel[faceDir](face, 0) = getAdvVel(xvalFace, faceDir);
                  a_macScalar[faceDir](face, 0) = getScalar(xvalFace);
                }
            }
        }
    }

  for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      //covered faces on the high side
      for(int ivof = 0; ivof < a_coveredFaceHi[faceDir].size(); ivof++)
        {
          const VolIndex& vof = a_coveredFaceHi[faceDir][ivof];
          RealVect xval;
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              if(idir != faceDir)
                {
                  xval[idir] = a_dx*(Real(vof.gridIndex()[idir]) + 0.5);
                }
              else
                {
                  //1 is because it is the high side
                  xval[idir] = a_dx*(Real(vof.gridIndex()[idir]) + 1.0);
                }
            }
          (*a_coveredAdvVelHi[faceDir])(vof, 0) = getAdvVel(xval, faceDir);
          (*a_coveredScalarHi[faceDir])(vof, 0) = getScalar(xval);
        }


      //covered faces on the lowside
      for(int ivof = 0; ivof < a_coveredFaceLo[faceDir].size(); ivof++)
        {
          const VolIndex& vof = a_coveredFaceLo[faceDir][ivof];
          RealVect xval;
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              if(idir != faceDir)
                {
                  xval[idir] = a_dx*(Real(vof.gridIndex()[idir]) + 0.5);
                }
              else
                {
                  //no 1 is because it is the low side
                  xval[idir] = a_dx*(Real(vof.gridIndex()[idir]));
                }
            }
          (*a_coveredAdvVelLo[faceDir])(vof, 0) = getAdvVel(xval, faceDir);
          (*a_coveredScalarLo[faceDir])(vof, 0) = getScalar(xval);
        }
    }
}
/****/
void
setExactStuff(Vector< LevelData<EBCellFAB>* >                     &  a_udelsExact,
              Vector< LevelData<EBFluxFAB>* >                     &  a_macAdvVel,
              Vector< LevelData<EBFluxFAB>* >                     &  a_macScalar,
              Vector< LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelLo,
              Vector< LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelHi,
              Vector< LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredScalarLo,
              Vector< LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredScalarHi,
              Vector< LayoutData< Vector< Vector<VolIndex> > >* > &  a_coveredFaceLo,
              Vector< LayoutData< Vector< Vector<VolIndex> > >* > &  a_coveredFaceHi,
              const Vector<DisjointBoxLayout>                     &  a_grids,
              const Vector<EBISLayout>                            &  a_ebisl,
              const Vector<ProblemDomain>                         &  a_domain,
              const Vector<Real>                                  &  a_dx,
              int a_nlevels)
{
  for(int ilev = 0; ilev < a_nlevels; ilev++)
    {
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          setExactStuff((*a_udelsExact[ilev])[dit()],
                        (*a_macAdvVel[ilev])[dit()],
                        (*a_macScalar[ilev])[dit()],
                        (*a_coveredAdvVelLo[ilev])[dit()],
                        (*a_coveredAdvVelHi[ilev])[dit()],
                        (*a_coveredScalarLo[ilev])[dit()],
                        (*a_coveredScalarHi[ilev])[dit()],
                        (*a_coveredFaceLo[ilev])[dit()],
                        (*a_coveredFaceHi[ilev])[dit()],
                        a_grids[ilev].get(dit()),
                        a_ebisl[ilev][dit()],
                        a_domain[ilev], a_dx[ilev]);
        }
    }
}

/****/
void getError(Vector< LevelData<EBCellFAB>* >&           a_error,
              Vector< DisjointBoxLayout >&               a_grids,
              Vector< EBISLayout >&                      a_ebisl,
              const ProblemDomain&                       a_level0Domain,
              const AMRParameters&                       a_params,
              const Vector<Vector<Box> > &               a_vectBoxes)
{
  int flowDir;
  ParmParse pp;
  pp.get("flow_dir", flowDir);
  Vector<Real> centerVect(SpaceDim);
  pp.getarr("spot_center",centerVect,0,SpaceDim);
  RealVect center;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      center[idir] = centerVect[idir];
    }

  Real coreRadius, inflowVel;
  pp.get("spot_radius", coreRadius);
  pp.get("inflow_vel", inflowVel);

  Real viscosity = 0.0;
  pout() << "inviscid test" << endl;
  int orderEBBC = 2; //not used --inviscid test
  bool doslip = true;
  InflowOutflowIBCFactory ibc(flowDir, inflowVel, orderEBBC, doslip);;

  EBAMRNoSubcycle amr(a_params, ibc, a_level0Domain, viscosity);

  amr.setupForFixedHierarchyRun(a_vectBoxes);

  a_grids = amr.getGrids();
  a_ebisl = amr.getEBISLayouts();

  int nlevels = a_params.m_maxLevel + 1;
  a_error.resize(nlevels);
  Vector< LevelData<EBCellFAB>* > udelsCalcu(nlevels, NULL);
  Vector< LevelData<EBCellFAB>* > udelsExact(nlevels, NULL);

  Vector< LevelData<EBFluxFAB>* > macAdvVel(nlevels, NULL);
  Vector< LevelData<EBFluxFAB>* > macScalar(nlevels, NULL);

  Vector< LayoutData< Vector< BaseIVFAB<Real> * > >* >  coveredAdvVelLo(nlevels, NULL);
  Vector< LayoutData< Vector< BaseIVFAB<Real> * > >* >  coveredAdvVelHi(nlevels, NULL);
  Vector< LayoutData< Vector< BaseIVFAB<Real> * > >* >  coveredScalarLo(nlevels, NULL);
  Vector< LayoutData< Vector< BaseIVFAB<Real> * > >* >  coveredScalarHi(nlevels, NULL);

  Vector< LayoutData< Vector< IntVectSet > >* >      & coveredSetsLo = amr.getCoveredSetsLo();
  Vector< LayoutData< Vector< IntVectSet > >* >      & coveredSetsHi = amr.getCoveredSetsHi();
  Vector< LayoutData< Vector< Vector<VolIndex> > >* >& coveredFaceLo = amr.getCoveredFaceLo();
  Vector< LayoutData< Vector< Vector<VolIndex> > >* >& coveredFaceHi = amr.getCoveredFaceHi();

  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      a_error[ilev]    = new LevelData<EBCellFAB>(a_grids[ilev], 1,  IntVect::Zero, ebcellfact);
      udelsCalcu[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], 1,  IntVect::Zero, ebcellfact);
      udelsExact[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], 1,  IntVect::Zero, ebcellfact);

      EBFluxFactory ebfluxfact(a_ebisl[ilev]);
      macAdvVel[ilev] = new LevelData<EBFluxFAB>(a_grids[ilev], 1,  IntVect::Zero, ebfluxfact);
      macScalar[ilev] = new LevelData<EBFluxFAB>(a_grids[ilev], 1,  IntVect::Zero, ebfluxfact);

      coveredAdvVelLo[ilev] = new LayoutData<Vector<BaseIVFAB<Real> * > >(a_grids[ilev]);
      coveredAdvVelHi[ilev] = new LayoutData<Vector<BaseIVFAB<Real> * > >(a_grids[ilev]);
      coveredScalarLo[ilev] = new LayoutData<Vector<BaseIVFAB<Real> * > >(a_grids[ilev]);
      coveredScalarHi[ilev] = new LayoutData<Vector<BaseIVFAB<Real> * > >(a_grids[ilev]);
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*coveredAdvVelLo[ilev])[dit()].resize(SpaceDim, NULL);
          (*coveredAdvVelHi[ilev])[dit()].resize(SpaceDim, NULL);
          (*coveredScalarLo[ilev])[dit()].resize(SpaceDim, NULL);
          (*coveredScalarHi[ilev])[dit()].resize(SpaceDim, NULL);

          const EBGraph& ebgraph = a_ebisl[ilev][dit()].getEBGraph();
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              (*coveredAdvVelLo[ilev])[dit()][idir] = new BaseIVFAB<Real>((*coveredSetsLo[ilev])[dit()][idir], ebgraph, 1);
              (*coveredAdvVelHi[ilev])[dit()][idir] = new BaseIVFAB<Real>((*coveredSetsHi[ilev])[dit()][idir], ebgraph, 1);
              (*coveredScalarLo[ilev])[dit()][idir] = new BaseIVFAB<Real>((*coveredSetsLo[ilev])[dit()][idir], ebgraph, 1);
              (*coveredScalarHi[ilev])[dit()][idir] = new BaseIVFAB<Real>((*coveredSetsHi[ilev])[dit()][idir], ebgraph, 1);
            }
        }

    }

  Vector<Real>& dx = amr.getDx();
  Vector<ProblemDomain>& domain = amr.getDomain();

  //set advective velocity and scalar  and exact udels
  setExactStuff(udelsExact,
                macAdvVel,
                macScalar,
                coveredAdvVelLo,
                coveredAdvVelHi,
                coveredScalarLo,
                coveredScalarHi,
                coveredFaceLo,
                coveredFaceHi,
                a_grids,
                a_ebisl,
                domain, dx, nlevels);

  //use amr to get approximate udels from exact velocities and scalars
  bool useNonConsOnly = true;
  bool useConsOnly = false;
  amr.computeAdvectiveDerivative(udelsCalcu,
                                 macAdvVel,
                                 macScalar,
                                 coveredAdvVelLo,
                                 coveredAdvVelHi,
                                 coveredScalarLo,
                                 coveredScalarHi,
                                 useNonConsOnly,
                                 useConsOnly);


  //error = udels calc - udelu exact
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBLevelDataOps::setVal(*a_error[ilev], 0.0);
      EBLevelDataOps::incr(*a_error[ilev], *udelsCalcu[ilev], 1.0);
      EBLevelDataOps::incr(*a_error[ilev], *udelsExact[ilev],-1.0);
    }


  //clean up temporary data holders. a_error  has to survive  past this
  //routine even though it was allocated here.
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete udelsCalcu[ilev];
      delete udelsExact[ilev];

      delete macAdvVel[ilev];
      delete macScalar[ilev];

      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              delete (*coveredAdvVelLo[ilev])[dit()][idir];
              delete (*coveredAdvVelHi[ilev])[dit()][idir];
              delete (*coveredScalarLo[ilev])[dit()][idir];
              delete (*coveredScalarHi[ilev])[dit()][idir];
            }
        }

      delete coveredAdvVelLo[ilev];
      delete coveredAdvVelHi[ilev];
      delete coveredScalarLo[ilev];
      delete coveredScalarHi[ilev];

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
    ParmParse pp(argc-2, argv+2, NULL, inFile);

    ProblemDomain domainCoar, domainFine;
    AMRParameters params;
    getAMRINSParameters(params, domainCoar);
    domainFine = refine(domainCoar, 2);

    Vector<DisjointBoxLayout> gridsFine, gridsCoar;
    Vector<EBISLayout>        ebislFine, ebislCoar;

    int nlevels = params.m_maxLevel + 1;
    Vector<LevelData<EBCellFAB>* > errorFine(nlevels, NULL);
    Vector<LevelData<EBCellFAB>* > errorCoar(nlevels, NULL);

    pout() << "generating geometry only on finest level" << endl;
    AMRINSGeometry(params, domainFine);

    Vector<Vector<Box> > fineBoxes, mediBoxes, coarBoxes;
    getFixedGrids(fineBoxes, mediBoxes, coarBoxes, params, domainCoar);


    pout() << "generating fine error " << endl;
    getError(errorFine,
             gridsFine, ebislFine, domainFine, params, mediBoxes);

    pout() << "generating coar error" << endl;

    getError(errorCoar,
             gridsCoar, ebislCoar, domainCoar, params, coarBoxes);

    string testName("advective derivative error");
    Vector<string> names(1, "udotdels_error");
    compareError(errorFine,   errorCoar,
                 gridsFine,   gridsCoar,
                 ebislFine,   ebislCoar,
                 domainFine,  domainCoar,
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
        outputError(errorFine,   errorCoar,
                    gridsFine,   gridsCoar,
                    domainFine,  domainCoar,
                    fileFine,     fileCoar,
                    names, params);
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

