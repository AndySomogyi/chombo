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
  //  retval = 1.0;
  return retval;
}
/****/
Real getAdvVel(const RealVect& a_xval, int a_velDir)
{
  int variableDir = Max(1 -a_velDir, 0);
  Real retval = sin(pi*a_xval[variableDir]);

//  if(a_velDir == 0)
//    retval = 1.0;
//  else
//    retval = 0.0;

  return retval;
}
/****/
void
setExactStuff( EBFluxFAB                  &  a_macAdvVel,
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
  //first the non-covered stuff
  IntVectSet ivsGrid(a_grid);
  for(VoFIterator vofit(ivsGrid, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
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
setExactStuff(Vector< LevelData<EBFluxFAB>* >                     &  a_macAdvVel,
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
          setExactStuff((*a_macAdvVel[ilev])[dit()],
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
/******/
/******/
void getError(Vector< LevelData<EBCellFAB>* >&           a_errorCoar,
              const Vector< LevelData<EBCellFAB>* >&     a_solnCoar,
              const Vector< DisjointBoxLayout >&         a_gridsCoar,
              const Vector< EBISLayout >&                a_ebislCoar,
              const ProblemDomain&                       a_level0DomainCoar,
              const Vector< LevelData<EBCellFAB>* >&     a_solnFine,
              const Vector< DisjointBoxLayout >&         a_gridsFine,
              const Vector< EBISLayout >&                a_ebislFine,
              const ProblemDomain&                       a_level0DomainFine,
              const AMRParameters&                       a_params)
{

  int nlevels = a_params.m_maxLevel + 1;
  a_errorCoar.resize(nlevels);
  int nref = 2;  //nothing to do with param refinement ratio. this is the refinement between the two solutions
  int nvar = 1;
  Interval interv(0, nvar-1);

  ProblemDomain domLevCoar = a_level0DomainCoar;
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBCellFactory ebcellfact(a_ebislCoar[ilev]);
      int nvar = a_solnCoar[ilev]->nComp();

      a_errorCoar[ilev] = new LevelData<EBCellFAB>(a_gridsCoar[ilev], nvar,  IntVect::Zero, ebcellfact);

      EBCoarseAverage averageOp(a_gridsFine[ilev], a_gridsCoar[ilev], a_ebislFine[ilev], a_ebislCoar[ilev], domLevCoar, nref, nvar, Chombo_EBIS::instance());
      //here make error = Ave(fine)
      averageOp.average(*a_errorCoar[ilev], *a_solnFine[ilev], interv);
      //now subtract off coarse so error= Ave(Fine) - coar
      for(DataIterator dit = a_gridsCoar[ilev].dataIterator(); dit.ok(); ++dit)
        {
          EBCellFAB&       errorFAB = (*a_errorCoar[ilev])[dit()];
          const EBCellFAB& solnFAB  = (*a_solnCoar[ilev])[dit()];

          errorFAB -= solnFAB;
        }
      domLevCoar.refine(a_params.m_refRatio[ilev]);
    }
}
/****/
/****/
void getSolution(Vector< LevelData<EBCellFAB>* >&           a_udotdels,
                 Vector< DisjointBoxLayout >&               a_grids,
                 Vector< EBISLayout >&                      a_ebisl,
                 const ProblemDomain&                       a_level0Domain,
                 const AMRParameters&                       a_params)
{
  int flowDir;
  Real inflowVel = 1.0;
  ParmParse pp;
  pp.get("flow_dir", flowDir);
  pp.get("inflow_vel", inflowVel);


  Real viscosity = 0.0;

  pout() << "inviscid test" << endl;
  bool doslip = true;//noop
  bool orderebbc = 2; //noop
  InflowOutflowIBCFactory ibc(flowDir, inflowVel, orderebbc, doslip);

  EBAMRNoSubcycle amr(a_params, ibc, a_level0Domain, viscosity);

  amr.setupForAMRRun();

  a_grids = amr.getGrids();
  a_ebisl = amr.getEBISLayouts();

  int nlevels = a_params.m_maxLevel + 1;
  a_udotdels.resize(nlevels);

  Vector< LevelData<EBFluxFAB>* > macAdvVel(nlevels, NULL);
  Vector< LevelData<EBFluxFAB>* > macGradient(nlevels, NULL);
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
      a_udotdels[ilev]    = new LevelData<EBCellFAB>(a_grids[ilev], 1,  IntVect::Zero, ebcellfact);

      EBFluxFactory ebfluxfact(a_ebisl[ilev]);
      macAdvVel[ilev] = new LevelData<EBFluxFAB>(a_grids[ilev], 1,  IntVect::Zero, ebfluxfact);
      macGradient[ilev] = new LevelData<EBFluxFAB>(a_grids[ilev], 1,  IntVect::Zero, ebfluxfact);
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
  setExactStuff(macAdvVel,
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

  EBCompositeMACProjector& macproj = amr.getMACProj();

  pout() << "projecting advection velocity" << endl;

  macproj.project(macAdvVel, macGradient);

  //correct covered advective velocities
  //faceDir,  velcomp are the same thing for advection velocities
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      macproj.correctVelocityComponent(coveredAdvVelLo,
                                       coveredAdvVelHi,
                                       coveredFaceLo,
                                       coveredFaceHi,
                                       coveredSetsLo,
                                       coveredSetsHi,
                                       macGradient, idir, idir);
    }

  //use amr to get approximate udels from exact velocities and scalars
  int iconv;
  pp.get("use_only_nonconservative_diff", iconv);
  bool useNonConvOnly = (iconv == 1);
  bool useConsOnly = false;
  amr.computeAdvectiveDerivative(a_udotdels,
                                 macAdvVel,
                                 macScalar,
                                 coveredAdvVelLo,
                                 coveredAdvVelHi,
                                 coveredScalarLo,
                                 coveredScalarHi,
                                 useNonConvOnly,
                                 useConsOnly);

  //clean up temporary data holders. a_udotdels  has to survive  past this
  //routine even though it was allocated here.
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
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
/***************/
/******************************/
/******************************/
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

    Vector<DisjointBoxLayout> gridsFine, gridsCoar, gridsMedi;
    Vector<EBISLayout>        ebislFine, ebislCoar, ebislMedi;

    int nlevels = params.m_maxLevel + 1;
    Vector<LevelData<EBCellFAB>* > solnFine( nlevels, NULL);
    Vector<LevelData<EBCellFAB>* > solnMedi( nlevels, NULL);
    Vector<LevelData<EBCellFAB>* > solnCoar( nlevels, NULL);
    Vector<LevelData<EBCellFAB>* > errorMedi(nlevels, NULL);
    Vector<LevelData<EBCellFAB>* > errorCoar(nlevels, NULL);


    pout() << "generating geometry only on finest level" << endl;
    AMRINSGeometry(params, domainFine);


    pout() << "generating fine solution " << endl;
    getSolution(solnFine, gridsFine, ebislFine, domainFine, params);

    pout() << "generating medi solution and error" << endl;
    //AMRINSGeometry(params, domainMedi);
    getSolution(solnMedi, gridsMedi, ebislMedi, domainMedi, params);

    getError(errorMedi,
             solnMedi, gridsMedi, ebislMedi, domainMedi,
             solnFine, gridsFine, ebislFine, domainFine,
             params);

    pout() << "generating coar solution and error" << endl;
    //AMRINSGeometry(params, domainCoar);
    getSolution(solnCoar, gridsCoar, ebislCoar, domainCoar, params);

    getError(errorCoar,
             solnCoar, gridsCoar, ebislCoar, domainCoar,
             solnMedi, gridsMedi, ebislMedi, domainMedi,
             params);

    string testName("Solution error");
    Vector<string>   names(1, string("udeluError"));
    compareError(errorMedi,   errorCoar,
                 gridsMedi,   gridsCoar,
                 ebislMedi,   ebislCoar,
                 domainMedi,  domainCoar,
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
        outputError(errorMedi,   errorCoar,
                    gridsMedi,   gridsCoar,
                    domainMedi,  domainCoar,
                    fileFine,     fileCoar,
                    names, params);
      }

    for(int ilev = 0; ilev < nlevels; ilev++)
      {
        delete solnFine[ilev];
        delete solnMedi[ilev];
        delete solnCoar[ilev];
        solnFine[ilev] = NULL;
        solnMedi[ilev] = NULL;
        solnCoar[ilev] = NULL;

        delete errorMedi[ilev];
        delete errorCoar[ilev];
        errorMedi[ilev] = NULL;
        errorCoar[ilev] = NULL;
      }

  }// End scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif
}
/**********************/
/**********************/
