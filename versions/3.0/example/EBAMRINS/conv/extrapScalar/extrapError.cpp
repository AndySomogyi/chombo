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
#include "EBAMRNoSubcycle.H"
#include "EBCoarseAverage.H"
#include "InflowOutflowIBC.H"
#include "EBPatchGodunov.H"
#include "UsingNamespace.H"

using std::cerr;

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

  ///retval = 1.0;
  return retval;
}
/****/
Real getAdvVel(const RealVect& a_xval, int a_velDir)
{
  int variableDir = Max(1 -a_velDir, 0);
  Real retval = sin(pi*a_xval[variableDir]);

  if(a_velDir == 0)
    retval = 1.0;
  else
    retval = 0.0;

  return retval;
}
/****/
void
setExactStuff( EBFluxFAB           &  a_macAdvVel,
               EBCellFAB           &  a_cellScala,
               EBCellFAB           &  a_cellVeloc,
               const Box           &  a_grid,
               const EBISBox       &  a_ebisBox,
               const ProblemDomain &  a_domain,
               const Real          &  a_dx)
{
  IntVectSet ivsGrid(a_grid);
  for(VoFIterator vofit(ivsGrid, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect xvalVoF;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          xvalVoF[idir] = a_dx*(Real(vof.gridIndex()[idir]) + 0.5);
        }

      a_cellScala(vof, 0) = getScalar(xvalVoF);
      for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          a_cellVeloc(vof, faceDir) = getAdvVel(xvalVoF, faceDir);
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
                }
            }
        }
    }
}

/****/
void
setExactStuff(Vector< LevelData<EBFluxFAB>* > &  a_macAdvVel,
              Vector< LevelData<EBCellFAB>* > &  a_cellScala,
              Vector< LevelData<EBCellFAB>* > &  a_cellVeloc,
              const Vector<DisjointBoxLayout> &  a_grids,
              const Vector<EBISLayout>        &  a_ebisl,
              const Vector<ProblemDomain>     &  a_domain,
              const Vector<Real>              &  a_dx,
              int a_nlevels)
{
  for(int ilev = 0; ilev < a_nlevels; ilev++)
    {
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          setExactStuff((*a_macAdvVel[ilev])[dit()],
                        (*a_cellScala[ilev])[dit()],
                        (*a_cellVeloc[ilev])[dit()],
                        a_grids[ilev].get(dit()),
                        a_ebisl[ilev][dit()],
                        a_domain[ilev], a_dx[ilev]);
        }
    }
}
/******/
Real
getMaxVel(const Vector< LevelData<EBCellFAB>* > &  a_cellVeloc,
          const Vector<DisjointBoxLayout> &  a_grids,
          const Vector<EBISLayout>        &  a_ebisl,
          int a_nlevels)

{
  //find the max on this proc
  Real maxVelLoc = 0.0;
  for(int ilev = 0; ilev < a_nlevels; ilev++)
    {
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          const EBCellFAB& velFAB = (*a_cellVeloc[ilev])[dit()];
          const Box& box =         a_grids[ilev].get(dit());
          const EBISBox& ebisBox =  a_ebisl[ilev][dit()];
          IntVectSet ivs(box);
          for(VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              for(int idir = 0; idir < SpaceDim; idir++)
                {
                  maxVelLoc = Max(maxVelLoc, Abs(velFAB(vofit(), idir)));
                }
            }
        }
    }

  //take the max over all procs
  int baseProc = 0;
  Vector<Real> maxVec;
  gather(maxVec,   maxVelLoc, baseProc);
  Real maxVel  = 0.0;
  if(procID() == baseProc)
    {
      CH_assert(maxVec.size() == numProc());
      for(int ivec = 0; ivec < numProc(); ivec++)
        {
          maxVel = max(Abs(maxVec[ivec]), maxVel);
        }
    }
  broadcast(maxVel,  baseProc);

  return maxVel;
}
/******/
void getError(Vector< LevelData<EBFluxFAB>* >&           a_errorCoar,
              const Vector< LevelData<EBFluxFAB>* >&     a_scalaCoar,
              const Vector< DisjointBoxLayout >&         a_gridsCoar,
              const Vector< EBISLayout >&                a_ebislCoar,
              const ProblemDomain&                       a_level0DomainCoar,
              const Vector< LevelData<EBFluxFAB>* >&     a_scalaFine,
              const Vector< DisjointBoxLayout >&         a_gridsFine,
              const Vector< EBISLayout >&                a_ebislFine,
              const ProblemDomain&                       a_level0DomainFine,
              const AMRParameters&                       a_params)
{

  int nlevels = a_params.m_maxLevel + 1;
  a_errorCoar.resize(nlevels);
  int nref = 2;  //nothing to do with param refinement ratio. this is the refinement between the two scalas
  int nvar = 1;
  Interval interv(0, nvar-1);

  ProblemDomain domLevCoar = a_level0DomainCoar;
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBFluxFactory ebfluxfact(a_ebislCoar[ilev]);
      int nvar = a_scalaCoar[ilev]->nComp();

      a_errorCoar[ilev] = new LevelData<EBFluxFAB>(a_gridsCoar[ilev], nvar,  IntVect::Zero, ebfluxfact);

      EBLevelDataOps::setVal(*a_errorCoar[ilev], 0.0);

      EBCoarseAverage averageOp(a_gridsFine[ilev], a_gridsCoar[ilev],
                                a_ebislFine[ilev], a_ebislCoar[ilev],
                                domLevCoar, nref, nvar, Chombo_EBIS::instance());

      //here make error = Ave(fine)
      averageOp.average(*a_errorCoar[ilev], *a_scalaFine[ilev], interv);
      //now subtract off coarse so error= Ave(Fine) - coar
      int ibox = 0;
      for(DataIterator dit = a_gridsCoar[ilev].dataIterator(); dit.ok(); ++dit)
        {
          for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
            {
              EBFaceFAB&       errorFAB = (*a_errorCoar[ilev])[dit()][faceDir];
              const EBFaceFAB& scalaFAB = (*a_scalaCoar[ilev])[dit()][faceDir];
              errorFAB -= scalaFAB;
              ibox++;
            }
        }
      domLevCoar.refine(a_params.m_refRatio[ilev]);
    }
}
/****/
/****/
void getExtrapScalar(Vector< LevelData<EBFluxFAB>* >&           a_extScal,
                     Real&                                      a_dt,
                     Vector< DisjointBoxLayout >&               a_grids,
                     Vector< EBISLayout >&                      a_ebisl,
                     const ProblemDomain&                       a_level0Domain,
                     const AMRParameters&                       a_params,
                     bool                                       a_atFinestLevel,
                     const Vector<Vector<Box> > &               a_vectBoxes)
{
  int flowDir;
  Real inflowVel = 1.0;
  ParmParse pp;
  pp.get("flow_dir", flowDir);
  pp.get("inflow_vel", inflowVel);

  Real viscosity = 0.0;
  pout() << "inviscid calc" << endl;
  bool doslip = true; //no op
  int orderebbc = 2;  //no op
  InflowOutflowIBCFactory ibcfact(flowDir, inflowVel, orderebbc, doslip);

  EBAMRNoSubcycle amr(a_params, ibcfact, a_level0Domain, viscosity);

  amr.setupForFixedHierarchyRun(a_vectBoxes);

  a_grids = amr.getGrids();
  a_ebisl = amr.getEBISLayouts();

  int nlevels = a_params.m_maxLevel + 1;
  a_extScal.resize(nlevels);

  Vector< LevelData<EBCellFAB>* > cellVeloc(nlevels, NULL);
  Vector< LevelData<EBCellFAB>* > cellGradi(nlevels, NULL);
  Vector< LevelData<EBCellFAB>* > cellScala(nlevels, NULL);
  Vector< LevelData<EBFluxFAB>* > macAdvVel(nlevels, NULL);
  Vector< LevelData<EBFluxFAB>* > macGradie(nlevels, NULL);

  Vector< LayoutData< Vector< BaseIVFAB<Real> * > >* >  coveredScalarLo(nlevels, NULL);
  Vector< LayoutData< Vector< BaseIVFAB<Real> * > >* >  coveredScalarHi(nlevels, NULL);

  Vector< LayoutData< Vector< IntVectSet > >* >      & coveredSetsLo = amr.getCoveredSetsLo();
  Vector< LayoutData< Vector< IntVectSet > >* >      & coveredSetsHi = amr.getCoveredSetsHi();

  IntVect ivghost = 4*IntVect::Unit;
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBFluxFactory ebfluxfact(a_ebisl[ilev]);
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      a_extScal[ilev] = new LevelData<EBFluxFAB>(a_grids[ilev],        1,  ivghost, ebfluxfact);
      macAdvVel[ilev] = new LevelData<EBFluxFAB>(a_grids[ilev],        1,  ivghost, ebfluxfact);
      macGradie[ilev] = new LevelData<EBFluxFAB>(a_grids[ilev],        1,  ivghost, ebfluxfact);
      cellScala[ilev] = new LevelData<EBCellFAB>(a_grids[ilev],        1,  ivghost, ebcellfact);
      cellVeloc[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,  ivghost, ebcellfact);
      cellGradi[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,  ivghost, ebcellfact);

      EBLevelDataOps::setVal(*cellScala[ilev], 0.0);
      EBLevelDataOps::setVal(*cellVeloc[ilev], 0.0);
      EBLevelDataOps::setVal(*cellGradi[ilev], 0.0);
      EBLevelDataOps::setVal(*macAdvVel[ilev], 0.0);
      EBLevelDataOps::setVal(*macGradie[ilev], 0.0);
      EBLevelDataOps::setVal(*a_extScal[ilev], 0.0);

      coveredScalarLo[ilev] = new LayoutData<Vector<BaseIVFAB<Real> * > >(a_grids[ilev]);
      coveredScalarHi[ilev] = new LayoutData<Vector<BaseIVFAB<Real> * > >(a_grids[ilev]);
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*coveredScalarLo[ilev])[dit()].resize(SpaceDim, NULL);
          (*coveredScalarHi[ilev])[dit()].resize(SpaceDim, NULL);

          const EBGraph& ebgraph = a_ebisl[ilev][dit()].getEBGraph();
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              (*coveredScalarLo[ilev])[dit()][idir] = new BaseIVFAB<Real>((*coveredSetsLo[ilev])[dit()][idir], ebgraph, 1);
              (*coveredScalarHi[ilev])[dit()][idir] = new BaseIVFAB<Real>((*coveredSetsHi[ilev])[dit()][idir], ebgraph, 1);
            }
        }

    }

  Vector<Real>& dx = amr.getDx();
  Vector<ProblemDomain>& domain = amr.getDomain();
  //set advective velocity and scalar  and exact udels
  setExactStuff(macAdvVel,
                cellScala,
                cellVeloc,
                a_grids,
                a_ebisl,
                domain, dx, nlevels);

  EBCompositeMACProjector& macproj = amr.getMACProj();
  EBCompositeCCProjector&   ccproj = amr.getCCProj();

  pout() << "projecting advection velocity" << endl;

  macproj.project(macAdvVel, macGradie);

  pout() << "projecting cell-centered velocity" << endl;

  ccproj.project(cellVeloc, cellGradi);

  pout() << "advecting scalar to faces" << endl;

  ///klugy mechanism to make sure everyone is using the same time step
  if(a_atFinestLevel)
    {
      Real maxVel = getMaxVel(cellVeloc, a_grids, a_ebisl, nlevels);
      const Vector<Real>& dx = amr.getDx();
      CH_assert(maxVel > 0.0);
      a_dt = 0.5*dx[nlevels-1]/maxVel;
    }
  amr.useFixedDt(a_dt);

  pout() << "inviscid test"  << endl;
  EBIBC* ibcPtr = ibcfact.create();
  RefCountedPtr<EBPhysIBCFactory> advectBC  = ibcPtr->getScalarAdvectBC(0);
  Vector<LevelData<EBCellFAB>* > * source = NULL;
  EBPatchGodunov::setCurComp(0);
  EBPatchGodunov::setDoingVel(0);
  amr.extrapolateScalarCol(a_extScal,
                           coveredScalarLo,
                           coveredScalarHi,
                           advectBC,
                           macAdvVel,
                           source,
                           cellScala,
                           cellVeloc);


  pout() << "cleaning up memory" << endl;
  delete ibcPtr;
  //clean up temporary data holders. a_extScal  has to survive  past this
  //routine even though it was allocated here.
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete macAdvVel[ilev];
      delete macGradie[ilev];
      delete cellScala[ilev];
      delete cellVeloc[ilev];

      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              delete (*coveredScalarLo[ilev])[dit()][idir];
              delete (*coveredScalarHi[ilev])[dit()][idir];
            }
        }

      delete coveredScalarLo[ilev];
      delete coveredScalarHi[ilev];

    }
}
/***************/
void
norm(Real& a_sumTot,
     Real& a_areaTot,
     const EBFaceFAB& a_error,
     const EBISBox& a_ebisBox,
     const IntVectSet& a_set,
     const int& a_idir,
     const int& a_comp,
     const int& a_normtype)
{
  CH_assert(a_normtype >= 0);
  CH_assert(a_normtype <= 2);
  FaceIterator faceit(a_set, a_ebisBox.getEBGraph(), a_idir,
                      FaceStop::SurroundingWithBoundary);
  a_areaTot = 0.0;
  a_sumTot = 0.0;
  for(faceit.reset(); faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      Real areaFrac = a_ebisBox.areaFrac(face);
      Real valFace = a_error(face, a_comp);
      a_areaTot += areaFrac;
      if(a_normtype == 0)
        a_sumTot = Max(Abs(valFace), a_sumTot);
      else if(a_normtype == 1)
        a_sumTot += Abs(valFace)*areaFrac;
      else
        a_sumTot += valFace*valFace*areaFrac;
    }
}
/***************/
Real
norm(const Vector<LevelData<EBFluxFAB> * >&    a_error,
     const Vector< DisjointBoxLayout >&        a_grids,
     const Vector< EBISLayout >&               a_ebisl,
     const AMRParameters&                      a_params,
     const int&                                a_idir,
     const int&                                a_comp,
     const int&                                a_normtype)
{

  int nlevels = a_params.m_maxLevel + 1;
  Real areaLoc = 0.0;
  Real sumLoc  = 0.0;

  //get the sums of area and integrand on this proc
  for(int ilev = 0; ilev < nlevels; ilev++)
    {

      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          const Box& box         = a_grids[ilev].get(dit());
          const EBISBox& ebisBox = a_ebisl[ilev][dit()];
          IntVectSet ivs(box);

          const EBFaceFAB& faceError = (*a_error[ilev])[dit()][a_idir];
          Real areaBox, sumBox;
          norm(sumBox, areaBox, faceError, ebisBox, ivs, a_idir, a_comp, a_normtype);
          areaLoc += areaBox;
          sumLoc  +=  sumBox;
        }
    }

  //do the broadcast-gather thing so we all know what the sums are,
  // now for the multi-processor fandango
  //gather all the rhodots onto a vector and add them up
  int baseProc = 0;
  Vector<Real> sumVec;
  Vector<Real> areaVec;
  gather(sumVec,   sumLoc, baseProc);
  gather(areaVec, areaLoc, baseProc);

  Real sumTot  = 0.0;
  Real areaTot = 0.0;
  if(procID() == baseProc)
    {
      CH_assert(sumVec.size() == numProc());
      CH_assert(areaVec.size() == numProc());
      for(int ivec = 0; ivec < numProc(); ivec++)
        {
          areaTot += areaVec[ivec];
          if(a_normtype == 0)
            {
              sumTot  = Max(Abs(sumVec[ivec]), sumTot);
            }
          else
            {
              sumTot += sumVec[ivec];
            }
        }
    }
  //broadcast the sum to all processors.
  broadcast(sumTot,  baseProc);
  broadcast(areaTot, baseProc);


  //figure out the norm based on sums of area and integrand
  Real normTot = sumTot;
  if((a_normtype > 0 ) && (areaTot > 0.0))
    normTot /= areaTot;
  if(a_normtype == 2)
    normTot = sqrt(normTot);

  return normTot;
}
/***************/
void compareError(const Vector< LevelData<EBFluxFAB>* >&   a_errorFine,
                  const Vector< LevelData<EBFluxFAB>* >&   a_errorCoar,
                  const Vector< DisjointBoxLayout >&       a_gridsFine,
                  const Vector< DisjointBoxLayout >&       a_gridsCoar,
                  const Vector< EBISLayout >&              a_ebislFine,
                  const Vector< EBISLayout >&              a_ebislCoar,
                  const ProblemDomain&                     a_level0DomainFine,
                  const ProblemDomain&                     a_level0DomainCoar,
                  const Vector<string>&                    a_names,
                  const string&                            a_testName,
                  const AMRParameters&                     a_params)
{
  int ncomp = a_errorFine[0]->nComp();
  pout() << "==============================================" << endl;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      pout() << "Direction = " << idir << endl;
      pout() << "==============================================" << endl;
      for(int comp = 0; comp < ncomp; comp++)
        {
          pout() << "component = " << comp << endl;
          pout() << "==============================================" << endl;
          for(int inorm = 0; inorm < 3; inorm++)
            {

              if(inorm == 0)
                {
                  pout() << endl << "Using max norm." << endl;
                }
              else
                {
                  pout() << endl << "Using L-" << inorm << "norm." << endl;
                }
              Real coarnorm = norm(a_errorCoar, a_gridsCoar, a_ebislCoar, a_params,  idir, comp, inorm);
              Real finenorm = norm(a_errorFine, a_gridsFine, a_ebislFine, a_params,  idir, comp, inorm);

              pout() << "Coarse Error Norm = " << coarnorm << endl;
              pout() << "Fine   Error Norm = " << finenorm << endl;
              if((Abs(finenorm) > 1.0e-8) && (Abs(coarnorm) > 1.0e-8))
                {
                  Real order = log(coarnorm/finenorm)/log(2.0);
                  pout() << "Order of scheme = " << order << endl;
                }
            }
        }
      pout() << "==============================================" << endl;
    }
}
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
    Vector<LevelData<EBFluxFAB>* > scalFine( nlevels, NULL);
    Vector<LevelData<EBFluxFAB>* > scalMedi( nlevels, NULL);
    Vector<LevelData<EBFluxFAB>* > scalCoar( nlevels, NULL);
    Vector<LevelData<EBFluxFAB>* > errorMedi(nlevels, NULL);
    Vector<LevelData<EBFluxFAB>* > errorCoar(nlevels, NULL);


    pout() << "generating geometry only on finest level" << endl;
    AMRINSGeometry(params, domainFine);

    Vector<Vector<Box> > fineBoxes, mediBoxes, coarBoxes;
    getFixedGrids(fineBoxes, mediBoxes, coarBoxes, params, domainCoar);

    pout() << "generating fine extrapolated scalar " << endl;
    Real dt;
    getExtrapScalar(scalFine, dt, gridsFine, ebislFine, domainFine, params, true, fineBoxes);

    pout() << "generating medium extrapolated scalar and error" << endl;
    getExtrapScalar(scalMedi, dt, gridsMedi, ebislMedi, domainMedi, params, false, mediBoxes);

    getError(errorMedi,
             scalMedi, gridsMedi, ebislMedi, domainMedi,
             scalFine, gridsFine, ebislFine, domainFine,
             params);

    pout() << "generating coarse  extrapolated scalar and error" << endl;
    getExtrapScalar(scalCoar, dt, gridsCoar, ebislCoar, domainCoar, params, false, coarBoxes);

    getError(errorCoar,
             scalCoar, gridsCoar, ebislCoar, domainCoar,
             scalMedi, gridsMedi, ebislMedi, domainMedi,
             params);

    string testName("Scalar extrap error");
    Vector<string>   names(1, string("extrapError"));
    compareError(errorMedi,   errorCoar,
                 gridsMedi,   gridsCoar,
                 ebislMedi,   ebislCoar,
                 domainMedi,  domainCoar,
                 names, testName, params);

  }// End scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif
}
/**********************/
/**********************/
