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
#include "EBAMRDataOps.H"



#define CHECK_FLOATING_PT 0

#if CHECK_FLOATING_PT==1
#include <fenv.h>
#endif
RealVect
getVelFromUTheta(const Real a_uTheta,
                 const RealVect& a_x,
                 const RealVect& a_vortCenter)
{//this function computes the velocity vector for radially symmetric cylindrical flow (u_r=u_z=0) about the prob_lo
  RealVect vel = RealVect::Zero;

  RealVect change = a_x - a_vortCenter;
  Real cylRadius = sqrt(pow(change[0],2)+pow(change[1],2));

  vel[0] = -change[1]*a_uTheta/cylRadius;
  vel[1] = change[0]*a_uTheta/cylRadius;

  return vel;
}
Real getVelExact(const VolIndex& a_vof, const Real& a_dx, const Real& a_domLen, const Real& a_freq,
                 const Real& a_magnitude, const int& a_idir, const int& a_velType)
{
  Real exactVal = 1.e99;
  IntVect iv = a_vof.gridIndex();
  RealVect x;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      x[idir] = a_dx*(0.5 + Real(iv[idir]));
    }
  const RealVect vortCenter = 0.5*RealVect::Unit;
  const RealVect change = x - vortCenter;
  const Real radius = sqrt(change.dotProduct(change));

  if(a_velType == 0)
    {
      Real pi = 4.*atan(1.0);
      exactVal = sin(a_freq*pi*x[a_idir]/a_domLen);
    }
  else if(a_velType == 1)
    {
      Real u_theta  = 0.5+0.5*tanh((radius-0.2)*50);
      u_theta      *= 0.5+0.5*tanh((-radius+0.3)*50);
      RealVect vel = getVelFromUTheta(u_theta,x,vortCenter);

      exactVal = vel[a_idir];
    }
  else if(a_velType == 2)
    {
      Real u_theta = exp(-pow(10*radius,2));
      u_theta *= radius;
      RealVect vel = getVelFromUTheta(u_theta,x,vortCenter);

      exactVal = vel[a_idir];
    }
  else if(a_velType == 3)
    {
      Real u_theta = exp(-pow(4*radius,2));
      u_theta *= radius;
      RealVect vel = getVelFromUTheta(u_theta,x,vortCenter);
      vel *= sin(200*x[a_idir]);
      exactVal = vel[a_idir];
    }
  else if(a_velType == 4)
    {
      if(a_idir==0)
        {
          //exactVal = x[0] - 0.5;
          exactVal = tanh(x[0]) - tanh(1.0)/2;
        }
      else
        {
          exactVal = 0.0;
        }
    }
  else
    {
      MayDay::Error("cellProjTest::bad velType");
    }
  return exactVal;
}
/******/
void setExactVeloc(LevelData<EBCellFAB>&                 a_veloc,
                   const DisjointBoxLayout&              a_grids,
                   const EBISLayout&                     a_ebisl,
                   const RealVect&                       a_dx,
                   const PoissonParameters&              a_params)
{
  ParmParse pp;
  Vector<Real> frequencies(SpaceDim, 1.0);
  Vector<Real> magnitudes(SpaceDim, 1.0);
  pp.getarr("velocity_frequencies", frequencies, 0, SpaceDim);
  pp.getarr("velocity_magnitudes", magnitudes, 0, SpaceDim);
  int velType;
  pp.get("velocity_type", velType);
  RealVect domainLength = a_params.domainLength;
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
              Real velexact = getVelExact(vofit(), a_dx[idir], domainLength[idir], frequencies[idir], magnitudes[idir], idir, velType);
              vel(vofit(), idir) = velexact;
            }
        }
    }
}
/******/
void getError(Vector< LevelData<EBCellFAB>* >&     a_error,
              const Vector< DisjointBoxLayout >&   a_grids,
              const Vector< EBISLayout >&          a_ebisl,
              const PoissonParameters&             a_params,
              const int&                           a_dofileout,
              const bool&                          a_isFine)
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
      setExactVeloc(*velo[ilev], a_grids[ilev], a_ebisl[ilev], dxLev, a_params);

      EBLevelDataOps::setToZero(*(a_error[ilev]));
      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);

    }

  Real domVal = 0.0;
  NeumannPoissonDomainBCFactory* domBCPhi = new NeumannPoissonDomainBCFactory();
  RefCountedPtr<BaseDomainBCFactory> baseDomainBCPhi = RefCountedPtr<BaseDomainBCFactory>(domBCPhi);
  domBCPhi->setValue(domVal);
  DirichletPoissonDomainBCFactory* domBCVel = new DirichletPoissonDomainBCFactory();
  RefCountedPtr<BaseDomainBCFactory> baseDomainBCVel = RefCountedPtr<BaseDomainBCFactory>(domBCVel);
  domBCVel->setValue(domVal);

  NeumannPoissonEBBCFactory*      ebBCPhi = new NeumannPoissonEBBCFactory();
  ebBCPhi->setValue(domVal);
  RefCountedPtr<BaseEBBCFactory>     baseEBBCPhi     = RefCountedPtr<BaseEBBCFactory>(ebBCPhi);


  const int bottomSolverType = 1;
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
                                       true, -1, 3 ,40,1.e99, 1,bottomSolverType);

  //get starting divergence
  Real norm[3];
  projectinator.kappaDivergence(a_error, velo);
  for(int inorm = 0; inorm < 3; inorm++)
    {
      norm[inorm] = EBArith::norm(a_error, a_grids, a_ebisl,a_params.refRatio, 0, inorm, EBNormType::OverBoth);
    }
  pout() << "div(vel) before cell project: " << "L_inf = " << norm[0]  << ", L_1 = " << norm[1] << ", L_2 = " << norm[2] << endl;

  pout() << "projecting velocity  "  << endl;
#ifdef CH_USE_HDF5
  if(a_dofileout == 1)
    {
      EBAMRDataOps::setCoveredVal(velo,0.0);
      if(a_isFine)
        {
          writeEBAMRname(   &velo,"veloBeforeFine.hdf5");
          writeEBAMRname(&a_error,"divuBeforeFine.hdf5");
        }
      else
        {
          writeEBAMRname(   &velo,"veloBeforeCoar.hdf5");
          writeEBAMRname(&a_error,"divuBeforeCoar.hdf5");
        }
    }
#endif
  projectinator.project(velo, gphi);
#ifdef CH_USE_HDF5
  if(a_dofileout == 1)
    {
      EBAMRDataOps::setCoveredVal(velo,0.0);
      if(a_isFine)
        {
          writeEBAMRname(&velo,"veloAfterFine.hdf5");
          writeEBAMRname(&gphi,"gphiAfterFine.hdf5");
        }
      else
        {
          writeEBAMRname(&velo,"veloAfterCoar.hdf5");
          writeEBAMRname(&gphi,"gphiAfterCoar.hdf5");
        }
    }
#endif
  projectinator.kappaDivergence(a_error, velo);
  for(int inorm = 0; inorm < 3; inorm++)
    {
      norm[inorm] = EBArith::norm(a_error, a_grids, a_ebisl,a_params.refRatio, 0, inorm, EBNormType::OverBoth);
    }
  pout() << "div(vel) after cell project: " << "L_inf = " << norm[0]  << ", L_1 = " << norm[1] << ", L_2 = " << norm[2] << endl;


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

    int dofileout;
    pp.get("do_error_output", dofileout);

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
    getError(errorFine, gridsFine, ebislFine, paramFine, dofileout, true);

    pout() << "defining coarse geometry" << endl;
    definePoissonGeometry(paramCoar);

    getCoarseLayoutsFromFine(gridsCoar, ebislCoar, gridsFine, paramCoar);

    pout() << "generating coarse error" << endl;
    getError(errorCoar, gridsCoar, ebislCoar, paramCoar, dofileout, false);

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
