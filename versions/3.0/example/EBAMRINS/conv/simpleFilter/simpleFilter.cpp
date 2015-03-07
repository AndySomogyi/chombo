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
#include "EBGradDivFilter.H"
#include "EBCoarseAverage.H"
#include "NeumannPoissonDomainBC.H"
#include "DirichletPoissonDomainBC.H"
#include "EBArith.H"
#include "UsingNamespace.H"


#define CHECK_FLOATING_PT 0

#if CHECK_FLOATING_PT==1
#include <fenv.h>
#endif


Real getVelExact(const VolIndex& a_vof, const RealVect& a_dx, const Real& a_freq,
                 const Real& a_magnitude, int a_idir)
{
  Real exactVal = 1.e99;
  const int probNum = 3;
  RealVect x = EBArith::getVofLocation(a_vof,a_dx,RealVect::Zero);
  if(probNum==1)
    {//sin function: u=f(y);v=f(x);w=f(x)
      int dDir = 0;
      if(a_idir  == 0)
        {
          dDir = 1;
        }
      exactVal = sin(a_freq*M_PI*x[dDir]);
    }
  else if(probNum==2)
    {//sin function: u=f(x);v=f(y);w=f(z)
      exactVal = sin(a_freq*M_PI*x[a_idir]);
    }
  else if(probNum==3)
    {//checkerboard
      int ivof = a_vof.gridIndex()[0];
      int jvof = a_vof.gridIndex()[1];
      if((ivof + jvof)%2 == 0)
        {
          exactVal =1.0;
        }
      else
        {
          exactVal =-1.0;
        }
    }
  else if(probNum==4)
    {//linear function: u=f(x);v=f(y);w=f(z)
      exactVal = x[a_idir];
    }
  else
    {
      MayDay::Error("simpleFilter::bad probNum");
    }
  return exactVal;
}
/******/
void setExactVeloc(LevelData<EBCellFAB>&                 a_veloc,
                   const DisjointBoxLayout&              a_grids,
                   const EBISLayout&                     a_ebisl,
                   const RealVect&                       a_dx,
                   const AMRParameters&                  a_params)
{
  IntVect ivdebug(D_DECL(46, 45, 0));
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
              Real velexact  = getVelExact(vofit(), a_dx,  frequencies[idir], magnitudes[idir], idir);
              vel(vofit(), idir) = velexact;
            }
        }
    }
}
/***/
void getFilteredVel(Vector< LevelData<EBCellFAB>* >&           a_veloc,
                    const Vector< DisjointBoxLayout >&         a_grids,
                    const Vector< EBISLayout >&                a_ebisl,
                    const ProblemDomain&                       a_level0Domain,
                    const RealVect&                            a_level0Dx,
                    const AMRParameters&                       a_params)
{
  int nlevels = a_grids.size();
  RealVect dxLev = a_level0Dx;
  Vector< LevelData<EBFluxFAB>* > fluxVel(nlevels, NULL);
  //set to exact velocity
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      EBFluxFactory ebfluxfact(a_ebisl[ilev]);
      a_veloc[ilev]      = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,  4*IntVect::Unit, ebcellfact);
      fluxVel[ilev]      = new LevelData<EBFluxFAB>(a_grids[ilev],        1,  4*IntVect::Unit, ebfluxfact);
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*fluxVel[ilev])[dit()].setVal(0.);
        }
      setExactVeloc(*a_veloc[ilev], a_grids[ilev], a_ebisl[ilev], dxLev, a_params);
      dxLev /= a_params.m_refRatio[ilev];
    }

  ParmParse pp;
  int numfilt;
  pp.get("num_filter_iterations", numfilt);
  for(int ifilt = 0; ifilt < numfilt; ifilt++)
    {
      pout() << "filter iteration = " << ifilt << endl;
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
              coarVelocPtr =  a_veloc[ilev-1];
            }
          int refRat = a_params.m_refRatio[ilev];
          EBGradDivFilter gdFilt(a_grids[ilev], coarGridsPtr, a_ebisl[ilev], coarEBISLPtr, domLev, dxLev, refRat);

          bool lowOrderOneSided, noExtrapToCovered;
          ParmParse pp;
          pp.get("lowOrderOneSided", lowOrderOneSided);
          pp.get("noExtrapToCovered", noExtrapToCovered);
          gdFilt.filter(*a_veloc[ilev], *fluxVel[ilev], coarVelocPtr, lowOrderOneSided, noExtrapToCovered);

          dxLev /= refRat;
          domLev.refine(refRat);
        }
    }
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete fluxVel[ilev];
    }

}
/***************/
void outputErrorLoc(const Vector< LevelData<EBCellFAB>* >&   a_errorFine,
                    const Vector< DisjointBoxLayout >&       a_gridsFine,
                    const ProblemDomain&                     a_level0DomainFine,
                    const string&                            a_fileFine,
                    const AMRParameters&                     a_params)

{
#ifdef CH_USE_HDF5
  bool replaceCovered = true;
  int nvar = SpaceDim;
  Vector<Real> coveredValues(nvar, 0.0);
  //values that don't matter in output file
  Real dxFine = 1.0;
  Real time   = 1.0;
  Real dt     = 1.0;

  ParmParse pp;

  Vector<int> refRatio = a_params.m_refRatio;
  int numlevels = a_params.m_maxLevel + 1;
  Box domainFine = a_level0DomainFine.domainBox();
  Vector<string> names(SpaceDim);
  for(int icomp = 0; icomp < SpaceDim; icomp++)
    {
      char charstr[90];
      sprintf(charstr, "velcomp%d",icomp);
      names[icomp] = string(charstr);
    }
  writeEBHDF5(a_fileFine, a_gridsFine, a_errorFine,  names,
              domainFine, dxFine, dt, time, refRatio, numlevels,
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

    ProblemDomain domainCoar, domainFine, domainMedi;
    AMRParameters params;
    getAMRINSParameters(params, domainCoar);
    domainMedi = refine(domainCoar, 2);
    domainFine = refine(domainMedi, 2);
    RealVect dxFine = RealVect::Unit;
    for(int idir = 0; idir < SpaceDim; idir++)
      {
        dxFine[idir] /= (Real(domainFine.size(idir)));
      }

    Vector<DisjointBoxLayout> gridsFine, gridsCoar, gridsMedi;
    Vector<EBISLayout>        ebislFine, ebislCoar, ebislMedi;

    int nlevels = params.m_maxLevel + 1;
    Vector<LevelData<EBCellFAB>* > filteredvelFine( nlevels, NULL);

    pout() << "generating geometry only on finest level" << endl;
    AMRINSGeometry(params, domainFine);

    Vector<Vector<Box> > fineBoxes, mediBoxes, coarBoxes;
    getFixedLayouts(gridsFine,   gridsMedi,  gridsCoar,
                    ebislFine,   ebislMedi,  ebislCoar,
                    domainFine, domainMedi, domainCoar, params);


    pout() << "generating fine filteredvel " << endl;
    getFilteredVel(filteredvelFine, gridsFine, ebislFine, domainFine, dxFine, params);


#if CH_SPACEDIM==2
    string fileFine("fineFilteredVel.2d.hdf5");
#else
    string fileFine("fineFilteredVel.3d.hdf5");
#endif
    int dofileout;
    pp.get("do_error_output", dofileout);
    if(dofileout == 1)
      {
        pout() << "outputting filtered velocity" << endl;
        outputErrorLoc(filteredvelFine,
                       gridsFine,
                       domainFine,
                       fileFine,
                       params);
      }

    for(int ilev = 0; ilev < nlevels; ilev++)
      {
        delete filteredvelFine[ilev];
      }

  }// End scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif
}

