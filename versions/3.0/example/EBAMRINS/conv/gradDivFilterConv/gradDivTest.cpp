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
#include "EBArith.H"
#include "EBAMRDataOps.H"
#include "NeumannPoissonEBBC.H"
#include "UsingNamespace.H"

using std::cerr;

#define CHECK_FLOATING_PT 0

#if CHECK_FLOATING_PT==1
#include <fenv.h>
#endif

//-1 =>  no divergence (exact solution = 0)
//0  =>  vel = sine functiton
//1  =>  vel = cubic
const int  g_exactProb = 1;
const Real g_pi = 4.*atan(1.0);
/******/
Real getVelExact(const RealVect& a_xval, const RealVect& a_freq,
                 const RealVect& a_magnitude, int a_idir, const Real& a_dx)
{
  int difDir = 0;
  if(a_idir  == 0)
    {
      difDir = 1;
    }
  //-1 =>  no divergence (exact solution = 0)
  //0  =>  vel = sine functiton
  //1  =>  vel = cubic
  Real valNoDiv   = sin(a_freq[difDir]*g_pi*a_xval[difDir]);
  Real valWithDiv;
  if(g_exactProb == 0)
    {
      valWithDiv = sin(a_freq[a_idir]*g_pi*a_xval[a_idir]);
    }
  else if(g_exactProb == 1)
    {
      valWithDiv = a_xval[a_idir]*a_xval[a_idir]*a_xval[a_idir];
    }
  else if(g_exactProb == -1)
    {
      valWithDiv = 0.;
    }
  else
    {
      MayDay::Error("bogus problem flag in getVelExact");
    }

  Real exactVal = valNoDiv + valWithDiv;

  return exactVal;
}

/******/
Real getGDUExact(const RealVect& a_xval, const RealVect& a_freq,
                 const RealVect& a_magnitude, int a_idir, const Real& a_dx)
{
  Real coeff = a_freq[a_idir]*g_pi;
  Real exactVal;

  if(g_exactProb == 0)
    {
      //corresponds to 2 derivatives of
      //Real valWithDiv = sin(a_freq[a_idir]*g_pi*a_xval[a_idir]);
      exactVal = -coeff*coeff*sin(coeff*a_xval[a_idir]);
    }
  else if(g_exactProb == 1)
    {
      //corresponds to 2 derivatives of
      //valWithDiv = a_xval[a_idir]*a_xval[a_idir]*a_xval[a_idir];
      exactVal = 6.*a_xval[a_idir];
    }
  else if(g_exactProb == -1)
    {
      exactVal = 0.;
    }
  else
    {
      MayDay::Error("bogus problem flag in getGDUExact");
    }

  return exactVal;
}
/******/
Real getVelExact(const VolIndex& a_vof, const Real& a_dx,
                 const RealVect& a_freq,
                 const RealVect& a_magnitude, int a_idir)
{
  RealVect vectDx = a_dx*RealVect::Unit;
  RealVect xval = EBArith::getVofLocation(a_vof, vectDx, RealVect::Zero);
  return getVelExact(xval, a_freq, a_magnitude, a_idir, a_dx);
}
/******/
Real getGDUExact(const VolIndex& a_vof, const Real& a_dx,
                 const RealVect& a_freq,
                 const RealVect& a_magnitude, int a_idir,
                 const RealVect& a_centroid)
{
  RealVect vectDx = a_dx*RealVect::Unit;
  RealVect offset = a_centroid;
  offset *= a_dx;
  RealVect xval = EBArith::getVofLocation(a_vof, vectDx, offset);
  return getGDUExact(xval, a_freq, a_magnitude, a_idir, a_dx);
}
/******/
Real getVelExact(const FaceIndex& a_face, const Real& a_dx,
                 const RealVect& a_freq,
                 const RealVect& a_magnitude)
{
  RealVect vectDx = a_dx*RealVect::Unit;
  RealVect xval = EBArith::getFaceLocation(a_face, vectDx, RealVect::Zero);
  return getVelExact(xval, a_freq, a_magnitude, a_face.direction(), a_dx);
}
/******/
void setExactStuff(LevelData<EBCellFAB>&                 a_cellVeloc,
                   LevelData<EBCellFAB>&                 a_gradDiv,
                   LevelData<EBFluxFAB>&                 a_fluxVeloc,
                   const DisjointBoxLayout&              a_grids,
                   const EBISLayout&                     a_ebisl,
                   const Real&                           a_dx,
                   const AMRParameters&                  a_params)
{
  ParmParse pp;
  Vector<Real> frequenciesVec(SpaceDim, 1.0);
  Vector<Real> magnitudesVec(SpaceDim, 1.0);
  pp.getarr("velocity_frequencies", frequenciesVec, 0, SpaceDim);
  pp.getarr("velocity_magnitudes",   magnitudesVec, 0, SpaceDim);
  RealVect frequencies, magnitudes;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      frequencies[idir]= frequenciesVec[idir];
      magnitudes [idir]= magnitudesVec [idir];
    }
  for(DataIterator dit = a_grids.dataIterator(); dit.ok();++dit)
    {
      EBCellFAB& cellVel = a_cellVeloc[dit()];
      EBCellFAB& cellGDU = a_gradDiv[dit()];
      EBFluxFAB& fluxVel = a_fluxVeloc[dit()];
      cellVel.setVal(0.);
      cellGDU.setVal(0.);
      fluxVel.setVal(0.);
      IntVectSet ivsBox = IntVectSet(cellVel.getRegion());
      const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          for(VoFIterator vofit(ivsBox, ebgraph); vofit.ok(); ++vofit)
            {
              Real velexact  = getVelExact(vofit(), a_dx,  frequencies, magnitudes, idir);
              RealVect centroid = a_ebisl[dit()].centroid(vofit());
              cellVel(vofit(), idir) = velexact;
              Real gduexact  = getGDUExact(vofit(), a_dx,  frequencies, magnitudes, idir, centroid);
              cellGDU(vofit(), idir) = gduexact;
            }
          FaceIterator faceit(ivsBox, ebgraph, idir, FaceStop::SurroundingWithBoundary);
          for(faceit.reset(); faceit.ok(); ++faceit)
            {
              fluxVel[idir](faceit(), 0) = getVelExact(faceit(), a_dx,  frequencies, magnitudes);
            }
        }
    }
}
/***/
void getError(Vector< LevelData<EBCellFAB>* >&     a_error,
              const Vector< DisjointBoxLayout >&         a_grids,
              const Vector< EBISLayout >&                a_ebisl,
              const ProblemDomain&                       a_level0Domain,
              const Real&                                a_level0Dx,
              const AMRParameters&                       a_params)
{
  int nlevels = a_grids.size();
  Vector<LevelData<EBCellFAB>* > veloc(nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > gradDiv(nlevels, NULL);
  Vector<LevelData<EBFluxFAB>* > fluxVel(nlevels, NULL);
  a_error.resize(nlevels, NULL);

  Real dxLev = a_level0Dx;
  //set to exact velocity
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBFluxFactory ebfluxfact(a_ebisl[ilev]);
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      veloc[ilev]      = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim, 4*IntVect::Unit, ebcellfact);
      gradDiv[ilev]    = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,   IntVect::Zero, ebcellfact);
      a_error[ilev]    = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,   IntVect::Zero, ebcellfact);
      fluxVel[ilev]    = new LevelData<EBFluxFAB>(a_grids[ilev],        1, 4*IntVect::Unit, ebfluxfact);

      //set veloc to exact velocity, gradDiv exact grad div
      setExactStuff(*veloc[ilev], *gradDiv[ilev], *fluxVel[ilev], a_grids[ilev], a_ebisl[ilev], dxLev, a_params);
      dxLev /= a_params.m_refRatio[ilev];
    }

  EBAMRDataOps::setToZero(a_error);

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
          coarVelocPtr =    veloc[ilev-1];
        }
      int refRat = a_params.m_refRatio[ilev];
      EBGradDivFilter gdFilt(a_grids[ilev], coarGridsPtr, a_ebisl[ilev], coarEBISLPtr, domLev, dxLev*RealVect::Unit, refRat);

      //this puts calculated grad div into error
      bool lowOrderOneSided, noExtrapToCovered;
      ParmParse pp;
      pp.get("lowOrderOneSided", lowOrderOneSided);
      pp.get("noExtrapToCovered", noExtrapToCovered);
      gdFilt.gradDiv(*a_error[ilev], *veloc[ilev], *fluxVel[ilev], coarVelocPtr, lowOrderOneSided, noExtrapToCovered);
      //gradDiv already holds exact grad div, error holds calculated gradDiv
      //this subracts off the exact grad div so error will now hold error = caculated gradDiv-exact gradDiv
      EBLevelDataOps::incr(*a_error[ilev], *gradDiv[ilev], -1.0);
      dxLev /= refRat;
      domLev.refine(refRat);
    }

  //need to keep error around so it can be used
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete veloc[ilev];
      delete gradDiv[ilev];
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
    Vector<LevelData<EBCellFAB>* > errorFine( nlevels, NULL);
    Vector<LevelData<EBCellFAB>* > errorMedi( nlevels, NULL);


    pout() << "generating geometry only on finest level" << endl;
    AMRINSGeometry(params, domainFine);

    Vector<Vector<Box> > fineBoxes, mediBoxes, coarBoxes;
    getFixedLayouts(gridsFine,   gridsMedi,  gridsCoar,
                    ebislFine,   ebislMedi,  ebislCoar,
                    domainFine, domainMedi, domainCoar, params);


    pout() << "generating fine error" << endl;
    getError(errorFine, gridsFine, ebislFine, domainFine, dxFine, params);

    pout() << "generating medium error" << endl;
    getError(errorMedi, gridsMedi, ebislMedi, domainMedi, dxMedi, params);


    Vector<string> names(SpaceDim, string("velcomp"));
    string testName("Filter error");
    compareError(errorFine,   errorMedi,
                 gridsFine,             gridsMedi,
                 ebislFine,             ebislMedi,
                 domainFine,           domainMedi,
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
        Vector<string> names(errorFine[0]->nComp(), string("var"));
        outputError(errorFine,    errorMedi,
                    gridsFine,   gridsMedi,
                    domainFine,  domainMedi,
                    fileFine,     fileCoar,
                    names,
                    params);
      }

    for(int ilev = 0; ilev < nlevels; ilev++)
      {
        delete errorFine[ilev];
        delete errorMedi[ilev];
      }

  }// End scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif
}

