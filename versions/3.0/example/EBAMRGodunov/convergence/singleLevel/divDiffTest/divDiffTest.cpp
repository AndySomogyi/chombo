#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <cmath>
#include <iostream>

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"

#include "EBExplosionIBCFactory.H"
#include "ModianoIBCFactory.H"
#include "EBPatchPolytropicFactory.H"

#include "EBPatchPolytropic.H"

#include "BaseIVFactory.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "VoFIterator.H"
#include "EBArith.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBLGIntegrator.H"
#include "EBAMRIO.H"
#include "PlaneIF.H"
#include "IntersectionIF.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#if CH_SPACEDIM==2
IntVect ivdebug(7,1);
#elif CH_SPACEDIM==3
IntVect ivdebug(1, 1, 0);
#endif

/***************/
/***************/
void
makeGeometry(Box& a_domain,
             Real& a_dx)
{
  //parse input file.  single level
  ParmParse pp;
  RealVect origin = RealVect::Zero;
  Vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);

 CH_assert(n_cell.size() == SpaceDim);
  IntVect lo = IntVect::Zero;
  IntVect hi;
  for(int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if(n_cell[ivec] <= 0)
        {
          pout() << " bogus number of cells input = " << n_cell[ivec];
          MayDay::Error();
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  a_domain = Box(lo, hi);

  Real prob_hi;
  int numOpen = n_cell[0];
  pp.get("domain_length",prob_hi);
  a_dx = prob_hi/numOpen;

  pout() << "channel geometry" << endl;
  RealVect channelNormal;
  vector<Real>  channelNormalVect(SpaceDim);
  pp.getarr("channel_normal",channelNormalVect, 0, SpaceDim);
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      channelNormal[idir] = channelNormalVect[idir];
    }

  Real channelRadius;
  pp.get("channel_radius", channelRadius);

  Real norm = 0.0;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      norm += channelNormal[idir] * channelNormal[idir];
    }
  norm = sqrt(norm);
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      channelNormal[idir] /= norm;
    }

  RealVect botPoint = RealVect::Zero;
  botPoint -= channelNormal;
  botPoint *= channelRadius;

  RealVect topPoint = RealVect::Zero;
  topPoint += channelNormal;
  topPoint *= channelRadius;

  PlaneIF bot(channelNormal,botPoint,true);
  PlaneIF top(channelNormal,topPoint,false);
  IntersectionIF channel(top,bot);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(channel,0,vectDx);
  //this generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  int bigboxsize=2048;
  ebisPtr->define(a_domain, origin, a_dx, workshop, bigboxsize);

  Vector<Box> vbox(1, a_domain);
  Vector<int> proc(1, 0);
  DisjointBoxLayout dbl(vbox, proc);
  EBISLayout ebisl;
  ebisPtr->fillEBISLayout(ebisl, dbl, a_domain, 0);
}

/***************/
/***************/
Real
norm(const BaseIVFAB<Real>& a_error,
     const IntVectSet& a_set,
     const EBISBox& a_ebisBox,
     const int& a_comp,
     const int& a_normtype)
{
 CH_assert(a_normtype >= 0);
 CH_assert(a_normtype <= 2);
  Real retval = 0.;
  VoFIterator vofit(a_set, a_ebisBox.getEBGraph());
  if(a_normtype == 0)
    {
      //max norm
      for(vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real val = a_error(vof, a_comp);
          retval = Max(Abs(val), retval);
        }
    }
  else
    {
      //integral norm
      Real volTot = 0.0;
      Real normTot = 0.0;
      for(vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real volFrac = a_ebisBox.volFrac(vof);
          Real val = a_error(vof, a_comp);
          volTot += volFrac;
          if(a_normtype == 1)
            normTot += Abs(val)*volFrac;
          else
            normTot += val*val*volFrac;
        }
      if(volTot > 1.0e-8)
        normTot /= volTot;
      if(a_normtype == 2)
        normTot = sqrt(normTot);

      retval = normTot;
    }
  return retval;
}
/***************/
/***************/
void
compareError(const BaseIVFAB<Real>& a_errorFine,
             const EBISBox& a_ebisBoxFine,
             const BaseIVFAB<Real>& a_errorCoar,
             const EBISBox& a_ebisBoxCoar)
{
  EBPatchPolytropic patchInt;
  IntVectSet ivsFine = a_errorFine.getIVS();
  IntVectSet ivsCoar = a_errorCoar.getIVS();
  for(int comp = 0; comp < a_errorFine.nComp(); comp++)
    {
      pout() << "==============================================" << endl;
      pout() << "Comparing error in variable  " << comp << endl;
      for(int inorm = 0; inorm < 2; inorm++)
        {

          if(inorm == 0)
            {
              pout() << endl << "Using max norm." << endl;
            }
          else
            {
              pout() << endl << "Using L-" << inorm << "norm." << endl;
            }
          Real coarnorm = norm(a_errorCoar,
                               ivsCoar, a_ebisBoxCoar,
                               comp, inorm);
          Real finenorm = norm(a_errorFine,
                               ivsFine, a_ebisBoxFine,
                               comp, inorm);

          pout() << "Coarse Error Norm = " << coarnorm << endl;
          pout() << "Fine   Error Norm = " << finenorm << endl;
          if((Abs(finenorm) > 1.0e-8) && (Abs(coarnorm) > 1.0e-8))
            {
              Real order = log(coarnorm/finenorm)/log(2.0);
              pout() << "Order of scheme = " << order << endl;
            }
        }
      pout() << "==============================================" << endl;
    }
}
/************/
/************/
void
diffDivUAutopsy()
{
  //make layouts == domain
  Box domainBoxFine, domainBoxCoar;
  Real dxFine, dxCoar;
  bool is_periodic[SpaceDim];
  for (int i = 0; i < SpaceDim; ++i)
    {
      is_periodic[i] = false;
    }

  makeGeometry(domainBoxFine, dxFine);
  domainBoxCoar = coarsen(domainBoxFine, 2);
  ProblemDomain domainFine(domainBoxFine.smallEnd(),
                           domainBoxFine.bigEnd(),
                           is_periodic);
  ProblemDomain domainCoar(domainBoxCoar.smallEnd(),
                           domainBoxCoar.bigEnd(),
                           is_periodic);

  dxCoar = 2.*dxFine;
  Vector<Box> boxFine(1, domainBoxFine);
  Vector<Box> boxCoar(1, domainBoxCoar);
  Vector<int> proc(1, 0);
  DisjointBoxLayout dblFine(boxFine, proc);
  DisjointBoxLayout dblCoar;
  coarsen(dblCoar, dblFine, 2);
  //source is a placeholder.  The others
  //get defined and used for both coarse and
  //fine but i don't use them after
  EBCellFAB source;

  RealVect modianoAxis, center;
  Real waveAmp, waveWidth, gamma;

  ParmParse pp;
  pp.get("wave_amplitude", waveAmp);
  pp.get("wave_width", waveWidth);
  pp.get("gamma",gamma);

  int idoNegWave;
  pp.get("use_negative_wave",idoNegWave);
  bool useNegativeWave = (idoNegWave==1);

  if(useNegativeWave)
    {
      pout() << "u - c wave " << endl;
    }
  else
    {
      pout() << "u + c wave " << endl;
    }
  RealVect channelNormal;
  vector<Real>  channelNormalVect(SpaceDim);
  vector<Real> centerpp(SpaceDim,0.5);
  pp.getarr("channel_normal",channelNormalVect, 0, SpaceDim);
  pp.getarr("wave_center",centerpp,0,SpaceDim);

  for(int idir = 0; idir < SpaceDim; idir++)
    {
      channelNormal[idir] = channelNormalVect[idir];
      center[idir] = centerpp[idir];
    }
  modianoAxis = mibcGetNormalVector(channelNormal);;

  ModianoIBC ibcFine(gamma, waveAmp, waveWidth, center, modianoAxis, false, useNegativeWave);
  ModianoIBC ibcCoar(gamma, waveAmp, waveWidth, center, modianoAxis, false, useNegativeWave);
  ModianoIBCFactory ibcFactFine(gamma, waveAmp, waveWidth, center, modianoAxis, false, useNegativeWave);
  ModianoIBCFactory ibcFactCoar(gamma, waveAmp, waveWidth, center, modianoAxis, false, useNegativeWave);

  //make ebislayouts
  //no need for ghost since we are using the domain
  //for the box
  Real cfl;
  pp.get("cfl", cfl);

  EBISLayout ebislFine, ebislCoar;
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  ebisPtr->fillEBISLayout(ebislFine, dblFine, domainBoxFine, 0);
  ebisPtr->fillEBISLayout(ebislCoar, dblCoar, domainBoxCoar, 0);
  int nCons = CNUM;
  int ifourth, iflatten, iartvisc;
  pp.get("use_fourth_order_slopes", ifourth);
  pp.get("use_flattening"         , iflatten);
  pp.get("use_art_visc"           , iartvisc);
  bool useFourthOrderSlopes = (ifourth  ==1);
  bool useFlattening        = false;
  bool useArtificialVisc    = (iartvisc ==1);

  for(DataIterator dit = dblFine.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBoxFine = ebislFine[dit()];
      const EBISBox& ebisBoxCoar = ebislCoar[dit()];

      //these get filled by modianoibc
      EBCellFAB consStateFine(ebisBoxFine, domainBoxFine, nCons);
      EBCellFAB consStateCoar(ebisBoxCoar, domainBoxCoar, nCons);
      ibcFine.define(domainFine, dxFine);
      ibcCoar.define(domainCoar, dxCoar);
      ibcFine.initialize(consStateFine, ebisBoxFine);
      ibcCoar.initialize(consStateCoar, ebisBoxCoar);

      EBPatchPolytropic ebppFine, ebppCoar;
      ebppFine.define(domainFine, dxFine);
      ebppCoar.define(domainCoar, dxCoar);

      //get the max wave speed and set the time step
      IntVectSet cfivs;
      ebppFine.setValidBox(domainBoxFine, ebisBoxFine, cfivs,0., 0.);
      ebppCoar.setValidBox(domainBoxCoar, ebisBoxCoar, cfivs,0., 0.);
      ebppFine.setGamma(gamma);
      ebppCoar.setGamma(gamma);
      bool useLimiting = false;
      ebppFine.setSlopeParameters(useFourthOrderSlopes, useFlattening, useLimiting);
      ebppCoar.setSlopeParameters(useFourthOrderSlopes, useFlattening, useLimiting);
      ebppFine.artificialViscosity(useArtificialVisc);
      ebppCoar.artificialViscosity(useArtificialVisc);
      ebppFine.setEBPhysIBC(ibcFactFine);
      ebppCoar.setEBPhysIBC(ibcFactCoar);
      Real maxWaveSpeedFine=
        ebppFine.getMaxWaveSpeed(consStateFine, domainBoxFine);
      Real maxWaveSpeedCoar=
        ebppCoar.getMaxWaveSpeed(consStateCoar, domainBoxCoar);
      Real dtFine = cfl*dxFine/maxWaveSpeedFine;
      Real dtCoar = cfl*dxCoar/maxWaveSpeedCoar;

      ebppFine.setValidBox(domainBoxFine, ebisBoxFine, cfivs,0., dtFine);
      ebppCoar.setValidBox(domainBoxCoar, ebisBoxCoar, cfivs,0., dtCoar);


      EBFluxFAB fluxFine(ebisBoxFine, domainBoxFine, nCons);
      EBFluxFAB fluxCoar(ebisBoxCoar, domainBoxCoar, nCons);
      //these get defined in computeFluxes
      EBCellFAB  slopePrimFine[SpaceDim], slopeNLimFine[SpaceDim];
      EBCellFAB  slopePrimCoar[SpaceDim], slopeNLimCoar[SpaceDim];
      EBCellFAB& primStateFine = ebppFine.getPrimState();
      EBCellFAB& primStateCoar = ebppCoar.getPrimState();
      //computed covered fluxes
      BaseIVFAB<Real>*  coveredFluxMinuFine= ebppFine.getCoveredFluxMinu();
      BaseIVFAB<Real>*  coveredFluxPlusFine= ebppFine.getCoveredFluxPlus();
      BaseIVFAB<Real>*  coveredFluxMinuCoar= ebppCoar.getCoveredFluxMinu();
      BaseIVFAB<Real>*  coveredFluxPlusCoar= ebppCoar.getCoveredFluxPlus();

      //vofs over which covered vofs live
      Vector<VolIndex>* coveredFaceMinuFine= ebppFine.getCoveredFaceMinu();
      Vector<VolIndex>* coveredFacePlusFine= ebppFine.getCoveredFacePlus();
      Vector<VolIndex>* coveredFaceMinuCoar= ebppCoar.getCoveredFaceMinu();
      Vector<VolIndex>* coveredFacePlusCoar= ebppCoar.getCoveredFacePlus();

      pout() << "using computed fluxes" << endl;
      EBCellFAB flat;
      bool verbose = false;
      ebppFine.computeFluxes(fluxFine,
                             coveredFluxMinuFine, coveredFluxPlusFine,
                             coveredFaceMinuFine, coveredFacePlusFine,
                             primStateFine, slopePrimFine, slopeNLimFine, flat,
                             consStateFine, source, domainBoxFine, dit(),verbose);

      ebppCoar.computeFluxes(fluxCoar,
                             coveredFluxMinuCoar, coveredFluxPlusCoar,
                             coveredFaceMinuCoar, coveredFacePlusCoar,
                             primStateCoar, slopePrimCoar, slopeNLimCoar, flat,
                             consStateCoar, source, domainBoxCoar, dit(),verbose);


      EBCellFAB divFNCFine(ebisBoxFine, domainBoxFine, nCons);
      EBCellFAB divFNCCoar(ebisBoxCoar, domainBoxCoar, nCons);
      divFNCFine.setVal(0.);
      divFNCCoar.setVal(0.);

      //this is the stable, non-conservative estimate of the solution update
      ebppFine.nonconservativeDivergence(divFNCFine, fluxFine,
                                         coveredFluxMinuFine, coveredFluxPlusFine,
                                         coveredFaceMinuFine, coveredFacePlusFine,
                                         domainBoxFine);

      ebppCoar.nonconservativeDivergence(divFNCCoar, fluxCoar,
                                         coveredFluxMinuCoar, coveredFluxPlusCoar,
                                         coveredFaceMinuCoar, coveredFacePlusCoar,
                                         domainBoxCoar);

      IntVectSet ivsIrregFine = ebisBoxFine.getIrregIVS(domainBoxFine);
      IntVectSet ivsIrregCoar = ebisBoxCoar.getIrregIVS(domainBoxCoar);
      IntVectSet ivsGrownFine = grow(ivsIrregFine, 1);
      IntVectSet ivsGrownCoar = grow(ivsIrregCoar, 1);
      ivsGrownFine &= domainBoxFine;
      ivsGrownCoar &= domainBoxCoar;

      // create sparse fluxes
      BaseIFFAB<Real>  centroidFluxFine[SpaceDim];
      BaseIFFAB<Real>  centroidFluxCoar[SpaceDim];
      BaseIFFAB<Real>* fluxInterpolantCoar[SpaceDim];
      BaseIFFAB<Real>* fluxInterpolantFine[SpaceDim];
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          centroidFluxFine[idir].define(ivsIrregFine, ebisBoxFine.getEBGraph(), idir, nCons);
          centroidFluxCoar[idir].define(ivsIrregCoar, ebisBoxCoar.getEBGraph(), idir, nCons);

          fluxInterpolantFine[idir] = new BaseIFFAB<Real>(ivsGrownFine, ebisBoxFine.getEBGraph(), idir, nCons);
          fluxInterpolantCoar[idir] = new BaseIFFAB<Real>(ivsGrownCoar, ebisBoxCoar.getEBGraph(), idir, nCons);
        }

      //copy fluxes into sparse fluxes
      FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          BaseIFFAB<Real>& interpolFine = *fluxInterpolantFine[idir];
          BaseIFFAB<Real>& interpolCoar = *fluxInterpolantCoar[idir];
          EBFaceFAB& fluxDirFine = fluxFine[idir];
          EBFaceFAB& fluxDirCoar = fluxCoar[idir];

          for(FaceIterator faceit(ivsGrownFine, ebisBoxFine.getEBGraph(), idir, stopCrit);
              faceit.ok(); ++faceit)
            {
              for(int ivar = 0; ivar < nCons; ivar++)
                {
                  interpolFine(faceit(), ivar) = fluxDirFine(faceit(), ivar);
                }
            }

          for(FaceIterator faceit(ivsGrownCoar, ebisBoxCoar.getEBGraph(), idir, stopCrit);
              faceit.ok(); ++faceit)
            {
              for(int ivar = 0; ivar < nCons; ivar++)
                {
                  interpolCoar(faceit(), ivar) = fluxDirCoar(faceit(), ivar);
                }
            }
        }
      //interpolate flux to centroids
      ebppFine.interpolateFluxToCentroids(centroidFluxFine,
                                          fluxInterpolantFine,
                                          ivsIrregFine);
      ebppCoar.interpolateFluxToCentroids(centroidFluxCoar,
                                          fluxInterpolantCoar,
                                          ivsIrregCoar);

      //clean up interpolant
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          delete fluxInterpolantFine[idir];
          delete fluxInterpolantCoar[idir];
          fluxInterpolantFine[idir] = NULL;
          fluxInterpolantCoar[idir] = NULL;
        }

      //compute flux at irregular boundaries.
      BaseIVFAB<Real> ebIrregFluxFine(ivsIrregFine, ebisBoxFine.getEBGraph(), nCons);
      BaseIVFAB<Real> ebIrregFluxCoar(ivsIrregCoar, ebisBoxCoar.getEBGraph(), nCons);

      EBCellFAB source;
      ebppFine.computeEBIrregFlux(ebIrregFluxFine, primStateFine,
                                  slopePrimFine, ivsIrregFine, source);
      ebppCoar.computeEBIrregFlux(ebIrregFluxCoar, primStateCoar,
                                  slopePrimCoar, ivsIrregCoar, source);

      //compute conservative divergence
      BaseIVFAB<Real> kDivFCFine(ivsIrregFine, ebisBoxFine.getEBGraph(), nCons);
      BaseIVFAB<Real> kDivFCCoar(ivsIrregCoar, ebisBoxCoar.getEBGraph(), nCons);
      kDivFCFine.setVal(0.);
      kDivFCCoar.setVal(0.);
      ebppFine.consUndividedDivergence(kDivFCFine,
                                       centroidFluxFine, ebIrregFluxFine, ivsIrregFine);
      ebppCoar.consUndividedDivergence(kDivFCCoar,
                                       centroidFluxCoar, ebIrregFluxCoar, ivsIrregCoar);

      BaseIVFAB<Real> errorFine(ivsIrregFine, ebisBoxFine.getEBGraph(), nCons);
      BaseIVFAB<Real> errorCoar(ivsIrregCoar, ebisBoxCoar.getEBGraph(), nCons);
      errorFine.setVal(0.);
      errorCoar.setVal(0.);
      //estimate of truncation error = kappa*(divfnc - divfc)
      for(VoFIterator vofit(ivsIrregFine, ebisBoxFine.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          for(int ivar = 0; ivar < nCons; ivar++)
            {
              Real volFrac = ebisBoxFine.volFrac(vof);
              Real kdivfc = kDivFCFine(vof, ivar);
              Real divfnc = divFNCFine(vof, ivar);
              Real divDiff = volFrac*divfnc - kdivfc;
              errorFine(vof, ivar) = divDiff;
            }
        }
      for(VoFIterator vofit(ivsIrregCoar, ebisBoxCoar.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          for(int ivar = 0; ivar < nCons; ivar++)
            {
              Real volFrac = ebisBoxCoar.volFrac(vof);
              Real kdivfc = kDivFCCoar(vof, ivar);
              Real divfnc = divFNCCoar(vof, ivar);
              Real divDiff = volFrac*divfnc - kdivfc;
              errorCoar(vof, ivar) = divDiff;
            }
        }

      //do the whole convergence test thing.
      compareError(errorFine, ebisBoxFine,
                   errorCoar, ebisBoxCoar);

    }

}
/************/
/************/
int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
  {
    // setChomboMPIErrorHandler();
#endif
    // Check for an input file
    char* inFile = NULL;


    if (a_argc > 1)
      {
        inFile = a_argv[1];
      }
    else
      {
        pout() << "Usage: <executable name> <inputfile>" << endl;
        pout() << "No input file specified" << endl;
        return -1;
      }
    ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

    diffDivUAutopsy();

#ifdef CH_MPI
    dumpmemoryatexit();
  }
  MPI_Finalize();
#endif
}

