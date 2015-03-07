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
#include "AllRegularService.H"
#include "PlaneIF.H"
#include "IntersectionIF.H"
#include "TiltedCylinderIF.H"

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

  RealVect channelNormal;
  vector<Real>  channelNormalVect(SpaceDim);
  pp.getarr("channel_normal",channelNormalVect, 0, SpaceDim);
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      channelNormal[idir] = channelNormalVect[idir];
    }

  Real channelRadius;
  pp.get("channel_radius", channelRadius);

  int whichGeom;
  pp.get("which_geom", whichGeom);
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  int biggridsize = 2048;
  if(whichGeom == 0)
    {
      pout() << "using all reg geom" << endl;
      AllRegularService allreg;
      ebisPtr->define(a_domain, origin, a_dx, allreg, biggridsize);
    }
  else if(whichGeom == 1)
    {
      Real alpha;
      pp.get("alpha", alpha);

      RealVect channelPoint = RealVect::Zero;
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (channelNormal[idir] != 0.0)
        {
          channelPoint[idir] = alpha / channelNormal[idir];
          break;
        }
      }

      bool inside = false;

      PlaneIF ramp(channelNormal,channelPoint,inside);

      pout() << "using a ramp" << endl;

      RealVect vectDx = RealVect::Unit;
      vectDx *= a_dx;

      GeometryShop workshop(ramp,0,vectDx);
      ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize);
    }
  else if(whichGeom == 2)
    {
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

      pout() << "using a channel" << endl;

      RealVect vectDx = RealVect::Unit;
      vectDx *= a_dx;

      GeometryShop workshop(channel,0,vectDx);
      ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize);
    }
  else if(whichGeom == 3)
    {
      pout() << "using a tilted cylinder" << endl;
      vector<Real>  cylinderAxisVect(SpaceDim);
      pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
      RealVect cylinderAxis;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          cylinderAxis[idir] = cylinderAxisVect[idir];
        }
      Real sum;
      PolyGeom::unifyVector(cylinderAxis, sum);

      pout() << "using a tilted cylinder implicit function" << endl;
      RealVect corner = RealVect::Zero;
      bool negativeInside = true;
      TiltedCylinderIF tunnel(channelRadius, cylinderAxis, corner, negativeInside);

      RealVect vectDx = RealVect::Unit;
      vectDx *= a_dx;

      GeometryShop workshop(tunnel,0,vectDx);
      ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize);

    }
  else
    {
      MayDay::Error("bogus geometry flag");
    }

  Vector<Box> vbox(1, a_domain);
  Vector<int> proc(1, 0);
  DisjointBoxLayout dbl(vbox, proc);
  EBISLayout ebisl;
  ebisPtr->fillEBISLayout(ebisl, dbl, a_domain, 0);
}
/***************/
/***************/
void
compareError(const EBCellFAB& a_errorFine,
             const EBISBox& a_ebisBoxFine,
             const EBCellFAB& a_errorCoar,
             const EBISBox& a_ebisBoxCoar)
{
  EBPatchPolytropic patchInt;
  Box gridFine = a_errorFine.getRegion();
  Box gridCoar = a_errorCoar.getRegion();
  pout() << "==============================================" << endl;
  for(int comp = 0; comp < a_errorFine.nComp(); comp++)
    {
      pout() << "Comparing error in variable  " << comp << endl;
      pout() << "==============================================" << endl;
      for(int itype = 0; itype < 3; itype++)
        {
          EBNormType::NormMode normtype;
          if(itype == 0)
            {
              normtype = EBNormType::OverBoth;
              pout() << endl << "Using all uncovered cells." << endl;
            }
          else if(itype == 1)
            {
              normtype = EBNormType::OverOnlyRegular;
              pout() << endl << "Using only regular cells." << endl;
            }
          else
            {
              normtype = EBNormType::OverOnlyIrregular;
              pout() << endl << "Using only irregular cells." << endl;
            }
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
              Real coarnorm = EBArith::norm(a_errorCoar,
                                            gridCoar, a_ebisBoxCoar,
                                            comp, inorm, normtype);
              Real finenorm = EBArith::norm(a_errorFine,
                                            gridFine, a_ebisBoxFine,
                                            comp, inorm, normtype);

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
/************/
/************/
void
divUAutopsy()
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
  EBCellFAB source, slopePrim[SpaceDim], slopeNLim[SpaceDim];

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
  vector<Real> centerpp(SpaceDim,0.5);
  pp.getarr("wave_center",centerpp,0,SpaceDim);
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      center[idir] = centerpp[idir];
    }

  int whichGeom;
  pp.get("which_geom", whichGeom);
  if(whichGeom == 3)
    {
      pout() << "uses cylinder axis for direction" << endl;
      vector<Real>  cylinderAxisVect(SpaceDim);
      pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
      RealVect cylinderAxis;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          cylinderAxis[idir] = cylinderAxisVect[idir];
        }
      Real sum;
      PolyGeom::unifyVector(cylinderAxis, sum);
      modianoAxis = cylinderAxis;
    }
  else
    {
      pout() << "uses channel slope for direction" << endl;
      RealVect channelNormal;
      vector<Real>  channelNormalVect(SpaceDim);
      pp.getarr("channel_normal",channelNormalVect, 0, SpaceDim);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          channelNormal[idir] = channelNormalVect[idir];
        }
      modianoAxis = mibcGetNormalVector(channelNormal);;
    }


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
      ebppFine.setValidBox(domainBoxFine, ebisBoxFine, cfivs, 0., 0.);
      ebppCoar.setValidBox(domainBoxCoar, ebisBoxCoar, cfivs, 0., 0.);
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

      ebppFine.setValidBox(domainBoxFine, ebisBoxFine, cfivs, 0., dtFine);
      ebppCoar.setValidBox(domainBoxCoar, ebisBoxCoar, cfivs, 0., dtCoar);


      EBCellFAB  errorFine(ebisBoxFine, domainBoxFine, nCons);
      EBCellFAB  errorCoar(ebisBoxCoar, domainBoxCoar, nCons);

      //exact solution on grid
      EBCellFAB  exactFineOld(ebisBoxFine, domainBoxFine, nCons);
      EBCellFAB  exactCoarOld(ebisBoxCoar, domainBoxCoar, nCons);
      EBCellFAB  exactFineNew(ebisBoxFine, domainBoxFine, nCons);
      EBCellFAB  exactCoarNew(ebisBoxCoar, domainBoxCoar, nCons);

      EBCellFAB nonConsDivFine(ebisBoxFine, domainBoxFine, nCons);
      EBCellFAB nonConsDivCoar(ebisBoxCoar, domainBoxCoar, nCons);

      nonConsDivFine.setVal(0.);
      nonConsDivCoar.setVal(0.);
      exactFineOld.setVal(0.);
      exactFineNew.setVal(0.);
      exactCoarOld.setVal(0.);
      exactCoarNew.setVal(0.);
      errorFine.setVal(0.);
      errorCoar.setVal(0.);

      //set exact values of state
      Real timeOldFine = 0.0;
      Real timeOldCoar = 0.0;

      Real timeNewFine = dtFine;
      Real timeNewCoar = dtCoar;

      ibcFine.setToExact(exactFineOld, ebisBoxFine, timeOldFine);
      ibcFine.setToExact(exactFineNew, ebisBoxFine, timeNewFine);
      ibcCoar.setToExact(exactCoarOld, ebisBoxCoar, timeOldCoar);
      ibcCoar.setToExact(exactCoarNew, ebisBoxCoar, timeNewCoar);

      EBFluxFAB fluxFine(ebisBoxFine, domainBoxFine, nCons);
      EBFluxFAB fluxCoar(ebisBoxCoar, domainBoxCoar, nCons);
      nonConsDivFine.setVal(0.);
      nonConsDivCoar.setVal(0.);

      EBCellFAB flat;
      pout() << "using computed fluxes" << endl;
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
      ebppFine.computeFluxes(fluxFine,
                             coveredFluxMinuFine, coveredFluxPlusFine,
                             coveredFaceMinuFine, coveredFacePlusFine,
                             primStateFine, slopePrim, slopeNLim, flat, consStateFine,
                             source, domainBoxFine, dit(),false);

      ebppCoar.computeFluxes(fluxCoar,
                             coveredFluxMinuCoar, coveredFluxPlusCoar,
                             coveredFaceMinuCoar, coveredFacePlusCoar,
                             primStateCoar, slopePrim, slopeNLim, flat, consStateCoar,
                             source, domainBoxCoar, dit(),false);


      //this is the stable, non-conservative estimate of the solution update
      ebppFine.nonconservativeDivergence(nonConsDivFine, fluxFine,
                                         coveredFluxMinuFine, coveredFluxPlusFine,
                                         coveredFaceMinuFine, coveredFacePlusFine,
                                         domainBoxFine);

      ebppCoar.nonconservativeDivergence(nonConsDivCoar, fluxCoar,
                                         coveredFluxMinuCoar, coveredFluxPlusCoar,
                                         coveredFaceMinuCoar, coveredFacePlusCoar,
                                         domainBoxCoar);

      //estimate of truncation error = 1/dt*(unew-uold) + divF
      Interval interv(0, nCons-1);

      errorFine.copy(domainBoxFine, interv, domainBoxFine, exactFineNew, interv);
      errorFine -=   exactFineOld;
      errorFine  /=  dtFine;
      errorFine  += nonConsDivFine;

      errorCoar.copy(domainBoxCoar, interv, domainBoxCoar, exactCoarNew, interv);
      errorCoar -=   exactCoarOld;
      errorCoar  /=  dtCoar;
      errorCoar  += nonConsDivCoar;


      //do the whole convergence test thing.
      compareError(errorFine, ebisBoxFine,
                   errorCoar, ebisBoxCoar);

      //output stuff to files

      Vector<Real> covValues(nCons, 0.0);
      string errorFileFine("errorFine.hdf5");
      string errorFileCoar("errorCoar.hdf5");
      Vector<string> names(nCons);
      names = ebppFine.stateNames();

#ifdef CH_USE_HDF5
      if(numProc() == 1)
        {
          writeEBHDF5(errorFileFine, domainBoxFine, errorFine, names,
                      domainBoxFine, dxFine, dtFine, timeNewFine, true, covValues);
          writeEBHDF5(errorFileCoar, domainBoxCoar, errorCoar, names,
                      domainBoxCoar, dxCoar, dtCoar, timeNewCoar, true, covValues);
        }
#endif
    }

}
/************/
/************/
int
main(int a_argc, char* a_argv[])
{
  int eekflag = 0;
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

    divUAutopsy();
#ifdef CH_MPI

  }
  dumpmemoryatexit();
  MPI_Finalize();
#endif
  return eekflag;
}

