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
#include "AllRegularService.H"
#include "TiltedCylinderIF.H"
#include "SquareCylinderIF.H"
#include "TransformIF.H"
#include "PlaneIF.H"
#include "IntersectionIF.H"
#include "EBAMRIO.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

//vofs to print out in detail
void getDebugIVS(IntVectSet& a_ivs,
                 const bool& a_thisCalcCoar)
{
  a_ivs = IntVectSet();
  if(!a_thisCalcCoar)
    {
      IntVect ivdebugloFine(D_DECL(66,10,10));
      IntVect ivdebughiFine(D_DECL(68,11,10));
      a_ivs |= ivdebugloFine;
      a_ivs |= ivdebughiFine;
    }
}
/***************/
/***************/
void
makeFinestDomain(Box& a_domain,
                 Real& a_dx)
{
  //parse input file.  single level
  ParmParse pp;
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
}
void
makeGeometry(const Box& a_domain,
             Real& a_dx)
{
  RealVect origin = RealVect::Zero;
  ParmParse  pp;
  int whichGeom;
  pp.get("which_geom", whichGeom);
  Real channelRadius;
  pp.get("channel_radius", channelRadius);
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  int biggridsize;
  pp.get("max_grid_size", biggridsize);
  int maxcoarsen = 0;
  if(whichGeom == 0)
    {
      pout() << "using all reg geom" << endl;
      AllRegularService allreg;
      ebisPtr->define(a_domain, origin, a_dx, allreg, biggridsize, maxcoarsen);
    }
  else if(whichGeom == 1)
    {
      RealVect channelNormal;
      vector<Real>  channelNormalVect(SpaceDim);
      pp.getarr("channel_normal",channelNormalVect, 0, SpaceDim);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          channelNormal[idir] = channelNormalVect[idir];
        }
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
      ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize, maxcoarsen);
    }
  else if(whichGeom == 2)
    {
      RealVect channelNormal;
      vector<Real>  channelNormalVect(SpaceDim);
      pp.getarr("channel_normal",channelNormalVect, 0, SpaceDim);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          channelNormal[idir] = channelNormalVect[idir];
        }

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
      ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize, maxcoarsen);
    }
  else if(whichGeom == 3)
    {
      vector<Real>  cylinderAxisVect(SpaceDim);
      pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
      RealVect cylinderAxis;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          cylinderAxis[idir] = cylinderAxisVect[idir];
        }
      Real sum;
      PolyGeom::unifyVector(cylinderAxis, sum);

      int cylinderType;
      pp.get("cylinder_type", cylinderType);
      RealVect corner = RealVect::Zero;
      bool negativeInside = true;
      if(cylinderType == 0)
        {
          pout() << "using a tilted circular cylinder" << endl;

          TiltedCylinderIF tunnel(channelRadius, cylinderAxis, corner, negativeInside);

          RealVect vectDx = RealVect::Unit;
          vectDx *= a_dx;

          GeometryShop workshop(tunnel,0,vectDx);
          ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize, 0);
        }
      else if(cylinderType == 1)
        {
          pout() << "using a tilted square cylinder" << endl;
          //make square cylinder around x axis
          SquareCylinderIF sqcyl(channelRadius, negativeInside);
          TransformIF transformif(sqcyl);
          RealVect xAxis= BASISREALV(0);
          transformif.rotate(xAxis, cylinderAxis);
          transformif.translate(corner);

          RealVect vectDx = RealVect::Unit;
          vectDx *= a_dx;

          GeometryShop workshop(transformif,0,vectDx);
          ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize, 0);
        }
      else
        {
          MayDay::Error("bogus cylinder type");
        }
    }
  else
    {
      MayDay::Error("bogus geometry flag");
    }

  //Vector<Box> vbox(1, a_domain);
  //Vector<int> proc(1, 0);
  //DisjointBoxLayout dbl(vbox, proc);
  //EBISLayout ebisl;
  //ebisPtr->fillEBISLayout(ebisl, dbl, a_domain, 0);
}
/************/
/************/
void
compareError(const BaseIVFAB<Real>& a_ebErrorFine,
             const BaseIVFAB<Real>& a_ebErrorCoar,
             const IntVectSet&      a_ivsIrregFine,
             const IntVectSet&      a_ivsIrregCoar,
             const EBISBox&         a_ebisBoxFine,
             const EBISBox&         a_ebisBoxCoar,
             const string&  a_testname)
{
  pout() << "==============================================" << endl;
  pout() << a_testname << " test  " << endl;
  for(int inorm = 0; inorm <= 2; inorm++)
    {
      if(inorm == 0)
        {
          pout()  << "Using max norm." << endl;
        }
      else
        {
          pout()  << "Using L-" << inorm << " norm." << endl;
        }

      Real ebIrregNormCoar, ebIrregNormFine;

      int comp = 0;
      EBArith::irregNorm(ebIrregNormCoar,
                         a_ebErrorCoar,
                         a_ivsIrregCoar,
                         a_ebisBoxCoar,
                         comp,  inorm);

      EBArith::irregNorm(ebIrregNormFine,
                         a_ebErrorFine,
                         a_ivsIrregFine,
                         a_ebisBoxFine,
                         comp,  inorm);


      if(a_ivsIrregCoar.isEmpty())
        {
          pout() << "no irregular fine vofs" << endl;
        }
      else
        {
          pout() << "Coar  " << a_testname << " Error Norm = " << ebIrregNormCoar << endl;
        }
      if(a_ivsIrregFine.isEmpty())
        {
          pout() << "no irregular fine vofs" << endl;
        }
      else
        {
          pout() << "Fine " << a_testname  << " Error Norm = " << ebIrregNormFine << endl;
        }
      if((Abs(ebIrregNormCoar) > 1.0e-12) && (Abs(ebIrregNormFine) > 1.0e-12))
        {
          Real order = log(ebIrregNormCoar/ebIrregNormFine)/log(2.0);
          pout() << "Order of " << a_testname <<  " at irregular boundary = " << order << endl << endl;
        }

    }
}
/********/
void
outputError(const BaseIVFAB<Real>& a_ebPressError,
            const BaseIVFAB<Real>& a_ebVNormError,
            const BaseIVFAB<Real>& a_ebNormDotAxisError,
            const BaseIVFAB<Real>& a_normError,
            const BaseIVFAB<Real>& a_pstarNoVNormError,
            const BaseIVFAB<Real>& a_vDotTrueNorm,
            const IntVectSet&      a_ivsIrreg,
            const EBISBox&         a_ebisBox,
            bool                   a_isCoarse)
{
#ifdef USE_HDF5
  char filename[100];
  if(a_isCoarse)
    {
      sprintf(filename,"errorCoar.%dd.hdf5",SpaceDim);
    }
  else
    {
      sprintf(filename,"errorFine.%dd.hdf5",SpaceDim);
    }

  Box domain = a_ebisBox.getDomain();
  const int nvar = 11 + 2*SpaceDim;

  Vector<string> names(nvar);
  names[0]  = string("pstarError");
  names[1]  = string("pstarErrorXbndryArea");
  names[2]  = string("velNormError");
  names[3]  = string("velNormErrorXbndryArea");
  names[4]  = string("normDotAxisError");
  names[5]  = string("normDotAxisXbndryArea");
  names[6]  = string("pstarNoVNormError");
  names[7]  = string("pstarNoVNormErrorXBndryArea");
  names[8]  = string("velDotTrueNormError");
  names[9]  = string("velDotTrueNormErrorXbndryArea");
  for(int idir = 0;  idir < SpaceDim; idir++)
    {
      char nameNorm[80];
      char nameNormBA[80];
      sprintf(nameNorm, "normError%d", idir);
      sprintf(nameNormBA, "normErrorXBndryArea%d", idir);
      names[10+idir]  = string(nameNorm);
      names[10+idir+SpaceDim] = string(nameNormBA);
    }
  names[nvar-1] = string("BoundaryArea");

  EBCellFAB fabData(a_ebisBox, domain, nvar);
  fabData.setVal(0.);
  for(VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const Real& pressError = a_ebPressError(vof, 0);
      const Real& vnormError = a_ebVNormError(vof, 0);
      const Real& normDotAxisError = a_ebNormDotAxisError(vof, 0);
      const Real& pstarNoVNormError = a_pstarNoVNormError(vof, 0);
      const Real& vDotTrueNormError = a_vDotTrueNorm(vof,0);
      const Real& bndryArea = a_ebisBox.bndryArea(vof);

      fabData(vof, 0)  = pressError;
      fabData(vof, 1)  = pressError*bndryArea;
      fabData(vof, 2)  = vnormError;
      fabData(vof, 3)  = vnormError*bndryArea;
      fabData(vof, 4)  = normDotAxisError;
      fabData(vof, 5)  = normDotAxisError*bndryArea;
      fabData(vof, 6)  = pstarNoVNormError;
      fabData(vof, 7)  = pstarNoVNormError*bndryArea;
      fabData(vof, 8)  = vDotTrueNormError;
      fabData(vof, 9)  = vDotTrueNormError*bndryArea;
      for(int idir = 0;  idir < SpaceDim; idir++)
        {
          const Real& normError = a_normError(vof, idir);
          fabData(vof, 10 + idir) = normError;
          fabData(vof, 10 + SpaceDim + idir) = normError*bndryArea;
        }
      fabData(vof, nvar-1) = bndryArea;
    }


  Real dx=1.0;
  Real dt=1.0;
  Real time=1.0;
  Vector<Real> covval(nvar, 0.0);
  bool repcov = false;
  writeEBHDF5(filename, domain, fabData, names, domain, dx, dt, time, repcov, covval);
#endif
}
/************/
/************/
void
printOpenSelect(const BaseIVFAB<Real>&  a_ffError,
                const EBISBox&          a_ebisBox,
                const IntVectSet&       a_ivs)
{
  IntVectSet ivs = a_ivs;
  ivs &= a_ebisBox.getRegion();
  int ncomp = a_ffError.nComp();

  for(VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      pout() << setw(12)
             << setprecision(6)
             << setiosflags(ios::showpoint)
             << setiosflags(ios::scientific) ;
      const VolIndex& vof = vofit();
      pout() << vof << "   " ;
      for(int icomp = 0; icomp < ncomp; icomp++)
        {
          Real error = a_ffError(vof, icomp);
          pout() << error << "  ";
        }
      pout() << endl;
    }
}
/*************/
void
fluxAutopsy()
{
  //make layouts == domain
  Box domainBoxFine, domainBoxCoar;
  Real dxFine, dxCoar;
  bool is_periodic[SpaceDim];
  for (int i = 0; i < SpaceDim; ++i)
    {
      is_periodic[i] = false;
    }

  makeFinestDomain(domainBoxFine, dxFine);
  domainBoxCoar = coarsen(domainBoxFine, 2);
  dxCoar = 2.*dxFine;
  ProblemDomain domainFine(domainBoxFine.smallEnd(),
                           domainBoxFine.bigEnd(),
                           is_periodic);
  ProblemDomain domainCoar(domainBoxCoar.smallEnd(),
                           domainBoxCoar.bigEnd(),
                           is_periodic);

  Vector<Box> boxFine(1, domainBoxFine);
  Vector<Box> boxCoar(1, domainBoxCoar);
  Vector<int> proc(1, 0);
  DisjointBoxLayout dblFine(boxFine, proc);
  DisjointBoxLayout dblCoar;
  coarsen(dblCoar, dblFine, 2);
  //source is a placeholder.
  EBCellFAB source;
  EBCellFAB slopePrimCoar[SpaceDim], slopeNLimCoar[SpaceDim];
  EBCellFAB slopePrimFine[SpaceDim], slopeNLimFine[SpaceDim];

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
  int idoFreeStream;
  pp.get("free_stream_prob",idoFreeStream);
  bool doFreeStream = (idoFreeStream==1);
  if(doFreeStream)
    {
      pout() << "free stream modiano bump " << endl;
    }
  else
    {
      pout() << "simple wave prob " << endl;
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


  ModianoIBC ibcFine(gamma, waveAmp, waveWidth, center, modianoAxis, doFreeStream, useNegativeWave);
  ModianoIBC ibcCoar(gamma, waveAmp, waveWidth, center, modianoAxis, doFreeStream, useNegativeWave);
  ModianoIBCFactory ibcFactFine(gamma, waveAmp, waveWidth, center, modianoAxis, doFreeStream, useNegativeWave);
  ModianoIBCFactory ibcFactCoar(gamma, waveAmp, waveWidth, center, modianoAxis, doFreeStream, useNegativeWave);

  //make ebislayouts
  //no need for ghost since we are using the domain
  //for the box
  Real cfl;
  pp.get("cfl", cfl);

  EBISLayout ebislFine, ebislCoar;
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();

  makeGeometry(domainBoxFine, dxFine);
  ebisPtr->fillEBISLayout(ebislFine, dblFine, domainBoxFine, 0);
  makeGeometry(domainBoxCoar, dxCoar);
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
      //cfivs is empty because this is a single-level calc
      IntVectSet cfivs;
      ebppFine.setValidBox(domainBoxFine, ebisBoxFine, cfivs, 0., 0.);
      ebppCoar.setValidBox(domainBoxCoar, ebisBoxCoar, cfivs, 0., 0.);
      ebppFine.setGamma(gamma);
      ebppCoar.setGamma(gamma);
      bool useLimiting = false;
      ebppFine.setSlopeParameters(useFourthOrderSlopes, useFlattening, useLimiting);
      ebppCoar.setSlopeParameters(useFourthOrderSlopes, useFlattening, useLimiting);
      ebppFine.doRZCoords(false);
      ebppCoar.doRZCoords(false);
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

      //To make the exact stuff work the same as
      //the  computed stuff, define the EBFluxFAB data
      //before going into routines and define BaseIVFAB
      //data inside routines.
      //computed fluxes on grid
      EBFluxFAB  cfluxFine(ebisBoxFine, domainBoxFine, nCons);
      EBFluxFAB  cfluxCoar(ebisBoxCoar, domainBoxCoar, nCons);
      //exact fluxes on grid
      EBFluxFAB  exactFine(ebisBoxFine, domainBoxFine, nCons);
      EBFluxFAB  exactCoar(ebisBoxCoar, domainBoxCoar, nCons);
      cfluxFine.setVal(0.);
      cfluxCoar.setVal(0.);
      exactFine.setVal(0.);
      exactCoar.setVal(0.);

      //exact covered fluxes
      BaseIVFAB<Real>  coveredExactMinuFine[SpaceDim];
      BaseIVFAB<Real>  coveredExactPlusFine[SpaceDim];
      BaseIVFAB<Real>  coveredExactMinuCoar[SpaceDim];
      BaseIVFAB<Real>  coveredExactPlusCoar[SpaceDim];

      //computed covered fluxes
      BaseIVFAB<Real>*  coveredCfluxMinuFine= ebppFine.getCoveredFluxMinu();
      BaseIVFAB<Real>*  coveredCfluxPlusFine= ebppFine.getCoveredFluxPlus();
      BaseIVFAB<Real>*  coveredCfluxMinuCoar= ebppCoar.getCoveredFluxMinu();
      BaseIVFAB<Real>*  coveredCfluxPlusCoar= ebppCoar.getCoveredFluxPlus();

      //vofs over which covered vofs live
      Vector<VolIndex>* coveredFacesMinuFine= ebppFine.getCoveredFaceMinu();
      Vector<VolIndex>* coveredFacesPlusFine= ebppFine.getCoveredFacePlus();
      Vector<VolIndex>* coveredFacesMinuCoar= ebppCoar.getCoveredFaceMinu();
      Vector<VolIndex>* coveredFacesPlusCoar= ebppCoar.getCoveredFacePlus();

      //irregular sets for eb fluxes
      IntVectSet ivsIrregFine = ebisBoxFine.getIrregIVS(domainBoxFine);
      IntVectSet ivsIrregCoar = ebisBoxCoar.getIrregIVS(domainBoxCoar);
      // debug vofs== vofs to dump open flux errors of
      IntVectSet ivsDebugFine, ivsDebugCoar;
      getDebugIVS(ivsDebugFine, false);
      getDebugIVS(ivsDebugCoar, true);

      ivsDebugFine &= ivsIrregFine;
      ivsDebugCoar &= ivsIrregCoar;

      //eb irregular fluxes
      BaseIVFAB<Real> ebPressCalcuFine(ivsIrregFine, ebisBoxFine.getEBGraph(), 1);
      BaseIVFAB<Real> ebPressCalcuCoar(ivsIrregCoar, ebisBoxCoar.getEBGraph(), 1);
      BaseIVFAB<Real> ebPressExactFine(ivsIrregFine, ebisBoxFine.getEBGraph(), 1);
      BaseIVFAB<Real> ebPressExactCoar(ivsIrregCoar, ebisBoxCoar.getEBGraph(), 1);
      BaseIVFAB<Real> ebPressErrorFine(ivsIrregFine, ebisBoxFine.getEBGraph(), 1);
      BaseIVFAB<Real> ebPressErrorCoar(ivsIrregCoar, ebisBoxCoar.getEBGraph(), 1);
      BaseIVFAB<Real> ebPstarNoVNormFine(ivsIrregFine, ebisBoxFine.getEBGraph(), 1);
      BaseIVFAB<Real> ebPstarNoVNormCoar(ivsIrregCoar, ebisBoxCoar.getEBGraph(), 1);
      BaseIVFAB<Real> ebVNormErrorFine(ivsIrregFine, ebisBoxFine.getEBGraph(), 1);
      BaseIVFAB<Real> ebVNormErrorCoar(ivsIrregCoar, ebisBoxCoar.getEBGraph(), 1);
      BaseIVFAB<Real> ebNormDotAxisErrorFine(ivsIrregFine, ebisBoxFine.getEBGraph(), 1);
      BaseIVFAB<Real> ebNormDotAxisErrorCoar(ivsIrregCoar, ebisBoxCoar.getEBGraph(), 1);
      BaseIVFAB<Real> ebPstarNoVNormErrorCoar(ivsIrregCoar, ebisBoxCoar.getEBGraph(), 1);
      BaseIVFAB<Real> ebPstarNoVNormErrorFine(ivsIrregFine, ebisBoxFine.getEBGraph(), 1);
      BaseIVFAB<Real> ebVDotTrueNormErrorCoar(ivsIrregCoar, ebisBoxCoar.getEBGraph(), 1);
      BaseIVFAB<Real> ebVDotTrueNormErrorFine(ivsIrregFine, ebisBoxFine.getEBGraph(), 1);

      BaseIVFAB<Real> ebNormErrorCoar(ivsIrregCoar, ebisBoxCoar.getEBGraph(), SpaceDim);
      BaseIVFAB<Real> ebNormErrorFine(ivsIrregFine, ebisBoxFine.getEBGraph(), SpaceDim);

      // nonconservative divergence necessary for irreg fluxes
      EBCellFAB nonConsDivFFine(ebisBoxFine, domainBoxFine, nCons);
      EBCellFAB nonConsDivFCoar(ebisBoxCoar, domainBoxCoar, nCons);
      EBCellFAB flat;
      bool verbose = false;
      pout() << "computing fine fluxes " << endl;
      EBCellFAB& primStateFine = ebppFine.getPrimState();
      EBCellFAB& primStateCoar = ebppCoar.getPrimState();
      ebppFine.computeFluxes(cfluxFine,
                             coveredCfluxMinuFine, coveredCfluxPlusFine,
                             coveredFacesMinuFine, coveredFacesPlusFine,
                             primStateFine, slopePrimFine, slopeNLimFine,
                             flat, consStateFine, source, domainBoxFine,
                             dit(),verbose);

      pout() << "taking fine divergence " << endl;
      ebppFine.nonconservativeDivergence(nonConsDivFFine, cfluxFine,
                                         coveredCfluxMinuFine, coveredCfluxPlusFine,
                                         coveredFacesMinuFine, coveredFacesPlusFine,
                                         domainBoxFine);

      EBCellFAB source;
      RealVect modianoCorner = RealVect::Zero;
      pout() << "computing fine boundary fields " << endl;
      ebppFine.computeBoundaryPress(ebPressCalcuFine, ebVNormErrorFine,
                                    ebNormDotAxisErrorFine, ebNormErrorFine,
                                    ebPstarNoVNormFine, ebVDotTrueNormErrorFine,
                                    modianoAxis, modianoCorner,
                                    primStateFine, nonConsDivFFine,
                                    slopePrimFine, ivsIrregFine, source);

      pout() << "computing coarse fluxes " << endl;
      ebppCoar.computeFluxes(cfluxCoar,
                             coveredCfluxMinuCoar, coveredCfluxPlusCoar,
                             coveredFacesMinuCoar, coveredFacesPlusCoar,
                             primStateCoar, slopePrimCoar, slopeNLimCoar,
                             flat, consStateCoar, source, domainBoxCoar,
                             dit(),verbose);



      pout() << "taking coarse divergence " << endl;
      ebppCoar.nonconservativeDivergence(nonConsDivFCoar, cfluxCoar,
                                         coveredCfluxMinuCoar, coveredCfluxPlusCoar,
                                         coveredFacesMinuCoar, coveredFacesPlusCoar,
                                         domainBoxCoar);

      pout() << "computing fine coarse fields " << endl;
      ebppCoar.computeBoundaryPress(ebPressCalcuCoar, ebVNormErrorCoar,
                                    ebNormDotAxisErrorCoar, ebNormErrorCoar,
                                    ebPstarNoVNormCoar, ebVDotTrueNormErrorCoar,
                                    modianoAxis, modianoCorner,
                                    primStateCoar, nonConsDivFCoar,
                                    slopePrimCoar, ivsIrregCoar, source);

      IntVectSet coveredSetsMinuFine[SpaceDim];
      IntVectSet coveredSetsPlusFine[SpaceDim];
      IntVectSet coveredSetsMinuCoar[SpaceDim];
      IntVectSet coveredSetsPlusCoar[SpaceDim];
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          coveredSetsMinuFine[idir]=coveredCfluxMinuFine[idir].getIVS();
          coveredSetsMinuCoar[idir]=coveredCfluxMinuCoar[idir].getIVS();
          coveredSetsPlusFine[idir]=coveredCfluxPlusFine[idir].getIVS();
          coveredSetsPlusCoar[idir]=coveredCfluxPlusCoar[idir].getIVS();
        }
      pout() << "computing exact fluxes fields " << endl;
      //get the exact fluxes
      ibcFine.computeExactFluxes(exactFine,
                                 coveredExactMinuFine, coveredExactPlusFine,
                                 coveredFacesMinuFine, coveredFacesPlusFine,
                                 coveredSetsMinuFine, coveredSetsPlusFine,
                                 ebPressExactFine, ivsIrregFine,
                                 ebisBoxFine, domainBoxFine,
                                 domainBoxFine, 0.5*dtFine);

      ibcCoar.computeExactFluxes(exactCoar,
                                 coveredExactMinuCoar, coveredExactPlusCoar,
                                 coveredFacesMinuCoar, coveredFacesPlusCoar,
                                 coveredSetsMinuCoar, coveredSetsPlusCoar,
                                 ebPressExactCoar, ivsIrregCoar,
                                 ebisBoxCoar, domainBoxCoar,
                                 domainBoxCoar,0.5*dtCoar);


      for(VoFIterator vofit(ivsIrregFine, ebisBoxFine.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          ebPressErrorFine(vof, 0) =  ebPressCalcuFine(vof, 0) - ebPressExactFine(vof, 0);
          ebPstarNoVNormErrorFine(vof, 0) =  ebPstarNoVNormFine(vof, 0) - ebPressExactFine(vof, 0);
        }
      for(VoFIterator vofit(ivsIrregCoar, ebisBoxCoar.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          ebPressErrorCoar(vof, 0) =  ebPressCalcuCoar(vof, 0) - ebPressExactCoar(vof, 0);
          ebPstarNoVNormErrorCoar(vof, 0) =  ebPstarNoVNormCoar(vof, 0) - ebPressExactCoar(vof, 0);
        }
      //vnorm compared with 0 so no subtraction needed
      //same goes for modiano norm

      pout() << "comparing error " << endl;
      //do the whole convergence test thing.

      string testname("boundary pressure");
      compareError(ebPressErrorFine, ebPressErrorCoar,
                   ivsIrregFine, ivsIrregCoar,
                   ebisBoxFine, ebisBoxCoar,
                   testname);

      testname = string("normal velocity");
      compareError(ebVNormErrorFine, ebVNormErrorCoar,
                   ivsIrregFine, ivsIrregCoar,
                   ebisBoxFine, ebisBoxCoar,
                   testname);


      testname = string("eb normal dot axis");
      compareError(ebNormDotAxisErrorFine, ebNormDotAxisErrorCoar,
                   ivsIrregFine, ivsIrregCoar,
                   ebisBoxFine, ebisBoxCoar,
                   testname);


      testname = string("eb normal component ");
      compareError(ebNormErrorFine, ebNormErrorCoar,
                   ivsIrregFine, ivsIrregCoar,
                   ebisBoxFine, ebisBoxCoar,
                   testname);


      testname = string("pstar (no vnorm)");
      compareError(ebPstarNoVNormErrorFine, ebPstarNoVNormErrorCoar,
                   ivsIrregFine, ivsIrregCoar,
                   ebisBoxFine, ebisBoxCoar,
                   testname);

      testname = string("vel dot true normal");
      compareError(ebVDotTrueNormErrorFine, ebVDotTrueNormErrorCoar,
                   ivsIrregFine, ivsIrregCoar,
                   ebisBoxFine, ebisBoxCoar,
                   testname);


      pout() << "outputting debug vofs of ebPressErrorCoar " << endl;
      printOpenSelect(ebPressErrorCoar, ebisBoxCoar, ivsDebugCoar);
      pout() << "outputting debug vofs of ebPressErrorFine " << endl;
      printOpenSelect(ebPressErrorFine, ebisBoxFine, ivsDebugFine);


      pout() << "outputting debug vofs of ebVNormErrorCoar " << endl;
      printOpenSelect(ebVNormErrorCoar, ebisBoxCoar, ivsDebugCoar);
      pout() << "outputting debug vofs of ebVNormErrorFine " << endl;
      printOpenSelect(ebVNormErrorFine, ebisBoxFine, ivsDebugFine);

      pout() << "outputting debug vofs of ebNormDotAxisErrorCoar " << endl;
      printOpenSelect(ebNormDotAxisErrorCoar, ebisBoxCoar, ivsDebugCoar);
      pout() << "outputting debug vofs of ebNormDotAxisErrorFine " << endl;
      printOpenSelect(ebNormDotAxisErrorFine, ebisBoxFine, ivsDebugFine);


      pout() << "outputting debug vofs of ebNormErrorCoar " << endl;
      printOpenSelect(ebNormErrorCoar, ebisBoxCoar, ivsDebugCoar);
      pout() << "outputting debug vofs of ebNormErrorFine " << endl;
      printOpenSelect(ebNormErrorFine, ebisBoxFine, ivsDebugFine);

      pout() << "outputting debug vofs of ebPstarNoVNormErrorCoar " << endl;
      printOpenSelect(ebPstarNoVNormErrorCoar, ebisBoxCoar, ivsDebugCoar);
      pout() << "outputting debug vofs of ebPstarNoVNormErrorFine " << endl;
      printOpenSelect(ebPstarNoVNormErrorFine, ebisBoxFine, ivsDebugFine);

      pout() << "outputting debug vofs of ebVDotTrueNormErrorCoar " << endl;
      printOpenSelect(ebVDotTrueNormErrorCoar, ebisBoxCoar, ivsDebugCoar);
      pout() << "outputting debug vofs of ebebVDotTrueNormErrorFine " << endl;
      printOpenSelect(ebVDotTrueNormErrorFine, ebisBoxFine, ivsDebugFine);

      pout() << "outputting error " << endl;
      outputError(ebPressErrorCoar,
                  ebVNormErrorCoar,
                  ebNormDotAxisErrorCoar,
                  ebNormErrorCoar,
                  ebPstarNoVNormErrorCoar,
                  ebVDotTrueNormErrorCoar,
                  ivsIrregCoar,
                  ebisBoxCoar, true);

      outputError(ebPressErrorFine,
                  ebVNormErrorFine,
                  ebNormDotAxisErrorFine,
                  ebNormErrorFine,
                  ebPstarNoVNormErrorFine,
                  ebVDotTrueNormErrorFine,
                  ivsIrregFine,
                  ebisBoxFine, false);
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

    fluxAutopsy();

#ifdef CH_MPI
    dumpmemoryatexit();
  }
  MPI_Finalize();
#endif
}

