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
#include <cstdio>
#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
#include <iomanip>
#include <string>

#include "REAL.H"
#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "BRMeshRefine.H"
#include "DisjointBoxLayout.H"
#include "CH_HDF5.H"
#include "AMRLevel.H"
#include "AMR.H"
#include "parstream.H"
#include "Tuple.H"
#include "BoxIterator.H"

#include "ParmParse.H"
#include "parstream.H"

#include "EBExplosionIBCFactory.H"
#include "ModianoIBCFactory.H"
#include "EBPatchPolytropicFactory.H"

#include "EBPatchPolytropic.H"

#include "EBAMRGodunovFactory.H"
#include "EBAMRGodunov.H"
#include "AMRLevel.H"
#include "AMR.H"
#include "BaseIVFactory.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "VoFIterator.H"
#include "EBArith.H"
#include "AllRegularService.H"
#include "EBLevelRedist.H"
#include "RedistStencil.H"
#include "SlabService.H"
#include "SquareCylinderIF.H"
#include "TiltedCylinderIF.H"
#include "PlaneIF.H"
#include "IntersectionIF.H"
#include "TransformIF.H"
#include "EBAMRIO.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "GodunovGeom.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

int verbosity;

IntVect g_ivdebug(D_DECL(72,11,23));
int g_stenWidth = 1;
/************/
extern "C"
{
  void ebcfArea(const BaseEBCellFAB<Real>* const a_flux);
  void ebffArea(const BaseEBFaceFAB<Real>* const a_flux);
  void bivfPt(const BaseIVFAB<Real>* const a_flux);
  void ebcfPt(const BaseEBCellFAB<Real>* const a_flux);
  void ebffPt(const BaseEBFaceFAB<Real>* const a_flux);
  void stupidFunctionToFoolGDB()
  {
    bivfPt(NULL);
    ebcfPt(NULL);
    ebffPt(NULL);
    ebcfArea(NULL);
    ebffArea(NULL);
  }
}
/************/
void bivfPt(const BaseIVFAB<Real>* const a_flux)
{
  const BaseIVFAB<Real>& fab = *a_flux;
  int ncomp = fab.nComp();
  const EBGraph& ebgraph = fab.getEBGraph();
  const IntVectSet& ivs = fab.getIVS();
  if(ivs.contains(g_ivdebug))
    {
      Vector<VolIndex> vofs = ebgraph.getVoFs(g_ivdebug);
      for(int ivof = 0; ivof < vofs.size(); ivof++)
        {
          const VolIndex& vof = vofs[ivof];
          pout() << vof ;
          for(int ivar = 0; ivar < ncomp; ivar++ )
            {
              pout() << "  " << fab(vof, ivar);
            }
          pout() << endl;
        }
    }
  else
    {
      pout() << g_ivdebug << " not in the fab " << endl;
    }
}
/************/
void ebcfPt(const BaseEBCellFAB<Real>* const a_flux)
{
  const BaseEBCellFAB<Real>& fab = *a_flux;
  int ncomp = fab.nComp();
  const EBISBox& ebgraph = fab.getEBISBox();
  const Box& region = fab.getRegion();
  if(region.contains(g_ivdebug))
    {
      Vector<VolIndex> vofs = ebgraph.getVoFs(g_ivdebug);
      for(int ivof = 0; ivof < vofs.size(); ivof++)
        {
          const VolIndex& vof = vofs[ivof];
          pout() << vof ;
          for(int ivar = 0; ivar < ncomp; ivar++ )
            {
              pout() << "  " << fab(vof, ivar);
            }
          pout() << endl;
        }
    }
  else
    {
      pout() << g_ivdebug << " not in the fab " << endl;
    }
}

/************/
void ebcfArea(const BaseEBCellFAB<Real>* const a_flux)
{
  const BaseEBCellFAB<Real>& fab = *a_flux;
  int ncomp = fab.nComp();
  const EBISBox& ebisBox = fab.getEBISBox();
  const Box& region = fab.getRegion();
  if(region.contains(g_ivdebug))
    {
      Box barea(g_ivdebug, g_ivdebug);
      barea.grow(g_stenWidth);
      barea &= region;
      IntVectSet ivs(barea);
      for(VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          pout() << vof ;
          for(int ivar = 0; ivar < ncomp; ivar++ )
            {
              pout() << "  " << fab(vof, ivar);
            }
          pout() << endl;
        }
    }
  else
    {
      pout() << g_ivdebug << " not in the fab " << endl;
    }
}
/************/
void ebffArea(const BaseEBFaceFAB<Real>* const a_flux)
{
  const BaseEBFaceFAB<Real>& fab = *a_flux;
  int ncomp = fab.nComp();
  const EBISBox& ebisBox = fab.getEBISBox();
  const Box& region = fab.getCellRegion();
  int direction = fab.direction();
  if(region.contains(g_ivdebug))
    {
      Box barea(g_ivdebug, g_ivdebug);
      barea.grow(g_stenWidth);
      barea &= region;
      IntVectSet ivs(barea);
      for(FaceIterator faceit(ivs, ebisBox.getEBGraph(), direction,
                              FaceStop::SurroundingNoBoundary);
          faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          pout() << face.gridIndex(Side::Lo) << " "
                 << face.gridIndex(Side::Hi) << " "  ;
          for(int ivar = 0; ivar < ncomp; ivar++ )
            {
              pout() << "  " << fab(face, ivar);
            }
          pout() << endl;
        }
    }
  else
    {
      pout() << g_ivdebug << " not in the fab " << endl;
    }
}
/************/
void ebffPt(const BaseEBFaceFAB<Real>* const a_flux)
{
  const BaseEBFaceFAB<Real>& fab = *a_flux;
  int ncomp = fab.nComp();
  const EBISBox& ebgraph = fab.getEBISBox();
  const Box& region = fab.getRegion();
  int direction = fab.direction();
  if(region.contains(g_ivdebug))
    {
      Vector<VolIndex> vofs = ebgraph.getVoFs(g_ivdebug);
      for(int ivof = 0; ivof < vofs.size(); ivof++)
        {
          const VolIndex& vof = vofs[ivof];
          for(SideIterator sd; sd.ok(); ++sd)
            {
              Vector<FaceIndex> faces = ebgraph.getFaces(vof, direction, sd());
              for(int iface = 0; iface < faces.size(); iface++)
                {
                  const FaceIndex& face = faces[iface];
                  pout() << face.gridIndex(Side::Lo) << " "
                         << face.gridIndex(Side::Hi) << " "  ;
                  for(int ivar = 0; ivar < ncomp; ivar++ )
                    {
                      pout() << "  " << fab(face, ivar);
                    }
                  pout() << endl;
                }
            }
        }
    }
  else
    {
      pout() << g_ivdebug << " not in the fab " << endl;
    }
}
/************/
void
generateCalculatedData(LevelData<EBCellFAB>& a_calcData,
                       Real& a_finalTime,
                       DisjointBoxLayout& a_grids,
                       EBISLayout&        a_ebisl,
                       const Box&  a_domain,
                       const Real& a_dx,
                       const bool& a_thisCalcCoar)
{
  AMR amr;
  // read inputs
  ParmParse pp;

  Real gamma = 1.4;
  pp.get("gamma",gamma);


  int redistRad = 0;
  pp.get("redist_radius",redistRad);

  int nstop = 1;
  //only running one step
  Real stopTime = 1.0;

  bool is_periodic[SpaceDim];
  for (int i = 0; i < SpaceDim; ++i)
    {
      is_periodic[i] = false;
    }

  //single level calc
  int max_level = 0;
  // note this requires a ref_ratio to be defined for the
  // finest level (even though it will never be used)
  std::vector<int> ref_ratios(1,1);

  Real cfl = 0.8;
  pp.get("cfl",cfl);

  //only running one time step
  Real initialCFL = cfl;

  Real fixed_dt = -1;
  pp.get("fixed_dt",fixed_dt);
  if(a_thisCalcCoar && fixed_dt > 0)
    {
      fixed_dt *= 2.0;
    }

  //not relevant because only running one step
  Real max_dt_growth = 1.1;
  Real dt_tolerance_factor = 1.1;

  vector<Real>  length(SpaceDim);
  pp.getarr("domain_length",length, 0, SpaceDim);
  RealVect domainLength;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      domainLength[idir] = length[idir];
    }

  int ifourth, iflatten, iartvisc;
  pp.get("use_fourth_order_slopes", ifourth);
  pp.get("use_flattening"         , iflatten);
  pp.get("use_art_visc"           , iartvisc);
  bool useFourthOrderSlopes = (ifourth  ==1);
  bool useFlattening        = (iflatten ==1);
  bool useArtificialVisc    = (iartvisc ==1);

  ProblemDomain prob_domain(a_domain.smallEnd(),
                            a_domain.bigEnd(),
                            is_periodic);

  // create initial/boundary condition object

  if(verbosity > 1)
    pout() << "Modiano initial and boundary conditions" << endl;

  RealVect modianoAxis;
  // check to see if we are using the cylinder

  int whichGeom;
  pp.get("which_geom", whichGeom);
  if(whichGeom != 2)
    {
      if(verbosity > 1)
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
      RealVect channelNormal;
      vector<Real>  channelNormalVect(SpaceDim);
      pp.getarr("channel_normal",channelNormalVect, 0, SpaceDim);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          channelNormal[idir] = channelNormalVect[idir];
        }
      if(verbosity > 1)
        pout() << "uses channel slope for direction" << endl;
      modianoAxis = mibcGetNormalVector(channelNormal);;
      //debug
      //modianoAxis = BASISREALV(1);
    }
  Real waveAmp, waveWidth;
  pp.get("wave_amplitude", waveAmp);
  pp.get("wave_width", waveWidth);

  vector<Real> centerpp(SpaceDim,0.5);
  RealVect center;
  pp.getarr("wave_center",centerpp,0,SpaceDim);
  for (int i = 0; i < SpaceDim; i++)
    center[i] = centerpp[i];

  int idoNegWave;
  pp.get("use_negative_wave",idoNegWave);
  bool useNegativeWave = (idoNegWave==1);
  if(verbosity > 1)
    {
      if(useNegativeWave)
        {
          pout() << "u - c wave " << endl;
        }
      else
        {
          pout() << "u + c wave " << endl;
        }
    }

  int idoFreeStream;
  pp.get("free_stream_prob",idoFreeStream);
  bool doFreeStream = (idoFreeStream==1);
  if(verbosity > 1)
    {
      if(doFreeStream)
        {
          pout() << "free stream modiano bump " << endl;
        }
      else
        {
          pout() << "simple wave prob " << endl;
        }
    }

  ModianoIBCFactory bcfactory(gamma, waveAmp, waveWidth,
                              center, modianoAxis, doFreeStream, useNegativeWave);

  int uselim, usesmush;
  pp.get("use_limiting", uselim);
  bool useLimiting = (uselim==1);
  bool doRZCoords = false;
  bool hasSourceTerm = false;
  pp.get("use_smushing", usesmush);

  bool doSmushing = usesmush==1;
  //create patch integrator
  EBPatchPolytropicFactory patchGamma(&bcfactory,
                                      gamma,
                                      useFourthOrderSlopes,
                                      useFlattening,
                                      useArtificialVisc,
                                      useLimiting,
                                      doRZCoords);


  //dummies
  Real refineThresh = 0.7;
  int tagBufferSize = 1;
  int iusemassredist;
  pp.get("use_mass_redist", iusemassredist);
  bool useMassRedist    = (iusemassredist ==1);
  int amrverbosity;
  pp.get("amrverbosity",amrverbosity);

  EBAMRGodunovFactory amrg_fact(initialCFL,
                                cfl,
                                redistRad,
                                domainLength,
                                refineThresh,
                                tagBufferSize,
                                amrverbosity,
                                useMassRedist,
                                doSmushing,
                                doRZCoords,
                                hasSourceTerm,
                                &patchGamma);

  int max_grid_size = 32;
  pp.get("max_grid_size",max_grid_size);

  amr.define(max_level,ref_ratios,prob_domain,&amrg_fact);
  amr.maxGridSize(max_grid_size);

  if (fixed_dt > 0)
    {
      amr.fixedDt(fixed_dt);
    }

  // the hyperbolic codes use a grid buffer of 1
  if(a_thisCalcCoar)
    {
      std::string prefix("pltCoar");
      amr.plotPrefix(prefix);
    }
  else
    {
      std::string prefix("pltFine");
      amr.plotPrefix(prefix);
    }
  amr.gridBufferSize(1);
  amr.checkpointInterval(-1);
  int dofileout;
  pp.get("do_fileoutput", dofileout);
  if(dofileout == 1)
    {
      amr.plotInterval(1);
    }
  else
    {
      amr.plotInterval(-1);
    }
  amr.maxDtGrow(max_dt_growth);
  amr.dtToleranceFactor(dt_tolerance_factor);
  amr.verbosity(amrverbosity);

  // initiailze hierarchy of levels
  if(verbosity > 1)
    {
      pout() << "starting amr setup " << endl;
    }
  amr.setupForNewAMRRun();
  // run
  if(verbosity > 1)
    {
      pout() << "starting amr run " << endl;
    }
#if 1
  //goes too slow
  amr.run(stopTime,nstop);
  // output last pltfile and statistics
  if(verbosity > 1)
    {
      pout() << "starting amr conclude " << endl;
    }
  amr.conclude();

  Vector<AMRLevel*> amrLevelData = amr.getAMRLevels();
  CH_assert(amrLevelData.size() == 1);

  EBAMRGodunov* amrGodunovData = dynamic_cast<EBAMRGodunov*>(amrLevelData[0]);

  LevelData<EBCellFAB>& amrData = amrGodunovData->getStateNew();

  amrGodunovData->fillConsAndPrim(a_calcData);
  a_grids = amrData.disjointBoxLayout();
  a_ebisl = amrGodunovData->getEBISLayout();
  a_finalTime = amr.getCurrentTime();
#endif
}
/***************/
void
generateExactData(LevelData<EBCellFAB>& a_exactData,
                  const DisjointBoxLayout& a_grids,
                  const EBISLayout&        a_ebisl,
                  const Box& a_domain,
                  const Real& a_dx,
                  const Real& a_finalTime,
                  const int&  a_nvar)
{
  EBCellFactory ebcf(a_ebisl);
  a_exactData.define(a_grids, a_nvar, IntVect::Zero, ebcf);

  RealVect modianoAxis, center;
  Real waveAmp, waveWidth, gamma;
  int whichGeom;
  ParmParse pp;
  pp.get("verbosity", verbosity);
  pp.get("which_geom", whichGeom);
  if(whichGeom == 3)
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
      modianoAxis = cylinderAxis;
      if(verbosity > 1)
        {
          pout() << " tilted cylinder " << endl;
          pout() << "uses cylinder axis for direction" << endl;
          pout() << "axis = " << modianoAxis << endl;
        }
    }
  else
    {
      RealVect channelNormal;
      vector<Real>  channelNormalVect(SpaceDim);
      pp.getarr("channel_normal",channelNormalVect, 0, SpaceDim);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          channelNormal[idir] = channelNormalVect[idir];
        }
      modianoAxis = mibcGetNormalVector(channelNormal);;
      if(verbosity > 1)
        {
          pout() << "channel geom" << endl;
          pout() << "uses channel slope for direction" << endl;
          pout() << " channel norm = " << channelNormal << endl;
          pout() << " flow dir = " << modianoAxis << endl;
        }
    }
  pp.get("wave_amplitude", waveAmp);
  pp.get("wave_width", waveWidth);
  pp.get("gamma",gamma);

  vector<Real> centerpp(SpaceDim,0.5);
  pp.getarr("wave_center",centerpp,0,SpaceDim);
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      center[idir] = centerpp[idir];
    }

  int idoNegWave;
  pp.get("use_negative_wave",idoNegWave);
  bool useNegativeWave = (idoNegWave==1);

  int idoFreeStream;
  pp.get("free_stream_prob",idoFreeStream);
  bool doFreeStream = (idoFreeStream==1);
  if(verbosity > 1)
    {
      pout() << "wave amp = " << waveAmp << " wave width = " << waveWidth << endl;

      pout() << " gamma = " << gamma << " center = " << center << endl;

      if(doFreeStream)
        pout() << "free stream problem " << endl;
      else
        pout() << "full Modiano problem" << endl;

      if(useNegativeWave)
        pout() << "with negative wave" << endl;
      else
        pout() << "with positive wave" << endl;
    }

  ModianoIBC mibc(gamma, waveAmp, waveWidth,
                  center, modianoAxis, doFreeStream, useNegativeWave);

  bool is_periodic[SpaceDim];
  for (int i = 0; i < SpaceDim; ++i)
    {
      is_periodic[i] = false;
    }
  ProblemDomain probdomain(a_domain.smallEnd(), a_domain.bigEnd(), is_periodic);
  RealVect dxvect = a_dx*RealVect::Unit;
  mibc.define(probdomain, dxvect);
  mibc.setToExactConsAndPrim(a_exactData, a_ebisl, a_finalTime);
}
/***************/
Vector<string> getNames()
{
  Vector<string> retval;
  EBPatchPolytropic patchInt;
  Vector<string> consNames  = patchInt.stateNames();
  Vector<string> primNames  = patchInt.primNamesNoLog();
  retval = consNames;
  retval.append(primNames);
  return retval;
}
/***************/
void
generateError(LevelData<EBCellFAB>& a_error,
              DisjointBoxLayout&    a_grids,
              EBISLayout&           a_ebisl,
              const Box&            a_domain,
              const Real&           a_dx,
              const bool&           a_thisCalcCoar)
{
  LevelData<EBCellFAB> calcData;
  Real finalTime;
  generateCalculatedData(calcData,  finalTime, a_grids, a_ebisl,
                         a_domain, a_dx, a_thisCalcCoar);

#if 1
  //goes too slow
  LevelData<EBCellFAB> exactData;

  int nvar = calcData.nComp();
  generateExactData(exactData, a_grids, a_ebisl,
                    a_domain, a_dx, finalTime, nvar);

  EBCellFactory ebcf(a_ebisl);
  a_error.define(a_grids, nvar, IntVect::Zero, ebcf);

  Interval interv(0, nvar-1);
  //truncation error is calc-exact/dt
  for(DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& errorFAB = a_error[dit()];
      EBCellFAB& calcFAB = calcData[dit()];
      EBCellFAB& exacFAB = exactData[dit()];

      //error is exact-calculated data
      errorFAB.copy(exacFAB);
      errorFAB -= calcFAB;
      CH_assert(finalTime > 0);
      errorFAB /= finalTime;
      for(int icomp = 0; icomp < errorFAB.nComp(); icomp++)
        {
          errorFAB.setCoveredCellVal(0.0, icomp);
        }
    }

#ifdef CH_USE_HDF5

#if CH_SPACEDIM==2
  string fileFine("pltFineExact.2d.hdf5");
  string fileCoar("pltCoarExact.2d.hdf5");
#else
  string fileFine("pltFineExact.3d.hdf5");
  string fileCoar("pltCoarExact.3d.hdf5");
#endif

  string filename;
  if(a_thisCalcCoar)
    filename = fileCoar;
  else
    filename = fileFine;
  Vector<string> names  = getNames();

  bool replaceCovered = true;
  Vector<Real> coveredValues(nvar, 0.0);
  //values that don't matter in output file
  Real dx = 1.0;
  Real time   = 1.0;
  Real dt     = 1.0;
  /**/
  ParmParse pp;
  int dofileout;
  pp.get("do_fileoutput", dofileout);
  if(dofileout == 1)
    {
      writeEBHDF5(filename, a_grids, exactData, names,
                  a_domain, dx, dt, time,
                  replaceCovered, coveredValues);
    }
#endif
  /**/
#endif
}
/***************/
void
outputError(const LevelData<EBCellFAB>& a_errorFine,
            const EBISLayout& a_ebislFine,
            const DisjointBoxLayout& a_gridsFine,
            const Box& a_domainFine,
            const LevelData<EBCellFAB>& a_errorCoar,
            const EBISLayout& a_ebislCoar,
            const DisjointBoxLayout& a_gridsCoar,
            const Box& a_domainCoar)
{
#ifdef CH_USE_HDF5
#if CH_SPACEDIM==2
  string fileFine("pltFineError.2d.hdf5");
  string fileCoar("pltCoarError.2d.hdf5");
#else
  string fileFine("pltFineError.3d.hdf5");
  string fileCoar("pltCoarError.3d.hdf5");
#endif
  EBPatchPolytropic patchInt;
  Vector<string> names  = getNames();
  bool replaceCovered = true;
  int nvar = a_errorFine.nComp();
  Vector<Real> coveredValues(nvar, 0.0);
  //values that don't matter in output file
  Real dxFine = 1.0;
  Real dxCoar = 1.0;
  Real time   = 1.0;
  Real dt     = 1.0;

  ParmParse pp;
  int dofileout;
  pp.get("do_erroroutput", dofileout);
  if(dofileout == 1)
    {
      writeEBHDF5(fileFine, a_gridsFine, a_errorFine, names,
                  a_domainFine, dxFine, dt, time,
                  replaceCovered, coveredValues);

      writeEBHDF5(fileCoar, a_gridsCoar, a_errorCoar, names,
                  a_domainCoar, dxCoar, dt, time,
                  replaceCovered, coveredValues);
    }
#endif
}
/***************/
int iscript(int icomp, int inorm, int ncomp)
{
  return icomp + inorm*ncomp;
}
/***************/
void
compareError(const LevelData<EBCellFAB>& a_errorFine,
             const EBISLayout& a_ebislFine,
             const DisjointBoxLayout& a_gridsFine,
             const Box& a_domainFine,
             const LevelData<EBCellFAB>& a_errorCoar,
             const EBISLayout& a_ebislCoar,
             const DisjointBoxLayout& a_gridsCoar,
             const Box& a_domainCoar)
{
  EBPatchPolytropic patchInt;
  Vector<string> names  = getNames();
  EBCellFactory factFine(a_ebislFine);
  EBCellFactory factCoar(a_ebislCoar);
  IntVect ivghost = IntVect::Zero;
  const int ncomp = names.size();
  const int nnorm = 3;
  Real* normsCoar = new Real[ncomp*nnorm];
  Real* normsFine = new Real[ncomp*nnorm];
  Real* orders    = new Real[ncomp*nnorm];
  for(int icomp = 0; icomp < ncomp; icomp++)
    {
      orders[icomp] = 0.0;
      for(int inorm = 0; inorm < nnorm; inorm++)
        {
          normsCoar[iscript(icomp, inorm, ncomp)] = 0;
          normsFine[iscript(icomp, inorm, ncomp)] = 0;
        }
    }

  if(verbosity > 1)
    pout() << "==============================================" << endl;
  for(int comp = 0; comp < names.size(); comp++)
    {
      if(verbosity > 1)
        {
          pout() << "Comparing error in variable  " << names[comp] << endl;
          pout() << "==============================================" << endl;
        }
      for(int itype = 0; itype < 3; itype++)
        {
          EBNormType::NormMode normtype;
          if(itype == 0)
            {
              normtype = EBNormType::OverBoth;
              if(verbosity > 1)
                pout() << endl << "Using all uncovered cells." << endl;
            }
          else if(itype == 1)
            {
              normtype = EBNormType::OverOnlyRegular;
              if(verbosity > 1)
                pout() << endl << "Using only regular cells." << endl;
            }
          else
            {
              normtype = EBNormType::OverOnlyIrregular;
              if(verbosity > 1)
                pout() << endl << "Using only irregular cells." << endl;
            }
          for(int inorm = 0; inorm <= 2; inorm++)
            {

              if(inorm == 0)
                {
                  if(verbosity > 1)
                    pout() << endl << "Using max norm." << endl;
                }
              else
                {
                  if(verbosity > 1)
                    pout() << endl << "Using L-" << inorm << "norm." << endl;
                }
              VolIndex& vofmax = EBArith::getVoFMax();
              Real&     valmax = EBArith::getValMax();

              bool dovmax = (comp==0) && (inorm==0) && (normtype == EBNormType::OverBoth);
              if(dovmax)
                {
                  valmax= 0.0;
                }
              Real coarnorm = EBArith::norm(a_errorCoar,
                                            a_gridsCoar, a_ebislCoar,
                                            comp, inorm, normtype);
              if(dovmax)
                {
                  pout() << vofmax << "is coar vof with greatest density error = " << valmax << endl;
                  valmax = 0.0;
                }

              Real finenorm = EBArith::norm(a_errorFine,
                                            a_gridsFine, a_ebislFine,
                                            comp, inorm, normtype);
              if(dovmax)
                {
                  pout() << vofmax << "is fine vof with greatest density error = " << valmax << endl;
                  valmax = 0.0;

                }
              if(verbosity > 1)
                {
                  pout() << "Coarse Error Norm = " << coarnorm << endl;
                  pout() << "Fine   Error Norm = " << finenorm << endl;
                }
              if(itype == 0)
                {
                  normsCoar[iscript(comp,inorm,ncomp)] = coarnorm;
                  normsFine[iscript(comp,inorm,ncomp)] = finenorm;
                }
              if((Abs(finenorm) > 1.0e-8) && (Abs(coarnorm) > 1.0e-8))
                {
                  Real order = log(coarnorm/finenorm)/log(2.0);
                  if(verbosity > 1)
                    pout() << "Order of scheme = " << order << endl;
                  if(itype == 0)
                    {
                      orders[iscript(comp,inorm,ncomp)] = order;
                    }
                }
            }
        }
      if(verbosity > 1)
        pout() << "==============================================" << endl;
    }
  //output in latex format to be safe

  int nfine = a_domainFine.size(0);
  if(verbosity > 0)
    {
      pout() << setw(12)
             << setprecision(6)
             << setiosflags(ios::showpoint)
             << setiosflags(ios::scientific) ;
      for (int inorm = 0; inorm <= 2; inorm++)
        {
          pout() << "\\begin{table}[p]" << endl;
          pout() << "\\begin{center}" << endl;
          pout() << "\\begin{tabular}{|c|c|c|c|} \\hline" << endl;
          pout() << "Variable & Coarse Error & Fine Error & Order\\\\" << endl;;
          pout() << "\\hline \\hline " << endl;
          for(int icomp = 0; icomp < ncomp; icomp++)
            {
              int iindex = iscript(icomp,inorm,ncomp);
              pout() << names[icomp] << " &    \t "
                     << setw(12)
                     << setprecision(6)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << normsCoar[iindex]  << " & "
                     << setw(12)
                     << setprecision(6)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << normsFine[iindex] << " & "
                     << setw(12)
                     << setprecision(6)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << orders[iindex];
              pout() << " \\\\ " << endl <<   "\\hline"  <<  endl;
            }
          pout() << "\\end{tabular}" << endl;
          pout() << "\\end{center}" << endl;
          pout() << "\\caption{Truncation error convergence rates using L-" << inorm << " norm. " << endl;
          pout() << "$h_f = \\frac{1}{" << nfine << "}$ and $h_c = 2 h_f$, $D = " << SpaceDim << "$ }"  << endl;
          pout() << "\\end{table}" << endl;
          pout() << endl << endl;
        }
    }
  delete[] normsCoar;
  delete[] normsFine;
  delete[] orders   ;

}
/***************/
void
multiplyByKappa(LevelData<EBCellFAB>& a_error,
                const EBISLayout& a_ebisl,
                const DisjointBoxLayout& a_grids)
{
 for(DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = a_grids.get(dit());
      const EBISBox& ebisBox = a_ebisl[dit()];
      IntVectSet ivs(box);
      for(VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real volFrac = ebisBox.volFrac(vof);
          for(int ivar = 0; ivar < a_error.nComp(); ivar++)
            {
              a_error[dit()](vof, ivar) *= volFrac;
            }
        }
    }
}
/***************/
void
compareKEMax(LevelData<EBCellFAB>& a_errorFine,
             const EBISLayout& a_ebislFine,
             const DisjointBoxLayout& a_gridsFine,
             const Box& a_domainFine,
             LevelData<EBCellFAB>& a_errorCoar,
             const EBISLayout& a_ebislCoar,
             const DisjointBoxLayout& a_gridsCoar,
             const Box& a_domainCoar)
{
  //first multiply by kappa
  multiplyByKappa(a_errorFine, a_ebislFine, a_gridsFine);
  multiplyByKappa(a_errorCoar, a_ebislCoar, a_gridsCoar);

  EBPatchPolytropic patchInt;
  Vector<string> names  = getNames();
  EBCellFactory factFine(a_ebislFine);
  EBCellFactory factCoar(a_ebislCoar);
  IntVect ivghost = IntVect::Zero;
  const int ncomp = names.size();
  Real* normsCoar = new Real[ncomp];
  Real* normsFine = new Real[ncomp];
  Real* orders    = new Real[ncomp];
  for(int icomp = 0; icomp < ncomp; icomp++)
    {
      orders[icomp] = 0.0;
      normsCoar[icomp] = 0;
      normsFine[icomp] = 0;
    }

  if(verbosity > 1)
    pout() << "==============================================" << endl;
  for(int comp = 0; comp < names.size(); comp++)
    {
      if(verbosity > 1)
        {
          pout() << "Comparing error in variable  " << names[comp] << " multiplied by vol frac" << endl;
          pout() << "==============================================" << endl;
        }
      EBNormType::NormMode normtype = EBNormType::OverBoth;
      if(verbosity > 1)
        pout() << endl << "Using all uncovered cells and max norm." << endl;

      int inorm = 0;
      Real coarnorm = EBArith::norm(a_errorCoar,
                                    a_gridsCoar, a_ebislCoar,
                                    comp, inorm, normtype);

      Real finenorm = EBArith::norm(a_errorFine,
                                    a_gridsFine, a_ebislFine,
                                    comp, inorm, normtype);
      if(verbosity > 1)
        {
          pout() << "Coarse Error Norm = " << coarnorm << endl;
          pout() << "Fine   Error Norm = " << finenorm << endl;
        }

      normsCoar[comp] = coarnorm;
      normsFine[comp] = finenorm;

      if((Abs(finenorm) > 1.0e-8) && (Abs(coarnorm) > 1.0e-8))
        {
          Real order = log(coarnorm/finenorm)/log(2.0);
          if(verbosity > 1)
            pout() << "Order of scheme = " << order << endl;
          orders[comp] = order;
        }
    }
  if(verbosity > 0)
    {
      pout() << "==============================================" << endl;
      int nfine = a_domainFine.size(0);
      pout() << setw(12)
             << setprecision(6)
             << setiosflags(ios::showpoint)
             << setiosflags(ios::scientific) ;
      pout() << "\\begin{table}[p]" << endl;
      pout() << "\\begin{center}" << endl;
      pout() << "\\begin{tabular}{|c|c|c|c|} \\hline" << endl;
      pout() << "Variable & Coarse Error & Fine Error & Order\\\\" << endl;;
      pout() << "\\hline \\hline " << endl;
      for(int icomp = 0; icomp < ncomp; icomp++)
        {
          int iindex = icomp;
          pout() << names[icomp] << " &    \t "
                 << setw(12)
                 << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << normsCoar[iindex]  << " & "
                 << setw(12)
                 << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << normsFine[iindex] << " & "
                 << setw(12)
                 << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << orders[iindex];
          pout() << " \\\\ " << endl <<   "\\hline"  <<  endl;
        }
      pout() << "\\end{tabular}" << endl;
      pout() << "\\end{center}" << endl;
      pout() << "\\caption{Truncation error weighted with volume fraction convergence rates using L-0  norm. " << endl;
      pout() << "$h_f = \\frac{1}{" << nfine << "}$ and $h_c = 2 h_f$, $D = " << SpaceDim << "$ }"  << endl;
      pout() << "\\end{table}" << endl;
      pout() << endl << endl;
    }
  delete[] normsCoar;
  delete[] normsFine;
  delete[] orders   ;

}
/***************/
void
getFinestDomain(Box& a_domain,
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
             const Real& a_dx)
{
  //parse input file.  single level
  ParmParse pp;
  RealVect origin = RealVect::Zero;
  Real channelRadius;
  pp.get("channel_radius", channelRadius);

  int whichGeom;
  pp.get("which_geom", whichGeom);
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  int biggridsize;
  //  int  fineLength = a_domain.size(0);
  pp.get("max_grid_size", biggridsize);
  if(whichGeom == 0)
    {
      if(verbosity > 1)
        pout() << "using all reg geom" << endl;
      AllRegularService allreg;
      ebisPtr->define(a_domain, origin, a_dx, allreg, biggridsize, 0);
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

      if(verbosity > 1)
        pout() << "using a ramp" << endl;

      RealVect vectDx = RealVect::Unit;
      vectDx *= a_dx;

      GeometryShop workshop(ramp,0,vectDx);
      ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize, 0);
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

      if(verbosity > 1)
        pout() << "using a channel" << endl;

      RealVect vectDx = RealVect::Unit;
      vectDx *= a_dx;

      GeometryShop workshop(channel,0,vectDx);
      ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize, 0);
    }
  else if(whichGeom == 3)
    {
      vector<Real>  cylinderAxisVect(SpaceDim);
      pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
      RealVect cylinderAxis;
      RealVect cylinderOrigin = RealVect::Zero;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          cylinderAxis[idir]   = cylinderAxisVect[idir];
        }
      Real sum;
      PolyGeom::unifyVector(cylinderAxis, sum);

      int cylinderType;
      pp.get("cylinder_type", cylinderType);
      bool negativeInside = true;
      if(cylinderType == 0)
        {
          if(verbosity > 1)
            pout() << "using a tilted circular cylinder" << endl;

          TiltedCylinderIF tunnel(channelRadius, cylinderAxis, cylinderOrigin, negativeInside);

          RealVect vectDx = RealVect::Unit;
          vectDx *= a_dx;

          GeometryShop workshop(tunnel,0,vectDx);
          ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize, 0);
        }
      else if(cylinderType == 1)
        {
          if(verbosity > 1)
            pout() << "using a tilted square  cylinder" << endl;

          //make square cylinder around x axis
          SquareCylinderIF sqcyl(channelRadius, negativeInside);
          TransformIF transformif(sqcyl);
          RealVect xAxis= BASISREALV(0);
          transformif.rotate(xAxis, cylinderAxis);
          transformif.translate(cylinderOrigin);

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
}
/***/
void
scatterPlot(const LevelData<EBCellFAB>&  a_solnError,
            const EBISLayout&            a_ebisl,
            const DisjointBoxLayout&     a_grids,
            const Box&                   a_domain,
            const Real&                  a_dx,
            const bool&                  a_thisCalcCoar,
            const int&                   a_ivarError,
            const string&                a_solnName,
            IrregErrorFcn                a_fcn,
            const string&                a_geomName,
            const string&                a_filename,
            ofstream&                    a_ossh,
            bool                         a_multByKappa)
{
  //make data of x vs y for all
  Real maxSoln = 0.0;
  Real maxGeom = 0.0;
  ParmParse pp;
  ofstream  osdat(a_filename.c_str());
  Real tolerance = 1.0e-2;

  vector<Real>  cylinderAxisVect(SpaceDim);
  pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
  RealVect cylinderAxis;
  //find maxes
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      cylinderAxis[idir] = cylinderAxisVect[idir];
    }

  for(DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = a_ebisl[dit()];
      const Box& box = a_grids.get(dit());
      IntVectSet ivs = ebisBox.getIrregIVS(box);
      BaseIVFAB<Real> geomError(ivs, ebisBox.getEBGraph(), 1);
      //generate geometry error
      a_fcn(geomError, ivs, ebisBox, cylinderAxis, a_dx);
      for(VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real geomerr = Abs(geomError(vof, 0));
          Real solnerr = Abs(a_solnError[dit()](vof, a_ivarError));
          Real volfrac = ebisBox.volFrac(vof);
          if(a_multByKappa)
            {
              solnerr *= volfrac;
            }
          if(geomerr > maxGeom) maxGeom = geomerr;
          if(solnerr > maxSoln) maxSoln = solnerr;
          //no point in outputting  a ton of zeros
          if(solnerr > tolerance*maxSoln)
            {
              osdat << geomerr << "   " << solnerr << endl;
            }
        }
    }
  osdat.close();
  char gnuplotname[200];
  char postscriptname[200];
  sprintf(gnuplotname, "%s.gnuplot",a_filename.c_str());
  sprintf(postscriptname,   "%s.ps",a_filename.c_str());
  ofstream osgnu(gnuplotname);

  maxSoln *= 1.05;
  maxGeom *= 1.05;

  osgnu << "set xzeroaxis"  << endl;
  osgnu << "set yzeroaxis"  << endl;
  osgnu << "set xrange [0:" << maxGeom << "]" << endl;
  osgnu << "set yrange [0:" << maxSoln << "]" << endl;
  osgnu << "set xlabel \"" << a_geomName << "\"" << endl;
  osgnu << "set ylabel \"" << a_solnName << "\"" << endl;
  osgnu << "set nokey " << endl;
  osgnu << "set term post eps " << endl;
  osgnu << "plot \"" << a_filename << "\"" << endl;

  osgnu.close();

  a_ossh << "gnuplot " << gnuplotname << " > " << postscriptname << endl;
  a_ossh << "set s=$status ; if ( $s != 0 ) exit $s " << endl;
}
/***/
void
uberScatterPlot(const LevelData<EBCellFAB>&  a_error,
                const EBISLayout&            a_ebisl,
                const DisjointBoxLayout&     a_grids,
                const Box&                   a_domain,
                const Real&                  a_dx,
                const bool&                  a_thisCalcCoar,
                ofstream&                    a_ossh)
{
  string geomName, errorName, filename;
  //if we cycle through vars, need to modify errorName, filename
  //lets just deal win conserved for now
  EBPatchPolytropic patchInt;
  int imultbykappa;
  ParmParse pp;
  pp.get("multbyvolfrac",imultbykappa);
  bool bmult = (imultbykappa==1);
  Vector<string> names  = getNames();
  CH_assert(names.size() <= a_error.nComp());
  char corf[10];
  if(a_thisCalcCoar)
    {
      sprintf(corf, "Coar");
    }
  else
    {
      sprintf(corf, "Fine");
    }
  for(int ivar = 0; ivar < names.size(); ivar++)
    {
      char scratch[200];
      if(bmult)
        {
          sprintf(scratch, "kappa(%s error)", names[ivar].c_str());
        }
      else
        {
          sprintf(scratch, "%s error", names[ivar].c_str());
        }
      errorName = string(scratch);

      sprintf(scratch, "%sErrorVsNormMinuTNorm%s.dat", names[ivar].c_str(), corf);
      filename = string(scratch);
      geomName = string("|norm - true norm|");
      scatterPlot(a_error, a_ebisl, a_grids, a_domain, a_dx, a_thisCalcCoar,
                  ivar, errorName, getNormMinTrueNorm, geomName, filename, a_ossh,bmult);

      sprintf(scratch, "%sErrorVsNormDotAxis%s.dat", names[ivar].c_str(), corf);
      filename = string(scratch);
      geomName = string("normal dot cylinder axis");
      scatterPlot(a_error, a_ebisl, a_grids, a_domain, a_dx, a_thisCalcCoar,
                  ivar, errorName, getNormalDotAxis, geomName, filename, a_ossh,bmult);


      sprintf(scratch, "%sErrorVsNormDotTrueNormM1%s.dat", names[ivar].c_str(), corf);
      filename = string(scratch);
      geomName = string("normal dot true normal - 1");
      scatterPlot(a_error, a_ebisl, a_grids, a_domain, a_dx, a_thisCalcCoar,
                  ivar, errorName, getNormalDotTrueNormM1, geomName, filename, a_ossh,bmult);

      sprintf(scratch, "%sErrorVsVolFrac%s.dat", names[ivar].c_str(), corf);
      filename = string(scratch);
      geomName = string("volume fraction");
      scatterPlot(a_error, a_ebisl, a_grids, a_domain, a_dx, a_thisCalcCoar,
                  ivar, errorName, getVolFrac, geomName, filename, a_ossh,bmult);

#if CH_SPACEDIM==3
      sprintf(scratch, "%sErrorVsNormDotCross%s.dat", names[ivar].c_str(), corf);
      filename = string(scratch);
      geomName = string("normal dot (axis cross truenorm)");
      scatterPlot(a_error, a_ebisl, a_grids, a_domain, a_dx, a_thisCalcCoar,
                  ivar, errorName, getNormalDotCrossVec, geomName, filename, a_ossh,bmult);
#endif
    }
}

/************/
int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
  {
    MPI_Barrier(Chombo_MPI::comm);
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

    Box domainFine, domainCoar;
    Real dxFine, dxCoar;
    LevelData<EBCellFAB>     errorFine, errorCoar;
    DisjointBoxLayout gridsFine, gridsCoar;
    EBISLayout        ebislFine, ebislCoar;

    //    UnfreedMemory();
    getFinestDomain(domainFine, dxFine);
    domainCoar = coarsen(domainFine, 2);
    dxCoar = 2.0*dxFine;

    pout() << "starting fine geom gen" << endl;
    makeGeometry(domainFine, dxFine);

    pout() << "generating fine error" << endl;
    generateError(errorFine, gridsFine, ebislFine, domainFine, dxFine, false);
    //    UnfreedMemory();


    pout() << "starting coar geom gen" << endl;
    makeGeometry(domainCoar, dxCoar);
    pout() << "generating coarse error" << endl;
    generateError(errorCoar, gridsCoar, ebislCoar, domainCoar, dxCoar, true);


    /* goes too slow. */
    //output both errors to file
    outputError(errorFine, ebislFine, gridsFine, domainFine,
                errorCoar, ebislCoar, gridsCoar, domainCoar);

    //compute convergence rates and output to stdout
    compareError(errorFine, ebislFine, gridsFine, domainFine,
                 errorCoar, ebislCoar, gridsCoar, domainCoar);

#ifndef CH_MPI
    int scatterout;
    pp.get("do_scatterplots", scatterout);
    if(scatterout == 1)
      {
        pout() << "scatter plotting truncation error vs geometric error" << endl;

        ofstream ossh("truncerr.sh");
        ossh << "#!/bin/csh" << endl;
        uberScatterPlot(errorCoar, ebislCoar, gridsCoar, domainCoar, dxCoar, true, ossh);
        uberScatterPlot(errorFine, ebislFine, gridsFine, domainFine, dxFine, false, ossh);
        ossh.close();
      }
#endif

    //error gets multiplied by kappa here.
    compareKEMax(errorFine, ebislFine, gridsFine, domainFine,
                 errorCoar, ebislCoar, gridsCoar, domainCoar);
#ifdef CH_MPI
  }
  MPI_Finalize();
#endif
}
/***************/
