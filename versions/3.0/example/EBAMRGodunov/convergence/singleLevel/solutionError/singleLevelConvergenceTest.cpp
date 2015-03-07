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

#include "ParmParse.H"
#include "CH_HDF5.H"
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
#include "EBAMRIO.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "AllRegularService.H"
#include "TiltedCylinderIF.H"
#include "GodunovGeom.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#if CH_SPACEDIM==2
IntVect ivdebug(3, 5);
#elif CH_SPACEDIM==3
IntVect ivdebug(1, 1, 0);
#endif
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
generateExactData(LevelData<EBCellFAB>& a_exactData,
                  const DisjointBoxLayout& a_grids,
                  const EBISLayout&        a_ebisl,
                  const Box& a_domain,
                  const Real& a_dx,
                  const Real& a_finalTime,
                  const Real& a_dt,
                  const int&  a_nvar)
{
  EBCellFactory ebcf(a_ebisl);
  a_exactData.define(a_grids, a_nvar, IntVect::Zero, ebcf);

  RealVect modianoAxis, center;
  Real waveAmp, waveWidth, gamma;

  ParmParse pp;
  pp.get("wave_amplitude", waveAmp);
  pp.get("wave_width", waveWidth);
  pp.get("gamma",gamma);
  int idoNegWave;
  pp.get("use_negative_wave",idoNegWave);
  bool useNegativeWave = (idoNegWave==1);

  vector<Real> centerpp(SpaceDim,0.5);
  pp.getarr("wave_center",centerpp,0,SpaceDim);

  for(int idir = 0; idir < SpaceDim; idir++)
    {
      center[idir] = centerpp[idir];
    }

  int testverbosity;
  pp.get("testverbosity", testverbosity);
  int whichGeom;
  pp.get("which_geom", whichGeom);
  if(whichGeom == 3)
    {
      if(testverbosity > 1)
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
      if(testverbosity > 1)
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


  int idoFreeStream;
  pp.get("free_stream_prob",idoFreeStream);
  bool doFreeStream = (idoFreeStream==1);

  ModianoIBC mibc(gamma, waveAmp, waveWidth,
                  center, modianoAxis, doFreeStream, useNegativeWave);

  bool is_periodic[SpaceDim];
  for (int i = 0; i < SpaceDim; ++i)
    {
      is_periodic[i] = false;
    }
  ProblemDomain probdomain(a_domain.smallEnd(), a_domain.bigEnd(), is_periodic);
  mibc.define(probdomain, a_dx);
  mibc.setToExactConsAndPrim(a_exactData, a_ebisl, a_finalTime);
}
/************/
/************/
void
generateCalculatedData(LevelData<EBCellFAB>& a_calcData,
                       Real& a_dt,
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

  int testverbosity;
  pp.get("testverbosity", testverbosity);
  int verbosity;
  pp.get("verbosity",verbosity);
 CH_assert(verbosity >= 0);

  Real gamma = 1.4;
  pp.get("gamma",gamma);


  int redistRad = 0;
  pp.get("redist_radius",redistRad);

  int nstop = 0;
  pp.get("max_step",nstop);
  if(a_thisCalcCoar)
    {
      nstop /= 2;
    }

  Real stopTime = 0.0;
  pp.get("max_time",stopTime);

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

  Real initialCFL = 0.1;
  pp.get("initial_cfl",initialCFL);

  Real fixed_dt = -1;
  pp.get("fixed_dt",fixed_dt);
  if(a_thisCalcCoar && fixed_dt > 0)
    {
      fixed_dt *= 2.0;
    }

  Real max_dt_growth = 1.1;
  pp.get("max_dt_growth",max_dt_growth);

  Real dt_tolerance_factor = 1.1;
  pp.get("dt_tolerance_factor",dt_tolerance_factor);
  Real domainLength;
  pp.get("domain_length",domainLength);

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
  if(testverbosity > 1)
    pout() << "Modiano initial and boundary conditions" << endl;

  RealVect modianoAxis;

  Real waveAmp, waveWidth;
  pp.get("wave_amplitude", waveAmp);
  pp.get("wave_width", waveWidth);
  int whichGeom;
  pp.get("which_geom", whichGeom);
  if(whichGeom == 3)
    {
      if(testverbosity > 1)
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
      if(testverbosity > 1)
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

  vector<Real> centerpp(SpaceDim,0.5);
  RealVect center;
  pp.getarr("wave_center",centerpp,0,SpaceDim);
  for (int i = 0; i < SpaceDim; i++)
    center[i] = centerpp[i];


  int idoNegWave;
  pp.get("use_negative_wave",idoNegWave);
  bool useNegativeWave = (idoNegWave==1);
  if(useNegativeWave)
    {
      if(testverbosity > 1)
        pout() << "u - c wave " << endl;
    }
  else
    {
      if(testverbosity > 1)
        pout() << "u + c wave " << endl;
    }

  int idoFreeStream;
  pp.get("free_stream_prob",idoFreeStream);
  bool doFreeStream = (idoFreeStream==1);
  if(doFreeStream)
    {
      if(testverbosity > 1)
        pout() << "free stream modiano bump " << endl;
    }
  else
    {
      if(testverbosity > 1)
        pout() << "simple wave prob " << endl;
    }

  ModianoIBCFactory bcfactory(gamma, waveAmp, waveWidth,
                              center, modianoAxis, doFreeStream, useNegativeWave);

  int uselim;
  pp.get("use_limiting", uselim);
  bool useLimiting = (uselim==1);
  bool doRZCoords = false;
  bool hasSourceTerm = false;
  bool doSmushing = true;
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

  EBAMRGodunovFactory amrg_fact(initialCFL,
                                cfl,
                                redistRad,
                                domainLength*RealVect::Unit,
                                refineThresh,
                                tagBufferSize,
                                verbosity,
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
  int plotInterval;
  pp.get("plot_interval", plotInterval);
  amr.checkpointInterval(-1);
  amr.plotInterval(plotInterval);
  amr.maxDtGrow(max_dt_growth);
  amr.dtToleranceFactor(dt_tolerance_factor);
  amr.verbosity(verbosity);

  // initialize hierarchy of levels
  amr.setupForNewAMRRun();
  // run
  amr.run(stopTime,nstop);
  // output last pltfile and stapptistics
  amr.conclude();

  Vector<AMRLevel*> amrLevelData = amr.getAMRLevels();
 CH_assert(amrLevelData.size() == 1);

  EBAMRGodunov* amrGodunovData = dynamic_cast<EBAMRGodunov*>(amrLevelData[0]);
  LevelData<EBCellFAB>& amrData = amrGodunovData->getStateNew();

  a_grids = amrData.disjointBoxLayout();
  a_finalTime = amr.getCurrentTime();
  a_dt = amrGodunovData->getDt();
  a_ebisl = amrGodunovData->getEBISLayout();
  amrGodunovData->fillConsAndPrim(a_calcData);
}
/***************/
/***************/
void
generateError(LevelData<EBCellFAB>& a_error,
              DisjointBoxLayout&    a_grids,
              EBISLayout&           a_ebisl,
              const Box&            a_domain,
              const Real&           a_dx,
              const bool&           a_thisCalcCoar)
{
  LevelData<EBCellFAB> calc;
  LevelData<EBCellFAB> exac;

  Real finalTime, dt;
  generateCalculatedData(calc, dt,
                         finalTime, a_grids, a_ebisl,
                         a_domain, a_dx, a_thisCalcCoar);


  EBPatchPolytropic patchInt;
  bool is_periodic[SpaceDim];
  for (int i = 0; i < SpaceDim; ++i)
    {
      is_periodic[i] = false;
    }
  ProblemDomain probdomain(a_domain.smallEnd(), a_domain.bigEnd(), is_periodic);
  patchInt.define(probdomain, a_dx);

  int nCons = patchInt.numConserved();
  int nPrim = patchInt.numPrimitives();
  int consAndPrim = nCons+nPrim;
  generateExactData( exac, a_grids, a_ebisl,
                    a_domain, a_dx, finalTime, dt, consAndPrim);


  EBCellFactory ebcf(a_ebisl);

  a_error.define(a_grids, consAndPrim, IntVect::Zero, ebcf);

  Interval interv(0, consAndPrim-1);
  for(DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = a_grids.get(dit());
      //error is exact-calculated data
      a_error[dit()].copy(box, interv, box, exac[dit()], interv);
      a_error[dit()] -= calc[dit()];

      for(int icomp = 0; icomp < a_error[dit()].nComp(); icomp++)
        {
          a_error[dit()].setCoveredCellVal(0.0, icomp);
        }

    }

}

/***************/
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
  bool rCov = true;
  int nCons = patchInt.numConserved();
  int nPrim = patchInt.numPrimitives();
  int consAndPrim = nCons+nPrim;
  Vector<Real> cVal(consAndPrim,-10.0);
  //values that don't matter in output file
  Real dxFine = 1.0;
  Real dxCoar = 2.0;
  Real time   = 1.0;
  Real dt     = 1.0;

  writeEBHDF5(fileFine, a_gridsFine, a_errorFine, names,
              a_domainFine, dxFine, dt, time,
              rCov, cVal);

  writeEBHDF5(fileCoar, a_gridsCoar, a_errorCoar, names,
              a_domainCoar, dxCoar, dt, time,
              rCov, cVal);
#endif
}
/***************/
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

  EBCellFactory factFine(a_ebislFine);
  EBCellFactory factCoar(a_ebislCoar);
  IntVect ivghost = IntVect::Zero;
  ParmParse pp;
  int testverbosity;
  pp.get("testverbosity", testverbosity);
  if(testverbosity > 1)
    pout() << "==============================================" << endl;
  for(int comp = 0; comp < names.size(); comp++)
    {
      if(testverbosity > 1)
        pout() << "Comparing error in variable  " << names[comp] << endl;
      if(testverbosity > 1)
        pout() << "==============================================" << endl;
      for(int itype = 0; itype < 3; itype++)
        {
          EBNormType::NormMode normtype;
          if(itype == 0)
            {
              normtype = EBNormType::OverBoth;
              if(testverbosity > 1)
                pout() << endl << "Using all uncovered cells." << endl  ;
            }
          else if(itype == 1)
            {
              normtype = EBNormType::OverOnlyRegular;
              if(testverbosity > 1)
                pout() << endl << "Using only regular cells." << endl ;
            }
          else
            {
              normtype = EBNormType::OverOnlyIrregular;
              if(testverbosity > 1)
                pout() << endl << "Using only irregular cells." << endl;
            }

          for(int inorm = 0; inorm <= 2; inorm++)
            {

              if(inorm == 0)
                {
                  if(testverbosity > 1)
                    pout() << endl << "Using max norm." << endl;
                }
              else
                {
                  if(testverbosity > 1)
                    pout() << endl << "Using L-" << inorm << "norm." << endl;
                }
              Real coarnorm = EBArith::norm(a_errorCoar,
                                            a_gridsCoar, a_ebislCoar,
                                            comp, inorm, normtype);
              Real finenorm = EBArith::norm(a_errorFine,
                                            a_gridsFine, a_ebislFine,
                                            comp, inorm, normtype);
              if(testverbosity > 1)
                pout() << "Coarse Error Norm = " << coarnorm << endl;
              if(testverbosity > 1)
                pout() << "Fine   Error Norm = " << finenorm << endl;
              if(itype == 0)
                {
                  normsCoar[iscript(comp,inorm,ncomp)] = coarnorm;
                  normsFine[iscript(comp,inorm,ncomp)] = finenorm;
                }
              if((Abs(finenorm) > 1.0e-12) && (Abs(coarnorm) > 1.0e-12))
                {
                  Real order = log(Abs(coarnorm/finenorm))/log(2.0);
                  if(testverbosity > 1)
                    pout() << "Order of scheme = " << order << endl;
                  if(itype == 0)
                    {
                      orders[iscript(comp,inorm,ncomp)] = order;
                    }
                }
            }
        }
      if(testverbosity > 1)
        pout() << "==============================================" << endl ;;
    }


  //output in latex format to be safe
  int nfine = a_domainFine.size(0);
  if(testverbosity > 0)
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
          pout() << "\\caption{Solution error convergence rates using L-" << inorm << " norm. " << endl;
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
/***************/
void
makeFinestDomain(Box& a_domain,
                 Real& a_dx)
{
  //parse input file.  single level
  ParmParse pp;
  Vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);

  int testverbosity;
  pp.get("testverbosity", testverbosity);
 CH_assert(n_cell.size() == SpaceDim);
  IntVect lo = IntVect::Zero;
  IntVect hi;
  for(int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if(n_cell[ivec] <= 0)
        {
          if(testverbosity > 1)
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
            ofstream&                    a_ossh)
{
  //make data of x vs y for all
  Real maxSoln = 0.0;
  Real maxGeom = 0.0;
  ParmParse pp;
  vector<Real>  cylinderAxisVect(SpaceDim);
  pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
  RealVect cylinderAxis;
  //find maxes
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      cylinderAxis[idir] = cylinderAxisVect[idir];
    }
  ofstream  osdat(a_filename.c_str());
  Real tolerance = 1.0e-3;
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
          if(geomerr > maxGeom) maxGeom = geomerr;
          if(solnerr > maxSoln) maxSoln = solnerr;
        }
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
      sprintf(scratch, "%s error", names[ivar].c_str());
      errorName = string(scratch);

      sprintf(scratch, "%sErrorVsNormMinuTNorm%s.dat", names[ivar].c_str(), corf);
      filename = string(scratch);
      geomName = string("|norm - true norm|");
      scatterPlot(a_error, a_ebisl, a_grids, a_domain, a_dx, a_thisCalcCoar,
                  ivar, errorName, getNormMinTrueNorm, geomName, filename, a_ossh);

      sprintf(scratch, "%sErrorVsNormDotAxis%s.dat", names[ivar].c_str(), corf);
      filename = string(scratch);
      geomName = string("normal dot cylinder axis");
      scatterPlot(a_error, a_ebisl, a_grids, a_domain, a_dx, a_thisCalcCoar,
                  ivar, errorName, getNormalDotAxis, geomName, filename, a_ossh);


      sprintf(scratch, "%sErrorVsNormDotTrueNormM1%s.dat", names[ivar].c_str(), corf);
      filename = string(scratch);
      geomName = string("normal dot true normal - 1");
      scatterPlot(a_error, a_ebisl, a_grids, a_domain, a_dx, a_thisCalcCoar,
                  ivar, errorName, getNormalDotTrueNormM1, geomName, filename, a_ossh);

#if CH_SPACEDIM==3
      sprintf(scratch, "%sErrorVsNormDotCross%s.dat", names[ivar].c_str(), corf);
      filename = string(scratch);
      geomName = string("normal dot (axis cross truenorm)");
      scatterPlot(a_error, a_ebisl, a_grids, a_domain, a_dx, a_thisCalcCoar,
                  ivar, errorName, getNormalDotCrossVec, geomName, filename, a_ossh);
#endif
    }
}

/***/
void
makeGeometry(const Box& a_domain,
             const Real& a_dx)
{

  RealVect origin = RealVect::Zero;
  ParmParse pp;
  Real channelRadius;
  pp.get("channel_radius", channelRadius);

  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  int biggridsize = 2048;
  char ebfilename[80];
  int  fineLength = a_domain.size(0);
  RealVect cylinderOrigin = RealVect::Zero;

  pp.get("max_grid_size", biggridsize);
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

  bool negativeInside = true;
  TiltedCylinderIF tunnel(channelRadius, cylinderAxis, cylinderOrigin, negativeInside);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(tunnel,0,vectDx);
  ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize);
  sprintf(ebfilename, "cyl%d.%dd.hdf5",fineLength, SpaceDim);
}
/***************/
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
  int testverbosity;
  ParmParse pp;
  pp.get("testverbosity", testverbosity);
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

  if(testverbosity > 1)
    pout() << "==============================================" << endl;
  for(int comp = 0; comp < names.size(); comp++)
    {
      if(testverbosity > 1)
        {
          pout() << "Comparing error in variable  " << names[comp] << " multiplied by vol frac" << endl;
          pout() << "==============================================" << endl;
        }
      EBNormType::NormMode normtype = EBNormType::OverBoth;
      if(testverbosity > 1)
        pout() << endl << "Using all uncovered cells and max norm." << endl;

      int inorm = 0;
      Real coarnorm = EBArith::norm(a_errorCoar,
                                    a_gridsCoar, a_ebislCoar,
                                    comp, inorm, normtype);

      Real finenorm = EBArith::norm(a_errorFine,
                                    a_gridsFine, a_ebislFine,
                                    comp, inorm, normtype);
      if(testverbosity > 1)
        {
          pout() << "Coarse Error Norm = " << coarnorm << endl;
          pout() << "Fine   Error Norm = " << finenorm << endl;
        }

      normsCoar[comp] = coarnorm;
      normsFine[comp] = finenorm;

      if((Abs(finenorm) > 1.0e-8) && (Abs(coarnorm) > 1.0e-8))
        {
          Real order = log(coarnorm/finenorm)/log(2.0);
          if(testverbosity > 1)
            pout() << "Order of scheme = " << order << endl;
          orders[comp] = order;
        }
    }
  if(testverbosity > 0)
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
      pout() << "\\caption{Solution error weighted with volume fraction convergence rates using L-0  norm. " << endl;
      pout() << "$h_f = \\frac{1}{" << nfine << "}$ and $h_c = 2 h_f$, $D = " << SpaceDim << "$ }"  << endl;
      pout() << "\\end{table}" << endl;
      pout() << endl << endl;
    }
  delete[] normsCoar;
  delete[] normsFine;
  delete[] orders   ;

}
/***************/
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

    Box domainFine, domainCoar;
    Real dxFine, dxCoar;
    makeFinestDomain(domainFine, dxFine);
    domainCoar = coarsen(domainFine, 2);
    dxCoar = 2.0*dxFine;

    LevelData<EBCellFAB>     errorFine, errorCoar;
    DisjointBoxLayout gridsFine, gridsCoar;
    EBISLayout        ebislFine, ebislCoar;
    int testverbosity;
    pp.get("testverbosity", testverbosity);
    if(testverbosity >= 1)
      pout() << "generating fine domain" << endl;

    makeGeometry(domainFine, dxFine);

    if(testverbosity >= 1)
      pout() << "generating fine error" << endl;

    generateError(errorFine,
                   gridsFine, ebislFine, domainFine, dxFine, false);

    if(testverbosity >= 1)
      pout() << "generating coarse domain" << endl;

    makeGeometry(domainCoar, dxCoar);

    if(testverbosity >= 1)
      pout() << "generating coarse error" << endl;

    generateError(errorCoar,
                  gridsCoar, ebislCoar, domainCoar, dxCoar, true);

    int fileout;
    pp.get("do_erroroutput", fileout);

    if(fileout == 1)
      {
        if(testverbosity > 1)
          pout() << "Outputting error to file" << endl;

        outputError(errorFine,
                    ebislFine, gridsFine, domainFine,
                    errorCoar,
                    ebislCoar, gridsCoar, domainCoar);
      }


    if(testverbosity > 1)
      pout() << "Comparing Error in Solution" << endl;

    //compute convergence rates and output to stdout
    compareError(errorFine, ebislFine, gridsFine, domainFine,
                 errorCoar, ebislCoar, gridsCoar, domainCoar);

#ifndef CH_MPI
    int scatterout;
    pp.get("do_scatterplots", scatterout);
    if(scatterout == 1)
      {
        if(testverbosity > 1)
          pout() << "scatter plotting solution error vs geometric error" << endl;

        ofstream ossh("solnerr.sh");
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
    dumpmemoryatexit();
  }
  MPI_Finalize();
#endif
}

