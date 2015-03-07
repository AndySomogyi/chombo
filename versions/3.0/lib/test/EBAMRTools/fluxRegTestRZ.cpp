#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// dtgraves Mon Feb 6 2006

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
#include "PlaneIF.H"
#include "DebugDump.H"
#include "SlabService.H"
#include "EBDebugDump.H"
#include "EBCoarseAverage.H"
#include "EBFluxRegister.H"
#include "EBFABView.H"

#include "UsingNamespace.H"

/******/
void getAllIrregRefinedLayouts(Vector<DisjointBoxLayout>& a_grids,
                               Vector<EBISLayout>&        a_ebisl,
                               const Vector<int>&         a_refRatio,
                               const Box&                 a_coarsestDomain,
                               const int&                 a_numLevels)
{
  ParmParse pp;
  int maxGridSize, blockFactor, bufferSize;
  Real fillRatio;
  pp.get("max_grid_size",maxGridSize);
  pp.get("block_factor", blockFactor);
  pp.get("tag_buffer", bufferSize);
  pp.get("fill_ratio", fillRatio);

  //split up coarsest domain by max box size and
  //make a dbl at coarsest level
  Vector<Box> boxesCoarsest;
  Vector<int> procsCoarsest;

  domainSplit(a_coarsestDomain, boxesCoarsest,
              maxGridSize, blockFactor);

  LoadBalance(procsCoarsest, boxesCoarsest);
  DisjointBoxLayout dblCoarsest(boxesCoarsest, procsCoarsest);

  //make a ebislayout at coarsest level.
  EBISLayout ebislCoarsest;
  const EBIndexSpace* const  ebisPtr = Chombo_EBIS::instance();
  ebisPtr->fillEBISLayout(ebislCoarsest,dblCoarsest,
                          a_coarsestDomain, 0);

  //make the tags
  IntVectSet tagsCoarsestLocal;
  //tag all irregular coarse iv's
  for(DataIterator dit = dblCoarsest.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = ebislCoarsest[dit()];
      const Box&       box = dblCoarsest.get(dit());
      tagsCoarsestLocal |= ebisBox.getIrregIVS(box);
    }

  //need to do parallel stuff here to make sure we all agree about tags
  int baseProc = 0;
  Vector<IntVectSet> allTags;
  gather(allTags, tagsCoarsestLocal, baseProc);
  IntVectSet tagsCoarsest;
  if(procID() == baseProc)
    {
      for(int p=0; p< allTags.size(); p++)
        {
          tagsCoarsest |= allTags[p];
        }
    }
  broadcast(tagsCoarsest, baseProc);

  //generate vector of grids.
  BRMeshRefine gridder(a_coarsestDomain, a_refRatio,
                       fillRatio,        blockFactor,
                       bufferSize,       maxGridSize);

  Vector<Vector<Box> > newMeshes(a_numLevels);
  Vector<Vector<Box> > oldMeshes(a_numLevels);
  oldMeshes[0]= boxesCoarsest;
  for(int ilev = 1; ilev < a_numLevels; ilev++)
    {
      oldMeshes[ilev] = Vector<Box>(1, refine(oldMeshes[ilev-1][0], a_refRatio[ilev-1]));
    }
  int baseLevel = 0;
  int maxLevel = a_numLevels-1;
  gridder.regrid(newMeshes, tagsCoarsest, baseLevel, maxLevel, oldMeshes);

  Vector<Vector<int> > newProcs(a_numLevels);
  for(int ilev = 0; ilev < a_numLevels; ilev++)
    {
      LoadBalance(newProcs[ilev], newMeshes[ilev]);
    }

  a_grids.resize(a_numLevels);
  a_ebisl.resize(a_numLevels);
  int numGhost = 4;
  Box domainLev = a_coarsestDomain;
  for(int ilev = 0; ilev < a_numLevels; ilev++)
    {
      a_grids[ilev] = DisjointBoxLayout(newMeshes[ilev], newProcs[ilev]);
      //generate ebislayout
      ebisPtr->fillEBISLayout(a_ebisl[ilev],a_grids[ilev],
                              domainLev, numGhost);
      domainLev.refine(a_refRatio[ilev]);
    }
}
/**********/
int makeGeometry(Box& a_domain, Real& a_dx)
{
  int eekflag =  0;
  //parse input file
  ParmParse pp;
  RealVect origin = RealVect::Zero;
  Vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);
  Real probhi;
  pp.get("prob_hi", probhi);
  Real dx = probhi/n_cell[0];
  a_dx = dx;
  CH_assert(n_cell.size() == SpaceDim);
  IntVect lo = IntVect::Zero;
  IntVect hi;
  for(int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if(n_cell[ivec] <= 0)
        {
          pout() << " bogus number of cells input = " << n_cell[ivec];
          return(-1);
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  a_domain.setSmall(lo);
  a_domain.setBig(hi);

  int iverbose;
  pp.get("verbose", iverbose);
  bool verbose = (iverbose==1);

  if(verbose)
    pout() << "ramp geometry" << endl;
  int upDir;
  int indepVar;
  Real startPt;
  Real slope;
  pp.get("up_dir",upDir);
  pp.get("indep_var",indepVar);
  pp.get("start_pt", startPt);
  pp.get("ramp_slope", slope);

  RealVect normal = RealVect::Zero;
  normal[upDir] = 1.0;
  normal[indepVar] = -slope;

  RealVect point = RealVect::Zero;
  point[upDir] = -slope*startPt;

  bool normalInside = true;

  PlaneIF ramp(normal,point,normalInside);

  RealVect vectDx = RealVect::Unit;
  vectDx *= dx;

  GeometryShop workshop(ramp,0,vectDx);
  //this generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain, origin, dx, workshop, 2048);
  return eekflag;
}
/*****************/

int fluxRegTest()
{
  int eekflag = 0;
  Interval interv(0,0);
  ParmParse pp;
  int iverbose;
  pp.get("verbose", iverbose);
  bool verbose = (iverbose==1);

  Box domainCoar, domainFine;
  Real dxFine, dxCoar;
  eekflag = makeGeometry(domainFine, dxFine);
  if(eekflag != 0)
    return eekflag;

  int numLevels = 2;
  Vector<int> refRatio;
  pp.getarr("ref_ratio", refRatio, 0, numLevels);

  int nref = refRatio[0];
  domainCoar = coarsen(domainFine, nref);
  dxCoar = dxFine*nref;

  Vector<DisjointBoxLayout> grids;
  Vector<EBISLayout>        ebisl;

  getAllIrregRefinedLayouts(grids,
                            ebisl,
                            refRatio,
                            domainCoar,
                            numLevels);

  DisjointBoxLayout gridsFine = grids[1];
  DisjointBoxLayout gridsCoar = grids[0];
  EBISLayout        ebislFine = ebisl[1];
  EBISLayout        ebislCoar = ebisl[0];

  EBCellFactory factFine(ebislFine);
  EBCellFactory factCoar(ebislCoar);
  IntVect ivghost = IntVect::Unit;

  LevelData<EBCellFAB> fineData(gridsFine, 1,ivghost, factFine);
  LevelData<EBCellFAB> coarData(gridsCoar, 1,ivghost, factCoar);

  //set data and flux registers to zero
  for(DataIterator coarIt = coarData.dataIterator();
      coarIt.ok(); ++coarIt)
      coarData[coarIt()].setVal(0.);
  for(DataIterator fineIt = fineData.dataIterator();
      fineIt.ok(); ++fineIt)
    fineData[fineIt()].setVal(0.);

  EBFluxRegister fluxReg(gridsFine,
                         gridsCoar,
                         ebislFine,
                         ebislCoar,
                         domainCoar,
                         nref, 1, Chombo_EBIS::instance());

  fluxReg.setToZero();
  //loop through directions
  Real scale = 1.0;
  Real fluxVal = 4.77;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      if(verbose)
        {
          pout() << "After setting to zero, the fluxReg in dir " << idir << " is" << endl;
          for(SideIterator sit; sit.ok(); ++sit)
            {
              fluxReg.dumpCoar(idir, sit());
              fluxReg.dumpFine(idir, sit());
            }
        }

      //increment and decrement
      //flux registers with equal size fluxes
      for(DataIterator coarIt = coarData.dataIterator();
          coarIt.ok(); ++coarIt)
        {
          const Box&  boxCoar = gridsCoar.get(coarIt());
          const EBISBox& ebisBox = ebislCoar[coarIt()];
          IntVectSet ivsBC(boxCoar);
          EBFaceFAB edgeFluxReg(ebisBox, boxCoar, idir, 1);
          BaseIFFAB<Real> edgeFluxIrr(ivsBC, ebisBox.getEBGraph(), idir, 1);
          edgeFluxReg.setVal(fluxVal);
          edgeFluxIrr.setVal(fluxVal);
          if(SpaceDim==2)
            {
              fluxReg.incrementCoarseRegulRZ(edgeFluxReg, scale,
                                             coarIt(), interv, idir,   dxCoar);
              fluxReg.incrementCoarseIrregulRZ(edgeFluxIrr, scale,
                                               coarIt(), interv, idir, dxCoar);
            }
          else
            {
              fluxReg.incrementCoarseRegular(edgeFluxReg, scale,
                                             coarIt(), interv, idir);
              fluxReg.incrementCoarseIrregular(edgeFluxIrr, scale,
                                               coarIt(), interv, idir);
            }
        }
      if(verbose)
        {
          pout() << "After coarse increment, the fluxReg in dir " << idir << " is" << endl;
          for(SideIterator sit; sit.ok(); ++sit)
            {
              fluxReg.dumpCoar(idir, sit());
              fluxReg.dumpFine(idir, sit());
            }
        }

      for(DataIterator fineIt = fineData.dataIterator();
          fineIt.ok(); ++fineIt)
        {
          const Box&  boxFine = gridsFine.get(fineIt());
          const EBISBox& ebisBox = ebislFine[fineIt()];
          IntVectSet ivsBF(boxFine);
          EBFaceFAB edgeFluxReg(ebisBox, boxFine, idir, 1);
          BaseIFFAB<Real> edgeFluxIrr(ivsBF, ebisBox.getEBGraph(), idir, 1);
          edgeFluxReg.setVal(fluxVal);
          edgeFluxIrr.setVal(fluxVal);

          for(SideIterator sit; sit.ok(); ++sit)
            {
              if(SpaceDim==2)
                {
                  fluxReg.incrementFineRegulRZ(edgeFluxReg, scale,
                                               fineIt(),  interv, idir, sit(), dxFine);

                  fluxReg.incrementFineIrregulRZ(edgeFluxIrr, scale,
                                                 fineIt(),  interv, idir, sit(), dxFine);
                }
              else
                {
                  fluxReg.incrementFineRegular(edgeFluxReg, scale,
                                               fineIt(),  interv, idir, sit());

                  fluxReg.incrementFineIrregular(edgeFluxIrr, scale,
                                                 fineIt(),  interv, idir, sit());
                }
              if(verbose)
                fluxReg.dumpFine(idir, sit());
            }
        }
      if(verbose)
        {
          pout() << "After fine increment, the fluxReg in dir " << idir << " is" << endl;
          for(SideIterator sit; sit.ok(); ++sit)
            {
              fluxReg.dumpCoar(idir, sit());
              fluxReg.dumpFine(idir, sit());
            }
        }
    }

  //reflux what ought to be zero into zero and the result should be zero
  //except where the coarse-fine boundary gets crossed by the embedded
  //boundary.  That should get fixed by the extramass thing.
  //there are no EBs crossing CF here because we are tagging
  //all irregular cells so we can just say it ought to be zero here
  fluxReg.reflux(coarData, interv, scale);

  for(int idir = 0; idir < SpaceDim; idir++)
    {
      if(verbose)
        {
          pout() << "After reflux, the fluxReg in dir " << idir << " is " << endl;
          for(SideIterator sit; sit.ok(); ++sit)
            {
              fluxReg.dumpCoar(idir, sit());
              fluxReg.dumpFine(idir, sit());
            }
        }
    }

  DataIterator datIt = coarData.dataIterator();
  for(datIt.reset(); datIt.ok(); ++datIt)
    {
      const EBCellFAB& data = coarData[datIt()];
      IntVectSet ivsBox(gridsCoar.get(datIt()));
      const EBISBox& ebisBox = ebislCoar[datIt()];
      Real rmax = 0.;
      Real rmin = 0.;
      for(VoFIterator vofit(ivsBox, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          rmax = Max(rmax, Abs(data(vof, 0)));
          rmin = Min(rmax, Abs(data(vof, 0)));
          if((rmax > 1.0e-10)||(rmin > 1.0e-10))
            {
              pout() << "EBFluxRegister failed the test at  "
                     << " vof = " << vof <<endl;
              eekflag = 42;
              return eekflag;
            }
        }
    }

  return eekflag;
}

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int eekflag=0;
  {//begin forever present scoping trick

    const char* in_file = "fluxregrz.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);

    eekflag = fluxRegTest();
    if(eekflag == 0)
      pout() << "ebfluxregister test passed." << endl;
    else
      pout() << "ebfluxregister test FAILED with error code "
           << eekflag << endl;

  }//end omnipresent scoping trick
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return eekflag;
}
/*****************/
