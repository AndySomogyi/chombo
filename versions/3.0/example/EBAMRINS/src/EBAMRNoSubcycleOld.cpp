#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "ParmParse.H"
#include "EBAMRNoSubcycleOld.H"
#include "EBAMRNoSubcycleF_F.H"
#include "EBAMRPoissonOpF_F.H"
#include "MeshRefine.H"
#include "BRMeshRefine.H"
#include "EBEllipticLoadBalance.H"
#include "EBArith.H"
#include "EBPWLFineInterp.H"
#include "EBPWQuadFineInterp.H"
#include "EBCoarseAverage.H"
#include "EBFluxFactory.H"
#include "EBCellFactory.H"
#include "EBLevelAdvect.H"
#include "EBGradDivFilter.H"
#include "EBPatchAdvect.H"
#include "REAL.H"
#include "EBPhysIBCFactory.H"
#include "EBAMRIO.H"
#include "BaseIFFactory.H"
#include "EBLevelRedist.H"
#include "BaseIVFactory.H"
#include "EBConstantCFInterp.H"
#include "EBArith.H"
#include "EBAMRDataOps.H"
#include "NeumannPoissonEBBC.H"
#include "DirichletPoissonEBBC.H"
#include <iomanip>
#include <cmath>
#include <cstdio>

#include "memusage.H"
#include "memtrack.H"

#include "UsingNamespace.H"

void setCoveredStuffToZero(LevelData<EBCellFAB>& a_vort);

/**********************/
EBAMRNoSubcycleOld::
EBAMRNoSubcycleOld(const AMRParameters&    a_params,
                const EBIBCFactory&     a_ibcfact,
                const ProblemDomain&    a_coarsestDomain,
                Real                    a_viscosity)
{
  CH_TIME("EBAMRNoSubcycleOld::EBAMRNoSubcycleOld");
  if (a_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld: constructor" << endl;
    }
  //set parameters of the run
  m_params    = a_params;
  m_viscosity = a_viscosity;
  m_viscousCalc = (m_viscosity > 0);

  //create initial and boundary condition object
  m_ibc    =   a_ibcfact.create();

  //resize vectors and set them where we can
  Real coarsestDx = m_params.m_domainLength/Real(a_coarsestDomain.size(0));
  int nlevels = m_params.m_maxLevel + 1;
  m_domain.resize(nlevels);
  m_dx.resize(nlevels);
  m_grids.resize(nlevels);
  m_ebisl.resize(nlevels);
  m_eblg.resize(nlevels);

  m_quadCFI.resize(nlevels);
  m_aveOper.resize(nlevels);
  m_aveSpac.resize(nlevels);
  m_ebLevAd.resize(nlevels);
  m_veloOld.resize(nlevels, NULL);
  m_veloNew.resize(nlevels, NULL);
  m_presOld.resize(nlevels, NULL);
  m_presNew.resize(nlevels, NULL);
  m_gphiOld.resize(nlevels, NULL);
  m_gphiNew.resize(nlevels, NULL);
  m_coveredFaceLitLo.resize(nlevels, NULL);
  m_coveredFaceLitHi.resize(nlevels, NULL);
  m_coveredSetsLitLo.resize(nlevels, NULL);
  m_coveredSetsLitHi.resize(nlevels, NULL);

  allocateDataHolders();

  m_domain[0] = a_coarsestDomain;
  m_dx[0]     =   coarsestDx;
  for(int ilev = 1; ilev < nlevels; ilev++)
    {
      CH_assert(m_params.m_refRatio[ilev-1] > 0);
      m_domain[ilev] = refine(m_domain[ilev-1], m_params.m_refRatio[ilev-1]);
      m_dx[ilev] = m_dx[ilev-1]/Real(m_params.m_refRatio[ilev-1]);
    }
  m_prescribedDt = -1.0;
  m_useFixedDt = false;

  m_ccProjector  = NULL;
  m_macProjector = NULL;
  m_time = 0.0;
  m_curStep = 0;
  m_dt = -1.0;
  //setup still needs to get called
  m_isSetup  = false;

  m_pointsUpdated = 0;
}
/**********/
void
EBAMRNoSubcycleOld::
allocateDataHolders()
{
  CH_TIME("EBAMRNoSubcycleOld::allocateDataHolders");
  for(int ilev = 0; ilev <= m_params.m_maxLevel; ilev++)
    {
      m_veloOld[ilev] = new LevelData<EBCellFAB>();
      m_veloNew[ilev] = new LevelData<EBCellFAB>();
      m_presOld[ilev] = new LevelData<EBCellFAB>();
      m_presNew[ilev] = new LevelData<EBCellFAB>();
      m_gphiOld[ilev] = new LevelData<EBCellFAB>();
      m_gphiNew[ilev] = new LevelData<EBCellFAB>();


      m_coveredFaceLitLo[ilev] = new LayoutData< Vector< Vector<VolIndex> > >();
      m_coveredFaceLitHi[ilev] = new LayoutData< Vector< Vector<VolIndex> > >();
      m_coveredSetsLitLo[ilev] = new LayoutData< Vector< IntVectSet > >      ();
      m_coveredSetsLitHi[ilev] = new LayoutData< Vector< IntVectSet > >      ();

    }
}
/**********/
EBAMRNoSubcycleOld::
~EBAMRNoSubcycleOld()
{
  CH_TIME("EBAMRNoSubcycleOld::~EBAMRNoSubcycleOld");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld:destructor" << endl;
    }
  delete m_ibc;
  for(int ilev = 0; ilev <= m_params.m_maxLevel; ilev++)
    {
      delete m_veloOld[ilev];
      delete m_veloNew[ilev];
      delete m_presOld[ilev];
      delete m_presNew[ilev];
      delete m_gphiOld[ilev];
      delete m_gphiNew[ilev];

      delete m_coveredFaceLitLo[ilev];
      delete m_coveredFaceLitHi[ilev];
      delete m_coveredSetsLitLo[ilev];
      delete m_coveredSetsLitHi[ilev];

    }
  if(m_ccProjector !=  NULL)
    {
      delete m_ccProjector;
    }
  if(m_macProjector !=  NULL)
    {
      delete m_macProjector;
    }
}
/**********/
void
EBAMRNoSubcycleOld::
setupForAMRRun()
{
  CH_TIME("EBAMRNoSubcycleOld::setupForAMRMRun");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::setupForAMRRun" << endl;
    }
  m_isSetup= true;
  m_doRestart    = false;
  // we're generating the hiearachy dyamically
  // modelled on AMR::initialGrid()

  Vector<Vector<Box> > old_boxes(1);
  Vector<Vector<Box> > new_boxes;

  // also keep old tags around
  Vector<IntVectSet> oldTags;

  // define base mesh
  // chop up base level into grids to satisfy box size requirements
  domainSplit(m_domain[0], old_boxes[0], m_params.m_maxBoxSize,
              m_params.m_blockFactor);

  if(m_params.m_maxLevel == 0)
    {
      EBLevelDataOps::pruneCoveredBoxes(old_boxes[0], m_domain[0], Chombo_EBIS::instance());
    }

  m_finestLevel = 0;

  // now generate more levels if necessary
  int top_level = 0;

  bool moreLevels = (m_params.m_maxLevel > 0);

  // create grid generation object
  BRMeshRefine meshrefine;

  if (moreLevels)
    {
      meshrefine.define(m_domain[0],              m_params.m_refRatio,
                        m_params.m_fillRatio,     m_params.m_blockFactor,
                        m_params.m_nestingRadius, m_params.m_maxBoxSize);
    }

  // now intialize data for existing hierarchy
  initialGrid(old_boxes);

  //print_memory_line("before initialData");
  //UnfreedMemory();
  initialData();
  //print_memory_line("after initialData");
  //UnfreedMemory();

  while (moreLevels)
    {
      // default is moreLevels = false
      // (only repeat loop in the case where a new level
      // is generated which is still coarser than maxLevel)
      moreLevels = false;

      int base_level = 0;
      int old_top_level = top_level;

      Vector<IntVectSet> tagsVect(top_level+1);
      tagCells(tagsVect);

      //print_memory_line("before meshRefine regrid");
      //UnfreedMemory();
      int new_finest = meshrefine.regrid(new_boxes, tagsVect,
                                         base_level, top_level,
                                         old_boxes);
      //print_memory_line("after meshRefine regrid");
      //UnfreedMemory();

      if (new_finest > top_level) top_level++;

      old_boxes = new_boxes;

      // now see if we need another pass through grid generation
      if ((top_level<m_params.m_maxLevel) && (top_level > old_top_level))
        moreLevels = true;


      // if we added another level, reinintialize everything again
      if(top_level > old_top_level)
        {
          initialGrid(new_boxes);
          initialData();
        }
    } // end loop over regridding passes

  defineIrregularData();
  // finally, call post-initialization
  postInitialize();

}
void
EBAMRNoSubcycleOld::
filter( Vector<LevelData<EBCellFAB>* >&   a_veloc)
{
  //finestlevel vs. maxlevel fix
  CH_TIME("EBAMRNoSubcycleOld::filter");
  for(int ifilter = 0; ifilter < m_params.m_numFilterIterations; ifilter++)
    {
      averageDown(a_veloc);
      for(int ilev = 0; ilev <= m_finestLevel; ilev++)
        {
          const DisjointBoxLayout*           coarGridsPtr = NULL;
          const EBISLayout*                  coarEBISLPtr = NULL;
          const LevelData<EBCellFAB>*        coarVelocPtr = NULL;
          int refRat = -1;
          if(ilev > 0)
            {
              coarGridsPtr = &m_grids[ilev-1];
              coarEBISLPtr = &m_ebisl[ilev-1];
              coarVelocPtr =  a_veloc[ilev-1];
              refRat = m_params.m_refRatio[ilev-1];
            }
          EBGradDivFilter gdFilt(m_grids[ilev], coarGridsPtr, m_ebisl[ilev], coarEBISLPtr, m_domain[ilev], m_dx[ilev]*RealVect::Unit, refRat);
          //flux velocity used for boundary conditions.  this assumes (correctly) that the
          //projection has just happened.
          Vector< LevelData<EBFluxFAB>* >& fluxVel = m_ccProjector->getMacVelocity();
          gdFilt.filter(*a_veloc[ilev], *fluxVel[ilev], coarVelocPtr);

        }
    }
}
/**************************/
Real
EBAMRNoSubcycleOld::
run(Real a_maxTime, int a_maxStep)
{
  CH_TIME("EBAMRNoSubcycleOld::run");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::run" << endl;
    }
  CH_assert(m_isSetup);

  // only call computeInitialDt if we're not doing a restart
  if (!m_doRestart)
    {
      // do initial velocity projection here instead of
      //in postInitialize to save time if we are not using the object for
      //a run
      if (m_params.m_verbosity > 0)
        {
          pout() << "EBAMRNoSubcycleOld:cc projecting initial vel" << endl;
        }
      m_ccProjector->project(m_veloOld, m_gphiOld);

      filter(m_veloOld);

      Interval interv(0, SpaceDim-1);
      for(int ilev = 0; ilev <= m_finestLevel; ilev++)
        {
          //set gradient back to zero
          EBLevelDataOps::setVal(   *m_gphiOld[ilev], 0.0);
          EBLevelDataOps::setVal(   *m_gphiNew[ilev], 0.0);
          m_veloOld[ilev]->copyTo(interv, *m_veloNew[ilev], interv);
        }

      m_dt = computeInitialDt();
      m_curStep = 0;
      m_time = 0.0;
    }
  else
    {
      m_dt = computeDt();
    }

  for(int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      setCoveredStuffToZero(*m_veloNew[ilev]);
      setCoveredStuffToZero(*m_veloOld[ilev]);
      setCoveredStuffToZero(*m_presOld[ilev]);
      setCoveredStuffToZero(*m_presNew[ilev]);
    }

  if (m_params.m_verbosity > 0)
    {
      pout() << "EBAMRNoSubcycleOld:starting run " << endl;
    }
  /// advance solution until done

  if (!m_doRestart)
    {
      pout() << "EBAMRNoSubcycleOld: iterating get gphi " << endl;
      m_advanceGphiOnly = true;
      for(int iter = 0; iter < m_params.m_initIterations; iter++)
        {
          pout() << "######################################" <<  endl;
          pout() << "EBAMRNoSubcycleOld: gphi iteration " << iter  << ", dt =" <<  m_dt  << endl;
          advance();
          postTimeStep();
        }
    }
  m_advanceGphiOnly = false;

  pout() << "EBAMRNoSubcycleOld: starting time advance  " << endl;
  while((a_maxTime > m_time) && (m_curStep < a_maxStep))
    {
      pout() << "EBAMRNoSubcycleOld: nstep = " << m_curStep << ", time = " << m_time << ", dt =" <<  m_dt  << endl;

      // dump plotfile and checkpointfile before regridding
#ifdef CH_USE_HDF5
      if ((m_curStep%m_params.m_plotInterval == 0) && (m_params.m_plotInterval > 0))
        {
          pout() << "writing plot file" << endl;
          writePlotFile();
        }

      if ((m_curStep%m_params.m_checkpointInterval == 0) && (m_params.m_checkpointInterval > 0))
        {
          pout() << "writing checkpoint file" << endl;
          writeCheckpointFile();
        }
#endif

      // do regridding if appropriate
      if ((m_curStep != 0) && (m_curStep%m_params.m_regridInterval == 0)
          && (m_params.m_regridInterval > 0)
          && (m_params.m_maxLevel > 0)
          )
        {
          regrid();
        }

      // compute new dt
      m_dt = computeDt();

      // do timestep
      advance();

      postTimeStep();

      //advance time and step number
      ++m_curStep;
      m_time += m_dt;

      pout() << "########## End of time step ##########" <<  endl;

    }  // end loop over timesteps

  pout() << "total number of points updated = " << m_pointsUpdated << endl;
  pout() << "" << endl;

  conclude();

  return m_time;
}
/******************/
void
EBAMRNoSubcycleOld::
conclude()
{
  CH_TIME("EBAMRNoSubcycleOld::conclude");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::conclude" << endl;
    }
#ifdef CH_USE_HDF5
  if (m_params.m_plotInterval >= 0)
    {
      writePlotFile();
    }

  if (m_params.m_checkpointInterval >= 0)
    {
      writeCheckpointFile();
    }
#endif
}
/*****************/
void
EBAMRNoSubcycleOld::tagCells(Vector<IntVectSet>& a_tags)
{
  CH_TIME("EBAMRNoSubcycleOld::tagCells");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::tagCells" << endl;
    }
  int numlevels = a_tags.size();
  if (numlevels > m_finestLevel+1) numlevels = m_finestLevel+1;


  for (int lev=0; lev<numlevels; lev++)
    {
      IntVectSet& levelTags = a_tags[lev];

      tagCellsLevel(levelTags, lev);
    }

}
/*****************/
void
EBAMRNoSubcycleOld::regrid()
{
  CH_TIME("EBAMRNoSubcycleOld::regrid");
  //smoothing to make regridding less heinous
  preRegrid();

  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::regrid" << endl;
    }
  // don't regrid base grid
  int lbase = 0;

  if (m_params.m_maxLevel > 0)
    {
      // first, construct tags
      int top_level = Min(m_finestLevel, m_params.m_maxLevel-1);
      Vector<IntVectSet> tagsVect(top_level+1);
      Vector<Vector<Box> > new_grids;
      Vector<Vector<Box> > vectBoxes(top_level+1);
      for(int ilev = 0; ilev <= top_level; ilev++)
        {
//           vectBoxes[ilev] = Vector<Box>(1, m_domain[ilev].domainBox());
          domainSplit(m_domain[ilev], vectBoxes[ilev], m_params.m_maxBoxSize,
                      m_params.m_blockFactor);
        }

      // do tagging
      tagCells(tagsVect);


      int new_finest_level;

      BRMeshRefine meshrefine(m_domain[0], m_params.m_refRatio,
                              m_params.m_fillRatio, m_params.m_blockFactor,
                              m_params.m_nestingRadius, m_params.m_maxBoxSize);




      new_finest_level = meshrefine.regrid(new_grids,
                                           tagsVect,
                                           lbase,
                                           top_level,
                                           vectBoxes);



      // can only add one level at a time
      new_finest_level = Min(m_finestLevel+1, new_finest_level);

      if ((m_finestLevel != new_finest_level) && (m_params.m_verbosity >= 2))
        {
          pout() << "finest level changes here from "
                 << m_finestLevel << " to "
                 << new_finest_level << endl;
        }

      //allow for levels to change
      m_finestLevel = Min(new_finest_level, m_params.m_maxLevel);


      // now redefine grid hierarchy
      regrid(new_grids);

      defineIrregularData();

      // finish up
      postRegrid();
    } // end if max level > 0
}
/*********************/
void
EBAMRNoSubcycleOld::preRegrid()
{
  CH_TIME("EBAMRNoSubcycleOld::preRegrid");
  if(m_viscousCalc && m_params.m_doRegridSmoothing)
    {
      if (m_params.m_verbosity > 0)
        {
          pout() << " smoothing velocity before regridding " << endl;
        }
      //allocate bunch of scratch stuff
      allocateTemporaries();
      //make cellscratch2 == 0 for residual calc
      EBAMRDataOps::setVal(m_cellScratc2, 0.0);

      //make vel = (I-mu lapl)vel
      Real alpha = 1.0;
//       Real beta = -4.*m_dt*m_viscosity;
      Real beta = -4.0*m_dt;
//       pout() <<"preRegrid:  alpha = "<<alpha<<", beta= "<<beta<< endl;
      int coarsestLevel = 0;
      for(int velComp = 0; velComp < SpaceDim; velComp++)
        {
          DirichletPoissonEBBC::s_velComp = velComp;
          //copy velo into cellscratc1, make cellscratch2 == 0
          for(int ilev = 0; ilev <= m_finestLevel; ilev++)
            {
              Interval srcInterv(velComp, velComp);
              Interval dstInterv(0 , 0);
              m_veloOld[ilev]->copyTo(srcInterv, *m_cellScratc1[ilev], dstInterv);
            }
          Real junkNorm = 0;
          //apply the operator (by computing the residual with rhs = 0, and * -1
          m_tgaSolver[velComp]->resetAlphaAndBeta(alpha, beta);
          junkNorm = m_solver[velComp]->computeAMRResidual(m_cellScratch,  //comes out holding -(alpha + beta lapl) vel
                                                           m_cellScratc1, //holds velocity component
                                                           m_cellScratc2, //holds zero
                                                           m_finestLevel,
                                                           coarsestLevel,
                                                           false);  //not homogeneous bcs


          //make cellscratch hold (alpha + beta lapl) vel
          EBAMRDataOps::scale(m_cellScratch,-1.0);
          averageDown(m_cellScratch);

          //copy into velocity containers (makes vel := (alpha + beta lapl)vel)
          for(int ilev = 0; ilev <= m_finestLevel; ilev++)
            {
              Interval dstInterv(velComp, velComp);
              Interval srcInterv(0 , 0);
              m_cellScratch[ilev]->copyTo(srcInterv, *m_veloOld[ilev], dstInterv);
              m_cellScratch[ilev]->copyTo(srcInterv, *m_veloNew[ilev], dstInterv);
            }
        }
      //remove all that scratch space
      deleteTemporaries();
    }// end loop over velocity compoents
}
/*********************/
void
EBAMRNoSubcycleOld::postRegrid()
{
  CH_TIME("EBAMRNoSubcycleOld::postRegrid");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::postRegrid" << endl;
    }
  if(m_viscousCalc && m_params.m_doRegridSmoothing)
    {
      if (m_params.m_verbosity > 0)
        {
          pout() << " de smoothing velocity after regridding " << endl;
        }
      //allocate bunch of scratch stuff
      allocateTemporaries();

      //make vel = (I-mu lapl)^-1 vel
      Real alpha = 1.0;
//       Real beta = -4.0*m_dt*m_viscosity;
      Real beta = -4.0*m_dt;
//       pout() <<"postRegrid:  alpha = "<<alpha<<", beta= "<<beta<< endl;
      for(int velComp = 0; velComp < SpaceDim; velComp++)
        {
          DirichletPoissonEBBC::s_velComp = velComp;
          //copy velo into cellscratc1, make cellscratch2 == 0
          for(int ilev = 0; ilev <= m_finestLevel; ilev++)
            {
              Interval srcInterv(velComp, velComp);
              Interval dstInterv(0 , 0);
              m_veloOld[ilev]->copyTo(srcInterv, *m_cellScratc1[ilev], dstInterv);
            }

          EBAMRDataOps::setToZero(m_cellScratch);
          //solve (I - mu lapl)velnew = velold
          m_tgaSolver[velComp]->resetAlphaAndBeta(alpha, beta);
          m_solver[velComp]->solve(m_cellScratch, //comes out holding desmoothed vel
                                   m_cellScratc1, //rhs = vel
                                   m_finestLevel, 0);

          averageDown(m_cellScratch);

          //copy into velocity containers (makes vel := (alpha + beta lapl)^-1 vel)
          for(int ilev = 0; ilev <= m_finestLevel; ilev++)
            {
              Interval dstInterv(velComp, velComp);
              Interval srcInterv(0 , 0);
              m_cellScratch[ilev]->copyTo(srcInterv, *m_veloOld[ilev], dstInterv);
              m_cellScratch[ilev]->copyTo(srcInterv, *m_veloNew[ilev], dstInterv);
            }
        }
      //remove all that scratch space
      deleteTemporaries();
    }// end loop over velocity compoents

  if (m_params.m_verbosity > 0)
    {
      pout() << "EBAMRNoSubcycleOld:cc projecting after regrid" << endl;
    }
  m_ccProjector->project(m_veloNew, m_gphiNew);
  //put gphi back to gphi old. the gradient coming out of this projection is not
  //really correlated to the pressure gradient.
  Interval interv(0, SpaceDim-1);
  for(int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      m_gphiOld[ilev]->copyTo(interv, *m_gphiNew[ilev], interv);
    }

  filter(m_veloNew);
  for(int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      m_veloNew[ilev]->copyTo(interv, *m_veloOld[ilev], interv);
    }

  // re-initialize pressure here
  pout() << "EBAMRNoSubcycleOld: iterating get gphi after regrid" << endl;
  m_advanceGphiOnly = true;;
  if(m_params.m_gphiIterations > 0)
    {
      EBAMRDataOps::setToZero(m_gphiOld);
      EBAMRDataOps::setToZero(m_gphiNew);
      EBAMRDataOps::setToZero(m_presOld);
      EBAMRDataOps::setToZero(m_presNew);
    }
  for(int iter = 0; iter < m_params.m_gphiIterations; iter++)
    {
      pout() << "EBAMRNoSubcycleOld: gphi iteration " << iter  << endl;
      advance();
    }
  m_advanceGphiOnly = false;
}
/*********************/
void EBAMRNoSubcycleOld::
averageDown(Vector<LevelData<EBFluxFAB>* >&  a_data)
{
  CH_TIME("EBAMRNoSubcycleOld::averageDown(fluxData)");
  // do average down here
  for (int ilev = m_finestLevel; ilev > 0; ilev--)
    {
      int ncomp = a_data[ilev]->nComp();
      Interval interval(0, ncomp-1);
      RefCountedPtr<EBCoarseAverage> avePtr = m_aveOper[ilev];
      if(ncomp == SpaceDim)
        {
          avePtr = m_aveSpac[ilev];
        }
      EBCoarseAverage& ebaverage = *avePtr;

      ebaverage.average(*a_data[ilev-1],
                        *a_data[ilev  ],
                        interval);
    } // end loop over levels
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      int ncomp = a_data[ilev]->nComp();
      Interval interval(0, ncomp-1);
      a_data[ilev]->exchange(interval);
    }
}
/**************************/
void
EBAMRNoSubcycleOld::
averageDown(Vector<LevelData<EBCellFAB>* >&  a_data)
{
  CH_TIME("EBAMRNoSubcycleOld::averageDown(celldata)");
  // do average down here
  for (int ilev = m_finestLevel; ilev > 0; ilev--)
    {
      int ncomp = a_data[ilev]->nComp();
      Interval interval(0, ncomp-1);
      RefCountedPtr<EBCoarseAverage> avePtr = m_aveOper[ilev];
      if(ncomp == SpaceDim)
        {
          avePtr = m_aveSpac[ilev];
        }
      EBCoarseAverage& ebaverage = *avePtr;

      ebaverage.average(*a_data[ilev-1],
                        *a_data[ilev  ],
                        interval);
    } // end loop over levels
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      int ncomp = a_data[ilev]->nComp();
      Interval interval(0, ncomp-1);
      a_data[ilev]->exchange(interval);
    }
}
/*********************/
void
EBAMRNoSubcycleOld::defineGrids(const Vector<Vector<Box> >& a_vectBoxes)
{
  CH_TIME("EBAMRNoSubcycleOld::defineGrids");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::defineGrids" << endl;
    }
  m_finestLevel = 0;
  // now define all of the storage we need
  for (int ilev = 0; ilev< a_vectBoxes.size(); ilev++)
    {
      if (a_vectBoxes[ilev].size() > 0)
        {
          m_finestLevel = ilev;
          // first do load balance
          Vector<int> procAssign;
          //print_memory_line("before defineGrids mortonOrdering");
          //UnfreedMemory();
          mortonOrdering((Vector<Box>&)(a_vectBoxes[ilev]));
          //print_memory_line("after defineGrids mortonOrdering");
          //UnfreedMemory();
          EBEllipticLoadBalance(procAssign,  a_vectBoxes[ilev], m_domain[ilev]);
          //print_memory_line("after defineGrids EBEllipticLoadBalance");
          //UnfreedMemory();
          m_grids[ilev] = DisjointBoxLayout();
          m_grids[ilev].define(a_vectBoxes[ilev], procAssign);
        }
    }
}
/*********************/
void
EBAMRNoSubcycleOld::defineEBISLs()
{
  CH_TIME("EBAMRNoSubcycleOld::defineEBISLs");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::defineEBISLs" << endl;
    }
  int numEBG =   4;
  m_eblg.resize(m_finestLevel+1);
  // now define all of the storage we need
  RefCountedPtr<EBPhysIBCFactory> advectBC = m_ibc->getVelAdvectBC(0); //this gets reset
  RefCountedPtr<EBPatchAdvectFactory> fact = RefCountedPtr<EBPatchAdvectFactory> (new EBPatchAdvectFactory(advectBC, m_params.m_useLimiting));

  for (int ilev = 0; ilev<= m_finestLevel; ilev++)
    {
      const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
      m_eblg[ilev] = EBLevelGrid(m_grids[ilev], m_domain[ilev], numEBG, ebisPtr);
      m_ebisl[ilev] = m_eblg[ilev].getEBISL();
      DisjointBoxLayout coarDBL;
      EBISLayout        coarEBISL;
      int refRat = 2;
       if(ilev > 0)
         {
           coarDBL =  m_grids[ilev-1];
           coarEBISL= m_ebisl[ilev-1];
           refRat = m_params.m_refRatio[ilev-1];
         }
       bool hasCoarser = (ilev > 0);
       bool hasFiner   = (ilev < m_finestLevel);
       m_ebLevAd[ilev]  = RefCountedPtr<EBLevelAdvect>(new EBLevelAdvect(m_grids[ilev],
                                                                          coarDBL,
                                                                          m_ebisl[ilev],
                                                                          coarEBISL,
                                                                          ProblemDomain(m_domain[ilev]),
                                                                          refRat,
                                                                          m_dx[ilev]*RealVect::Unit,
                                                                          hasCoarser,
                                                                          hasFiner,
                                                                          &(*fact)));
      if(ilev > 0)
        {
          //always one component for quadcfi---only way to get reuse
          int nvarquad = 1;
          m_quadCFI[ilev]  = RefCountedPtr<EBQuadCFInterp>(new  EBQuadCFInterp(m_grids[ilev], m_grids[ilev-1],
                                                                               m_ebisl[ilev], m_ebisl[ilev-1],
                                                                               m_domain[ilev-1],
                                                                               m_params.m_refRatio[ilev-1], nvarquad,
                                                                               *m_eblg[ilev].getCFIVS()));

          m_aveOper[ilev]  = RefCountedPtr<EBCoarseAverage>(new EBCoarseAverage(m_grids[ilev], m_grids[ilev-1],
                                                                                m_ebisl[ilev], m_ebisl[ilev-1],
                                                                                m_domain[ilev-1],
                                                                                m_params.m_refRatio[ilev-1], nvarquad,
                                                                                Chombo_EBIS::instance()));

          m_aveSpac[ilev]  = RefCountedPtr<EBCoarseAverage>(new EBCoarseAverage(m_grids[ilev], m_grids[ilev-1],
                                                                                m_ebisl[ilev], m_ebisl[ilev-1],
                                                                                m_domain[ilev-1],
                                                                                m_params.m_refRatio[ilev-1], SpaceDim,
                                                                                Chombo_EBIS::instance()));
        }
    }
  long long totalPoints = 0;
  long long totalBoxes  = 0;
  int numLevels = m_finestLevel + 1;
  for(int ilev = 0; ilev < numLevels; ilev++)
    {
      long long pointsThisLevel = 0;
      for(LayoutIterator lit = m_grids[ilev].layoutIterator(); lit.ok(); ++lit)
        {
          pointsThisLevel += m_grids[ilev][lit()].numPts();
        }
      totalPoints += pointsThisLevel;
      totalBoxes += m_grids[ilev].size();
      pout() << "getAllIrregRefineLayouts:level[" << ilev
             << "], number of boxes = " << m_grids[ilev].size()
             << ", number of points = " << pointsThisLevel << endl;
    }
  pout() << "getAllIrregRefineLayouts:"
         <<  "   total boxes = " << totalBoxes
         <<  ", total points = " << totalPoints <<  endl;
}
/*********************/
void
EBAMRNoSubcycleOld::initialGrid(const Vector<Vector<Box> >& a_vectBoxes)
{
  CH_TIME("EBAMRNoSubcycleOld::initialGrid");
  if (m_params.m_verbosity >= 3)
    {
      pout () << "EBAMRNoSubcycleOld::initialGrid "  << endl;
    }
  //define grids and finest level
  //print_memory_line("before initialGrid defineGrids");
  //UnfreedMemory();
  defineGrids(a_vectBoxes);
  //print_memory_line("before initialGrid defineEBISLs");
  //UnfreedMemory();
  defineEBISLs();
  //print_memory_line("before initialGrid defineOld");
  //UnfreedMemory();
  defineOld();
  //print_memory_line("before initialGrid defineNew");
  //UnfreedMemory();
  defineNew();
  //print_memory_line("before initialGrid defineProjections");
  //UnfreedMemory();
  defineProjections();
  //print_memory_line("after initialGrid defineProjections");
  //UnfreedMemory();
}
/*********************/
void
EBAMRNoSubcycleOld::
defineOld()
{
  CH_TIME("EBAMRNoSubcycleOld::defineOld");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::defineOld" << endl;
    }
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      EBCellFactory ebcellfact(m_ebisl[ilev]);
      m_veloOld[ilev]->define(m_grids[ilev], SpaceDim,  4*IntVect::Unit, ebcellfact);
      m_gphiOld[ilev]->define(m_grids[ilev], SpaceDim,    IntVect::Unit, ebcellfact);
      m_presOld[ilev]->define(m_grids[ilev],        1,    IntVect::Unit, ebcellfact);
    }
  EBAMRDataOps::setToZero(m_veloOld);
  EBAMRDataOps::setToZero(m_presOld);
  EBAMRDataOps::setToZero(m_gphiOld);
}
/*********************/
void
EBAMRNoSubcycleOld::
defineNew()
{
  CH_TIME("EBAMRNoSubcycleOld::defineNew");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::defineNew" << endl;
    }
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      EBCellFactory ebcellfact(m_ebisl[ilev]);
      m_veloNew[ilev]->define(m_grids[ilev], SpaceDim,  4*IntVect::Unit, ebcellfact);
      m_gphiNew[ilev]->define(m_grids[ilev], SpaceDim,    IntVect::Unit, ebcellfact);
      m_presNew[ilev]->define(m_grids[ilev],        1,    IntVect::Unit, ebcellfact);
    }
  EBAMRDataOps::setToZero(m_veloNew);
  EBAMRDataOps::setToZero(m_gphiNew);
}
/*********************/
void
EBAMRNoSubcycleOld::
defineProjections()
{
  CH_TIME("EBAMRNoSubcycleOld::defineProjections");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::defineProjections" << endl;
    }
  if(m_ccProjector != NULL)
    {
      delete m_ccProjector;
      m_ccProjector = NULL;
    }
  if(m_macProjector != NULL)
    {
      delete m_macProjector;
      m_macProjector = NULL;
    }
  RefCountedPtr<BaseDomainBCFactory> macBCVel = m_ibc->getMACVelBC();
  RefCountedPtr<BaseDomainBCFactory> celBCPhi = m_ibc->getPressBC();
  RefCountedPtr<BaseEBBCFactory>     ebbcVelo = m_ibc->getVelocityEBBC(0);
  RefCountedPtr<BaseEBBCFactory>     ebbcPhi  = m_ibc->getPressureEBBC();

  int numLevels = m_finestLevel + 1;
  RealVect coarDxVect = m_dx[0]*RealVect::Unit;

  ParmParse pp;

  pp.query("mg_num_smooths", m_params.m_numSmooth);

  pp.query("mg_relax_type", m_params.m_relaxType);
  bool lazy = false;
  pp.query("mg_relax_lazy", lazy);

  int bottomSolverType = 0;
  pp.query("mg_bottom_solver", bottomSolverType);

  int numPreCond = 4;
  pp.query("mg_num_precond", numPreCond);

  pp.query("mg_hang", m_params.m_hang);
  pp.query("mg_tolerance", m_params.m_tolerance);

  EBAMRPoissonOp::doLazyRelax(lazy);

  m_macProjector =   new EBCompositeMACProjector(m_eblg, m_params.m_refRatio, m_quadCFI,
                                                 coarDxVect, RealVect::Zero,
                                                 macBCVel, celBCPhi, ebbcPhi,
                                                 m_params.m_subtractOffMean, numLevels,
                                                 m_params.m_verbosity,numPreCond,m_time,m_params.m_relaxType,bottomSolverType);


  m_macProjector->setSolverParams(m_params.m_numSmooth, m_params.m_iterMax, m_params.m_mgCycle,
                                  m_params.m_hang, m_params.m_tolerance,
                                  m_params.m_verbosity);

  m_ccProjector  =   new EBCompositeCCProjector( m_eblg, m_params.m_refRatio, m_quadCFI,
                                                 coarDxVect, RealVect::Zero,
                                                 macBCVel, celBCPhi, ebbcPhi,
                                                 m_params.m_subtractOffMean, numLevels,
                                                 m_params.m_verbosity, numPreCond, m_time, m_params.m_relaxType,
                                                 bottomSolverType,m_macProjector);
  if(m_viscousCalc)
    {
      pout() << "using TGA for viscous solver" << endl;

      int lbase = 0;
      int lmax = m_finestLevel;

      LinearSolver<LevelData<EBCellFAB> >* bottomSolverPtr = NULL;
      if(bottomSolverType == 0)//BiCGStab
        {
          bottomSolverPtr = &m_bottomSolver;
        }
      else if(bottomSolverType == 1)//simple relaxation
        {
          bottomSolverPtr = &m_bottomSolverSimp;
        }
      else
        {
          MayDay::Error("EBAMRNoSubcycleOld::defineProjections -- bad bottom solver type");
        }

      Vector<LevelData<EBCellFAB>* > phi(numLevels);
      Vector<LevelData<EBCellFAB>* > rhs(numLevels);
      for(int ilev = 0; ilev < numLevels; ilev++)
        {
          EBCellFactory ebcellfact(m_ebisl[ilev]);
          phi[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], 1, 4*IntVect::Unit, ebcellfact);
          rhs[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], 1, 4*IntVect::Unit, ebcellfact);
        }

      m_bottomSolver.m_verbosity = 0;
      m_bottomSolverSimp.setNumSmooths(64*m_params.m_numSmooth);
      pp.query("mg_num_cycles", m_params.m_mgCycle);
      pp.query("mg_iter_max", m_params.m_iterMax);
      pp.query("mg_norm_thresh", m_params.m_normThresh);
      int numLevels = m_finestLevel + 1;

      Real alpha = 1.0;
      Real beta = m_viscosity;

      ProblemDomain level0Dom = m_eblg[0].getDomain();

      for(int idir = 0;  idir < SpaceDim; idir++)
        {
          DirichletPoissonEBBC::s_velComp = idir;
          RefCountedPtr<BaseDomainBCFactory> celBCVel = m_ibc->getVelBC(idir);

          EBAMRPoissonOpFactory opFactory(m_eblg,
                                          m_params.m_refRatio,
                                          m_quadCFI,
                                          coarDxVect,
                                          RealVect::Zero,
                                          numPreCond,
                                          m_params.m_relaxType,
                                          celBCVel,
                                          ebbcVelo,
                                          alpha,
                                          beta,
                                          m_time,
                                          4*IntVect::Unit,
                                          4*IntVect::Unit);

          m_solver[idir] = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > (new AMRMultiGrid<LevelData<EBCellFAB> >);
          m_solver[idir]->define(level0Dom,
                                 opFactory,
                                 bottomSolverPtr,
                                 numLevels);

          m_solver[idir]->m_pre        =  m_params.m_numSmooth;
          m_solver[idir]->m_post       =  m_params.m_numSmooth;
          m_solver[idir]->m_bottom     =  m_params.m_numSmooth;
          m_solver[idir]->m_hang       =  m_params.m_hang;
          m_solver[idir]->m_eps        =  m_params.m_tolerance;
          m_solver[idir]->m_verbosity  =  m_params.m_verbosity;
          m_solver[idir]->m_iterMax    =  m_params.m_iterMax;
          m_solver[idir]->m_normThresh =  m_params.m_normThresh;
          m_solver[idir]->setMGCycle(m_params.m_mgCycle);

          m_solver[idir]->init(phi,
                               rhs,
                               lmax,
                               lbase);

          m_tgaSolver[idir] = RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > >
            (new AMRTGA<LevelData<EBCellFAB> >(m_solver[idir],
                                               opFactory,
                                               level0Dom,
                                               m_params.m_refRatio,
                                               numLevels,
                                               m_params.m_verbosity));
        }
      for(int ilev = 0; ilev < numLevels; ilev++)
        {
          delete phi[ilev];
          delete rhs[ilev];
        }
    }
}
/*********************/
void
EBAMRNoSubcycleOld::initialData()
{
  CH_TIME("EBAMRNoSubcycleOld::initialData");
  if (m_params.m_verbosity >= 3)
    {
      pout () << "EBAMRNoSubcycleOld::initialData "  << endl;
    }

  for (int ilev = 0; ilev<= m_finestLevel; ilev++)
    {
      RealVect dxLev = m_dx[ilev]*RealVect::Unit;
      m_ibc->initializeVelocity(*m_veloOld[ilev], m_grids[ilev], m_ebisl[ilev], m_domain[ilev], RealVect::Zero, m_time, dxLev);
      m_ibc->initializeVelocity(*m_veloNew[ilev], m_grids[ilev], m_ebisl[ilev], m_domain[ilev], RealVect::Zero, m_time, dxLev);
      EBLevelDataOps::setVal(   *m_gphiOld[ilev], 0.0);
      EBLevelDataOps::setVal(   *m_gphiNew[ilev], 0.0);
      EBLevelDataOps::setVal(   *m_presOld[ilev], 0.0);
      EBLevelDataOps::setVal(   *m_presNew[ilev], 0.0);
    }
}
/*********************/
void
EBAMRNoSubcycleOld::postInitialize()
{
  CH_TIME("EBAMRNoSubcycleOld::postInitialize");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::postInitialize" << endl;
    }
}
/*********************/
void
EBAMRNoSubcycleOld::useFixedDt(Real a_dt)
{
  CH_TIME("EBAMRNoSubcycleOld::useFixedDt");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::useFixedDt" << endl;
    }
  m_dt = a_dt;
  m_prescribedDt = a_dt;
  m_useFixedDt = true;
}
/*********************/
Real
EBAMRNoSubcycleOld::computeDt()
{

  CH_TIMERS("EBAMRNoSubcycleOld::computeDt");
  CH_TIMER("computeDt_fillMask",t1);
  CH_TIMER("computeDt_allReduce",t2);
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRNoSubcycleOld::computeDt " << endl;
    }

  Real dt = m_params.m_maxDt;

  // if prescribedDt > 0 then use that
  // will take care of that
  if (m_useFixedDt)
    {
      dt = m_prescribedDt;
    }
  else
    {
      int numLevels = m_finestLevel+1;
      Real maxVel = 0.0;
      for (int ilev=0; ilev < numLevels; ilev++)
        {
          LevelData<EBCellFAB>& velNewLD = *m_veloNew[ilev];

          for(DataIterator dit = velNewLD.dataIterator(); dit.ok(); ++dit)
            {
              BaseFab<Real>& velFAB = (velNewLD[dit()]).getSingleValuedFAB();
              int iRegIrregCovered;
              CH_START(t1);
              const BaseFab<int>& maskFAB = m_ebisl[ilev][dit()].getEBGraph().getMask(iRegIrregCovered);
              CH_STOP(t1);
              if(iRegIrregCovered != -1)//not all covered
                {
                  if(iRegIrregCovered == 0)//has irreg
                    {
                      Box box = m_grids[ilev].get(dit());
                      for(int idir = 0; idir < SpaceDim; idir++)
                        {
                          FORT_MAXNORMMASK(CHF_REAL(maxVel),
                                           CHF_CONST_FRA1(velFAB,idir),
                                           CHF_BOX(box),
                                           CHF_CONST_FRA1(maskFAB,0));
                        }
                      IntVectSet ivs = m_ebisl[ilev][dit()].getMultiCells(box);
                      for(VoFIterator vofit(ivs, m_ebisl[ilev][dit()].getEBGraph()); vofit.ok(); ++vofit)
                        {
                          for(int idir = 0; idir < SpaceDim; idir++)
                            {
                              if(Abs(velNewLD[dit()](vofit(), idir)) > maxVel) {maxVel = Abs(velNewLD[dit()](vofit(), idir));}
                            }
                        }
                    }
                  else//all reg
                    {
                      Box box = m_grids[ilev].get(dit());
                      for(int idir = 0; idir < SpaceDim; idir++)
                        {
                          FORT_MAXNORM(CHF_REAL(maxVel),
                                       CHF_CONST_FRA1(velFAB,idir),
                                       CHF_BOX(box));
                        }
                    }
                }
            }
        }
      CH_START(t2);
#ifdef CH_MPI
      Real tmp = 1.;
      int result = MPI_Allreduce(&maxVel, &tmp, 1, MPI_CH_REAL,
                                 MPI_MAX, Chombo_MPI::comm);
      if (result != MPI_SUCCESS)
        { //bark!!!
          MayDay::Error("sorry, but I had a communcation error on norm");
        }
      maxVel = tmp;
#endif
      //        Real volume=1.;
      //        EBLevelDataOps::gatherBroadCast(maxVel, volume, 0);
      CH_STOP(t2);

      if(maxVel > 1.0e-10)
        {
          dt  = m_params.m_cfl*m_dx[m_finestLevel]/maxVel;
        }
      else
        {
          dt = m_params.m_maxDt;
        }
    } // end if no prescribed dt

  // finally, enforce max dt grow factor
  if ((m_dt > 0) && (m_params.m_maxDtGrow > 0) &&
       (dt > m_params.m_maxDtGrow*m_dt))
    {
      dt = m_params.m_maxDtGrow*m_dt;
    }

  return dt;
}
/*********************/
Real
EBAMRNoSubcycleOld::
computeInitialDt()
{
  CH_TIME("EBAMRNoSubcycleOld::computeInitialDt");
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRNoSubcycleOld::computeInitialDt " << endl;
    }
  Real retval;
  if(m_useFixedDt)
    {
      retval = m_prescribedDt;
    }
  else
    {
      retval =  (m_params.m_initCFL/m_params.m_cfl)*computeDt();
    }
  return retval;
}
/*********************/
void
EBAMRNoSubcycleOld::copyOldToNew(int a_ilev)
{
  CH_TIME("EBAMRNoSubcycleOld::copyOldToNew");
  Interval vecInterv(0, SpaceDim-1);
  Interval scaInterv(0, 0);
  m_veloOld[a_ilev]->copyTo(vecInterv, *m_veloNew[a_ilev], vecInterv);
  m_presOld[a_ilev]->copyTo(scaInterv, *m_presNew[a_ilev], scaInterv);
  m_gphiOld[a_ilev]->copyTo(vecInterv, *m_gphiNew[a_ilev], vecInterv);
}
/*********************/
void
EBAMRNoSubcycleOld::copyNewToOld(int a_ilev)
{
  CH_TIME("EBAMRNoSubcycleOld::copyNewToOld");
  Interval vecInterv(0, SpaceDim-1);
  Interval scaInterv(0, 0);
  m_veloNew[a_ilev]->copyTo(vecInterv, *m_veloOld[a_ilev], vecInterv);
  m_presNew[a_ilev]->copyTo(scaInterv, *m_presOld[a_ilev], scaInterv);
  m_gphiNew[a_ilev]->copyTo(vecInterv, *m_gphiOld[a_ilev], vecInterv);
}
/*********************/
void
EBAMRNoSubcycleOld::
regrid(const Vector<Vector<Box> >& a_newGrids)
{
  CH_TIME("EBAMRNoSubcycleOld::regrid");

  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRNoSubcycleOld::regrid " << endl;
    }
  //save state into old state
  swapOldAndNewStates();

  //this changes m_grids and m_ebisl
  defineGrids(a_newGrids);
  defineEBISLs();

  //this redefines new data with new set of grids
  //this can also change m_finestLevel
  defineNew();
  defineProjections();
  Interval vecInterv(0, SpaceDim-1);
  Interval scaInterv(0, 0);
  copyOldToNew(0);
  // now fill new data holders
  for (int ilev=1; ilev<= m_finestLevel; ilev++)
    {
      //interpolate everywhere
      EBPWLFineInterp ebInterpVec(m_grids[ ilev  ],
                                  m_grids[ ilev-1],
                                  m_ebisl[ ilev  ],
                                  m_ebisl[ ilev-1],
                                  m_domain[ilev-1],
                                  m_params.m_refRatio[ilev-1],
                                  SpaceDim);

      ebInterpVec.interpolate(*m_veloNew[ilev  ],
                              *m_veloNew[ilev-1],
                              vecInterv);

      //copy from old layout's data
      copyOldToNew(ilev);
    } // end loop over levels

  //no average down here because it was expected to have been done previously

  //now redefine old data and copy the new data over
  defineOld();
  for (int ilev=0; ilev<= m_finestLevel; ilev++)
    {
      copyNewToOld(ilev);
    } // end loop over levels

  //swap the pointers back so that confusion does not erupt
  swapOldAndNewStates();
}
/**********/
void
EBAMRNoSubcycleOld::swapOldAndNewStates()
{
  CH_TIME("EBAMRNoSubcycleOld::swapOldAndNewStates");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::swapOldAndNewStates" << endl;
    }
  LevelData<EBCellFAB>* tempLDPtr;

  for (int ilev=0; ilev<= m_finestLevel; ilev++)
    {
      tempLDPtr = m_veloOld[ilev];
      m_veloOld[ilev] = m_veloNew[ilev];
      m_veloNew[ilev] = tempLDPtr;


      tempLDPtr = m_gphiOld[ilev];
      m_gphiOld[ilev] = m_gphiNew[ilev];
      m_gphiNew[ilev] = tempLDPtr;


      tempLDPtr = m_presOld[ilev];
      m_presOld[ilev] = m_presNew[ilev];
      m_presNew[ilev] = tempLDPtr;

    } // end loop over levels
}
/*****************/
void
EBAMRNoSubcycleOld::postTimeStep()
{
  CH_TIME("EBAMRNoSubcycleOld::postTimeStep");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::postTimeStep" << endl;
    }
  // do average down here
  averageDown(m_veloNew);
  averageDown(m_presNew);
  averageDown(m_gphiNew);
  averageDown(m_veloOld);
  averageDown(m_presOld);
  averageDown(m_gphiOld);

}
/*****************/
void
EBAMRNoSubcycleOld::
computeVorticity(LevelData<EBCellFAB>& a_vort, int a_level)
{
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRStokes::computeVorticity" << endl;
    }
  CH_assert(m_isSetup);
  EBCellFactory ebcellfact(m_ebisl[a_level]);
  //define the data holder
  //vorticity is a vector in 3d, scalar in 2d
  int ncomp = 1;
  if(SpaceDim==3)
    {
      ncomp = SpaceDim;
    }
  a_vort.define(m_grids[a_level], ncomp, IntVect::Zero, ebcellfact);

  //interpolate velocity at coarse-fine interfaces
  //need some scratch space to do quadcfi
  for (int ilev=0; ilev < m_cellScratch.size(); ilev++)
    {
      if(m_cellScratch[ilev] != NULL)
        {
          delete m_cellScratch[ilev];
          m_cellScratch[ilev] = NULL;
        }
    }

  int numLevels = m_finestLevel+1;
  m_cellScratch.resize(numLevels, NULL);
  for (int ilev=0; ilev < numLevels; ilev++)
    {
      EBCellFactory ebcellfact(m_ebisl[ilev]);
      m_cellScratch[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], 1, 4*IntVect::Unit, ebcellfact);
    }

  if(a_level > 0)
    {

      LevelData<EBCellFAB>& phiCoar = *m_cellScratch[a_level-1];
      LevelData<EBCellFAB>& phiFine = *m_cellScratch[a_level  ];
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          Interval phiInterv(0, 0);
          Interval velInterv(idir, idir);
          m_veloNew[a_level-1]->copyTo(velInterv, phiCoar, phiInterv);
          m_veloNew[a_level  ]->copyTo(velInterv, phiFine, phiInterv);
          m_quadCFI[a_level]->interpolate(phiFine,
                                          phiCoar,
                                          phiInterv);

          //on copy back, we need ghost cells, so do the data iterator loop
          for(DataIterator dit = phiFine.dataIterator(); dit.ok(); ++dit)
            {
              Box region = phiFine[dit()].getRegion(); //includes ghost cells
              (*m_veloNew[a_level])[dit()].copy(region, velInterv, region, phiFine[dit()], phiInterv);
            }

          EBLevelDataOps::setVal(phiFine, 0.0);
          EBLevelDataOps::setVal(phiCoar, 0.0);
        }

    }
  m_veloNew[a_level]->exchange(Interval(0, SpaceDim-1));

  //clean up temp space
  for (int ilev=0; ilev < m_cellScratch.size(); ilev++)
    {
      if(m_cellScratch[ilev] != NULL)
        {
          delete m_cellScratch[ilev];
          m_cellScratch[ilev] = NULL;
        }
    }

  Real dxLev = m_dx[a_level];
  //do actual computation
  for(DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB&  vortFAB =                a_vort[dit()];
      EBCellFAB&  veloFAB = (*m_veloNew[a_level])[dit()];
      const EBGraph& ebgraph =   m_ebisl[a_level][dit()].getEBGraph();
      vortFAB.setVal(0.);

      BaseFab<Real>& regVort = vortFAB.getSingleValuedFAB();
      BaseFab<Real>& regVelo = veloFAB.getSingleValuedFAB();
      Box interiorBox = m_grids[a_level].get(dit());
      interiorBox.grow(1);
      interiorBox &= m_domain[a_level];
      interiorBox.grow(-1);
      int vortIndex;
      //some shucking and jiving to get around the vector/scalar
      //vorticity thing
#if CH_SPACEDIM==3
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          vortIndex = idir;
#else
          vortIndex = 0;
          int idir = 3;
#endif

          //compute on all cells as regular
          FORT_COMPUTEVORT(CHF_FRA1(regVort, vortIndex),
                           CHF_CONST_FRA(regVelo),
                           CHF_BOX(interiorBox),
                           CHF_CONST_REAL(dxLev),
                           CHF_CONST_INT(idir));;

          //do cells on domain boundary as if they were irregular
          IntVectSet ivsIrreg(m_grids[a_level].get(dit()));
          ivsIrreg -= interiorBox;
          ivsIrreg |= ebgraph.getIrregCells(interiorBox);
          int diffDirVec[2];
          if(SpaceDim==2 || idir ==2)
            {
              diffDirVec[0] = 0;
              diffDirVec[1] = 1;
            }
          else if(idir == 0)
            {
              diffDirVec[0] = 1;
              diffDirVec[1] = 2;
            }
          else if(idir==1)
            {
              diffDirVec[0] = 2;
              diffDirVec[1] = 0;
            }
          else
            {
              MayDay::Error("missed a case");
            }

          for(VoFIterator vofit(ivsIrreg, ebgraph); vofit.ok(); ++vofit)
            {
              const VolIndex vof = vofit();
              Real vortValue = 0.0;

              //taking the derivative d(udvdir)/d(xdiffdir)
              for(int idiff = 0; idiff < 2; idiff++)
                {

                  Real signDiff = 0;
                  int diffDir= -1;
                  int velComp = -1;
                  if(idiff == 0)
                    {
                      signDiff = 1.0;
                      diffDir  = diffDirVec[0];
                      velComp  = diffDirVec[1];
                    }
                  else if(idiff == 1)
                    {
                      signDiff = -1.0;
                      diffDir  = diffDirVec[1];
                      velComp  = diffDirVec[0];
                    }
                  else
                    {
                      MayDay::Error("missed a case");
                    }

                  Vector<FaceIndex> hiFaces =  ebgraph.getFaces(vof, diffDir, Side::Hi);
                  Vector<FaceIndex> loFaces =  ebgraph.getFaces(vof, diffDir, Side::Lo);

                  bool hasHi = (hiFaces.size() == 1) && (!hiFaces[0].isBoundary());
                  bool hasLo = (loFaces.size() == 1) && (!loFaces[0].isBoundary());
                  Real diffValue = 0.0;
                  if(hasHi && hasLo)
                    {
                      const VolIndex& vofHi = hiFaces[0].getVoF(Side::Hi);
                      const VolIndex& vofLo = loFaces[0].getVoF(Side::Lo);
                      Real hiValue = veloFAB(vofHi, velComp);
                      Real loValue = veloFAB(vofLo, velComp);
                      diffValue =0.5*(hiValue - loValue)/dxLev;
                    }
                  else if(hasHi)
                    {
                      const VolIndex& vofHi = hiFaces[0].getVoF(Side::Hi);
                      Real hiValue = veloFAB(vofHi, velComp);
                      Real loValue = veloFAB(vof,   velComp);
                      diffValue =(hiValue - loValue)/dxLev;
                    }
                  else if(hasLo)
                    {
                      const VolIndex& vofLo = loFaces[0].getVoF(Side::Lo);
                      Real hiValue = veloFAB(vof  , velComp);
                      Real loValue = veloFAB(vofLo, velComp);
                      diffValue =  (hiValue - loValue)/dxLev;
                    }
                  else
                    {
                      diffValue = 0.0;
                    }
                  vortValue += signDiff*diffValue;
                } //end loop over idiff

              vortFAB(vof, vortIndex) = vortValue;

            } //end loop over irregular vofs
#if CH_SPACEDIM==3
        } //end loop over vort components in 3d
#endif

    }
}
/**********/
void
EBAMRNoSubcycleOld::
tagCellsLevel(IntVectSet& a_tags, int a_level)
{
  CH_TIME("EBAMRNoSubcycleOld::tagCellsLevel");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::tagCellsLevel" << endl;
    }
  CH_assert(m_isSetup);
  LevelData<EBCellFAB> vort;
  computeVorticity(vort, a_level);
  a_tags.makeEmpty();

  for(DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& vortFAB = vort[dit()];
      const Box& grid =        m_grids[a_level].get(dit());
      const EBGraph& ebgraph = m_ebisl[a_level][dit()].getEBGraph();
      IntVectSet ivsTot(grid);
      Box shrunkDom = m_domain[a_level].domainBox();
      int shrinkNumCells = m_params.m_tagShrinkDomain;
      int shrinkRefRatio = 1;
      for(int ilev = 1; ilev <= a_level; ilev++)
        {
          shrinkRefRatio *= m_params.m_refRatio[ilev-1];
        }
      for(int idir = 0;  idir < SpaceDim; idir++)
        {
          if(idir != m_params.m_flowDir)
            {
              shrunkDom.grow(idir, -shrinkNumCells*shrinkRefRatio);
            }
        }
      ivsTot &= shrunkDom;

      for(VoFIterator vofit(ivsTot, ebgraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect& iv = vof.gridIndex();
          Real vortmag = 0.0;

          for(int icomp = 0; icomp < vortFAB.nComp(); icomp++)
            {
              Real vortDirVal = vortFAB(vof, icomp);
              vortmag += vortDirVal*vortDirVal;
            }
          vortmag = sqrt(vortmag);
          if(vortmag >= m_params.m_refineThreshold)
            {
              a_tags |= iv;
            }
        } //end loop over vofs

      //refine all irregular cells
      IntVectSet irregIVS = ebgraph.getIrregCells(grid);
      a_tags |= irregIVS;
    }

  a_tags.grow(m_params.m_tagBuffer);

}
/*****************/
void
EBAMRNoSubcycleOld::
defineIrregularData()
{
  CH_TIME("EBAMRNoSubcycleOld::defineIrregularData");
  for(int ilev = 0; ilev <=  m_finestLevel; ilev++)
    {
      m_coveredFaceLitLo[ilev]->define(m_grids[ilev]);
      m_coveredFaceLitHi[ilev]->define(m_grids[ilev]);
      m_coveredSetsLitLo[ilev]->define(m_grids[ilev]);
      m_coveredSetsLitHi[ilev]->define(m_grids[ilev]);
      for(DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          Box litBox = m_grids[ilev].get(dit());
          litBox.grow(1);
          litBox &= m_domain[ilev];
          (*m_coveredFaceLitLo[ilev])[dit()].resize(SpaceDim);
          (*m_coveredFaceLitHi[ilev])[dit()].resize(SpaceDim);
          (*m_coveredSetsLitLo[ilev])[dit()].resize(SpaceDim);
          (*m_coveredSetsLitHi[ilev])[dit()].resize(SpaceDim);

          for(int idir = 0; idir < SpaceDim; idir++)
            {

              IntVectSet irregIVSPlus, irregIVSMinu;
              //get the covered sets and faces
              EBArith::computeCoveredFaces((*m_coveredFaceLitHi[ilev])[dit()][idir],
                                           (*m_coveredSetsLitHi[ilev])[dit()][idir],
                                           irregIVSPlus,idir, Side::Hi, m_ebisl[ilev][dit()], litBox);
              EBArith::computeCoveredFaces((*m_coveredFaceLitLo[ilev])[dit()][idir],
                                           (*m_coveredSetsLitLo[ilev])[dit()][idir],
                                           irregIVSMinu, idir, Side::Lo,  m_ebisl[ilev][dit()], litBox);

            }
        }

    }//end loop over levels
}
/*****************/
void
EBAMRNoSubcycleOld::
allocateTemporaries()
{
  CH_TIME("EBAMRNoSubcycleOld::allocateTemporaries");
  // number of active levels
  int numLevels = m_finestLevel+1;
  m_advVel.resize(numLevels, NULL);
  m_coveredAdvVelLo.resize(numLevels, NULL);
  m_coveredAdvVelHi.resize(numLevels, NULL);
  m_uDotDelU.resize(numLevels, NULL);
  m_uDotDelS.resize(numLevels, NULL);

  m_cellScratch.resize(numLevels, NULL);
  m_cellScratc1.resize(numLevels, NULL);
  m_cellScratc2.resize(numLevels, NULL);
  m_cellScratc3.resize(numLevels, NULL);
  m_macGradient.resize(numLevels, NULL);
  m_macScratch1.resize(numLevels, NULL);
  m_macScratch2.resize(numLevels, NULL);
  m_coveredScratchLo.resize(numLevels, NULL);
  m_coveredScratchHi.resize(numLevels, NULL);

  // allocate storage
  for (int ilev=0; ilev < numLevels; ilev++)
    {
      EBFluxFactory ebfluxfact(m_ebisl[ilev]);
      EBCellFactory ebcellfact(m_ebisl[ilev]);

      m_uDotDelU[ilev]    = new LevelData<EBCellFAB>(m_grids[ilev], SpaceDim, 4*IntVect::Unit, ebcellfact);
      m_uDotDelS[ilev]    = new LevelData<EBCellFAB>(m_grids[ilev],        1, 4*IntVect::Unit, ebcellfact);
      m_cellScratch[ilev] = new LevelData<EBCellFAB>(m_grids[ilev],        1, 4*IntVect::Unit, ebcellfact);
      m_cellScratc1[ilev] = new LevelData<EBCellFAB>(m_grids[ilev],        1, 4*IntVect::Unit, ebcellfact);
      m_cellScratc2[ilev] = new LevelData<EBCellFAB>(m_grids[ilev],        1, 4*IntVect::Unit, ebcellfact);
      m_cellScratc3[ilev] = new LevelData<EBCellFAB>(m_grids[ilev],        1, 4*IntVect::Unit, ebcellfact);


      m_advVel[ilev]      = new LevelData<EBFluxFAB>(m_grids[ilev],        1, 4*IntVect::Unit, ebfluxfact);
      m_macGradient[ilev] = new LevelData<EBFluxFAB>(m_grids[ilev],        1, 4*IntVect::Unit, ebfluxfact);
      m_macScratch1[ilev] = new LevelData<EBFluxFAB>(m_grids[ilev],        1, 4*IntVect::Unit, ebfluxfact);
      m_macScratch2[ilev] = new LevelData<EBFluxFAB>(m_grids[ilev],        1, 4*IntVect::Unit, ebfluxfact);

      m_coveredAdvVelLo[ilev]  = new LayoutData<Vector<BaseIVFAB<Real> * > >(m_grids[ilev]);
      m_coveredAdvVelHi[ilev]  = new LayoutData<Vector<BaseIVFAB<Real> * > >(m_grids[ilev]);
      m_coveredScratchLo[ilev] = new LayoutData<Vector<BaseIVFAB<Real> * > >(m_grids[ilev]);
      m_coveredScratchHi[ilev] = new LayoutData<Vector<BaseIVFAB<Real> * > >(m_grids[ilev]);

      for(DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*m_coveredAdvVelLo[ilev])[dit()].resize(SpaceDim, NULL);
          (*m_coveredAdvVelHi[ilev])[dit()].resize(SpaceDim, NULL);
          (*m_coveredScratchLo[ilev])[dit()].resize(SpaceDim, NULL);
          (*m_coveredScratchHi[ilev])[dit()].resize(SpaceDim, NULL);

          const EBGraph& ebgraph = m_ebisl[ilev][dit()].getEBGraph();
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              (*m_coveredAdvVelLo[ilev])[dit()][idir]  = new BaseIVFAB<Real>((*m_coveredSetsLitLo[ilev])[dit()][idir], ebgraph, 1);
              (*m_coveredAdvVelHi[ilev])[dit()][idir]  = new BaseIVFAB<Real>((*m_coveredSetsLitHi[ilev])[dit()][idir], ebgraph, 1);
              (*m_coveredScratchLo[ilev])[dit()][idir] = new BaseIVFAB<Real>((*m_coveredSetsLitLo[ilev])[dit()][idir], ebgraph, 1);
              (*m_coveredScratchHi[ilev])[dit()][idir] = new BaseIVFAB<Real>((*m_coveredSetsLitHi[ilev])[dit()][idir], ebgraph, 1);
            }
        }
    }
  EBAMRDataOps::setToZero(m_uDotDelU  );
  EBAMRDataOps::setToZero(m_uDotDelS  );
  EBAMRDataOps::setToZero(m_cellScratch);
  EBAMRDataOps::setToZero(m_cellScratc1);
  EBAMRDataOps::setToZero(m_cellScratc2);
  EBAMRDataOps::setToZero(m_cellScratc3);

}
/*****************/
void
EBAMRNoSubcycleOld::
deleteTemporaries()
{
  CH_TIME("EBAMRNoSubcycleOld::deleteTemporaries");
  // number of active levels
  int numLevels = m_finestLevel+1;
  // clean up storage
  for (int ilev=0; ilev<numLevels; ilev++)
    {
      delete m_advVel[ilev];
      delete m_uDotDelU[ilev];
      delete m_uDotDelS[ilev];
      delete m_cellScratch[ilev];
      delete m_cellScratc1[ilev];
      delete m_cellScratc2[ilev];
      delete m_cellScratc3[ilev];
      delete m_macGradient[ilev];
      delete m_macScratch1[ilev];
      delete m_macScratch2[ilev];
      for(DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              delete (*m_coveredAdvVelLo[ilev])[dit()][idir];
              delete (*m_coveredAdvVelHi[ilev])[dit()][idir];
              delete (*m_coveredScratchLo[ilev])[dit()][idir];
              delete (*m_coveredScratchHi[ilev])[dit()][idir];
            }
        }
      delete m_coveredAdvVelLo[ilev];
      delete m_coveredAdvVelHi[ilev];
      delete m_coveredScratchLo[ilev];
      delete m_coveredScratchHi[ilev];
    }

  for (int ilev=0; ilev<numLevels; ilev++)
    {
      m_advVel[ilev]           = NULL;
      m_uDotDelU[ilev]         = NULL;
      m_uDotDelS[ilev]         = NULL;
      m_cellScratch[ilev]      = NULL;
      m_cellScratc1[ilev]      = NULL;
      m_cellScratc2[ilev]      = NULL;
      m_cellScratc3[ilev]      = NULL;
      m_macGradient[ilev]      = NULL;
      m_macScratch1[ilev]      = NULL;
      m_macScratch2[ilev]      = NULL;
      m_coveredAdvVelLo[ilev]  = NULL;
      m_coveredAdvVelHi[ilev]  = NULL;
      m_coveredScratchLo[ilev] = NULL;
      m_coveredScratchHi[ilev] = NULL;
    }
}
/*****************/
void
EBAMRNoSubcycleOld::
advance()
{

  CH_TIME("EBAMRNoSubcycleOld::advance");
  if (m_params.m_verbosity > 1)
    {
      pout() << "EBAMRNoSubcycleOld::advance, step " << m_curStep
             << ", starting time = "  << m_time
             << ", dt = " << m_dt << endl;
    }

  //allocate space for integration
  allocateTemporaries();

  if(!m_params.m_stokesFlow)
    {
      pout() << "computing advection velocities" << endl;
      computeAdvectionVelocities(m_advVel, m_coveredAdvVelLo, m_coveredAdvVelHi, m_veloOld);
      pout() << "computing udelu" << endl;
      computeAdvectiveDerivative(m_uDotDelU, m_veloOld, true);
      pout() << "computing udels" << endl; //left this in so the printouts will not change
    }
  else
    {
      pout() << "EBAMRNoSubcycleOld::advance: <<<Stokes flow---> udelu==udels == 0>>>" << endl;
      EBAMRDataOps::setToZero(m_uDotDelU);
      EBAMRDataOps::setToZero(m_uDotDelS);
    }

  pout() << "advancing solution" << endl;
  advanceVelocity();

  EBAMRDataOps::setCoveredVal(m_veloNew, 0.0);
  EBAMRDataOps::setCoveredVal(m_gphiNew, 0.0);
  EBAMRDataOps::setCoveredVal(m_presNew, 0.0);
  //copy new stuff to old stuff --- false for truncation error tests
  if(m_params.m_copyOverOld  || m_advanceGphiOnly)
    {
      for(int ilev = 0; ilev <= m_finestLevel; ilev++)
        {
          copyNewToOld(ilev);
        }
    }

  pointsUpdated();

  //remove all unnecessary data
  deleteTemporaries();

}

/******************/
void
EBAMRNoSubcycleOld::
pointsUpdated()
{
  CH_TIME("EBAMRNoSubcycleOld::pointsUpdated");
  long long totalPoints = 0;
  long long totalBoxes  = 0;
  int numLevels = m_finestLevel + 1;
  for(int ilev = 0; ilev < numLevels; ilev++)
    {
      long long pointsThisLevel = 0;
      for(LayoutIterator lit = m_grids[ilev].layoutIterator(); lit.ok(); ++lit)
        {
          pointsThisLevel += m_grids[ilev][lit()].numPts();
        }
      totalPoints += pointsThisLevel;
      totalBoxes += m_grids[ilev].size();
    }

  m_pointsUpdated += totalPoints;
}

/*****/
void
EBAMRNoSubcycleOld::
extrapolateScalarCol(Vector<LevelData<EBFluxFAB>* >                     &  a_macScalar,
                     Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacLo,
                     Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacHi,
                     const RefCountedPtr<EBPhysIBCFactory>              &  a_advectBC,
                     const Vector<LevelData<EBFluxFAB>* >               &  a_advectiveVel,
                     const Vector<LevelData<EBCellFAB>*>*                  a_sourceTerm,
                     const Vector<LevelData<EBCellFAB>* >               &  a_cellScalar,
                     const Vector<LevelData<EBCellFAB>* >               &  a_cellVelocity)
{
  CH_TIME("EBAMRNoSubcycleOld::extrapolateScalarCol");

  //extrapolate velocity component to faces
  for(int ilev = 0; ilev <= m_finestLevel; ilev++)
    {

      EBPatchGodunov::setCurLevel(ilev);
      LevelData<EBCellFAB>* coarDataOld = NULL;
      LevelData<EBCellFAB>* coarDataNew = NULL;

      LevelData<EBCellFAB>* coarVeloOld = NULL;
      LevelData<EBCellFAB>* coarVeloNew = NULL;

      DisjointBoxLayout coarDBL;
      EBISLayout        coarEBISL;
      int refRat = 2;
      if(ilev > 0)
        {
          coarDBL =  m_grids[ilev-1];
          coarEBISL= m_ebisl[ilev-1];
          coarDataOld = a_cellScalar[ilev-1];
          coarDataNew = a_cellScalar[ilev-1];
          coarVeloOld = a_cellVelocity[ilev-1];
          coarVeloNew = a_cellVelocity[ilev-1];
          refRat = m_params.m_refRatio[ilev-1];
        }

      EBLevelAdvect& ebLevelAdvect = *m_ebLevAd[ilev];
      //need to reset boundary conditions because this object was defined with dummy bcs
      //because it is reused over several different variables

      ebLevelAdvect.resetBCs(a_advectBC);

      //this advects all SpaceDim components of velocity to
      //each face.  put into macscratch.
      LevelData<EBCellFAB>* source = NULL;
      if(a_sourceTerm != NULL)
        {
          source = (*a_sourceTerm)[ilev];
        }
      ebLevelAdvect.advectToFacesCol(*a_macScalar[ilev],  //extrapolated componenent idir
                                     *a_coveredMacLo[ilev],
                                     *a_coveredMacHi[ilev],
                                     *m_coveredFaceLitLo[ilev],
                                     *m_coveredFaceLitHi[ilev],
                                     *m_coveredSetsLitLo[ilev],
                                     *m_coveredSetsLitHi[ilev],
                                     *a_cellScalar[ilev],    //consstate
                                     *a_cellVelocity[ilev],        //use veloOld as normal velocity
                                     *a_advectiveVel[ilev],   //contains initial advective velocity
                                     coarDataOld,
                                     coarDataNew,
                                     coarVeloOld,
                                     coarVeloNew,
                                     m_time, m_time, m_time, m_dt, source);

      a_macScalar[ilev]->exchange(Interval(0,0));
    }

}

/*****************/
void
EBAMRNoSubcycleOld::
extrapolateToCoveredFaces(Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacLo,
                          Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacHi,
                          Vector<LevelData<EBFluxFAB>* >&                       a_macOpen,
                          Vector<LevelData<EBCellFAB>* >&                       a_cellOpen,
                          int                                                   a_idir)
{
  CH_TIME("EBAMRNoSubcycleOld::extrapolateToCoveredFaces");
  //extrapolate the projected velocity to covered faces
  for(int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      ProblemDomain curDomain(m_domain[ilev]);

      for(DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {

          EBPatchAdvect& patcher  = m_ebLevAd[ilev]->getPatchAdvect(dit());
          patcher.extrapToCoveredFaces(*(*a_coveredMacLo[ilev])[dit()][a_idir],
                                       (*a_macOpen[ilev])[dit()][a_idir],
                                       (*a_cellOpen[ilev])[dit()],
                                       (*m_coveredFaceLitLo[ilev])[dit()][a_idir],
                                       a_idir, Side::Lo, m_grids[ilev].get(dit()));

          patcher.extrapToCoveredFaces(*(*a_coveredMacHi[ilev])[dit()][a_idir],
                                       (*a_macOpen[ilev])[dit()][a_idir],
                                       (*a_cellOpen[ilev])[dit()],
                                       (*m_coveredFaceLitHi[ilev])[dit()][a_idir],
                                       a_idir, Side::Hi, m_grids[ilev].get(dit()));
        }
    }
}
/*****************/
void
EBAMRNoSubcycleOld::
viscousSourceForAdvect(Vector<LevelData<EBCellFAB>* >&       a_source,
                       Vector<LevelData<EBCellFAB>* >&       a_velComp,
                       Vector<LevelData<EBCellFAB>* >&       a_zero,
                       int                                   a_icomp)
{
  CH_TIME("EBAMRNoSubcycleOld::viscousSourceForAdvect");
  EBAMRDataOps::setToZero(a_source);
  EBAMRDataOps::setToZero(a_zero);
  CH_assert(a_velComp[0]->nComp() == 1);

  applyEBAMROp(a_source,   //returns holding lapl(vel comp)
               a_velComp,   //holds cell centered vel comp
               a_zero, //holds the zero for the residual calc
               a_icomp);  //velocity component

  //make cellscratch2 hold nu*lapl.
// viscosity already lives in operator
//   EBAMRDataOps::scale(a_source, m_viscosity);

  //fill ghost cells over coarse fine interface

  //fill fine-fine ghost cells
  EBAMRDataOps::exchangeAll(a_source);

  //now fill the ghost cells of the laplacian over coarse-fine interfaces with
  //constant extrapolation from neighboring valid values
  //this includes an exchange
  for(int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      IntVect ivGhost = a_source[ilev]->ghostVect();
      EBConstantCFInterp interpolator(m_grids[ilev], m_ebisl[ilev], m_domain[ilev], ivGhost);
      interpolator.interpolate(*a_source[ilev]);
    }
}
/****/
void
EBAMRNoSubcycleOld::
applyEBAMROp(Vector<LevelData<EBCellFAB>* >&       a_lap,
             Vector<LevelData<EBCellFAB>* >&       a_phi,
             Vector<LevelData<EBCellFAB>* >&       a_zero,
             int                                   a_velComp)
{
  CH_TIME("EBAMRNoSubcycleOld::applyEBAMROp");
  Real alpha = 0.0;
  Real beta = 1;
  m_tgaSolver[a_velComp]->resetAlphaAndBeta(alpha, beta);
  int coarsestLevel = 0;
  bool homogeneousBC = false;
  //apply the operator (by computing the residual with rhs = 0, and * -1
  m_solver[a_velComp]->computeAMRResidual(a_lap,
                                          a_phi,
                                          a_zero,
                                          m_finestLevel,
                                          coarsestLevel,
                                          homogeneousBC);

  EBAMRDataOps::scale(a_lap,-1.0);
  EBAMRDataOps::setCoveredVal(a_lap,0.0);
  EBAMRDataOps::setCoveredAMRVal(a_lap,m_ebisl,m_params.m_refRatio,0.0);

}
/*****************/
void
EBAMRNoSubcycleOld::
computeAdvectionVelocities(Vector<LevelData<EBFluxFAB> *>&                       a_advVel,
                           Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelLo,
                           Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelHi,
                           const Vector<LevelData<EBCellFAB> *>&                 a_veloOld)
{
  CH_TIMERS("EBAMRNoSubcycleOld::computeAdvectionVelocites");
  CH_TIMER("extrapolation_to_faces", t1);
  CH_TIMER("mac_projection", t2);
  CH_TIMER("post_projection", t3);

  CH_START(t1);
  for(int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      //initally fill advective velocity with average to faces of
      //cell-centered solution
      RealVect dxLev = m_dx[ilev]*RealVect::Unit;
      LayoutData<IntVectSet> cfivs;
      ccpAverageVelocityToFaces(*m_macScratch1[ilev], *a_veloOld[ilev],
                                m_grids[ilev], m_ebisl[ilev], m_domain[ilev], dxLev,
                                cfivs);
    }
  //extrapolate to  get advective velocity at covered faces
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      for(int ilev = 0; ilev <= m_finestLevel; ilev++)
        {
          //copy cell centered velocity into scratch space for each direction
          Interval srcInterv(idir, idir);
          Interval dstInterv(0 , 0);
          a_veloOld[ilev]->copyTo(srcInterv, *m_cellScratch[ilev], dstInterv);
        }
      extrapolateToCoveredFaces(m_coveredScratchLo,
                                m_coveredScratchHi,
                                m_macScratch1,
                                m_cellScratch, idir);
    }


  //extrapolate velocities to faces
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      for(int ilev = 0; ilev <= m_finestLevel; ilev++)
        {
          //copy cell centered velocity into scratch space for each direction
          Interval srcInterv(idir, idir);
          Interval dstInterv(0 , 0);
          EBLevelDataOps::setVal(*m_cellScratch[ilev], 0.0);
          a_veloOld[ilev]->copyTo(srcInterv, *m_cellScratch[ilev], dstInterv);
        }

      //fill in viscous source term where necessary
      Vector<LevelData<EBCellFAB>*>* source = NULL;
      if(m_viscousCalc)
        {
          //cellscratch already holds velocity component
          //make cellscratch3 = zero

          DirichletPoissonEBBC::s_velComp = idir;
          viscousSourceForAdvect(m_cellScratc2,   //returns holding source term = nu*lapl - gphiold
                                 m_cellScratch,   //holds cell centered vel comp n
                                 m_cellScratc3,   //will hold the zero for the residual calc (zeroed inside routine)
                                 idir);  //velocity component

          source = &m_cellScratc2;
        }

      pout() << "computeAdvectionVelocities: doing velocity component " << idir << endl;
      EBPatchGodunov::setCurComp(idir);
      EBPatchGodunov::setDoingVel(1);
      EBPatchGodunov::setDoingAdvVel(1);

      RefCountedPtr<EBPhysIBCFactory> advectBC = m_ibc->getVelAdvectBC(idir);
      //cellscratchsca used for consstate
      //extrapolate velocity component to faces
      extrapolateScalarCol(m_macScratch2,
                           m_coveredScratchLo,
                           m_coveredScratchHi,
                           advectBC,
                           m_macScratch1, //contains initial advective vel
                           source,
                           m_cellScratch,
                           a_veloOld);

      //now copy the result to the appropriate faces in advVel. (advVel only has normal velocities)
      for(int ilev = 0; ilev <= m_finestLevel; ilev++)
        {
          for(DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
            {
              Interval interv(0, 0);
              EBFaceFAB&  extrapFAB = (*m_macScratch2[ilev])[dit()][idir];
              EBFaceFAB&  advVelFAB =  (*a_advVel[ilev])[dit()][idir];
              Box region = extrapFAB.getCellRegion();
              advVelFAB.copy(region, interv, region, extrapFAB, interv);

              //same goes for covered velocites.  only normal ones get into covered adv vel.
              (*a_coveredAdvVelLo[ilev])[dit()][idir]->copy(region, interv, region, *(*m_coveredScratchLo[ilev])[dit()][idir], interv);
              (*a_coveredAdvVelHi[ilev])[dit()][idir]->copy(region, interv, region, *(*m_coveredScratchHi[ilev])[dit()][idir], interv);
            }

        }
    } //end loop over velocity directions

  CH_STOP(t1);

  //MAC project the velocity
  if (m_params.m_verbosity > 0)
    {
      pout() << "EBAMRNoSubcycleOld:mac projecting advection vel" << endl;
    }


  CH_START(t2);

  m_macProjector->project(a_advVel, m_macGradient);

  CH_STOP(t2);

  EBLevelMACProjector::setVerbose(false);
  if (m_params.m_verbosity > 3)
    {

      m_macProjector->kappaDivergence(m_cellScratch, a_advVel);
      averageDown(m_cellScratch);

      Real norm[3];
      for(int inorm = 0; inorm < 3; inorm++)
        {
          norm[inorm] = EBArith::norm(m_cellScratch, m_grids, m_ebisl,
                                      m_params.m_refRatio, 0, inorm, EBNormType::OverBoth);
        }
      pout() << setprecision(8)
             << setiosflags(ios::showpoint)
             << setiosflags(ios::scientific);
      pout() << "div(adv vel) after mac projection: " <<
        "L_inf = " << norm[0]  <<
        ", L_1 = " << norm[1] <<
        ", L_2 = " << norm[2] << endl;

    }

  CH_START(t3);
  //average down advection velocity so that it makes sense at coarse-fine interfaces
  averageDown(a_advVel     );
  averageDown(m_macGradient);

  for(int idir = 0; idir < SpaceDim; idir++)
    {
      //faceDir,  velcomp are the same thing for advection velocities
      int faceDir = idir;
      int velComp = idir;
      //correct covered velocity with extrapolation of pressure gradient
      m_macProjector->correctVelocityComponent(a_coveredAdvVelLo,
                                               a_coveredAdvVelHi,
                                               m_coveredFaceLitLo,
                                               m_coveredFaceLitHi,
                                               m_coveredSetsLitLo,
                                               m_coveredSetsLitHi,
                                               m_macGradient, faceDir, velComp);
    }
  CH_STOP(t3);
}
/*****************/
void
EBAMRNoSubcycleOld::
computeAdvectiveDerivative(Vector<LevelData<EBCellFAB>* >&    a_uDotDelS,
                           Vector<LevelData<EBCellFAB>* >&    a_scalOld,
                           bool                               a_reallyVelocity)
{
  CH_TIME("EBAMRNoSubcycleOld::computeAdvectiveDerivative");
  int ncomp = a_uDotDelS[0]->nComp();
  int nlevels = m_finestLevel+1;

  //make temporaries with right number of variables
  Vector<LevelData<EBFluxFAB>* >                       macScratchVec(nlevels, NULL);
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  coveredScratchVecLo(nlevels, NULL);
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  coveredScratchVecHi(nlevels, NULL);

  for(int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      EBFluxFactory ebfluxfact(m_ebisl[ilev]);
      macScratchVec[ilev]       = new LevelData<EBFluxFAB>(m_grids[ilev], ncomp, 4*IntVect::Unit, ebfluxfact);

      coveredScratchVecLo[ilev] = new LayoutData<Vector<BaseIVFAB<Real> * > >(m_grids[ilev]);
      coveredScratchVecHi[ilev] = new LayoutData<Vector<BaseIVFAB<Real> * > >(m_grids[ilev]);
      for(DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {

          (*coveredScratchVecLo[ilev])[dit()].resize(SpaceDim, NULL);
          (*coveredScratchVecHi[ilev])[dit()].resize(SpaceDim, NULL);

          const EBGraph& ebgraph = m_ebisl[ilev][dit()].getEBGraph();
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              (*coveredScratchVecLo[ilev])[dit()][idir] = new BaseIVFAB<Real>((*m_coveredSetsLitLo[ilev])[dit()][idir], ebgraph, ncomp);
              (*coveredScratchVecHi[ilev])[dit()][idir] = new BaseIVFAB<Real>((*m_coveredSetsLitHi[ilev])[dit()][idir], ebgraph, ncomp);
            }
        }
    }
  for(int icomp = 0; icomp < ncomp; icomp++)
    {
      EBPatchGodunov::setCurComp(icomp);
      EBPatchGodunov::setDoingAdvVel(0);

      RefCountedPtr<EBPhysIBCFactory> advectBC;
      if(a_reallyVelocity)
        {
          advectBC  = m_ibc->getVelAdvectBC(icomp);
          EBPatchGodunov::setDoingVel(1);
        }
      else
        {
          advectBC  = m_ibc->getScalarAdvectBC(icomp);
          EBPatchGodunov::setDoingVel(0);
        }

      for (int ilev=0; ilev <= m_finestLevel; ilev++)
        {
          Interval srcInterv(icomp, icomp);
          Interval dstInterv(0, 0);
          a_scalOld[ilev]->copyTo(srcInterv, *m_cellScratch[ilev], dstInterv);
        }

      //fill in viscous source term where necessary
      Vector<LevelData<EBCellFAB>* > * source = NULL;
      if(a_reallyVelocity && m_viscousCalc)
        {
          //cellscratch already holds velocity component
          //make cellscratch3 = zero

          viscousSourceForAdvect(m_cellScratc2,   //returns holding source term = nu*lapl - gphiold
                                 m_cellScratch,   //holds cell centered vel comp n
                                 m_cellScratc3,   //will hold the zero for the residual calc (zeroed inside routine)
                                 icomp);  //velocity component

          source = &m_cellScratc2;
        }

      //cellscratch used for consstate
      extrapolateScalarCol(m_macScratch1,
                           m_coveredScratchLo,
                           m_coveredScratchHi,
                           advectBC,
                           m_advVel,
                           source,
                           m_cellScratch,
                           m_veloOld);


      if(a_reallyVelocity)
        {
          //correct with previously stored pressure gradient
          m_macProjector->correctVelocityComponent(m_macScratch1, m_macGradient, icomp);
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              int faceDir = idir;
              int velComp = icomp;
              //correct covered velocity with extrapolation of pressure gradient
              m_macProjector->correctVelocityComponent(m_coveredScratchLo,
                                                       m_coveredScratchHi,
                                                       m_coveredFaceLitLo,
                                                       m_coveredFaceLitHi,
                                                       m_coveredSetsLitLo,
                                                       m_coveredSetsLitHi,
                                                       m_macGradient, faceDir, velComp);
            }

          //overwrite extrapolated with advective velocity if appropriate
          //(the normal velocity here had the wrong boundary conditions)
          for(int ilev = 0; ilev <= m_finestLevel; ilev++)
            {
              for(DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
                {
                  Box box = grow(m_grids[ilev].get(dit()), 1);
                  Interval interv(0,0);
                  //copy over covered advective velocity
                  (*m_coveredScratchLo[ilev])[dit()][icomp]->copy(box, interv, box, *(*m_coveredAdvVelLo[ilev])[dit()][icomp], interv);
                  (*m_coveredScratchHi[ilev])[dit()][icomp]->copy(box, interv, box, *(*m_coveredAdvVelHi[ilev])[dit()][icomp], interv);

                  EBFluxFAB& macExtrapFAB    = (*m_macScratch1[ilev])[dit()];
                  const EBFluxFAB& advVelFAB = (*m_advVel[ilev])[dit()];
                  //icomp is the the component of the velocity and
                  //therefore also the face for which it is the normal component
                  macExtrapFAB[icomp].copy(advVelFAB[icomp]);
                }
            }
        }

      //copy extrapolated values into vector holders
      for(int ilev = 0; ilev <= m_finestLevel; ilev++)
        {
          Interval srcInterv(0, 0);
          Interval dstInterv(icomp, icomp);
          m_macScratch1[ilev]->copyTo(srcInterv, *macScratchVec[ilev],  dstInterv);
          int ibox = 0;
          for(DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
            {
              Box box = grow(m_grids[ilev].get(dit()), 1);
              box &= m_domain[ilev];
              for(int idir = 0; idir < SpaceDim; idir++)
                {
                  (*coveredScratchVecLo[ilev])[dit()][idir]->copy(box, dstInterv, box, *(*m_coveredScratchLo[ilev])[dit()][idir], srcInterv);
                  (*coveredScratchVecHi[ilev])[dit()][idir]->copy(box, dstInterv, box, *(*m_coveredScratchHi[ilev])[dit()][idir], srcInterv);
                }

              ibox++;
            }
        }

    } //end loop over components

  //average down face centered stuff so that it makes sense at coarse-fine interfaces
  averageDown(macScratchVec);

  //compute the actual advective derivative
  //split this out to make it separately testable
  computeAdvectiveDerivative(a_uDotDelS, m_advVel, macScratchVec,
                             m_coveredAdvVelLo, m_coveredAdvVelHi,
                             coveredScratchVecLo , coveredScratchVecHi);

  for(int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      delete macScratchVec[ilev];
      for(DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              delete (*coveredScratchVecLo[ilev])[dit()][idir];
              delete (*coveredScratchVecHi[ilev])[dit()][idir];
            }
        }
      delete coveredScratchVecLo[ilev];
      delete coveredScratchVecHi[ilev];
    }
}
/*******************/
/*******************/
void
EBAMRNoSubcycleOld::
computeAdvectiveDerivative(Vector<LevelData<EBCellFAB>* >                     &  a_uDotDelS,
                           Vector<LevelData<EBFluxFAB>* >                     &  a_macAdvVel,
                           Vector<LevelData<EBFluxFAB>* >                     &  a_macScalar,
                           Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelLo,
                           Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelHi,
                           Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredScalarLo,
                           Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredScalarHi,
                           bool a_nonConsOnly,
                           bool a_consOnly)
{
  CH_TIME("EBAMRNoSubcycleOld::computeAdvectiveDerivative2");
  int ncomp = a_uDotDelS[0]->nComp();
  //compute the actual advective derivative
  for(int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      a_macScalar[ilev]->exchange(Interval(0,ncomp-1));
      a_macAdvVel[ilev]->exchange(Interval(0,0));
      //compute scalar advective derivative
      //then copy into the appropriate
      //component of the input data holder
      RealVect dxLev = m_dx[ilev]*RealVect::Unit;
      ProblemDomain domLev = m_domain[ilev];
      IntVectSet cfivs;
      for(DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          EBPatchAdvect& patcher = m_ebLevAd[ilev]->getPatchAdvect(dit());
          //EBCellFAB& udeluFAB = (*a_uDotDelS[ilev])[dit()];
          //this takes the non-conservative advective derivative = udotdels
          patcher.advectiveDerivative((*a_uDotDelS[ilev])[dit()],
                                      (*a_macScalar[ilev])[dit()],      //contains extrapolated vel component (facerho)
                                      (*a_macAdvVel[ilev])[dit()],           //facevel
                                      (*a_coveredScalarLo[ilev])[dit()], //coveredRhoMinu
                                      (*a_coveredScalarHi[ilev])[dit()], //coveredRhoPlus
                                      (*a_coveredAdvVelLo[ilev])[dit()],  //coveredVelMinu
                                      (*a_coveredAdvVelHi[ilev])[dit()],  //coveredVelPlus
                                      (*m_coveredFaceLitLo[ilev])[dit()], //coveredFaceMinu
                                      (*m_coveredFaceLitHi[ilev])[dit()], //coveredFacePlus
                                      m_grids[ilev].get(dit()) );


        }

      //these are not grown by one.
      LayoutData<IntVectSet> irregSetsSmall;
      //these are grown by one in the directions != idir
      LayoutData<IntVectSet> irregSetsGrown[SpaceDim];
      LevelData<BaseIFFAB<Real> > fluxInterpolants[SpaceDim];

      irregSetsSmall.define(m_grids[ilev]);
      for(DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          const EBISBox& ebisBox = m_ebisl[ilev][dit()];
          if(!ebisBox.isAllCovered())
            {
              const Box&  thisBox = m_grids[ilev].get(dit());
              irregSetsSmall[dit()] = ebisBox.getIrregIVS(thisBox);
            }
        }

      for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          EBArith::defineFluxInterpolant(fluxInterpolants[faceDir],
                                         irregSetsGrown  [faceDir],
                                         m_grids[ilev], m_ebisl[ilev],
                                         m_domain[ilev], ncomp, faceDir);
        }


      //set up flux = u*s on irregular sets
      //first put cell-face centered u*s into fluxInterpolant
      for(DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          //first put cell-face centered u*s into fluxInterpolant
          const EBISBox& ebisBox =  m_ebisl[ilev][dit()];
          const Box& cellBox = m_grids[ilev].get(dit());
          //put the interpolant = vel*s  into interpolantGrid
          for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
            {
              IntVectSet ivsIrregGrown = irregSetsGrown[faceDir][dit()];
              ivsIrregGrown &= cellBox;
              FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;

              BaseIFFAB<Real>& interpol = fluxInterpolants[faceDir][dit()];
              interpol.setVal(7.7777e7);
              EBFaceFAB& velDir = (*a_macAdvVel[ilev])[dit()][faceDir];
              EBFaceFAB& scaDir = (*a_macScalar[ilev])[dit()][faceDir];
              for(FaceIterator faceit(ivsIrregGrown, ebisBox.getEBGraph(),
                                      faceDir, stopCrit);
                  faceit.ok(); ++faceit)
                {
                  for(int ivar = 0; ivar < ncomp; ivar++)
                    {
                      Real vel = velDir(faceit(), 0);
                      Real sca = scaDir(faceit(), ivar);
                      interpol(faceit(), ivar) = vel*sca;
                    }
                }
            }
        }

      //exchange ghost cell data for flux interpolant
      for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          fluxInterpolants[faceDir].exchange(Interval(0, ncomp-1));
        }

      //just doing redistribution over a level here since the EB and CF do not intersect
      BaseIVFactory<Real> ivfact(m_ebisl[ilev], irregSetsSmall);
      LevelData<BaseIVFAB<Real> > massDiffLD(m_grids[ilev], ncomp, 2*IntVect::Unit, ivfact);
      Interval consInterv(0, ncomp-1);
      EBLevelRedist levelRedist(m_grids[ilev], m_ebisl[ilev], m_domain[ilev], ncomp);
      levelRedist.setToZero();

      for(DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          const EBISBox& ebisBox = m_ebisl[ilev][dit()];

          //define centroid flux and do interpolation
          BaseIFFAB<Real> centroidFlux[SpaceDim];
          const BaseIFFAB<Real>* interpolantGrid[SpaceDim];
          const IntVectSet& ivsIrregSmall = irregSetsSmall[dit()];
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              const BaseIFFAB<Real>& interpol = fluxInterpolants[idir][dit()];
              interpolantGrid[idir] = &interpol;
              BaseIFFAB<Real>& fluxDir= centroidFlux[idir];
              fluxDir.define(ivsIrregSmall, ebisBox.getEBGraph(), idir, ncomp);
            }

          EBPatchAdvect& patcher  = m_ebLevAd[ilev]->getPatchAdvect(dit());

          EBPatchGodunov& patchGod = patcher;
          patchGod.interpolateFluxToCentroids(centroidFlux,
                                              interpolantGrid,
                                              ivsIrregSmall);

          BaseIVFAB<Real>  ebFlux(ivsIrregSmall, ebisBox.getEBGraph(), ncomp);
          BaseIVFAB<Real> consDiv(ivsIrregSmall, ebisBox.getEBGraph(), ncomp);
          ebFlux.setVal(0.0);
          consDiv.setVal(0.0);

          //this fills  consDiv with kappa*div(US)
          patchGod.consUndividedDivergence(consDiv, centroidFlux, ebFlux, ivsIrregSmall);

          //a_udotdels holds udels(NC).  make it hold (kappa*div(us) + (1-kappa)udels(NC))
          EBCellFAB&        udelsFAB = (*a_uDotDelS[ilev])[dit()];
          BaseIVFAB<Real>&  massDiff = massDiffLD[dit()];
          massDiff.setVal(0.);
          for(VoFIterator vofit(ivsIrregSmall, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              Real kappa = ebisBox.volFrac(vof);
              for(int ivar = 0; ivar < ncomp; ivar++)
                {
                  Real nonConsDiv  = udelsFAB(vof, ivar);
                  Real kappaConsDiv = consDiv(vof, ivar);
                  udelsFAB(vof, ivar) = (1.0-kappa)*nonConsDiv + kappaConsDiv;
                  if(ilev == m_finestLevel)
                    {
                      massDiff(vof, ivar) = -(1.0-kappa)*(kappa*nonConsDiv - kappaConsDiv);
                    }
                  CH_assert(!a_nonConsOnly || !a_consOnly);
                  if(a_nonConsOnly)
                    {
                      udelsFAB(vof, ivar) = nonConsDiv;
                    }
                  if(a_consOnly)
                    {
                      udelsFAB(vof, ivar) = kappaConsDiv;
                    }
                }
            }
          if(ilev == m_finestLevel)
            {
              levelRedist.increment(massDiff, dit(), consInterv);
            }
        } //end loop over boxes
      //smush the mass differnce back in.
      if(!a_consOnly && !a_nonConsOnly && ilev== m_finestLevel)
        {
          levelRedist.redistribute(*a_uDotDelS[ilev], consInterv);
        }
    } //end loop over levels
}
/*****************/
void
EBAMRNoSubcycleOld::
advanceVelocity()
{
  CH_TIMERS("EBAMRNoSubcycleOld::advanceVelocity");
  CH_TIMER("invsicid_advance", t1);
  CH_TIMER("viscous_advance", t2);
  CH_TIMER("cell-centered_projection", t3);
  if(!m_viscousCalc)
    {
      CH_START(t1);
      //for inviscid calc, do not add pressure gradient into
      //the update so we don't have to iterate for the pressure gradient after regrid.
      //If you use the pressure grad in the update and do not iterate after regrid,
      //the solution can ring.
      for (int ilev=0; ilev<= m_finestLevel; ilev++)
        {
          for(DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
            {
              EBCellFAB& oldVel = (*m_veloOld[ilev])[dit()];
              EBCellFAB& newVel = (*m_veloNew[ilev])[dit()];
              EBCellFAB& uDotDelU = (*m_uDotDelU[ ilev])[dit()];
              EBCellFAB& gradPres = (*m_gphiOld [ ilev])[dit()];
              //make udotdelu = udotdelu + gradp
              uDotDelU += gradPres;

              //make udotdelu = -udotdelu - gradp
              uDotDelU *= -1.0;

              //now udotdelu = dt*(-udotdelu - gradp);
              uDotDelU *=  m_dt;

              //put change of velocity into newvel
              newVel.copy(uDotDelU);


              //finally make newvel = old vel - dt*(udotdelu + grad p)
              newVel += oldVel;
            }
        }
      CH_STOP(t1);
    }
  else
    {
      CH_START(t2);
      for(int ilev = 0; ilev <= m_finestLevel; ilev++)
        {
          setCoveredStuffToZero(*m_uDotDelU[ilev]);
          setCoveredStuffToZero(*m_gphiOld[ilev]);
        }
      //add gradient of pressure into udotdelu
      EBAMRDataOps::incr(m_uDotDelU, m_gphiOld, 1.0);
      //make udelu = -udelu-gradp == the source term of heat eqn
      EBAMRDataOps::scale(m_uDotDelU, -1.0);
      if (m_params.m_verbosity > 0)
        {
          pout() << "EBAMRNoSubcycleOld: using viscous solver for velocity components" << endl;
        }
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          //make cellscratch = udotdelu
          EBAMRDataOps::setToZero(m_cellScratch);
          //put component of rhs into cellscratch
          //put component of velo into cellscratch1
          //new velocity comes out in cellscratch2
          Interval vecInterv(idir, idir);
          Interval scaInterv(0, 0);
          for(int ilev = 0; ilev <= m_finestLevel; ilev++)
            {
              m_uDotDelU[ilev]->copyTo(vecInterv,  *m_cellScratch[ilev], scaInterv);
              m_veloOld[ ilev]->copyTo(vecInterv,  *m_cellScratc1[ilev], scaInterv);
              m_veloOld[ ilev]->copyTo(vecInterv,  *m_cellScratc2[ilev], scaInterv);
            }

          //tell EBBC which velocity component we are solving for
          DirichletPoissonEBBC::s_velComp = idir;

          //solve the stinking equation
          int lbase = 0;
          int lmax = m_finestLevel;
          m_tgaSolver[idir]->oneStep(m_cellScratc2,  //vel new
                                     m_cellScratc1, // vel old
                                     m_cellScratch, //source
                                     m_dt, lbase, lmax);

          //now copy the answer back from scrach into velo
          for(int ilev = 0; ilev <= m_finestLevel; ilev++)
            {
              m_cellScratc2[ilev]->copyTo(scaInterv, *m_veloNew[ilev],  vecInterv);
            }
        }
      CH_STOP(t2);
    }

  if (m_params.m_verbosity > 0)
    {
      pout() << "EBAMRNoSubcycleOld: cc projecting" << endl;
    }

  //make u* := u* + dt*gphiOld
  //this makes the output pressure (I-P)u*= gphiNew*dt
  EBAMRDataOps::incr(m_veloNew, m_gphiOld, m_dt);
  //this puts  into gphinew (new pressure gradient)*dt
  //since we have added only a pure gradient, P(u*) = u^n+1 still
  CH_START(t3);
  m_ccProjector->project(m_veloNew, m_gphiNew);
  const Vector<LevelData<EBCellFAB>* >& projPhi = m_ccProjector->getPhi();
  EBAMRDataOps::assign(m_presNew, projPhi);
  CH_STOP(t3);
  filter(m_veloNew);
  if (m_params.m_verbosity > 3)
    {
      m_ccProjector->kappaDivergence(m_cellScratch, m_veloNew);
      Real norm[3];
      for(int inorm = 0; inorm < 3; inorm++)
        {
          norm[inorm] = EBArith::norm(m_cellScratch, m_grids, m_ebisl,
                                      m_params.m_refRatio, 0, inorm, EBNormType::OverBoth);
        }
      pout() << "div(vel) after cell project and filter: " <<
        "L_inf = " << norm[0]  <<
        ", L_1 = " << norm[1] <<
        ", L_2 = " << norm[2] << endl;
    }


  //presnew now holds (new pressure)*dt
  //gphinew now holds (new pressure gradient)*dt
  //so divide out the dt
  EBAMRDataOps::scale(m_presNew, 1.0/m_dt);
  EBAMRDataOps::scale(m_gphiNew, 1.0/m_dt);

  if(m_advanceGphiOnly)
    {
      for (int ilev=0; ilev<= m_finestLevel; ilev++)
        {
          Interval vecInterv(0, SpaceDim-1);
          m_veloOld[ilev]->copyTo(vecInterv, *m_veloNew[ilev], vecInterv);
        }
    }
}
/*****************/
#ifdef CH_USE_HDF5
/*****************/
void
EBAMRNoSubcycleOld::writePlotFile()
{
  CH_TIME("EBAMRNoSubcycleOld::writePlotFile");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::writePlotFile" << endl;
    }
  int curNumLevels = m_finestLevel + 1;

  Vector<string> presnames(1, string("pressure"));
#if CH_SPACEDIM == 2
  Vector<string> vortnames(1, string("vorticity"));
#else
  Vector<string> vortnames(SpaceDim);
#endif

  Vector<string> velonames(SpaceDim);
  Vector<string> gphinames(SpaceDim);
  Vector<string> names;
  for(int idir = 0; idir < SpaceDim; idir++)
    {

      char velochar[100];
      char gphichar[100];
      sprintf(velochar, "velocity%d", idir);
      sprintf(gphichar, "gradPres%d", idir);
      velonames[idir] = string(velochar);
      gphinames[idir] = string(gphichar);

#if CH_SPACEDIM==3
      char vortchar[100];
      sprintf(vortchar, "vorticity%d", idir);
      vortnames[idir] = string(vortchar);
#endif

    }

  names = velonames;
  names.append(gphinames);
  names.append(presnames);
  names.append(vortnames);

  string numZeros;
  if(m_curStep < 10)
    {
      numZeros = string("0000");
    }
  else if(m_curStep < 100)
    {
      numZeros = string("000");
    }
  else if(m_curStep < 1000)
    {
      numZeros = string("00");
    }
  else if(m_curStep < 10000)
    {
      numZeros = string("0");
    }
  else
    {
      numZeros = string("");
    }

  char fileChar[100];
  int ncells = m_domain[0].size(0);
  sprintf(fileChar, "plot.nx%d.step%s%d.%dd.hdf5", ncells, numZeros.c_str(), m_curStep, SpaceDim);

  bool replaceCovered = false;
  Vector<Real> coveredValues;

  int nlev = m_finestLevel + 1;
  Vector<LevelData<EBCellFAB>* > vorticity(nlev, NULL);
  Vector<LevelData<EBCellFAB>* > outputData(nlev, NULL);

  for(int ilev = 0; ilev < nlev; ilev++)
    {
      vorticity[ilev] = new LevelData<EBCellFAB>();
      computeVorticity(*vorticity[ilev], ilev);

#if CH_SPACEDIM==2
      int nvar = 2*SpaceDim + 2;
#else
      int nvar = 3*SpaceDim + 1;
#endif
      EBCellFactory ebcellfact(m_ebisl[ilev]);
      outputData[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], nvar, IntVect::Zero, ebcellfact);

      Interval srcInterv, dstInterv;
      srcInterv = Interval(0, SpaceDim-1);
      dstInterv = Interval(0, SpaceDim-1);
      m_veloNew[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);

      srcInterv = Interval(0, SpaceDim-1);
      dstInterv = Interval(SpaceDim, 2*SpaceDim-1);
      m_gphiNew[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);

      srcInterv = Interval(0, 0);
      dstInterv = Interval(2*SpaceDim, 2*SpaceDim);
      m_presNew[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);

#if CH_SPACEDIM==2
      srcInterv = Interval(0, 0);
      dstInterv = Interval(2*SpaceDim+1, 2*SpaceDim+1);
      vorticity[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);
#else
      srcInterv = Interval(0, SpaceDim-1);
      dstInterv = Interval(2*SpaceDim+1, 3*SpaceDim);
      vorticity[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);
#endif
      setCoveredStuffToZero(*outputData[ilev]);
    }

  string filename(fileChar);
  writeEBHDF5(filename,
              m_grids,
              outputData,
              names,
              m_domain[0].domainBox(),
              m_dx[0],
              m_dt, m_time,
              m_params.m_refRatio,
              curNumLevels,
              replaceCovered,
              coveredValues);

  for(int ilev = 0; ilev <nlev; ilev++)
    {
      delete vorticity[ilev];
      delete outputData[ilev];
    }
}
/**********/
void
EBAMRNoSubcycleOld::writeCheckpointFile()
{
  CH_TIME("EBAMRNoSubcycleOld::writeCheckpointFile");
  CH_assert(m_isSetup);
  // Setup the level header information
  HDF5HeaderData header;

  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::writeCheckpointFile" << endl;
    }

  //all the stuff in m_params had to come in at define time
  //this is also true of viscosity, the domains, and the boundary conditions
  //viscous calc gets set there too

  //bool conversion to int
  int iuseFixedDt       = 0;
  if(m_useFixedDt) iuseFixedDt = 1;

  header.m_real["time"]                   = m_time;
  header.m_real["dt"]                     = m_dt;
  header.m_real["prescribed_dt"]          = m_prescribedDt;
  header.m_int ["cur_step"]               = m_curStep;
  header.m_int ["finest_level"]           = m_finestLevel;
  header.m_int ["use_fixed_dt"]           = iuseFixedDt;

  char iter_str[100];

  sprintf(iter_str, "check%d.%dd.hdf5",m_curStep, SpaceDim );

  HDF5Handle handleOut(iter_str, HDF5Handle::CREATE);
  // Write the header for this level
  header.writeToFile(handleOut);
  for(int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      handleOut.setGroupToLevel(ilev);
      write(handleOut,m_grids[ilev]);
      write(handleOut,*m_veloOld[ilev],"veloOld");
      write(handleOut,*m_veloNew[ilev],"veloNew");
      write(handleOut,*m_gphiOld[ilev],"gphiOld");
      write(handleOut,*m_gphiNew[ilev],"gphiNew");
      write(handleOut,*m_presOld[ilev],"presOld");
      write(handleOut,*m_presNew[ilev],"presNew");
    }
  handleOut.close();
}
/*****************/
void
EBAMRNoSubcycleOld::readCheckpointFile(const string& a_restartFile)
{
  CH_TIME("EBAMRNoSubcycleOld::readCheckpointFile");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::readCheckpointFile" << endl;
    }

  HDF5Handle handleIn(a_restartFile, HDF5Handle::OPEN_RDONLY);
  HDF5HeaderData header;
  header.readFromFile(handleIn);

  //all the stuff in m_params had to come in at define time
  //this is also true of viscosity, the domains, and the boundary conditions

  m_time          =   header.m_real["time"]                   ;
  m_dt            =   header.m_real["dt"]                     ;
  m_prescribedDt  =   header.m_real["prescribed_dt"]          ;
  m_curStep       =   header.m_int ["cur_step"]               ;
  m_finestLevel   =   header.m_int ["finest_level"]           ;
  int iuseFixedDt =   header.m_int["use_fixed_dt"];
  m_useFixedDt =  (iuseFixedDt == 1);

  //get all the grids
  for(int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      handleIn.setGroupToLevel(ilev);
      // Get the grids
      Vector<Box> vboxGrids;
      const int gridStatus = read(handleIn, vboxGrids);
      if (gridStatus != 0)
        {
          MayDay::Error("readCheckpointLevel: file has no grids");
        }

      Vector<int> proc_map;
      EBEllipticLoadBalance(proc_map,vboxGrids, m_domain[ilev]);

      m_grids[ilev]= DisjointBoxLayout(vboxGrids,proc_map);
    }

  //define stuff using grids
  defineEBISLs();
  defineOld();
  defineNew();
  //now input the actual data
  for(int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      handleIn.setGroupToLevel(ilev);
      read<EBCellFAB>(handleIn, *m_veloOld[ilev], "veloOld", m_grids[ilev], Interval(), false);
      read<EBCellFAB>(handleIn, *m_veloNew[ilev], "veloNew", m_grids[ilev], Interval(), false);
      read<EBCellFAB>(handleIn, *m_gphiOld[ilev], "gphiOld", m_grids[ilev], Interval(), false);
      read<EBCellFAB>(handleIn, *m_gphiNew[ilev], "gphiNew", m_grids[ilev], Interval(), false);
      read<EBCellFAB>(handleIn, *m_presOld[ilev], "presOld", m_grids[ilev], Interval(), false);
      read<EBCellFAB>(handleIn, *m_presNew[ilev], "presNew", m_grids[ilev], Interval(), false);
    }
  handleIn.close();
}
/*****************/
#endif // CH_USE_HDF5
/*****************/
void
EBAMRNoSubcycleOld::
setupForRestart(const string& a_restartFile)
{
  CH_TIME("EBAMRNoSubcycleOld::setupForRestart");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::setupForRestart" << endl;
    }
  m_isSetup = true;
  m_doRestart  = true;
#ifdef CH_USE_HDF5
  readCheckpointFile(a_restartFile);
#else
  MayDay::Error("cannot restart from checkpoint without hdf5");
#endif
  defineProjections();
  defineIrregularData();
  postInitialize();
}

/**********/
void
EBAMRNoSubcycleOld::
setupForFixedHierarchyRun(const Vector<Vector<Box> >& a_grids)
{
  CH_TIME("EBAMRNoSubcycleOld::setupForFixedHierarchyRun");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycleOld::setupForFixedHierarchyRun" << endl;
    }
  //turn off regridding
  m_params.m_regridInterval = -1;
  m_isSetup = true;

  m_finestLevel = a_grids.size() - 1;
  for(int ilev = 0; ilev < m_finestLevel; ilev++)
    {
      CH_assert(a_grids[ilev].size() > 0);
    }
  initialGrid(a_grids);

  initialData();
  defineIrregularData();
  // finally, call post-initialization
  postInitialize();
  m_doRestart  = false;
}
/**********/
