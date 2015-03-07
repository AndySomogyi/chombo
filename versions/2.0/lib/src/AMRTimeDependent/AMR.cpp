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
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>

#include "REAL.H"
#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "BRMeshRefine.H"
#include "DisjointBoxLayout.H"
#include "CH_HDF5.H"
#include "AMRLevel.H"
#include "parstream.H"
#include "Tuple.H"
#include "BoxIterator.H"
#include "AMR.H"
#include "NamespaceHeader.H"

#ifdef CH_USE_TIMER
using namespace Chombo;
#endif

#define SMALL_TIME 1.0e-8

void AMR::setDefaultValues()
{
  m_isDefined = false;
  m_isSetUp   = false;
  m_max_level = -1;
  m_finest_level = -1;
  m_checkpoint_interval= -1;
  m_plot_interval=-1;
  m_max_grid_size= 0;
  m_max_base_grid_size = m_max_grid_size;
  m_restart_step = -1;
  m_lastcheck_step = -1;
  m_cur_step =0;
  m_maxDtGrow = 1.1;
  m_time_eps = 1.0e-6;
  m_dt_base = -1;
  m_amrlevels.resize(0);
  m_use_meshrefine=false;
  m_plotfile_prefix = string("pltstate");
  m_checkpointfile_prefix = string("chk");
  m_verbosity = 0;
  m_cur_time = 0;
  m_dt_tolerance_factor = 1.1;
  m_fixedDt = -1;
  m_blockFactor = 4;
#ifdef CH_USE_TIMER
  m_timer = NULL ;  
#endif
}

AMR::AMR()
{
  setDefaultValues();
}

void AMR::clearMemory()
{
  // clean up memory
  if (m_isDefined)
    {
      for (int lev = 0; lev < m_amrlevels.size(); ++lev)
        {
          if (m_amrlevels[lev] != NULL)
            {
              delete m_amrlevels[lev];
              m_amrlevels[lev] = NULL;
            }
        }
      m_isDefined = false;
      m_isSetUp = false;
    }
}

AMR::~AMR()
{
  clearMemory();
}

void AMR::maxGridSize(int a_max_grid_size)
{
  CH_assert(a_max_grid_size >= 0);

  if (m_max_base_grid_size == 0)
    {
      m_max_base_grid_size = a_max_grid_size;
    }

  m_max_grid_size = a_max_grid_size;
  m_mesh_refine.maxSize(a_max_grid_size);
}

void AMR::maxBaseGridSize(int a_max_base_grid_size)
{
  CH_assert(a_max_base_grid_size > 0);

  m_max_base_grid_size = a_max_base_grid_size;
}

void AMR::dtToleranceFactor(Real a_dt_tolerance_factor)
{
  CH_assert(a_dt_tolerance_factor > 0);

  m_dt_tolerance_factor = a_dt_tolerance_factor;
}

void AMR::gridBufferSize(int a_grid_buffer_size)
{
  m_mesh_refine.bufferSize(a_grid_buffer_size);
}

int AMR::maxGridSize() const
{
  CH_assert(isDefined());

  return(m_max_grid_size);
}

int AMR::maxBaseGridSize() const
{
  CH_assert(isDefined());
  return(m_max_base_grid_size);
}

void AMR::plotPrefix(const std::string& a_plotfile_prefix)
{
  CH_assert(isDefined());

  m_plotfile_prefix = a_plotfile_prefix;
}

void AMR::blockFactor(int a_blockFactor)
{
  CH_assert(a_blockFactor >= 1);

  m_blockFactor = a_blockFactor;
  m_mesh_refine.blockFactor(a_blockFactor);
}

void AMR::fillRatio(Real a_fillRatio)
{
  CH_assert(isDefined());

  m_fillRatio = a_fillRatio;
  m_mesh_refine.fillRatio(a_fillRatio);
}

void AMR::checkpointPrefix(const std::string& a_checkpointfile_prefix)
{
  CH_assert(isDefined());

  m_checkpointfile_prefix = a_checkpointfile_prefix;
}

void AMR::regridIntervals(const Vector<int>& a_regridIntervals)
{
  CH_assert(isDefined());
  CH_assert(a_regridIntervals.size() >= m_max_level);

  m_regrid_intervals = a_regridIntervals;
}

void AMR::plotInterval(int a_plot_interval)
{
  m_plot_interval = a_plot_interval;
}

void AMR::checkpointInterval(int a_checkpoint_interval)
{
  m_checkpoint_interval = a_checkpoint_interval;
}

void AMR::define(int                          a_max_level,
                 const Vector<int>&           a_ref_ratios,
                 const Box&                   a_prob_domain,
                 const AMRLevelFactory* const a_amrLevelFact)
{
  ProblemDomain physdomain(a_prob_domain);

  define(a_max_level, a_ref_ratios, physdomain, a_amrLevelFact);
}

void AMR::define(int                          a_max_level,
                 const Vector<int>&           a_ref_ratios,
                 const ProblemDomain&         a_prob_domain,
                 const AMRLevelFactory* const a_amrLevelFact)
{
  clearMemory();

  // Function already called in AMR() ... so make sure you set any values
  // _after_ the define() function. (ndk)
  setDefaultValues();

  m_isDefined = true;
  CH_assert(a_ref_ratios.size() > a_max_level);

  if (m_verbosity >= 3)
    {
      pout() << "AMR::define(...)" << endl;
    }

  CH_assert(a_max_level >= 0);

  m_max_level = a_max_level;

  // set to zero to start
  m_finest_level=0;
  m_use_meshrefine = true;

  // import refinement ratios
  m_ref_ratios = a_ref_ratios;
  if (m_ref_ratios.size() == 0)
    {
      m_ref_ratios.resize(1,1);
    }

  // set the subcycling reduction factors to 1, i.e., no subcycling yet
  m_reduction_factor.resize(a_ref_ratios.size(),1);
  if (m_reduction_factor.size() == 0)
    {
      m_reduction_factor.resize(1,1);
    }

  // resize vectors
  m_dt_new.resize(m_max_level+1);
  m_dt_cur.resize(m_max_level+1);
  m_steps_since_regrid.resize(m_max_level+1,0);
  m_regrid_intervals.resize(m_max_level+1,0);
  m_cell_updates.resize(m_max_level+1,0);

  // create hierarchy of levels
  m_amrlevels.resize(m_max_level+1,NULL);

  // create base level
  m_amrlevels[0] = a_amrLevelFact->new_amrlevel();

  m_amrlevels[0]->define(NULL,
                         a_prob_domain,
                         0,
                         m_ref_ratios[0]);

  // create finer levels
  int ref_factor = 1;
  for (int level = 0; level < m_max_level; ++level)
    {
      ref_factor *= m_ref_ratios[level];

      const ProblemDomain level_prob_domain = refine(a_prob_domain, ref_factor);

      m_amrlevels[level+1] = a_amrLevelFact->new_amrlevel();
      m_amrlevels[level+1]->define(m_amrlevels[level],
                                   level_prob_domain,
                                   level+1,
                                   m_ref_ratios[level+1]);
      m_amrlevels[level]->finerLevelPtr(m_amrlevels[level+1]);
    }

  // build BRMeshRefine Object using some default values
  // these can be reset later, if desired.
  Real defaultFillRatio = 0.75;
  int defaultBlockFactor = 1;
  int defaultBufferSize = 1;

  m_mesh_refine.define(a_prob_domain, a_ref_ratios, defaultFillRatio,
                       defaultBlockFactor, defaultBufferSize,m_max_grid_size);
}

void AMR::setupForFixedHierarchyRun(const Vector<Vector<Box> >& a_amr_grids,
                                    int                         a_proper_nest)
{
  CH_assert(isDefined());

  m_isSetUp = true;

  if (m_verbosity >= 3)
    {
      pout() << "AMR::setFixedHierarchy(...)" << endl;
    }

  m_finest_level = a_amr_grids.size() -1;
  m_finest_level_old = m_finest_level;

  // set to zero to start
  m_finest_level=0;

  m_use_meshrefine = false;

  // set this to -1 initially
  m_dt_base = -1;

#ifndef NDEBUG
  // check to see that grids satisfy proper nesting
  // and that they're coarsenable
  int numLevels = a_amr_grids.size();
  if (numLevels > 1)
    {
      Box thisDomain;
      for (int lev = 1; lev < numLevels; lev++)
        {
          const Vector<Box>& crseGrids = a_amr_grids[lev-1];
          const Vector<Box>& fineGrids = a_amr_grids[lev];

          if (fineGrids.size() > 0)
            {
              // first build IntVectSet of coarse level
              IntVectSet PNdomainCrse;
              for (int icrse = 0; icrse < crseGrids.size(); icrse++)
                {
                  PNdomainCrse |= crseGrids[icrse];
                }

              // if base level, get base level domain
              if (lev==1)
                {
                  thisDomain = PNdomainCrse.minBox();
                }

              // get proper nesting domain
              PNdomainCrse.nestingRegion(a_proper_nest,thisDomain);

              // now loop through fine-level boxes, coarsen, and
              // make sure that they're contained in PN domain
              for (int ifine = 0; ifine < fineGrids.size(); ++ifine)
                {
                  Box crseBox(fineGrids[ifine]);
                  crseBox.coarsen(m_ref_ratios[lev-1]);
                  CH_assert(PNdomainCrse.contains(crseBox));

                  // also check to see if fine box is coarsenable
                  crseBox.refine(m_ref_ratios[lev-1]);
                  CH_assert(crseBox == fineGrids[ifine]);
                }

              // now refine physical domain to next level
              if (lev<numLevels-1)
                {
                  thisDomain.refine(m_ref_ratios[lev-1]);
                }
            } // end if fineGrids size isn't 0
        } // end loop over levels
    } // end if number of levels > 1; end proper nesting check
#endif

  m_amr_grids = a_amr_grids;
  m_regrid_intervals.resize(m_max_level, -1);
  m_cur_step = 0;
  m_finest_level = a_amr_grids.size() - 1;

  if (m_finest_level > m_max_level)
    {
      cerr << "AMR::setupForFixedHierarchy: too many levels in input grids" << endl;
      cerr << "max_level from define function = " << m_max_level << endl;;
      cerr << "finest level(top level of grids in setup) = = " << m_finest_level << endl;
      MayDay::Error("AMR::setupForFixedHierarchy: define call and setup call inconsistent in class AMR");
    }

  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_amrlevels[level]->initialGrid(m_amr_grids[level]);
      m_amrlevels[level]->initialData();
    }

  for (int level = m_finest_level + 1; level <= m_max_level; ++level)
    {
      m_amrlevels[level]->regrid(Vector<Box>());
    }

  // call post-initialize once all the levels have been defined
  for (int level = m_finest_level; level >= 0; --level)
    {
      m_amrlevels[level]->postInitialize();
    }

  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_dt_new[level] = m_amrlevels[level]->computeInitialDt();
      m_dt_cur[level] = m_dt_new[level];
    }

  assignDt();
}

void AMR::setupForNewAMRRun()
{
  CH_assert(m_isDefined);

  m_isSetUp = true;

  if (m_verbosity >= 3)
    {
      pout() << "AMR::setupForNewAMRRun" << endl;
    }

  m_cur_step = 0;

  initialGrid();

  m_finest_level_old = m_finest_level;
  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_dt_new[level] = m_amrlevels[level]->computeInitialDt();
      m_dt_cur[level] = m_dt_new[level];
    }

  assignDt();
}

#ifdef CH_USE_HDF5
// read checkpoint file
void AMR::setupForRestart(HDF5Handle& a_handle)
{
  CH_assert(m_isDefined);

  m_isSetUp = true;

  if (m_verbosity >= 3)
    {
      pout() << "AMR::restart" << endl;
    }

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (m_verbosity >= 3)
    {
      pout() << "hdf5 header data: " << endl;
      pout() << header << endl;
    }

  // read max level
  if (header.m_int.find("max_level") == header.m_int.end())
    {
      MayDay::Error("AMR::restart: checkpoint file does not contain max_level");
    }

  // note that this check should result in a warning rather than an error,
  // because you should be able to restart with a different number of levels
  // (DFM 2/27/02)
  int max_level_check = header.m_int ["max_level"];
  if (max_level_check != m_max_level)
    {
      pout() << "AMR::restart: checkpoint file inconsistent with inputs to define " << endl;
      pout() << "max level input to define = " << m_max_level << endl;
      pout() << "max level in checkpoint = " << max_level_check << endl;
      MayDay::Warning("AMR::restart: checkpoint file inconsistent with inputs to define ");
    }

  if (m_verbosity >= 2)
    {
      pout() << "read max_level = " << m_max_level << endl;
    }

  // read finest level
  if (header.m_int.find("num_levels") == header.m_int.end())
    {
      MayDay::Error("AMR::restart: checkpoint file does not contain num_levels");
    }
  int num_levels = header.m_int ["num_levels"];

  if (m_verbosity >= 2)
    {
      pout() << "read num_levels = " << num_levels << endl;
    }

  m_finest_level = num_levels - 1;
  m_finest_level_old = m_finest_level;

  if (m_finest_level > m_max_level)
    {
      pout() << "AMR::restart: checkpoint file inconsistent with inputs to define " << endl;
      pout() << "numlevels input to define = " << m_max_level + 1<< endl;
      pout() << "numlevels in checkpoint = " << num_levels << endl;
      MayDay::Error("AMR::restart: checkpoint file inconsistent with inputs to define ");
    }

  if (m_verbosity >= 2)
    {
      pout() << "set finest_level = " << m_finest_level << endl;
    }

  if (m_finest_level > m_max_level)
    {
      MayDay::Error("AMR::restart: finest_level > max_level");
    }

  if (header.m_int.find("iteration") == header.m_int.end())
    {
      MayDay::Error("AMR::restart: checkpoint file does not contain iteration");
    }

  m_cur_step = header.m_int ["iteration"];

  if (m_verbosity >= 2)
    {
      pout() << "read cur_step = " << m_cur_step << endl;
    }

  m_restart_step = m_cur_step;

  if (m_verbosity >= 2)
    {
      pout() << "set restart_step = " << m_restart_step << endl;
    }

  if (header.m_real.find("time") == header.m_real.end())
    {
      MayDay::Error("AMR::restart: checkpoint file does not contain time");
    }

  m_cur_time = header.m_real["time"];

  if (m_verbosity >= 2)
    {
      pout() << "read cur_time = " << m_cur_time << endl;
    }

  // read regrid intervals
  if (header.m_int.find("regrid_interval_0") == header.m_int.end())
    {
      MayDay::Error("AMR::restart: Error: regrid intervals are not in the checkpoint file");
    }
  else
    {
      //hey if max_level > 10^98, we have a problem
      char level_str[100];

      //don't go all the way to max level because you don't have to
      for (int level = 0; level < m_max_level; ++level)
        {
          sprintf(level_str, "%d", level);

          const std::string label = std::string("regrid_interval_") + level_str;
          if (header.m_int.find(label) == header.m_int.end())
            {
              cerr << "checkpoint file does not have " << label << endl;
              MayDay::Error("not enough regrid intervals in amr::restart");
            }
          else
            {
              m_regrid_intervals[level] = header.m_int [label];
            }
        }
    }

  // read physics class header data
  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_amrlevels[level]->readCheckpointHeader(a_handle);
      m_amrlevels[level]->readCheckpointLevel(a_handle);
    }

  for (int level = 0; level < m_finest_level; ++level)
    {
      int refratio_test = m_amrlevels[level]->refRatio();
      if (refratio_test != m_ref_ratios[level])
        {
          pout() << "AMR::restart: checkpoint file inconsistent with inputs to define " << endl;
          pout() << "for level " << level << endl;
          pout() << "refratio input to define = " << m_ref_ratios[level] << endl;
          pout() << "refratio in checkpoint = "  << refratio_test << endl;
          MayDay::Error("AMR::restart: checkpoint file inconsistent with inputs to define ");
        }
    }

  // maintain time steps
  m_dt_new.resize(m_max_level+1);
  m_dt_cur.resize(m_max_level+1);

  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_dt_new[level] = m_amrlevels[level]->dt();
      m_dt_cur[level] = m_dt_new[level];
    }

  assignDt();

  // maintain steps since regrid
  m_steps_since_regrid.resize(m_max_level+1,0);

  //restart cell updates(we could also output them to the chk file)
  m_cell_updates.resize(m_max_level+1,0);

  // final thing to do -- call initialGrid and initialData on undefined levels
  // (just in case there are setup things which need to be done there
  // (DFM 2/27/02)
  for (int level = m_finest_level+1; level <= m_max_level; ++level)
    {
      m_amrlevels[level]->initialGrid(Vector<Box>());
      m_amrlevels[level]->initialData();
    }
}
#endif

void AMR::conclude() const
{
  CH_assert(isDefined());
  CH_assert(isSetUp());

  if (m_verbosity >= 3)
    {
      pout() << "AMR::conclude" << endl;
    }

  if (m_plot_interval >= 0)
  {
    writePlotFile();
  }

  if ((m_checkpoint_interval >= 0)     &&
      (m_lastcheck_step != m_cur_step) &&
      (m_restart_step != m_cur_step))
    {
      writeCheckpointFile();
    }

  if (m_verbosity >= 2)
    {
      long long total_cell_updates = 0;
      for (int ll = 0; ll < m_max_level+1; ll++)
        {
          total_cell_updates += m_cell_updates[ll];

          pout() << "number of points updated at level " << ll <<" = " <<
              m_cell_updates[ll] << endl;
        }
#ifdef CH_OSF1
      pout() << "total number of points updated = " << (Real)total_cell_updates << endl;
#else
      pout() << "total number of points updated = " << total_cell_updates << endl;
#endif
    }
}

// go baby go
void AMR::run(Real a_max_time, int a_max_step)
{
  CH_assert(isDefined());
  CH_assert(isSetUp());

#ifdef CH_USE_TIMER
  double last_timestep_time = 0 ;
#endif
  if (m_verbosity >= 3)
    {
      pout() << "AMR::coarseTimeStep:" << endl;
      pout() << "max_time = " << a_max_time << endl;
      pout() << "max_step = " << a_max_step << endl;
    }

  Real old_dt_base = m_dt_base;

  for ( ; (m_cur_step < a_max_step) &&
          (a_max_time - m_cur_time > m_time_eps*m_dt_base);
        ++m_cur_step, m_cur_time += old_dt_base)
    {
#ifdef CH_USE_TIMER
      m_timer->start() ;
#endif
      old_dt_base = m_dt_base;
      for (int level = 0; level <= m_max_level; ++level)
        {
          m_amrlevels[level]->time(m_cur_time);
        }

      if ((m_checkpoint_interval > 0)      &&
          (m_lastcheck_step != m_cur_step) &&
          (m_restart_step != m_cur_step)   &&
          (m_cur_step % m_checkpoint_interval == 0))
        {
          writeCheckpointFile();
          m_lastcheck_step= m_cur_step;
        }

      if ((m_plot_interval > 0) &&
          (m_cur_step % m_plot_interval == 0))
        {
          writePlotFile();
        }

      int level = 0;
      int stepsLeft = 0;
      bool timeBoundary = true;
      (void)timeStep(level,stepsLeft,timeBoundary);

      assignDt();
#ifdef CH_USE_TIMER
      m_timer->stop();
#endif
      if(m_verbosity >= 1)
        {
          pout() << resetiosflags(ios::fixed)
                 << "coarse time step " << setw(3) << m_cur_step
                 << "  old time = " << setw(12)
                 << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << m_cur_time
                 << "  old dt = "   << setw(12)
                 << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << old_dt_base
#ifdef CH_USE_TIMER
                 << "  wallclocktime = "
                 << m_timer->wc_time() - last_timestep_time
#endif
                 << resetiosflags(ios::scientific)
                 << endl;
        }
#ifdef CH_USE_TIMER
      last_timestep_time = m_timer->wc_time() ;
#endif
    }
}

// add restriction on growth of new time step
// dtnew = min(dtold*factor, dtwant) default factor = 1.1
void AMR::assignDt()
{
  CH_assert(isDefined());

  if (m_verbosity >= 3)
    {
      pout() << "AMR::assignDt" << endl;
    }

  if (m_fixedDt > 0)
    {
      m_dt_base = m_fixedDt;
    }
  else
    {
      m_dt_base = m_dt_new[0];

      // Multiply time step for each level by refinement factor, then min
      // across all levels.  Only go to finest_level_old because the higher
      // level dt's have not been set
      for (int level = 1; level <= m_finest_level_old; ++level)
        {
          int ref_factor = 1;

          for (int ilev = 0; ilev < level; ++ilev)
            {
              ref_factor *= m_ref_ratios[ilev];
            }

          Real dt_base_equiv = m_dt_new[level] * ref_factor;
          m_dt_base = Min(m_dt_base, dt_base_equiv);
        }

      // Check all the actual dt's (scaled by the refinement factor) and only
      // allow a growth of "m_maxDtGrow".
      for (int level = 0; level <= m_finest_level_old; ++level)
        {
          int ref_factor = 1;

          for (int ilev = 0; ilev < level; ++ilev)
            {
              ref_factor *= m_ref_ratios[ilev];
            }

          Real dt_base_equiv = m_dt_cur[level] * ref_factor * m_maxDtGrow;
          m_dt_base = Min(m_dt_base, dt_base_equiv);
        }
    }

  // refine base time step for all levels
  for (int level = 0; level <= m_max_level; ++level)
    {
      int ref_factor = 1;

      // reset reduction factors as there is no subcycling going on yet
      m_reduction_factor[level] = 1;

      for (int ilev = 0; ilev < level; ++ilev)
        {
          ref_factor *= m_ref_ratios[ilev];
        }

      Real dt_level = m_dt_base / ref_factor;
      m_amrlevels[level]->dt(dt_level);
    }
}

// This function performs a time step of level "a_level" and all finer
// levels.  It does subcycling in time as necessary (see below).  For this
// reason the number of time steps left, "a_stepsLeft", for this level is
// passed in.  If no additional subcycling occurs then this is returned.
// If additional subcycling occurs then this is used to compute the new
// number of time steps remaining (for this level) and this is returned.
int AMR::timeStep(int a_level, int a_stepsLeft, bool a_coarseTimeBoundary)
{
  CH_assert(isDefined());
  CH_assert(isSetUp());

  if (m_verbosity >= 3)
    {
      pout() << "AMR::timeStep(" << a_level << ")" << endl;
    }

  if (a_level < m_max_level)
    {
      if (m_verbosity >= 4)
        {
          pout() << "Regrid (level " << a_level << ") needed - ";
        }

      // regrid if necessary
      if (needToRegrid(a_level,a_stepsLeft))
        {
          if (m_verbosity >= 4)
            {
              pout() << "yes" << endl;
            }

          regrid(a_level);
        }
      else
        {
          if (m_verbosity >= 4)
            {
              pout() << "no" << endl;
            }
        }
    }

  // If this wasn't just done by the next coarser level, check to see if
  // it is necessary to do additional subcycling in time.
  if (!a_coarseTimeBoundary)
    {
      // The factor by which the current time step at the current level
      // has been divided (so far) for subcycling.
      int maxFactor = m_reduction_factor[a_level];

      // Compute the new subcycling factor for this level and all finer
      // levels and find the maximum
      for (int i = a_level; i <= m_max_level; i++)
        {
          int factor;
          Real dtCur = m_amrlevels[i]->dt();
          Real dtNew = m_dt_new[i];

          // The current factor for level "i"
          factor = m_reduction_factor[i];

          // While the current dt exceeds the new (max) dt by a tolerance
          // double the subcycling factor and half the current dt
          while (dtCur > m_dt_tolerance_factor*dtNew)
            {
              factor *= 2;
              dtCur *= 0.5;
            }

          if (factor > maxFactor)
            {
              maxFactor = factor;
            }
        }

      // More subcycling is necessary
      if (maxFactor > m_reduction_factor[a_level])
        {
          if (m_verbosity >= 3)
            {
              pout() << "  Subcycling --- maxFactor: " << maxFactor << endl;
            }

          // Adjust the number of time steps left for the current level
          a_stepsLeft = (a_stepsLeft+1)*maxFactor/m_reduction_factor[a_level] - 1;

          // Adjust the dt's on this and all finer levels
          for (int i = a_level; i <= m_max_level; i++)
            {
              int factor;

              factor = maxFactor/m_reduction_factor[i];
              m_amrlevels[i]->dt(m_amrlevels[i]->dt()/factor);

              if (m_verbosity >= 4)
                {
                  pout() << "    Level " << i << ": factor: " << factor
                         << " (" << m_reduction_factor[i] << "), dt: "
                         << m_amrlevels[i]->dt() << endl;
                }

              m_reduction_factor[i] = maxFactor;
            }
        }
    }

  // advance this level
  m_amrlevels[a_level]->advance();

  // get the new dt
  Real dt_level = m_amrlevels[a_level]->computeDt();

  // Save the current dt and the new (max) dt.
  m_dt_cur[a_level] = m_amrlevels[a_level]->dt();
  m_dt_new[a_level] = dt_level;

  // increment counter that gives the number of cells updates.
  long numPts = 0;
  const Vector<Box >& levelBoxes = m_amrlevels[a_level]->boxes();

  for (int ll = 0;ll < levelBoxes.size();ll++)
    {
      numPts += levelBoxes[ll].numPts();
    }

  m_cell_updates[a_level] += numPts;

  if(a_level < m_max_level)
    {
      ++m_steps_since_regrid[a_level];
    }

  // advance the finer levels by subcycling
  if(a_level < m_finest_level)
    {
      int stepsLeft = m_ref_ratios[a_level];
      bool timeBoundary = true;

      while (stepsLeft > 0)
        {
          stepsLeft--;

          // Advance the finer level and take into account possible
          // subcycling by allowing for a change in "stepsLeft".
          //[NOTE: the if() test looks redundant with above, but it is not 
          //       because m_finest_level may change during a regrid();
          //       why you would regrid during a subcycle I don't know. <dbs>]
          if (a_level < m_finest_level)
            stepsLeft = timeStep(a_level+1,stepsLeft,timeBoundary);

          // The first time the next finer level time aligns with the current
          // level time.  After that this is not the case.
          //[NOTE: this if() test _is_ redundant. <dbs>]
          if (timeBoundary == true)
            {
              timeBoundary = false;
            }
        }
    }

  m_amrlevels[a_level]->postTimeStep();

  // Return the (possibly updated) number of time steps left on this level.
  return(a_stepsLeft);
}

bool AMR::needToRegrid(int a_level, int a_stepsLeft) const
{
  CH_assert(isDefined());
  CH_assert(isSetUp());

  bool regridThisLevel = false;
  bool regridNextCoarserLevel = false;

  if (m_verbosity >= 3)
    {
      pout() << "AMR::needToRegrid(" << a_level << ")" << endl;
    }

  regridThisLevel =
    (m_regrid_intervals[a_level] > 0) &&
    (m_steps_since_regrid[a_level] >= m_regrid_intervals[a_level]);

  int nextCoarserLevel = a_level-1;

  regridNextCoarserLevel =
    (a_stepsLeft == 0)                         &&
    (nextCoarserLevel >= 0)                    &&
    (m_regrid_intervals[nextCoarserLevel] > 0) &&
    (m_steps_since_regrid[nextCoarserLevel] >= m_regrid_intervals[nextCoarserLevel]);

  return(regridThisLevel && !regridNextCoarserLevel);
}

// generate new grid hierarchy
void AMR::regrid(int a_base_level)
{
  CH_assert(isDefined());
  CH_assert(isSetUp());

  if (m_verbosity >= 2)
    {
      pout() << "AMR::regrid(" << a_base_level << ")" << endl;
    }

  if (m_verbosity >= 2)
    {
      pout() << "AMR::regrid: base level = " << a_base_level << endl;
    }

  for (int level = a_base_level; level <= m_max_level; ++level)
    {
      m_steps_since_regrid[level] = 0;
    }

  m_finest_level_old = m_finest_level;
  int top_level = Min(m_finest_level, m_max_level - 1);

  Vector<Vector<Box> > old_grids(top_level+1);
  Vector<Vector<Box> > new_grids;
  Vector<IntVectSet> tags(top_level+1);
  Vector<ProblemDomain> problem_domains(top_level+1);

  if (m_use_meshrefine)
    {
      for (int level = a_base_level; level <= top_level; ++level)
        {
          m_amrlevels[level]->tagCells(tags[level]);
          old_grids[level] = m_amrlevels[level]->boxes();
          problem_domains[level] = m_amrlevels[level]->problemDomain();

          if (m_verbosity >= 2)
            {
              pout() << "AMR::regrid: problem domain[" << level << "]: " << problem_domains[level] << endl;
            }

          if (m_verbosity >= 4)
            {
              pout() << "AMR::regrid: old_grids[" << level << "]: " << endl;

              for (int i = 0; i< old_grids[level].size(); ++i)
                {
                  pout() << "  " << i << ": " << old_grids[level][i] << endl;
                }

              if (m_verbosity >= 5)
                {
                  pout() << "AMR::regrid: tags[" << level << "]: " << tags[level]
                         << endl;
                }
            }
        }

      int new_finest_level;
//    if (procID() == uniqueProc(SerialTask::compute))
//      {
      new_finest_level = m_mesh_refine.regrid(new_grids,
                                              tags,
                                              a_base_level,
                                              top_level,
                                              old_grids);
//      }

//    broadcast(new_finest_level, uniqueProc(SerialTask::compute));

      //can only add one level at a time
      new_finest_level = Min(m_finest_level+1, new_finest_level);

      if ((m_finest_level != new_finest_level) && (m_verbosity >= 2))
        {
          pout() << "finest level changes here from "
                 << m_finest_level << " to "
                 << new_finest_level << endl;
        }

      //allow for levels to change
      m_finest_level = Min(new_finest_level, m_max_level);

      //need to assign times if number of levels has grown
      if (m_finest_level > m_finest_level_old)
        {
          Real ratTime = m_amrlevels[m_finest_level_old]->time();

          for (int level = m_finest_level_old+1; level <= m_finest_level; ++level)
            {
              m_amrlevels[level]->time(ratTime);
            }
        }

//    broadcast(new_grids,uniqueProc(SerialTask::compute));

      if (m_verbosity >= 4)
        {
          if (new_grids.size() == 0)
            {
              pout() << "No new_grids" << endl;
              pout() << endl;
            }
          else
            {
              for (int level = a_base_level; level <= m_finest_level; ++level)
                {
                  pout() << "new_grids[" << level << "]: " << endl;

                  for (int i = 0; i< new_grids[level].size(); ++i)
                    {
                      pout() << "  " << i << ": " << new_grids[level][i] << endl;
                    }

                  pout() << endl;
                }
            }
        }

    }
  else
    {
      // use pre-defined grids
      new_grids = m_amr_grids;
    }

  // before regridding (but after tagging) allow for pre-regridding ops
  // (tjl 8/23/06 - this needed for mapped grid computations and a default
  // implementation is provided in AMRLevel to preserve existing codes)
  for (int level = m_finest_level; level >= a_base_level; --level)
    {
      m_amrlevels[level]->preRegrid(a_base_level);
    }

  for (int level = a_base_level + 1; level <= m_finest_level; ++level)
    {
      m_amrlevels[level]->regrid(new_grids[level]);
    }

  for (int level = m_finest_level + 1; level <= m_max_level; ++level)
    {
      m_amrlevels[level]->regrid(Vector<Box>());
    }

  // now that the new hierarchy is defined, do post-regridding ops
  // (dfm 8/26/05 -- call postRegrid on base_level as well, to 
  // cover the case where all levels finer than base_level are removed)
  for (int level = m_finest_level; level >= a_base_level; --level)
    {
      m_amrlevels[level]->postRegrid(a_base_level);
    }
}

Vector<AMRLevel*> AMR::getAMRLevels()
{
  return m_amrlevels;
}

void AMR::initialTime(Real a_initialTime)
{
  // this needs to be called before setup, but after define
  CH_assert(m_isDefined);
  CH_assert(!m_isSetUp);

  m_cur_time = a_initialTime;

  // propagate this time to all of the AMR levels
  for (int lev=0; lev<m_amrlevels.size(); lev++)
    {
      m_amrlevels[lev]->time(a_initialTime);
    }
}
  

Real AMR::getCurrentTime() const
{
  return m_cur_time;
}

void AMR::makeBaseLevelMesh(Vector<Box>& a_grids) const
{
  CH_assert(m_isDefined);

  if (m_max_base_grid_size == 0)
    {
      // define base level to be single grid
      a_grids.resize(1);
      a_grids[0] = m_amrlevels[0]->problemDomain().domainBox();
    }
  else
    {
      if (m_max_base_grid_size < m_blockFactor)
        {
          MayDay::Abort("Base grid size must be greater than blocking factor");
        }

      // chop base level up into grids of no more than m_max_grid_size on a side
      // first coarsen to enforce the blocking factor
      ProblemDomain problem_domain = coarsen(m_amrlevels[0]->problemDomain(),
                                             m_blockFactor);
      int max_grid_size = m_max_base_grid_size / m_blockFactor;
      Tuple<Vector<int>,SpaceDim> box_sizes;
      IntVect num_grids;
      IntVect base_size;

      for (int d = 0; d < SpaceDim; ++d)
        {
          int num_div = 1;
          int domain_size = problem_domain.domainBox().size(d);
          while (num_div * max_grid_size < domain_size)
            {
              ++num_div;
            }

          // int(x/y) +(x%y)?1:0 is integer division with rounding upwards, for x,y>0
          base_size[d] = int(domain_size/num_div) +((domain_size%num_div) ? 1 : 0);
          box_sizes[d].resize(num_div, base_size[d]);
          box_sizes[d][num_div-1] = domain_size -(num_div - 1) * base_size[d];
          num_grids[d] = num_div;
        }

      Box b(IntVect::Zero,num_grids - IntVect::Unit);
      const IntVect& domain_hi = problem_domain.domainBox().bigEnd();
      const IntVect& domain_lo = problem_domain.domainBox().smallEnd();

      BoxIterator bit(b);
      for (bit.begin(); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          IntVect lo = domain_lo + iv * base_size;
          IntVect hi = min(lo + base_size - IntVect::Unit, domain_hi);
          Box grid(lo,hi);

          a_grids.push_back(refine(grid, m_blockFactor) );
        }
    }
}

// generate initial grid hierarchy
void AMR::initialGrid()
{
  CH_assert(m_isDefined);

  if (m_verbosity >= 3)
    {
      pout() << "AMR::initialGrid" << endl;
    }

  Vector<Vector<Box> > old_grids(1);
  Vector<Vector<Box> > new_grids;

  makeBaseLevelMesh(old_grids[0]);

  // we will keep the old tags around while we generate new ones
  Vector<IntVectSet> old_tags(m_max_level);

  if (m_use_meshrefine)
    {
      for (int top_level = 0;
           top_level < Min(m_finest_level+1, m_max_level);
           ++top_level)
        {
          if (m_verbosity >= 2)
            {
              pout() << "AMR::initialGrid: top level = " << top_level << endl;
            }

          Vector<IntVectSet> tags(top_level+1);
          Vector<ProblemDomain> problem_domains(top_level+1);

          // loop over all levels and initialize grids,
          // then loop over all levels and initialize data
          // then loop over levels and generate tags
          // do this in three separate loops to handle
          // the case where initial data is generated
          // using a multilevel operation(for instance,
          // computing initial velocity from initial
          // vorticity through a multilevel elliptic solve)
          // DFM(11/28/2000)
          for (int level = 0; level <= top_level; ++level)
            {
              m_amrlevels[level]->initialGrid(old_grids[level]);
            }

          for (int level = 0; level <= top_level; ++level)
            {
              m_amrlevels[level]->initialData();
            }

          for (int level = 0; level <= top_level; ++level)
            {
              m_amrlevels[level]->tagCellsInit(tags[level]);

              // union old tags with current ones.  this prevents
              // you from unrefining a region you previously
              // decided you wanted refined
              tags[level] |= old_tags[level];

              problem_domains[level] = m_amrlevels[level]->problemDomain();

              if (m_verbosity >= 3)
                {
                  pout() << "AMR::initialGrid: problem domain["
                         << level << "]: " << problem_domains[level] << endl;
                }

              if (m_verbosity >= 4)
                {
                  pout() << "AMR::initialGrid: old_grids[" << level << "]: " << endl;

                  for (int i = 0; i< old_grids[level].size(); ++i)
                    {
                      pout() << "  " << i << ": " << old_grids[level][i] << endl;
                    }

                  if (m_verbosity >= 5)
                    {
                      pout() << "AMR::initialGrid: tags[" << level << "]: "
                             << tags[level] << endl;
                    }
                }
            }
 //       if(procID() == uniqueProc(SerialTask::compute))
//          {
          m_finest_level = m_mesh_refine.regrid(new_grids,
                                                tags,
                                                0,
                                                top_level,
                                                old_grids);
//          }
//        broadcast(m_finest_level, uniqueProc(SerialTask::compute));
//        broadcast(new_grids,      uniqueProc(SerialTask::compute));

          // do this only if a new level was generated
          if (m_finest_level > top_level)
            {
              old_grids = new_grids;

              // copy current tags to old_tags
              for (int lev = 0; lev <= top_level; lev++)
                {
                  old_tags[lev] = tags[lev];
                }
            }

          if (m_verbosity >= 4)
            {
              for (int level = 0; level < new_grids.size(); ++level)
                {
                  pout() << "new_grids[" << level << "]: " << endl;

                  for (int i = 0; i< new_grids[level].size(); ++i)
                    {
                      pout() << "  " << i << ": " << new_grids[level][i] << endl;
                    }

                  pout() << endl;
                }
            }
        }

      if (m_finest_level == 0)
        {
          new_grids.resize(1);
          new_grids[0] = old_grids[0];
        }
    }
  else
    {
      // use predefined grids
      CH_assert(m_amr_grids.size() > 0);

      new_grids = m_amr_grids;
      m_finest_level = m_amr_grids.size()-1;
    }

  // separate loops for initialGrid and initialData
  // to ensure that hierarchy is defined before we do initialization
  // (for case where initialization is a multilevel operatation)
  // (dfm 11/7/02)
  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_amrlevels[level]->initialGrid(new_grids[level]);
    }

  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_amrlevels[level]->initialData();
    }

  for (int level = m_finest_level + 1; level <= m_max_level; ++level)
    {
      // m_amrlevels[level]->regrid(Vector<Box>());
      m_amrlevels[level]->initialGrid(Vector<Box>());
      m_amrlevels[level]->initialData();
    }

  // call post-initialize once all the levels have been defined
  for (int level = m_finest_level; level >= 0; --level)
    {
      m_amrlevels[level]->postInitialize();
    }
}

void AMR::writePlotFile() const
{
  CH_assert(m_isDefined);

  if (m_verbosity >= 3)
    {
      pout() << "AMR::writePlotFile" << endl;
    }

#ifdef CH_USE_HDF5
  char iter_str[80];

  sprintf(iter_str,
          "%s%04d.%dd.hdf5",
          m_plotfile_prefix.c_str(), m_cur_step, SpaceDim);

  if (m_verbosity >= 2)
    {
      pout() << "plot file name = " << iter_str << endl;
    }

  HDF5Handle handle(iter_str, HDF5Handle::CREATE);

  // write amr data
  HDF5HeaderData header;
  header.m_int ["max_level"]  = m_max_level;
  header.m_int ["num_levels"] = m_finest_level + 1;
  header.m_int ["iteration"]  = m_cur_step;
  header.m_real["time"]       = m_cur_time;

  // should steps since regrid be in the checkpoint file?
  header.writeToFile(handle);

  if (m_verbosity >= 3)
    {
      pout() << header << endl;
    }

  // write physics class header data
  m_amrlevels[0]->writePlotHeader(handle);

  // write physics class per-level data
  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_amrlevels[level]->writePlotLevel(handle);
    }

  handle.close();
#endif
}

void AMR::writeCheckpointFile() const
{
  CH_assert(m_isDefined);

  if (m_verbosity >= 3)
    {
      pout() << "AMR::writeCheckpointFile" << endl;
    }

#ifdef CH_USE_HDF5
  char iter_str[100];

  sprintf(iter_str,
          "%s%d.%dd.hdf5",
          m_checkpointfile_prefix.c_str(), m_cur_step, SpaceDim );

  if (m_verbosity >= 2)
    {
      pout() << "checkpoint file name = " << iter_str << endl;
    }

  HDF5Handle handle(iter_str, HDF5Handle::CREATE);

  // write amr data
  HDF5HeaderData header;
  header.m_int ["max_level"]  = m_max_level;
  header.m_int ["num_levels"] = m_finest_level + 1;
  header.m_int ["iteration"]  = m_cur_step;
  header.m_real["time"]       = m_cur_time;

  for (int level = 0; level < m_regrid_intervals.size(); ++level)
    {
      char headername[100];
      sprintf(headername, "regrid_interval_%d", level);
      header.m_int[headername] = m_regrid_intervals[level];
    }

  // should steps since regrid be in the checkpoint file?
  header.writeToFile(handle);

  if (m_verbosity >= 3)
    {
      pout() << header << endl;
    }

  // write physics class data
  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_amrlevels[level]->writeCheckpointHeader(handle);
      m_amrlevels[level]->writeCheckpointLevel(handle);
    }

  handle.close();
#endif
}

void AMR::verbosity(int a_verbosity)
{
  m_verbosity = a_verbosity;

  if (m_verbosity >= 2)
    {
      pout() << " AMR::verbosity  m_verbosity=" << m_verbosity << endl;
    }
}

int AMR::verbosity() const
{
  return(m_verbosity);
}

void AMR::maxDtGrow(Real a_dtGrowFactor)
{
  m_maxDtGrow = a_dtGrowFactor;
}

Real AMR::maxDtGrow() const
{
  return m_maxDtGrow;
}

void AMR::timeEps(Real a_timeEps)
{
  m_time_eps = a_timeEps;
}

Real AMR::timeEps() const
{
  return m_time_eps;
}

void AMR::fixedDt(Real a_dt)
{
  m_fixedDt = a_dt;
}

Real AMR::fixedDt() const
{
  return m_fixedDt;
}

bool AMR::isDefined() const
{
  return m_isDefined;
}

bool AMR::isSetUp() const
{
  return m_isSetUp;
}

#ifdef CH_USE_TIMER
Timer * AMR::timer(Timer *a_timer )
{
  Timer * old_timer = m_timer ;
  if( a_timer != NULL ){
    m_timer = a_timer ;
  }
  return old_timer ;
}
#endif
#include "NamespaceFooter.H"
