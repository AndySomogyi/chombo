#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <iomanip>

#include "parstream.H"
#include "ParmParse.H"
#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "LevelFluxRegister.H"
#include "PiecewiseLinearFillPatch.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "FluxBox.H"

#include "AMRIO.H"

#include "AMRLevel.H"
#include "AMRLevelUpwind.H"
#include "UpwindF_F.H"

#include "computeSum.H"


#define TIME_EPS 1.0e-10


///
const LevelData<FArrayBox>&
AMRLevelUpwind::getStateNew() const
{
  return m_UNew;
}

///
const LevelData<FArrayBox>&
AMRLevelUpwind::getStateOld() const
{
  return m_UOld;
}
// Constructor
AMRLevelUpwind::AMRLevelUpwind()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelUpwind default constructor" << endl;
  }

  m_cfl = 0.8;
  m_domainLength = 1.0;
  m_refineThresh = 0.2;
  m_initial_dt_multiplier = 1.0;
}

// Destructor
AMRLevelUpwind::~AMRLevelUpwind()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelUpwind destructor" << endl;
  }

}

// This instance should never get called - historical
void AMRLevelUpwind::define(AMRLevel*  a_coarserLevelPtr,
                                   const Box& a_problemDomain,
                                   int        a_level,
                                   int        a_refRatio)
{
  ProblemDomain physdomain(a_problemDomain);

  MayDay::Error("AMRLevelUpwind::define -\n\tShould never be called with a Box for a problem domain");
}

// Define new AMR level
void AMRLevelUpwind::define(AMRLevel*            a_coarserLevelPtr,
                            const ProblemDomain& a_problemDomain,
                            int                  a_level,
                            int                  a_refRatio)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelUpwind::define " << a_level << endl;
    }
  
  // for now, hardwire this to one
  m_numStates = 1;
  m_stateNames.resize(m_numStates, "U");


  // Call inherited define
  AMRLevel::define(a_coarserLevelPtr,
                   a_problemDomain,
                   a_level,
                   a_refRatio);

  // Get setup information from the next coarser level
  if (a_coarserLevelPtr != NULL)
  {
    AMRLevelUpwind* amrGodPtr = dynamic_cast<AMRLevelUpwind*>(a_coarserLevelPtr);
    
    if (amrGodPtr != NULL)
      {
        m_cfl = amrGodPtr->m_cfl;
        m_domainLength = amrGodPtr->m_domainLength;
        m_refineThresh = amrGodPtr->m_refineThresh;
        m_tagBufferSize = amrGodPtr->m_tagBufferSize;
      }
    else
      {
        MayDay::Error("AMRLevelUpwind::define: a_coarserLevelPtr is not castable to AMRLevelUpwind*");
      }
  }
  
  // Compute the grid spacing
  m_dx = m_domainLength / a_problemDomain.domainBox().longside();
  
  // Nominally, one layer of ghost cells is maintained permanently and
  // individual computations may create local data with more
  m_numGhost = 1;

}

// Advance by one timestep
Real AMRLevelUpwind::advance()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelUpwind::advance level " 
           << m_level << " to time " << m_time + m_dt << endl;
  }
  
  // Copy the new to the old
  DataIterator dit = m_UNew.dataIterator();
  for ( ; dit.ok(); ++dit)
    {
      m_UOld[dit()].copy(m_UNew[dit()]);
    }
  
  // fill in ghost cells, if necessary
  
  // for now, only do periodic -- this is a simple example, so we'll
  // neglect physical BC's
  // assert that this is true
  for (int dir=0; dir<SpaceDim; dir++)
    {
      CH_assert(m_problem_domain.isPeriodic(dir));
    }
  
  AMRLevelUpwind* coarserLevelPtr = NULL;
  // interpolate from coarser level, if appropriate
  if (m_level > 0)
    {
      coarserLevelPtr = getCoarserLevel();
      
      // get old and new coarse-level data
      LevelData<FArrayBox>& crseDataOld = coarserLevelPtr->m_UOld;
      LevelData<FArrayBox>& crseDataNew = coarserLevelPtr->m_UNew;
      const DisjointBoxLayout& crseGrids = crseDataNew.getBoxes();

      Real newCrseTime = coarserLevelPtr->m_time;
      Real oldCrseTime = newCrseTime - coarserLevelPtr->m_dt;
      Real coeff = (m_time - oldCrseTime)/coarserLevelPtr->m_dt;

      if (abs(coeff) < TIME_EPS) coeff = 0;
      if (abs(coeff-1.0) < TIME_EPS) coeff = 1.0;

      const ProblemDomain& crseDomain = coarserLevelPtr->m_problem_domain;
      int nRefCrse = coarserLevelPtr->refRatio();
      int nGhost = 1;
      
      PiecewiseLinearFillPatch filpatcher(m_grids, crseGrids,
                                         m_UNew.nComp(), crseDomain,
                                         nRefCrse, nGhost);
                                         
      
      filpatcher.fillInterp(m_UOld, crseDataOld, crseDataNew, 
                            coeff, 0, 0, m_UNew.nComp());
                            
    }

  // exchange copies overlapping ghost cells on this level
  m_UOld.exchange();

  // now go patch-by-patch, compute upwind flux, and do update
  // iterator will only reference patches on this processor
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box gridBox = m_grids.get(dit());
      FArrayBox& thisOldSoln = m_UOld[dit];
      FArrayBox& thisNewSoln = m_UNew[dit];

      FluxBox fluxes(gridBox, thisOldSoln.nComp());

      // loop over directions
      for (int dir=0; dir<SpaceDim; dir++)
        {
          // note that gridbox will need to be the face-centered
          // grid box 
          Box faceBox = fluxes[dir].box();

          FORT_UPWIND(CHF_FRA(fluxes[dir]),
                      CHF_FRA(thisOldSoln),
                      CHF_REALVECT(m_advectionVel),
                      CHF_REAL(m_dt),
                      CHF_REAL(m_dx),
                      CHF_BOX(faceBox),
                      CHF_INT(dir));
          

                                             
          // increment flux registers with fluxes
          Interval UInterval = m_UNew.interval();
          
          if (m_hasFiner)
            {
              // this level's flux register goes between this level 
              // and the next finer level
              m_fluxRegister.incrementCoarse(fluxes[dir], m_dt, dit(),
                                             UInterval, UInterval, dir);
            }
          
          if (m_level > 0)
            {
              // 
              LevelFluxRegister& crseFluxReg = coarserLevelPtr->m_fluxRegister;

              crseFluxReg.incrementFine(fluxes[dir], m_dt, dit(), UInterval,
                                        UInterval, dir, Side::Lo);
              crseFluxReg.incrementFine(fluxes[dir], m_dt, dit(), UInterval,
                                        UInterval, dir, Side::Hi);
            }
                      
        } // end loop over directions

      // do flux difference to increment solution

      thisNewSoln.copy(thisOldSoln);

      for (int dir=0; dir<SpaceDim; dir++)
        {
          FORT_INCREMENTDIVDIR(CHF_FRA(thisNewSoln),
                               CHF_FRA(fluxes[dir]),
                               CHF_BOX(gridBox),
                               CHF_REAL(m_dx),
                               CHF_REAL(m_dt),
                               CHF_INT(dir));
        } 

    } // end loop over grid boxes



  // Update the time and store the new timestep
  m_time += m_dt;

  // kind of a hack, but since advection velocity is constant, 
  // no harm, no foul...
  return m_dt;
}

// Things to do after a timestep
void AMRLevelUpwind::postTimeStep()
{
  // Used for conservation tests
  static Real orig_integral = 0.0;
  static Real last_integral = 0.0;
  static bool first = true;

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelUpwind::postTimeStep " << m_level << endl;
  }

  if (m_hasFiner)
  {
    // Reflux
    Real scale = -1.0/m_dx;
    m_fluxRegister.reflux(m_UNew,scale);

    // Average from finer level data
    AMRLevelUpwind* finerLevelPtr = getFinerLevel();

    finerLevelPtr->m_coarseAverage.averageToCoarse(m_UNew,
                                                   finerLevelPtr->m_UNew);
  }

  if (s_verbosity >= 2 && m_level == 0)
  {
    int nRefFine = 1;

    pout() << "AMRLevelUpwind::postTimeStep:" << endl;
    pout() << "  Sums:" << endl;
    for (int comp = 0; comp < m_numStates; comp++)
      {
      Interval curComp(comp,comp);
      Real integral = computeSum(m_UNew,NULL,nRefFine,m_dx,curComp);

      pout() << "    " << setw(23)
                       << setprecision(16)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << integral
             << " --- " << m_stateNames[comp];

      if (comp == 0 && !first)
        {
        pout() << " (" << setw(23)
                       << setprecision(16)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << (integral-last_integral)
               << " " << setw(23)
                      << setprecision(16)
                      << setiosflags(ios::showpoint)
                      << setiosflags(ios::scientific)
                      << (integral-orig_integral)
               << ")";
      }

      pout() << endl;

      if (comp == 0)
      {
        if (first)
        {
          orig_integral = integral;
          first = false;
        }

        last_integral = integral;
      }
    }
  }

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelUpwind::postTimeStep " << m_level << " finished" << endl;
  }
}

// Create tags for regridding
void AMRLevelUpwind::tagCells(IntVectSet& a_tags)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelUpwind::tagCells " << m_level << endl;
    }
  
  // Create tags based on undivided gradient of density
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
  IntVectSet localTags;
  // If there is a coarser level interpolate undefined ghost cells
  if (m_hasCoarser)
    {
      const AMRLevelUpwind* amrGodCoarserPtr = getCoarserLevel();
      
      PiecewiseLinearFillPatch pwl(levelDomain,
                                   amrGodCoarserPtr->m_UNew.disjointBoxLayout(),
                                   m_numStates,
                                   amrGodCoarserPtr->m_problem_domain,
                                   amrGodCoarserPtr->m_ref_ratio,
                                   1);
      
      pwl.fillInterp(m_UNew,
                     amrGodCoarserPtr->m_UNew,
                     amrGodCoarserPtr->m_UNew,
                     1.0,
                     0,
                     0,
                     m_numStates);
    }
  m_UNew.exchange(Interval(0,m_numStates-1));
  
  // Compute relative gradient
  DataIterator dit = levelDomain.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& b = levelDomain[dit()];
      FArrayBox gradFab(b,SpaceDim);
      const FArrayBox& UFab = m_UNew[dit()];
      
      for (int dir = 0; dir < SpaceDim; ++dir)
        {
          const Box bCenter = b & grow(m_problem_domain,-BASISV(dir));
          const Box bLo     = b & adjCellLo(bCenter,dir);
          const int hasLo = ! bLo.isEmpty();
          const Box bHi     = b & adjCellHi(bCenter,dir);
          const int hasHi = ! bHi.isEmpty();
          FORT_GETRELGRADF(CHF_FRA1(gradFab,dir),
                           CHF_CONST_FRA1(UFab,0),
                           CHF_CONST_INT(dir),
                           CHF_BOX(bLo),
                           CHF_CONST_INT(hasLo),
                           CHF_BOX(bHi),
                           CHF_CONST_INT(hasHi),
                           CHF_BOX(bCenter));
        }
      
      FArrayBox gradMagFab(b,1);
      FORT_OLDMAGNITUDEF(CHF_FRA1(gradMagFab,0),
                         CHF_CONST_FRA(gradFab),
                         CHF_BOX(b));
      
      // Tag where gradient exceeds threshold
      BoxIterator bit(b);
      for (bit.begin(); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          
          if (gradMagFab(iv) >= m_refineThresh)
            {
              localTags |= iv;
            }
        }
      
    }
  
  localTags.grow(m_tagBufferSize);
  
  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = localTags.minBox();
  localTagsBox &= m_problem_domain;
  localTags &= localTagsBox;
  
  a_tags = localTags;
}

// Create tags at initialization
void AMRLevelUpwind::tagCellsInit(IntVectSet& a_tags)
{
  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelUpwind::tagCellsInit " << m_level << endl;
    }
  
  tagCells(a_tags);
}

// Set up data on this level after regridding
void AMRLevelUpwind::regrid(const Vector<Box>& a_newGrids)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelUpwind::regrid " << m_level << endl;
    }
  
  // Save original grids and load balance
  m_level_grids = a_newGrids;
  m_grids = loadBalance(a_newGrids);
  
  if (s_verbosity >= 4)
    {
      // Indicate/guarantee that the indexing below is only for reading
      // otherwise an error/assertion failure occurs
      const DisjointBoxLayout& constGrids = m_grids;
      
      pout() << "new grids: " << endl;
      
      for (LayoutIterator lit = constGrids.layoutIterator(); lit.ok(); ++lit)
        {
          pout() << constGrids[lit()] << endl;
        }
    }
  
  // Save data for later
  DataIterator dit = m_UNew.dataIterator();
  for ( ; dit.ok(); ++dit)
    {
      m_UOld[dit()].copy(m_UNew[dit()]);
    }
  
  // Reshape state with new grids
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_UNew.define(m_grids,m_numStates,ivGhost);
  
  // Set up data structures
  levelSetup();
  
  // Interpolate from coarser level
  if (m_hasCoarser)
    {
      AMRLevelUpwind* amrGodCoarserPtr = getCoarserLevel();
      m_fineInterp.interpToFine(m_UNew,amrGodCoarserPtr->m_UNew);
    }
  
  // Copy from old state
  m_UOld.copyTo(m_UOld.interval(),
                m_UNew,
                m_UNew.interval());
  
  m_UOld.define(m_grids,m_numStates,ivGhost);
}

// Initialize grids
void AMRLevelUpwind::initialGrid(const Vector<Box>& a_newGrids)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelUpwind::initialGrid " << m_level << endl;
    }
  
  // Save original grids and load balance
  m_level_grids = a_newGrids;
  m_grids = loadBalance(a_newGrids);
  
  if (s_verbosity >= 4)
    {
      // Indicate/guarantee that the indexing below is only for reading
      // otherwise an error/assertion failure occurs
      const DisjointBoxLayout& constGrids = m_grids;
      
      pout() << "new grids: " << endl;
      for (LayoutIterator lit = constGrids.layoutIterator(); lit.ok(); ++lit)
        {
          pout() << constGrids[lit()] << endl;
        }
    }
  
  // Define old and new state data structures
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_UNew.define(m_grids,m_numStates,ivGhost);
  m_UOld.define(m_grids,m_numStates,ivGhost);
  
  // Set up data structures
  levelSetup();
}

// Initialize data
void AMRLevelUpwind::initialData()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelUpwind::initialData " << m_level << endl;
    }
  
  DataIterator dit = m_UNew.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FORT_INITSCAL(CHF_FRA(m_UNew[dit]),
                    CHF_REAL(m_dx));
    }
  
}

// Things to do after initialization
void AMRLevelUpwind::postInitialize()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelUpwind::postInitialize " << m_level << endl;
    }
}

#ifdef CH_USE_HDF5

// Write checkpoint header
void AMRLevelUpwind::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelUpwind::writeCheckpointHeader" << endl;
  }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numStates;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = m_stateNames[comp];
  }

  // Write the header
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }
}

// Write checkpoint data for this level
void AMRLevelUpwind::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelUpwind::writeCheckpointLevel" << endl;
  }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]       = m_ref_ratio;
  header.m_int ["tag_buffer_size"] = m_tagBufferSize;
  header.m_real["dx"]              = m_dx;
  header.m_real["dt"]              = m_dt;
  header.m_real["time"]            = m_time;
  header.m_box ["prob_domain"]     = m_problem_domain.domainBox();

  // Setup the periodicity info
  D_TERM(
         if (m_problem_domain.isPeriodic(0))
           header.m_int ["is_periodic_0"] = 1;
         else
           header.m_int ["is_periodic_0"] = 0; ,

         if (m_problem_domain.isPeriodic(1))
           header.m_int ["is_periodic_1"] = 1;
         else
           header.m_int ["is_periodic_1"] = 0; ,

         if (m_problem_domain.isPeriodic(2))
           header.m_int ["is_periodic_2"] = 1;
         else
           header.m_int ["is_periodic_2"] = 0; );

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }

  // Write the data for this level
  write(a_handle,m_UNew.boxLayout());
  write(a_handle,m_UNew,"data");
}

// Read checkpoint header
void AMRLevelUpwind::readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelUpwind::readCheckpointHeader" << endl;
  }

  // Reader the header
  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << "hdf5 header data:" << endl;
    pout() << header << endl;
  }

  // Get the number of components
  if (header.m_int.find("num_components") == header.m_int.end())
  {
    MayDay::Error("AMRLevelUpwind::readCheckpointHeader: checkpoint file does not have num_components");
  }

  int numStates = header.m_int["num_components"];
  if (numStates != m_numStates)
  {
    MayDay::Error("AMRLevelUpwind::readCheckpointHeader: num_components in checkpoint file does not match solver");
  }

  // Get the component names
  std::string stateName;
  char compStr[60];
  for (int comp = 0; comp < m_numStates; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    if (header.m_string.find(compStr) == header.m_string.end())
    {
      MayDay::Error("AMRLevelUpwind::readCheckpointHeader: checkpoint file does not have enough component names");
    }

    stateName = header.m_string[compStr];
    if (stateName != m_stateNames[comp])
    {
      MayDay::Error("AMRLevelUpwind::readCheckpointHeader: state_name in checkpoint does not match solver");
    }
  }
}

// Read checkpoint data for this level
void AMRLevelUpwind::readCheckpointLevel(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelUpwind::readCheckpointLevel" << endl;
  }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  // Read the header for this level
  a_handle.setGroup(label);

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << "hdf5 header data:" << endl;
    pout() << header << endl;
  }

  // Get the refinement ratio
  if (header.m_int.find("ref_ratio") == header.m_int.end())
  {
    MayDay::Error("AMRLevelUpwind::readCheckpointLevel: file does not contain ref_ratio");
  }
  m_ref_ratio = header.m_int["ref_ratio"];

  if (s_verbosity >= 2)
  {
    pout() << "read ref_ratio = " << m_ref_ratio << endl;
  }

  // Get the tag buffer size
  if (header.m_int.find("tag_buffer_size") == header.m_int.end())
  {
    MayDay::Error("AMRLevelUpwind::readCheckpointLevel: file does not contain tag_buffer_size");
  }
  m_tagBufferSize=  header.m_int["tag_buffer_size"];

  if (s_verbosity >= 2)
  {
    pout() << "read tag_buffer_size = " << m_tagBufferSize << endl;
  }

  // Get dx
  if (header.m_real.find("dx") == header.m_real.end())
  {
    MayDay::Error("AMRLevelUpwind::readCheckpointLevel: file does not contain dx");
  }
  m_dx = header.m_real["dx"];

  if (s_verbosity >= 2)
  {
    pout() << "read dx = " << m_dx << endl;
  }

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
  {
    MayDay::Error("AMRLevelUpwind::readCheckpointLevel: file does not contain dt");
  }
  m_dt = header.m_real["dt"];

  if (s_verbosity >= 2)
  {
    pout() << "read dt = " << m_dt << endl;
  }

  // Get time
  if (header.m_real.find("time") == header.m_real.end())
  {
    MayDay::Error("AMRLevelUpwind::readCheckpointLevel: file does not contain time");
  }
  m_time = header.m_real["time"];

  if (s_verbosity >= 2)
  {
    pout() << "read time = " << m_time << endl;
  }

  // Get the problem domain
  if (header.m_box.find("prob_domain") == header.m_box.end())
  {
    MayDay::Error("AMRLevelUpwind::readCheckpointLevel: file does not contain prob_domain");
  }

  Box domainBox = header.m_box["prob_domain"];

  // Get the periodicity info -- this is more complicated than it really
  // needs to be in order to preserve backward compatibility
  bool isPeriodic[SpaceDim];
  D_TERM(if (!(header.m_int.find("is_periodic_0") == header.m_int.end()))
           isPeriodic[0] =  (header.m_int["is_periodic_0"] == 1);
         else
           isPeriodic[0] = false; ,

         if (!(header.m_int.find("is_periodic_1") == header.m_int.end()))
           isPeriodic[1] =  (header.m_int["is_periodic_1"] == 1);
         else
           isPeriodic[1] = false; ,

         if (!(header.m_int.find("is_periodic_2") == header.m_int.end()))
           isPeriodic[2] =  (header.m_int["is_periodic_2"] == 1);
         else
           isPeriodic[2] = false;);

  m_problem_domain = ProblemDomain(domainBox,isPeriodic);

  // Get the grids
  Vector<Box> grids;
  const int gridStatus = read(a_handle,grids);

  if (gridStatus != 0)
  {
    MayDay::Error("AMRLevelUpwind::readCheckpointLevel: file does not contain a Vector<Box>");
  }

  // Create level domain
  m_grids = loadBalance(grids);

  // Indicate/guarantee that the indexing below is only for reading
  // otherwise an error/assertion failure occurs
  const DisjointBoxLayout& constGrids = m_grids;

  LayoutIterator lit = constGrids.layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
  {
    const Box& b = constGrids[lit()];
    m_level_grids.push_back(b);
  }

  if (s_verbosity >= 4)
  {
    pout() << "read level domain: " << endl;
    LayoutIterator lit = m_grids.layoutIterator();
    for (lit.begin(); lit.ok(); ++lit)
    {
      const Box& b = m_grids[lit()];
      pout() << lit().intCode() << ": " << b << endl;
    }
    pout() << endl;
  }

  // Reshape state with new grids
  m_UNew.define(m_grids,m_numStates);
  const int dataStatus = read<FArrayBox>(a_handle,
                                         m_UNew,
                                         "data",
                                         m_grids);

  if (dataStatus != 0)
  {
    MayDay::Error("AMRLevelUpwind::readCheckpointLevel: file does not contain state data");
  }
  m_UOld.define(m_grids,m_numStates);

  // Set up data structures
  levelSetup();
}

// Write plotfile header
void AMRLevelUpwind::writePlotHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelUpwind::writePlotHeader" << endl;
  }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numStates;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = m_stateNames[comp];
  }

  // Write the header
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }
}

// Write plotfile data for this level
void AMRLevelUpwind::writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelUpwind::writePlotLevel" << endl;
  }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]   = m_ref_ratio;
  header.m_real["dx"]          = m_dx;
  header.m_real["dt"]          = m_dt;
  header.m_real["time"]        = m_time;
  header.m_box ["prob_domain"] = m_problem_domain.domainBox();

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }

  // Write the data for this level
  write(a_handle,m_UNew.boxLayout());
  write(a_handle,m_UNew,"data");
}

#endif

// Returns the dt computed earlier for this level
Real AMRLevelUpwind::computeDt()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelUpwind::computeDt " << m_level << endl;
  }


  Real maxVel = 0;
  for (int dir=0; dir<SpaceDim; dir++)
    {
      if (m_advectionVel[dir]> maxVel) maxVel = m_advectionVel[dir];
    }

  Real newDt;
  newDt = m_dx/maxVel;
  
  return newDt;
}

// Compute dt using initial data
Real AMRLevelUpwind::computeInitialDt()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelUpwind::computeInitialDt " << m_level << endl;
  }

  
  Real newDt = computeDt();
  
  newDt *= m_initial_dt_multiplier;
  
  return newDt;
}

// Set the CFL number
void AMRLevelUpwind::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
}

// Set the physical dimension of the longest side of the domain
void AMRLevelUpwind::domainLength(Real a_domainLength)
{
  m_domainLength = a_domainLength;
}

// Set the refinement threshold
void AMRLevelUpwind::refinementThreshold(Real a_refineThresh)
{
  m_refineThresh = a_refineThresh;
}

// Set the tag buffer size
void AMRLevelUpwind::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
}

// Create a load-balanced DisjointBoxLayout from a collection of Boxes
DisjointBoxLayout AMRLevelUpwind::loadBalance(const Vector<Box>& a_grids)
{
  // Load balance and create boxlayout
  Vector<int> procMap;
  //if (procID() == uniqueProc(SerialTask::compute))
  //{
  //  LoadBalance(procMap,a_grids);
  //}
  //broadcast(procMap,uniqueProc(SerialTask::compute));
  
  // appears to be faster for all procs to do the loadbalance (ndk)
  LoadBalance(procMap,a_grids);
  
  if (s_verbosity >= 4)
    {
      pout() << "AMRLevelUpwind::loadBalance: procesor map: " << endl;
      for (int igrid = 0; igrid < a_grids.size(); ++igrid)
        {
          pout() << igrid << ": " << procMap[igrid] << "  " << endl;
        }
      pout() << endl;
    }
  
  DisjointBoxLayout dbl(a_grids,procMap,m_problem_domain);
  dbl.close();
  
  return dbl;
}

// Setup menagerie of data structures
void AMRLevelUpwind::levelSetup()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelUpwind::levelSetup " << m_level << endl;
    }
  
  AMRLevelUpwind* amrGodCoarserPtr = getCoarserLevel();
  AMRLevelUpwind* amrGodFinerPtr   = getFinerLevel();
  
  m_hasCoarser = (amrGodCoarserPtr != NULL);
  m_hasFiner   = (amrGodFinerPtr   != NULL);
  
  if (m_hasCoarser)
    {
      int nRefCrse = m_coarser_level_ptr->refRatio();
      
      m_coarseAverage.define(m_grids,
                             m_numStates,
                             nRefCrse);
      
      m_fineInterp.define(m_grids,
                          m_numStates,
                          nRefCrse,
                          m_problem_domain);
      
      
      //const DisjointBoxLayout& coarserLevelDomain = amrGodCoarserPtr->m_grids;
      
      
      // This may look twisted but you have to do this this way because the
      // coarser levels get setup before the finer levels so, since a flux
      // register lives between this level and the next FINER level, the finer
      // level has to do the setup because it is the only one with the
      // information at the time of construction.
      
      // Maintain flux registers
      amrGodCoarserPtr->m_fluxRegister.define(m_grids,
                                              amrGodCoarserPtr->m_grids,
                                              m_problem_domain,
                                              amrGodCoarserPtr->m_ref_ratio,
                                              m_numStates);
      amrGodCoarserPtr->m_fluxRegister.setToZero();
    }
  
}

// Get the next coarser level
AMRLevelUpwind* AMRLevelUpwind::getCoarserLevel() const
{
  AMRLevelUpwind* amrGodCoarserPtr = NULL;
  
  if (m_coarser_level_ptr != NULL)
    {
      amrGodCoarserPtr = dynamic_cast<AMRLevelUpwind*>(m_coarser_level_ptr);
      
      if (amrGodCoarserPtr == NULL)
        {
          MayDay::Error("AMRLevelUpwind::getCoarserLevel: dynamic cast failed");
        }
    }
  
  return amrGodCoarserPtr;
}

// Get the next finer level
AMRLevelUpwind* AMRLevelUpwind::getFinerLevel() const
{
  AMRLevelUpwind* amrGodFinerPtr = NULL;
  
  if (m_finer_level_ptr != NULL)
    {
      amrGodFinerPtr = dynamic_cast<AMRLevelUpwind*>(m_finer_level_ptr);
      
      if (amrGodFinerPtr == NULL)
        {
          MayDay::Error("AMRLevelUpwind::getFinerLevel: dynamic cast failed");
        }
    }
  
  return amrGodFinerPtr;
}


/// Set the advection velocity
void 
AMRLevelUpwind::advectionVel(const RealVect& a_vel)
{
  m_advectionVel = a_vel;
}

/// accessor for the advection velocity

const RealVect& 
AMRLevelUpwind::advectionVel() const
{
  return m_advectionVel;
}
