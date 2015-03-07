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

#include <string>
#include "parstream.H"

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
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "computeSum.H"
#include "ComputeEnergy.H"
#include "PiecewiseLinearFillPatch.H"

#include "AMRIO.H"

#include "AMRLevel.H"
#include "AMRLevelWaveEquation.H"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Constructor
AMRLevelWaveEquation::AMRLevelWaveEquation()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation default constructor" << endl;
  }

  m_cfl = 0.8;
  m_domainLength = 1.0;
  m_refineThresh = 0.2;
  m_initial_dt_multiplier = 0.1;
}

////////////////////////////////////////////////////////////////////////////////

// Destructor
AMRLevelWaveEquation::~AMRLevelWaveEquation()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation destructor" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////

// Define new AMR level
void AMRLevelWaveEquation::define(AMRLevel*            a_coarserLevelPtr,
                                   const ProblemDomain& a_problemDomain,
                                   int                  a_level,
                                   int                  a_refRatio)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::define " << a_level << endl;
  }

  // Call inherited define
  AMRLevel::define(a_coarserLevelPtr,
                   a_problemDomain,
                   a_level,
                   a_refRatio);

  // Get setup information from the next coarser level
  if (a_coarserLevelPtr != NULL)
  {
    AMRLevelWaveEquation* amrWavePtr = dynamic_cast<AMRLevelWaveEquation*>(a_coarserLevelPtr);

    if (amrWavePtr != NULL)
    {
      m_cfl = amrWavePtr->m_cfl;
      m_domainLength = amrWavePtr->m_domainLength;
      m_refineThresh = amrWavePtr->m_refineThresh;
      m_tagBufferSize = amrWavePtr->m_tagBufferSize;
    }
    else
    {
      MayDay::Error("AMRLevelWaveEquation::define: a_coarserLevelPtr is not castable to AMRLevelWaveEquation*");
    }
  }

  // Compute the grid spacing
  m_dx = m_domainLength / a_problemDomain.domainBox().longside();

  // Nominally, one layer of ghost cells is maintained permanently and
  // individual computations may create local data with more

  // Begin application-dependent code - PC.

  m_numGhost = 1;
  m_numStates = 2;
  m_stateNames.push_back("phi_0");
  m_stateNames.push_back("phi_1");
  m_stateNames.push_back("pi_0");
  m_stateNames.push_back("pi_1");

  // End application-dependent code - PC.

}

////////////////////////////////////////////////////////////////////////////////

// Advance by one timestep
Real AMRLevelWaveEquation::advance()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::advance level " << m_level << " to time " << m_time << endl;
  }

  // Begin application-dependent code - PC.

  // Copy the new to the old
  m_phiNew.copyTo(m_phiNew.interval(),
                m_phiOld,
                m_phiOld.interval());
  m_piNew.copyTo(m_piNew.interval(),
                m_piOld,
                m_piOld.interval());

  // End application-dependent code - PC.

  //XXX -- unused
  //XXXReal newDt = 0.0;

  // Set up arguments to LevelGodunov::step based on whether there are
  // coarser and finer levels

  // Undefined flux register in case we need it
  LevelFluxRegister dummyFR;

  // Undefined leveldata in case we need it
  LevelData<FArrayBox> dummyData;

  // Set arguments to dummy values and then fix if real values are available
  LevelFluxRegister* coarserFR = &dummyFR;
  LevelFluxRegister* finerFR   = &dummyFR;

  Real tCoarserOld = 0.0;
  Real tCoarserNew = 0.0;

  LevelData<FArrayBox>* coarserPhiOld = &dummyData;
  LevelData<FArrayBox>* coarserPhiNew = &dummyData;
  // A coarser level exists
  if (m_hasCoarser)
  {
    AMRLevelWaveEquation* coarserPtr = getCoarserLevel();

    // Recall that my flux register goes between my level and the next
    // finer level
    coarserFR = &coarserPtr->m_fluxRegister;

    coarserPhiOld = &coarserPtr->m_phiOld;
    coarserPhiNew = &coarserPtr->m_phiNew;

    tCoarserNew = coarserPtr->m_time;
    tCoarserOld = tCoarserNew - coarserPtr->m_dt;
  }

  // A finer level exists
  if (m_hasFiner)
  {
    // Recall that my flux register goes between my level and the next
    // finer level
    finerFR = &m_fluxRegister;
    finerFR->setToZero();
  }

  // Advance wave equation by one time step using 4th-order
  // Runge-Kutta.

  Interval inter = m_phiOld.interval();
  LevelData<FArrayBox> LOfPhiTmp;
  LevelData<FArrayBox> phiTmp;
  LevelData<FArrayBox> piTmp;
  LOfPhiTmp.define(m_phiOld);
  phiTmp.define(m_phiOld);
  piTmp.define(m_piOld);
  m_phiOld.copyTo(inter,phiTmp,inter);
  m_piOld.copyTo(inter,piTmp,inter);

  // RK Stage 1: compute RHS.
  // LofPhiTmp gets laplacian(phiOld);
  // flux registers are also computed
  m_levelWaveOperator.eval(m_phiOld,LOfPhiTmp,
                           *finerFR,*coarserFR,
                           *coarserPhiOld,tCoarserOld,
                           *coarserPhiNew,tCoarserNew,
                           m_time,m_dt/6);

  // RK Stage 1: compute update {phi,pi}_1, partial update to new vals.

  updateODE(LOfPhiTmp,m_piNew,m_dt/6);
  updateODE(piTmp,m_phiNew,m_dt/6);

  updateODE(piTmp,phiTmp,m_dt/2);
  updateODE(LOfPhiTmp,piTmp,m_dt/2);

  // RK Stage 2: compute RHS.

  m_levelWaveOperator.eval(phiTmp,LOfPhiTmp,
                           *finerFR, *coarserFR,
                           *coarserPhiOld,tCoarserOld,
                           *coarserPhiNew,tCoarserNew,
                           m_time+m_dt/2,m_dt/3);


  // RK Stage 2: compute update {phi,pi}_2, partial update to new vals.

  updateODE(piTmp,m_phiNew,m_dt/3);
  updateODE(LOfPhiTmp,m_piNew,m_dt/3);

  m_phiOld.copyTo(inter,phiTmp,inter);
  updateODE(piTmp,phiTmp,m_dt/2);
  m_piOld.copyTo(inter,piTmp,inter);
  updateODE(LOfPhiTmp,piTmp,m_dt/2);

  // RK Stage 3: compute RHS.

  m_levelWaveOperator.eval(phiTmp,LOfPhiTmp,
                           *finerFR, *coarserFR,
                           *coarserPhiOld,tCoarserOld,
                           *coarserPhiNew,tCoarserNew,
                           m_time+m_dt/2,m_dt/3);


  // RK Stage 3: compute update {phi,pi}_3, partial update to new vals.

  updateODE(piTmp,m_phiNew,m_dt/3);
  updateODE(LOfPhiTmp,m_piNew,m_dt/3);

  m_phiOld.copyTo(inter,phiTmp,inter);
  updateODE(piTmp,phiTmp,m_dt);
  m_piOld.copyTo(inter,piTmp,inter);
  updateODE(LOfPhiTmp,piTmp,m_dt);

  // RK Stage 4: compute RHS.

  m_levelWaveOperator.eval(phiTmp,LOfPhiTmp,
                           *finerFR, *coarserFR,
                           *coarserPhiOld,tCoarserOld,
                           *coarserPhiNew,tCoarserNew,
                           m_time+m_dt,m_dt/6);

  // RK Stage 4: compute final update of solution.

  updateODE(piTmp,m_phiNew,m_dt/6);
  updateODE(LOfPhiTmp,m_piNew,m_dt/6);

  // End application-dependent code - PC.

  // Update the time and store the new timestep
  m_time += m_dt;
  Real returnDt = m_dt;

  m_dtNew = returnDt;

  return returnDt;
}

////////////////////////////////////////////////////////////////////////////////

// Level Update operator for ODE solver. [NOTE: application-dependent code. PC]
// solution[i] = solution[i] + rhs[i]*dt
void AMRLevelWaveEquation::updateODE(LevelData<FArrayBox>& a_rhs,
                                     LevelData<FArrayBox>& a_soln,
                                     Real a_dt)
{
  DataIterator dit = a_soln.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    FArrayBox& curSoln = a_soln[dit()];
    FArrayBox curRhs(a_rhs[dit()].box(),m_numStates);
    curRhs.copy(a_rhs[dit()]);
    curRhs *= a_dt;
    curSoln += curRhs;
  }
}

////////////////////////////////////////////////////////////////////////////////

// Data initializer for ODE solver. [NOTE: application-dependent code. PC]
// new[i] = old[i]
void AMRLevelWaveEquation::initializeODE(LevelData<FArrayBox>& a_oldSoln,
                                         LevelData<FArrayBox>& a_newSoln)
{
  DataIterator dit = a_oldSoln.dataIterator();
  for (dit.begin();dit.ok();++dit)
  {
    FArrayBox& curOldSoln = a_oldSoln[dit()];
    FArrayBox& curNewSoln = a_newSoln[dit()];
    curNewSoln.copy(curOldSoln);
  }
}

////////////////////////////////////////////////////////////////////////////////

// Things to do after a timestep
void AMRLevelWaveEquation::postTimeStep()
{
  // Used for conservation tests
  static Real orig_integral = 0.0;
  static Real last_integral = 0.0;
  static bool first = true;

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::postTimeStep " << m_level << endl;
  }

  // Begin application-dependent code - PC.

  if (m_hasFiner)
  {
    // Reflux
    Real scale = -1.0/m_dx;
    m_fluxRegister.reflux(m_piNew,scale);

    // Average from finer level data
    AMRLevelWaveEquation* amrWaveFinerPtr = getFinerLevel();

    amrWaveFinerPtr->m_coarseAverage.averageToCoarse(m_piNew,
                                                    amrWaveFinerPtr->m_piNew);


/*
    amrWaveFinerPtr->m_coarseAverage.averageToCoarse(m_phiNew,
                                                    amrWaveFinerPtr->m_phiNew);
*/
    amrWaveFinerPtr->m_levelWaveOperator.avgdown(
                                                    amrWaveFinerPtr->m_phiNew,
                                                    m_phiNew);
  }
  // End application-dependent code - PC.

  if (s_verbosity >= 2 && m_level == 0)
  {
    int nRefFine = 1;

    pout() << "AMRLevelWaveEquation::postTimeStep:" << endl;
    pout() << "  Sums:" << endl;
    for (int comp = 0; comp < m_numStates; comp++)
    {
      Interval curComp(comp,comp);
      // Begin application-dependent code - PC.
      Real integral = computeSum(m_piNew,NULL,nRefFine,m_dx,curComp);
      // End application-dependent code - PC.

      pout() << "    " << setw(23)
                       << setprecision(16)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << integral
             << " --- " << m_stateNames[comp];

      if (comp == 0 )
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
        if (first)
        {
          orig_integral = integral;
          first = false;
        }
        last_integral = integral;
      }
      pout() << endl;
    }
  }
  if (s_verbosity >= 1 && m_level == 0)
  {
    // Begin application-dependent code - <dbs>
    Real energy = computeEnergy();
    pout() << " integrated energy = " << energy << endl ;
    // End application-dependent code - <dbs>
  }
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::postTimeStep " << m_level << " finished" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////

// Create tags for regridding
void AMRLevelWaveEquation::tagCells(IntVectSet& a_tags)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::tagCells " << m_level << endl;
  }

  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::tagCellsInit " << m_level << endl;
  }

  // Begin application-dependent code - PC.

  // Create tags based on undivided second difference of phi.

  LevelData<FArrayBox> lOfPhi(m_phiNew.getBoxes(),m_numStates);

  // Set up arguments to LevelWaveOperator::eval based on whether there are
  // coarser and finer levels

  // Undefined flux register in case we need it
  LevelFluxRegister dummyFR;

  // Undefined leveldata in case we need it
  LevelData<FArrayBox> dummyData;

  // Set arguments to dummy values and then fix if real values are available
  LevelFluxRegister* coarserFR = &dummyFR;
  LevelFluxRegister* finerFR   = &dummyFR;

  Real tCoarserOld = 0.0;
  Real tCoarserNew = 0.0;

  LevelData<FArrayBox>* coarserPhiOld = &dummyData;
  LevelData<FArrayBox>* coarserPhiNew = &dummyData;
  // A coarser level exists
  if (m_hasCoarser)
  {
    AMRLevelWaveEquation* coarserPtr = getCoarserLevel();

    coarserPhiOld = &coarserPtr->m_phiOld;
    coarserPhiNew = &coarserPtr->m_phiNew;

    tCoarserNew = coarserPtr->m_time;
    tCoarserOld = tCoarserNew - coarserPtr->m_dt;
  }
  IntVectSet localTags;
  const DisjointBoxLayout& levelDomain = m_phiNew.disjointBoxLayout();
  LevelWaveOperator lwo;
  if (m_hasCoarser)
    {
        lwo.define(m_grids,
                   getCoarserLevel()->m_grids,
                   m_problem_domain,
                   getCoarserLevel()->m_ref_ratio,
                   m_numStates,
                   m_dx,
                   m_hasCoarser,
                   m_hasFiner);
    }
  else
    {
        lwo.define(m_grids,
                   DisjointBoxLayout(),
                   m_problem_domain,
                   m_ref_ratio,
                   m_numStates,
                   m_dx,
                   m_hasCoarser,
                   m_hasFiner);
    }
  lwo.eval(m_phiNew,lOfPhi,
           *finerFR,*coarserFR,
           *coarserPhiOld,tCoarserOld,
           *coarserPhiNew,tCoarserNew,
           m_time,m_dt);
  DataIterator dit = levelDomain.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& b = levelDomain[dit()];
    FArrayBox& lPhiLocal = lOfPhi[dit()];

    FArrayBox tagVals(b,1);
    FORT_MAGNITUDE(CHF_FRA1(tagVals,0),
                   CHF_CONST_FRA(lPhiLocal),
                   CHF_BOX(b));


  // Tag where gradient exceeds threshold

    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();

      if (tagVals(iv) >= m_refineThresh/m_dx/m_dx)
      {
        localTags |= iv;
      }
    }
  }

  // End application-dependent code - PC.

  localTags.grow(m_tagBufferSize);

  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = localTags.minBox();
  localTagsBox &= m_problem_domain;
  localTags &= localTagsBox;

  a_tags = localTags;
}

////////////////////////////////////////////////////////////////////////////////

// Create tags at initialization
void AMRLevelWaveEquation::tagCellsInit(IntVectSet& a_tags)
{
  tagCells(a_tags);
}

////////////////////////////////////////////////////////////////////////////////

// Set up data on this level after regridding
void AMRLevelWaveEquation::regrid(const Vector<Box>& a_newGrids)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::regrid " << m_level << endl;
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
  // Begin application-dependent code - PC.

  LevelData<FArrayBox> phiOld;
  LevelData<FArrayBox> piOld;
  phiOld.define(m_phiNew);
  piOld.define(m_piNew);
  m_phiNew.copyTo(m_phiNew.interval(),
                phiOld,
                phiOld.interval());

  m_piNew.copyTo(m_piNew.interval(),
                piOld,
                piOld.interval());

  // Reshape state with new grids
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_phiNew.define(m_grids,m_numStates,ivGhost);
  m_phiOld.define(m_grids,m_numStates,ivGhost);
  m_piNew.define(m_grids,m_numStates,ivGhost);
  m_piOld.define(m_grids,m_numStates,ivGhost);

  // Set up data structures
  levelSetup();

  // Interpolate from coarser level
  if (m_hasCoarser)
  {
    AMRLevelWaveEquation* amrWaveCoarserPtr = getCoarserLevel();
    m_fineInterp.interpToFine(m_phiNew,amrWaveCoarserPtr->m_phiNew);
    m_fineInterp.interpToFine(m_piNew,amrWaveCoarserPtr->m_piNew);

  // Begin application-dependent code - PC.

  }

  // Copy from old state
  phiOld.copyTo(phiOld.interval(),
                  m_phiNew,
                  m_phiNew.interval());
  piOld.copyTo(piOld.interval(),
                  m_piNew,
                  m_piNew.interval());
}

////////////////////////////////////////////////////////////////////////////////

// Initialize grids
void AMRLevelWaveEquation::initialGrid(const Vector<Box>& a_newGrids)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::initialGrid " << m_level << endl;
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
  // Begin application-dependent code - PC.

  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_phiNew.define(m_grids,m_numStates,ivGhost);
  m_piNew.define(m_grids,m_numStates,ivGhost);
  m_phiOld.define(m_grids,m_numStates,ivGhost);
  m_piOld.define(m_grids,m_numStates,ivGhost);
  // End application-dependent code - PC.

  // Set up data structures
  levelSetup();
}

////////////////////////////////////////////////////////////////////////////////

// Initialize data
void AMRLevelWaveEquation::initialData()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::initialData " << m_level << endl;
  }
// Begin application-dependent code - PC.

  m_physIBCPtr->initialize(m_phiNew,m_piNew,m_dx);

// End application-dependent code - PC.
}

////////////////////////////////////////////////////////////////////////////////

// Things to do after initialization
void AMRLevelWaveEquation::postInitialize()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::postInitialize " << m_level << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////

#ifdef CH_USE_HDF5

// Write checkpoint header
void AMRLevelWaveEquation::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::writeCheckpointHeader" << endl;
  }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numStates*2;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates*2; ++comp)
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
void AMRLevelWaveEquation::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::writeCheckpointLevel" << endl;
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
  LevelData<FArrayBox> outData(m_phiNew.getBoxes(),2*m_numStates);
  Interval interval0(0,m_numStates-1);
  Interval interval1(m_numStates,2*m_numStates-1);
  m_phiNew.copyTo(interval0,outData,interval0);
  m_piNew.copyTo(interval0,outData,interval1);
  write(a_handle,outData.boxLayout());
  write(a_handle,outData,"data");
}

// Read checkpoint header
void AMRLevelWaveEquation::readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::readCheckpointHeader" << endl;
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
    MayDay::Error("AMRLevelWaveEquation::readCheckpointHeader: checkpoint file does not have num_components");
  }

  int numStates = header.m_int["num_components"];
  if (numStates != m_numStates*2)
  {
    MayDay::Error("AMRLevelWaveEquation::readCheckpointHeader: num_components in checkpoint file does not match solver");
  }

  // Get the component names
  std::string stateName;
  char compStr[60];
  for (int comp = 0; comp < m_numStates*2; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    if (header.m_string.find(compStr) == header.m_string.end())
    {
      MayDay::Error("AMRLevelWaveEquation::readCheckpointHeader: checkpoint file does not have enough component names");
    }

    stateName = header.m_string[compStr];
    if (stateName != m_stateNames[comp])
    {
      MayDay::Error("AMRLevelWaveEquation::readCheckpointHeader: state_name in checkpoint does not match solver");
    }
  }
}

// Read checkpoint data for this level
void AMRLevelWaveEquation::readCheckpointLevel(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::readCheckpointLevel" << endl;
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
    MayDay::Error("AMRLevelWaveEquation::readCheckpointLevel: file does not contain ref_ratio");
  }
  m_ref_ratio = header.m_int["ref_ratio"];

  if (s_verbosity >= 2)
  {
    pout() << "read ref_ratio = " << m_ref_ratio << endl;
  }

  // Get the tag buffer size
  if (header.m_int.find("tag_buffer_size") == header.m_int.end())
  {
    MayDay::Error("AMRLevelWaveEquation::readCheckpointLevel: file does not contain tag_buffer_size");
  }
  m_tagBufferSize=  header.m_int["tag_buffer_size"];

  if (s_verbosity >= 2)
  {
    pout() << "read tag_buffer_size = " << m_tagBufferSize << endl;
  }

  // Get dx
  if (header.m_real.find("dx") == header.m_real.end())
  {
    MayDay::Error("AMRLevelWaveEquation::readCheckpointLevel: file does not contain dx");
  }
  m_dx = header.m_real["dx"];

  if (s_verbosity >= 2)
  {
    pout() << "read dx = " << m_dx << endl;
  }

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
  {
    MayDay::Error("AMRLevelWaveEquation::readCheckpointLevel: file does not contain dt");
  }
  m_dt = header.m_real["dt"];

  if (s_verbosity >= 2)
  {
    pout() << "read dt = " << m_dt << endl;
  }

  // Get time
  if (header.m_real.find("time") == header.m_real.end())
  {
    MayDay::Error("AMRLevelWaveEquation::readCheckpointLevel: file does not contain time");
  }
  m_time = header.m_real["time"];

  if (s_verbosity >= 2)
  {
    pout() << "read time = " << m_time << endl;
  }

  // Get the problem domain
  if (header.m_box.find("prob_domain") == header.m_box.end())
  {
    MayDay::Error("AMRLevelWaveEquation::readCheckpointLevel: file does not contain prob_domain");
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
    MayDay::Error("AMRLevelWaveEquation::readCheckpointLevel: file does not contain a Vector<Box>");
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
    LayoutIterator lit = constGrids.layoutIterator();
    for (lit.begin(); lit.ok(); ++lit)
    {
      const Box& b = constGrids[lit()];
      pout() << lit().intCode() << ": " << b << endl;
    }
    pout() << endl;
  }

  // Reshape state with new grids
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  LevelData<FArrayBox> inData ;
  inData.define(m_grids,m_numStates*2,ivGhost);
  int dataStatus = read<FArrayBox>(a_handle, inData, "data", m_grids);
  if (dataStatus != 0)
  {
    MayDay::Error("AMRLevelWaveEquation::readCheckpointLevel: file does not contain all state data");
  }
  m_phiNew.define(m_grids,m_numStates,ivGhost);
  m_phiOld.define(m_grids,m_numStates,ivGhost);
  // repeat for pi variable
  m_piNew.define(m_grids,m_numStates,ivGhost);
  m_piOld.define(m_grids,m_numStates,ivGhost);

  Interval interval0(0,m_numStates-1);
  Interval interval1(m_numStates,2*m_numStates-1);
  inData.copyTo(interval0,m_phiNew,interval0);
  inData.copyTo(interval1,m_piNew,interval0);

  // Set up data structures
  levelSetup();
}

// Write plotfile header
void AMRLevelWaveEquation::writePlotHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::writePlotHeader" << endl;
  }

  // Setup the number of components
  HDF5HeaderData header;
  char compStr[30];
// Begin application-dependent code - PC.
  header.m_int["num_components"] = 3*m_numStates;

  // Setup the component names
    sprintf(compStr,"component_%d",0);
    header.m_string[compStr] = "phi_0";
    sprintf(compStr,"component_%d",1);
    header.m_string[compStr] = "phi_1";

    sprintf(compStr,"component_%d",2);
    header.m_string[compStr] = "pi_0";
    sprintf(compStr,"component_%d",3);
    header.m_string[compStr] = "pi_1";

    sprintf(compStr,"component_%d",4);
    header.m_string[compStr] = "LOfPhi_0";
    sprintf(compStr,"component_%d",5);
    header.m_string[compStr] = "LOfPhi_1";
// End application-dependent code - PC.
  // Write the header
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }
}

// Write plotfile data for this level
void AMRLevelWaveEquation::writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::writePlotLevel" << endl;
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

  // Begin application-dependent code - PC.

  // Write the data for this level
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  LevelData<FArrayBox> lOfPhi(m_phiNew.getBoxes(),m_numStates);
  LevelData<FArrayBox> phiDummy(m_phiNew.getBoxes(),m_numStates,ivGhost);
  Interval interval0(0,m_numStates-1);
  Interval interval1(m_numStates,2*m_numStates-1);
  Interval interval2(2*m_numStates,3*m_numStates-1);

  // Set up arguments to LevelGodunov::step based on whether there are
  // coarser and finer levels
  // Make copy of phi.
  m_phiNew.copyTo(interval0,phiDummy,interval0);

  // Undefined flux register in case we need it
  LevelFluxRegister dummyFR;

  // Undefined leveldata in case we need it
  LevelData<FArrayBox> dummyData;

  // Set arguments to dummy values and then fix if real values are available
  LevelFluxRegister* coarserFR = &dummyFR;
  LevelFluxRegister* finerFR   = &dummyFR;

  Real tCoarserOld = 0.0;
  Real tCoarserNew = 0.0;

  LevelData<FArrayBox>* coarserPhiOld = &dummyData;
  LevelData<FArrayBox>* coarserPhiNew = &dummyData;
  // A coarser level exists
  if (m_hasCoarser)
  {
    AMRLevelWaveEquation* coarserPtr = getCoarserLevel();

    coarserPhiOld = &coarserPtr->m_phiOld;
    coarserPhiNew = &coarserPtr->m_phiNew;

    tCoarserNew = coarserPtr->m_time;
    tCoarserOld = tCoarserNew - coarserPtr->m_dt;
  }
  LevelWaveOperator lwo;
  if (m_hasCoarser)
    {
        lwo.define(m_grids,
                   getCoarserLevel()->m_grids,
                   m_problem_domain,
                   getCoarserLevel()->m_ref_ratio,
                   m_numStates,
                   m_dx,
                   m_hasCoarser,
                   m_hasFiner);
    }
  else
    {
        lwo.define(m_grids,
                   DisjointBoxLayout(),
                   m_problem_domain,
                   m_ref_ratio,
                   m_numStates,
                   m_dx,
                   m_hasCoarser,
                   m_hasFiner);
    }
  lwo.eval(phiDummy,lOfPhi,
           *finerFR,*coarserFR,
           *coarserPhiOld,tCoarserOld,
           *coarserPhiNew,tCoarserNew,
           m_time,m_dt);

  LevelData<FArrayBox> outData(m_phiNew.getBoxes(),3*m_numStates);
  lOfPhi.copyTo(interval0,outData,interval2);
  m_phiNew.copyTo(interval0,outData,interval0);
  m_piNew.copyTo(interval0,outData,interval1);
  write(a_handle,m_phiNew.boxLayout());
  write(a_handle,outData,"data");

  // End application-dependent code - PC.

}

#endif

////////////////////////////////////////////////////////////////////////////////

// Returns the dt computed earlier for this level
Real AMRLevelWaveEquation::computeDt()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::computeDt " << m_level << endl;
  }

  Real newDt;
  newDt = m_dtNew;

  return newDt;
}

////////////////////////////////////////////////////////////////////////////////

// Compute dt using initial data
Real AMRLevelWaveEquation::computeInitialDt()
{
  Real newDT = m_initial_dt_multiplier * m_dx;

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::computeInitialDt on level " << m_level << " = " << newDT << endl;
  }

  return newDT;
}

////////////////////////////////////////////////////////////////////////////////

// Compute integral of energy over this and all finer levels
// by chasing pointers up the level hierarchy
Real AMRLevelWaveEquation::computeEnergy()
{
  Real total_e = 0 ;
  AMRLevelWaveEquation * level_ptr = this ;
  while( level_ptr != NULL )
  {
    DisjointBoxLayout * fineDBL = NULL ;
    if( level_ptr->m_hasFiner ) fineDBL = &(level_ptr->getFinerLevel()->m_grids);
    Real lambda = 1 ;
    Real dx = level_ptr->m_dx ;
    int  ref_ratio = level_ptr->m_ref_ratio ;
    Real energy = ::computeEnergy(level_ptr->m_phiNew,level_ptr->m_piNew, fineDBL,
                                  ref_ratio, dx, lambda);
    total_e += energy ;
    // repeat with next finer level; use NULL to stop
    level_ptr = ( level_ptr->m_hasFiner
                  ? dynamic_cast<AMRLevelWaveEquation*>(level_ptr->m_finer_level_ptr)
                  : NULL ) ;
  }
  return total_e ;
}

////////////////////////////////////////////////////////////////////////////////

// Set the CFL number
void AMRLevelWaveEquation::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
}

////////////////////////////////////////////////////////////////////////////////

// Set the physical dimension of the longest side of the domain
void AMRLevelWaveEquation::domainLength(Real a_domainLength)
{
  m_domainLength = a_domainLength;
}

////////////////////////////////////////////////////////////////////////////////

// Set the refinement threshold
void AMRLevelWaveEquation::refinementThreshold(Real a_refineThresh)
{
  m_refineThresh = a_refineThresh;
}

////////////////////////////////////////////////////////////////////////////////

// Set the tag buffer size
void AMRLevelWaveEquation::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
}

////////////////////////////////////////////////////////////////////////////////

// Create a load-balanced DisjointBoxLayout from a collection of Boxes
DisjointBoxLayout AMRLevelWaveEquation::loadBalance(const Vector<Box>& a_grids)
{
  // Load balance and create boxlayout
  Vector<int> procMap;

  // appears to be faster for all procs to do the loadbalance (ndk)
  LoadBalance(procMap,a_grids);

  if (s_verbosity >= 4)
  {
    pout() << "AMRLevelWaveEquation::loadBalance: procesor map: " << endl;
    for (int igrid = 0; igrid < a_grids.size(); ++igrid)
    {
      pout() << igrid << ": " << procMap[igrid] << "  " << a_grids[igrid].volume() << endl;
    }
    pout() << endl;
  }

  DisjointBoxLayout dbl(a_grids,procMap,m_problem_domain);
  dbl.close();

  return dbl;
}

////////////////////////////////////////////////////////////////////////////////

// Setup menagerie of data structures
void AMRLevelWaveEquation::levelSetup()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelWaveEquation::levelSetup " << m_level << endl;
  }

  AMRLevelWaveEquation* amrWaveCoarserPtr = getCoarserLevel();
  AMRLevelWaveEquation* amrWaveFinerPtr   = getFinerLevel();

  m_hasCoarser = (amrWaveCoarserPtr != NULL);
  m_hasFiner   = (amrWaveFinerPtr   != NULL);

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

    const DisjointBoxLayout& coarserLevelDomain = amrWaveCoarserPtr->m_grids;

    // Maintain levelWaveOperator
    m_levelWaveOperator.define(m_grids,
                          coarserLevelDomain,
                          m_problem_domain,
                          nRefCrse,
                          m_numStates,
                          m_dx,
                          m_hasCoarser,
                          m_hasFiner);

    // This may look twisted but you have to do this this way because the
    // coarser levels get setup before the finer levels so, since a flux
    // register lives between this level and the next FINER level, the finer
    // level has to do the setup because it is the only one with the
    // information at the time of construction.

    // Maintain flux registers
    amrWaveCoarserPtr->m_fluxRegister.define(m_grids,
                                            amrWaveCoarserPtr->m_grids,
                                            m_problem_domain,
                                            amrWaveCoarserPtr->m_ref_ratio,
                                            m_numStates);
    amrWaveCoarserPtr->m_fluxRegister.setToZero();
  }
  else
  {
    m_levelWaveOperator.define(m_grids,
                          DisjointBoxLayout(),
                          m_problem_domain,
                          m_ref_ratio,
                          m_numStates,
                          m_dx,
                          m_hasCoarser,
                          m_hasFiner);
  }
}

////////////////////////////////////////////////////////////////////////////////

// Get the next coarser level
AMRLevelWaveEquation* AMRLevelWaveEquation::getCoarserLevel() const
{
  AMRLevelWaveEquation* amrWaveCoarserPtr = NULL;

  if (m_coarser_level_ptr != NULL)
  {
    amrWaveCoarserPtr = dynamic_cast<AMRLevelWaveEquation*>(m_coarser_level_ptr);

    if (amrWaveCoarserPtr == NULL)
    {
      MayDay::Error("AMRLevelWaveEquation::getCoarserLevel: dynamic cast failed");
    }
  }

  return amrWaveCoarserPtr;
}

////////////////////////////////////////////////////////////////////////////////

// Get the next finer level
AMRLevelWaveEquation* AMRLevelWaveEquation::getFinerLevel() const
{
  AMRLevelWaveEquation* amrWaveFinerPtr = NULL;

  if (m_finer_level_ptr != NULL)
  {
    amrWaveFinerPtr = dynamic_cast<AMRLevelWaveEquation*>(m_finer_level_ptr);

    if (amrWaveFinerPtr == NULL)
    {
      MayDay::Error("AMRLevelWaveEquation::getFinerLevel: dynamic cast failed");
    }
  }

  return amrWaveFinerPtr;
}
