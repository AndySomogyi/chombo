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

#include "Box.H"
#include "Vector.H"
#include "LayoutIterator.H"
#include "parstream.H"

#include "AMRLevel.H"
#include "NamespaceHeader.H"

int AMRLevel::s_verbosity = 0;

bool AMRLevel::isDefined() const
{
  return m_isDefined;
}

AMRLevel::~AMRLevel()
{
}

AMRLevel::AMRLevel()
{
  m_coarser_level_ptr = NULL;
  m_finer_level_ptr = NULL;
  m_isDefined = false;
  m_level = 0;
  m_time = 0;
  m_dt = 0;
  m_initial_dt_multiplier = 0.1;
}

void AMRLevel::define(AMRLevel*  a_coarser_level_ptr,
                      const Box& a_problem_domain,
                      int        a_level,
                      int        a_ref_ratio)
{
  ProblemDomain physDomain(a_problem_domain);

  define(a_coarser_level_ptr, physDomain, a_level, a_ref_ratio);
}

void AMRLevel::define(AMRLevel*            a_coarser_level_ptr,
                      const ProblemDomain& a_problem_domain,
                      int                  a_level,
                      int                  a_ref_ratio)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevel::define" << endl;
  }

  m_coarser_level_ptr = a_coarser_level_ptr;
  m_problem_domain = a_problem_domain;
  m_level = a_level;
  m_ref_ratio = a_ref_ratio;
  m_finer_level_ptr = NULL;
  m_isDefined = true;
}

void AMRLevel::finerLevelPtr(AMRLevel* a_finer_level_ptr)
{
  m_finer_level_ptr = a_finer_level_ptr;
}

void AMRLevel::dt(Real a_dt)
{
  m_dt = a_dt;
}

Real AMRLevel::dt() const
{
  return(m_dt);
}

const ProblemDomain& AMRLevel::problemDomain() const
{
  return(m_problem_domain);
}

Vector<Box> AMRLevel::boxes() const
{
  return(m_level_grids);
}

int AMRLevel::refRatio() const
{
  return(m_ref_ratio);
}

void AMRLevel::time(Real a_time)
{
  m_time = a_time;
}

Real AMRLevel::time() const
{
  return m_time;
}

void AMRLevel::initialDtMultiplier(Real a_initial_dt_multiplier)
{
  m_initial_dt_multiplier = a_initial_dt_multiplier;
}

Real AMRLevel::initialDtMultiplier() const
{
  return m_initial_dt_multiplier;
}

// static
void AMRLevel::verbosity(int a_verbosity)
{
  s_verbosity = a_verbosity;
}

// static
int AMRLevel::verbosity()
{
  return(s_verbosity);
}

void AMRLevel::preRegrid(int a_base_level)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevel::preRegrid" << endl;
  }
}

void AMRLevel::preRegrid(int a_base_level, const Vector<Vector<Box> >& a_new_grids)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevel::preRegrid" << endl;
  }
}

void AMRLevel::postRegrid(int a_base_level)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevel::postRegrid" << endl;
  }
}

Vector<AMRLevel*> AMRLevel::getAMRLevelHierarchy()
{
  Vector<AMRLevel*> retval;
  // First go to level 0
  AMRLevel* levelPtr = this;
  bool hasCoarser = ((levelPtr->m_coarser_level_ptr != NULL)
                  && (levelPtr->m_coarser_level_ptr->m_isDefined)
                  && (levelPtr->m_coarser_level_ptr->m_level_grids.size() > 0));
  while (hasCoarser)
    {
      levelPtr = levelPtr->m_coarser_level_ptr;
      hasCoarser = ((levelPtr->m_coarser_level_ptr != NULL)
                 && (levelPtr->m_coarser_level_ptr->m_isDefined)
                 && (levelPtr->m_coarser_level_ptr->m_level_grids.size() > 0));
    }

  // Now can accumulate the pointers by chasing finer level
  retval.push_back(levelPtr);

  bool hasFiner = ((levelPtr->m_finer_level_ptr != NULL)
                && (levelPtr->m_finer_level_ptr->m_isDefined)
                && (levelPtr->m_finer_level_ptr->m_level_grids.size() > 0));
  while (hasFiner)
    {
      levelPtr = levelPtr->m_finer_level_ptr;
      retval.push_back(levelPtr);

      hasFiner = ((levelPtr->m_finer_level_ptr != NULL)
               && (levelPtr->m_finer_level_ptr->m_isDefined)
               && (levelPtr->m_finer_level_ptr->m_level_grids.size() > 0));
    }

  return retval;
}

#include "NamespaceFooter.H"
