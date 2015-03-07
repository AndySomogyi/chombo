#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// NodeNoOpSmoother.cpp
// adapted from NoOpSmoother by DFMartin, Sun, May 5, 2002
// petermc, 5 June 2002

#include "NodeNoOpSmoother.H"

// ---------------------------------------------------------
NodeNoOpSmoother::NodeNoOpSmoother()
{
}


// ---------------------------------------------------------
NodeNoOpSmoother::~NodeNoOpSmoother()
{
}


// ---------------------------------------------------------
NodeBaseBottomSmoother*
NodeNoOpSmoother::new_bottomSmoother() const
{
  NodeNoOpSmoother* newsmoother = new NodeNoOpSmoother();
  if (newsmoother == NULL)
    {
      MayDay::Error("Out of Memory in NodeNoOpSmoother::new_bottomSmoother");
    }
  return static_cast<NodeBaseBottomSmoother*> (newsmoother);
}



// ---------------------------------------------------------
/***********************/
// Nothing gets done here.
/***********************/
void
NodeNoOpSmoother::doBottomSmooth(LevelData<NodeFArrayBox>& a_phi,
                                     const LevelData<NodeFArrayBox>& a_rhs,
                                     NodeLevelOp* a_levelop_ptr)
{
  // do nothing
}
