#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "DataIterator.H"
#include "NamespaceHeader.H"

#ifdef CH_MPI
DataIterator::DataIterator(const BoxLayout& plan,
                           const int* layoutID)
  :m_layout(plan), m_index(0), m_current(0, layoutID), m_procID(0)
{
  // in many programs a user will call MPI_Finalize, at the end of
  // the main routine.  But some objects, like LayoutData, have not
  // had their destructors called.
  m_procID = procID();
  begin(); // might result in ok()==false if there were no boxes on this processor
}

void DataIterator::end()
{
  m_index = m_layout.size();
}

#endif
#include "NamespaceFooter.H"
