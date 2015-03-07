#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "LayoutIterator.H"
#include "DataIterator.H"
#include "NamespaceHeader.H"

void LayoutIterator::reset()
{
  m_index = 0;

  if (ok())
  {
    m_current.m_index = (*(m_layout.m_boxes))[m_index].index;
  }
}

void LayoutIterator::end()
{
  m_index = m_layout.size();
}

LayoutIterator::LayoutIterator(const BoxLayout& a_boxlayout,
                               const int* a_layoutID)
  :
  m_layout(a_boxlayout),
  m_index(0),
  m_current(0, a_layoutID)
{
  if (ok())
  {
    m_current.m_index = (*(m_layout.m_boxes))[m_index].index;
  }
}
#include "NamespaceFooter.H"
