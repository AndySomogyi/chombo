#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <cstring>

#include "BaseFab.H"
#include "NamespaceHeader.H"

Real BaseFabRealSetVal = BASEFAB_REAL_SETVAL;

template < > void BaseFab<Real>::define()
{
 CH_assert(m_nvar > 0);
 CH_assert(m_dptr == 0);
 CH_assert(m_numpts > 0);
 CH_assert(!m_aliased);
  //CH_assert(!(The_FAB_Arena == 0)); // not a sufficient test!!!

#ifdef CH_USE_MEMORY_TRACKING
  if (s_Arena == NULL)
  {
    s_Arena = new BArena(name().c_str());
  }
#else
  if (s_Arena == NULL)
  {
    s_Arena = new BArena("");
  }
#endif

  // if (s_Arena == NULL)
  // {
  //   MayDay::Error("malloc in basefab failed");
  // }

  m_truesize = m_nvar * m_numpts;
  m_dptr     = static_cast<Real*>(s_Arena->alloc(m_truesize * sizeof(Real)));

#ifdef CH_USE_MEMORY_TRACKING
  s_Arena->bytes += m_truesize * sizeof(Real) + sizeof(BaseFab<Real>);

  if (s_Arena->bytes > s_Arena->peak)
  {
    s_Arena->peak = s_Arena->bytes;
  }
#endif

#ifdef CH_USE_SETVAL
  setVal(BaseFabRealSetVal);
#endif
}

template < > void BaseFab<int>::define()
{
 CH_assert(m_nvar > 0);
 CH_assert(m_dptr == 0);
 CH_assert(m_numpts > 0);
 CH_assert(!m_aliased);
  //CH_assert(!(The_FAB_Arena == 0)); // not a sufficient test!!!

#ifdef CH_USE_MEMORY_TRACKING
  if (s_Arena == NULL)
  {
    s_Arena = new BArena(name().c_str());
  }
#else
  if (s_Arena == NULL)
  {
    s_Arena = new BArena("");
  }
#endif

  // if(s_Arena == NULL)
  // {
  //   MayDay::Error("malloc in basefab failed");
  // }

  m_truesize = m_nvar * m_numpts;
  m_dptr     = static_cast<int*>(s_Arena->alloc(m_truesize * sizeof(int)));

#ifdef CH_USE_MEMORY_TRACKING
  s_Arena->bytes += m_truesize * sizeof(int) + sizeof(BaseFab<int>);

  if (s_Arena->bytes > s_Arena->peak)
  {
    s_Arena->peak = s_Arena->bytes;
  }
#endif
}

template < > void BaseFab<Real>::undefine()
{
  if (m_aliased)
  {
    m_dptr = 0;
    return;
  }

  if (m_dptr == 0)
  {
    return;
  }

  s_Arena->free(m_dptr);

#ifdef CH_USE_MEMORY_TRACKING
  s_Arena->bytes -= m_truesize * sizeof(Real) + sizeof(BaseFab<Real>);
#endif

  m_dptr = 0;
}

template < > void BaseFab<int>::undefine()
{
  if (m_aliased)
  {
    m_dptr = 0;
    return;
  }

  if (m_dptr == 0)
  {
    return;
  }

  s_Arena->free(m_dptr);

#ifdef CH_USE_MEMORY_TRACKING
  s_Arena->bytes -= m_truesize * sizeof(int) + sizeof(BaseFab<int>);
#endif

  m_dptr = 0;
}

template < > void BaseFab<Real>::setVal(Real a_val)
{
  if (a_val == 0)
  {
    memset(m_dptr, 0, m_truesize*sizeof(Real));
  }
  else
  {
    Real* end = m_dptr + m_truesize;
    for (Real* v = m_dptr; v<end; v++)
    {
      *v = a_val;
    }
  }
}
#include "NamespaceFooter.H"
