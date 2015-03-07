#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _DATAINDEX_H_
#define _DATAINDEX_H_

#include <cstdlib>
#include "NamespaceHeader.H"

#ifndef WRAPPER
class BoxLayout;

class DataIterator;
class LayoutIterator;

#endif /* WRAPPER */


///An index for LayoutIterator
class LayoutIndex
{
public:
  ///
  LayoutIndex();

  ~LayoutIndex()
  {}

  // default copy and assign should be fine.

  ///
  bool operator == (const LayoutIndex& a_rhs) const;

  bool eq(const LayoutIndex& a_rhs) const
  {
    return *this == a_rhs;
  }

  ///
  bool operator != (const LayoutIndex& a_rhs) const
  {
    return !(*this == a_rhs);
  }

  /// returns 'true' if this DataIndex has been null constructed.
  bool isNull() const;

  // not user function. breaks design encapsulation. method not intended for
  // users. The alternative is to make this class a friend of BoxLayoutData<T>,
  // which implicitly makes this class templated, which makes no bloody sense.
  // The behaviour of this method is undefined for users.  bvs
  int intCode() const;


private:
  friend class DataIterator;
  friend class LayoutIterator;
  friend class BoxLayout;
  friend class DataIndex;

  int        m_index;
  const int* m_layoutIntPtr;

  LayoutIndex(int        a_index,
              const int* a_layoutID)
    :
    m_index(a_index),
    m_layoutIntPtr(a_layoutID)
  {}
};

class DataIndex : public LayoutIndex
{
public:
  explicit DataIndex(const LayoutIndex& a_promotion)
    :
    LayoutIndex(a_promotion)
  {}

  DataIndex()
    :
    LayoutIndex()
  {}

  ~DataIndex()
  {}


private:
  friend class LayoutIterator;
  friend class DataIterator;
  friend class BoxLayout;
  friend class Copier;

  DataIndex(int        a_index,
            const int* a_layoutID)
    :
    LayoutIndex(a_index, a_layoutID)
  {}
};

#ifndef WRAPPER
inline LayoutIndex::LayoutIndex()
  :
  m_index(0),
  m_layoutIntPtr(NULL)
{}

inline bool LayoutIndex::operator == (const LayoutIndex& a_rhs) const
{
  return (m_index == a_rhs.m_index);
}

inline int LayoutIndex::intCode() const
{
  return m_index;
}

inline bool LayoutIndex::isNull() const
{
  return m_layoutIntPtr == NULL;
}
#endif

#include "NamespaceFooter.H"
#endif
