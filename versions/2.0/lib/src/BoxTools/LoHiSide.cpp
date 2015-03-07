#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "SPACE.H"

#include "LoHiSide.H"
#include "NamespaceHeader.H"

Side::LoHiSide flip(Side::LoHiSide a_side)
{
    CH_assert((a_side == Side::Lo) || (a_side == Side::Hi));

    return (a_side == Side::Lo) ? Side::Hi : Side::Lo;
}

Side::LoHiSide Side::flip(Side::LoHiSide a_side)
{
    CH_assert((a_side == Side::Lo) || (a_side == Side::Hi));

    return (a_side == Side::Lo) ? Side::Hi : Side::Lo;
}

int sign(Side::LoHiSide a_side)
{
    CH_assert((a_side == Side::Lo) || (a_side == Side::Hi));

    return (a_side == Side::Lo) ? -1 : 1;
}

SideIterator::SideIterator()
  :
  m_current(-1)
{
  reset();
}

void SideIterator::begin()
{
    m_current = 0;
}

void SideIterator::next()
{
    ++m_current;
}

void SideIterator::operator ++ ()
{
    ++m_current;
}

Side::LoHiSide SideIterator::operator () () const
{
    switch (m_current)
    {
      case 0:
          return Side::Lo;
          //break;

      case 1:
          return Side::Hi;
          //break;

      default:
          return Side::Invalid;
          //break;
    }
}

bool SideIterator::ok() const
{
    return ((m_current > -1) && (m_current < Side::NUMSIDES));
}
#include "NamespaceFooter.H"
