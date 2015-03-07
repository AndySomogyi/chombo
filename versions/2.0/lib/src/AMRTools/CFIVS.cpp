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

#include "CFIVS.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

bool CFIVS::isDefined() const
{
  return m_isDefined;
}

void CFIVS::coarsen(int ref)
{
  CH_assert(m_isDefined);
  m_fiinterpIVS.coarsen(ref);
  m_packedBox.coarsen(ref);
}
const IntVectSet& CFIVS::getFineIVS() const
{
 CH_assert(m_isDefined);

  return m_fiinterpIVS;
}

void CFIVS::setDefaultValues()
{
  m_fiinterpIVS.define();
  m_packed = false;

  m_isDefined = false;
}

CFIVS::CFIVS()
{
  setDefaultValues();
}

CFIVS::~CFIVS()
{
  setDefaultValues();
}

CFIVS::CFIVS(const Box&               a_domain,
             const Box&               a_boxIn,
             const DisjointBoxLayout& a_fineBoxes,
             int                      a_direction,
             Side::LoHiSide           a_hiorlo)
{
  setDefaultValues();
  ProblemDomain probdomain(a_domain);

  define(probdomain, a_boxIn, a_fineBoxes, a_direction, a_hiorlo);

}

CFIVS::CFIVS(const ProblemDomain&     a_domain,
             const Box&               a_boxIn,
             const DisjointBoxLayout& a_fineBoxes,
             int                      a_direction,
             Side::LoHiSide           a_hiorlo)
{
  setDefaultValues();

  define(a_domain, a_boxIn, a_fineBoxes, a_direction, a_hiorlo);

}

void CFIVS::define(const Box&               a_domain,
                   const Box&               a_boxIn,
                   const DisjointBoxLayout& a_fineBoxes,
                   int                      a_direction,
                   Side::LoHiSide           a_hiorlo)
{
  ProblemDomain probdomain(a_domain);

  define(probdomain, a_boxIn, a_fineBoxes, a_direction, a_hiorlo);
}

void CFIVS::define(const ProblemDomain&     a_domain,
                   const Box&               a_boxIn,
                   const DisjointBoxLayout& a_fineBoxes,
                   int                      a_direction,
                   Side::LoHiSide           a_hiorlo)
{
  CH_TIME("CFIVS::define(slow)");
  m_isDefined = true;

  CH_assert(a_direction >= 0);
  CH_assert(a_direction < SpaceDim);
  CH_assert(!a_domain.isEmpty());
  CH_assert(a_domain.contains(a_boxIn));
  CH_assert(a_fineBoxes.checkPeriodic(a_domain));
  CH_assert((a_hiorlo == Side::Lo) || (a_hiorlo == Side::Hi));

  // create fine stencil
  Box finebox = a_boxIn;
  Box edgebox;

  if (a_hiorlo == Side::Lo)
    {
      edgebox = adjCellLo(finebox,a_direction,1);
    }
  else
    {
      edgebox = adjCellHi(finebox,a_direction,1);
    }

  edgebox &= a_domain;

  if (!edgebox.isEmpty())
    {
      m_fiinterpIVS.define(edgebox);

      LayoutIterator lit = a_fineBoxes.layoutIterator();
      Box periodicTestBox(a_domain.domainBox());

      if (a_domain.isPeriodic())
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (a_domain.isPeriodic(idir))
                periodicTestBox.grow(idir,-1);
            }
        }

      for (lit.reset(); lit.ok(); ++lit)
        {
          m_fiinterpIVS -= a_fineBoxes[lit()];

          // only do this IF we're periodic _and_ both boxes
          // adjoin the domain box boundary somewhere
          if (a_domain.isPeriodic() && !periodicTestBox.contains(edgebox)
              && !periodicTestBox.contains(a_fineBoxes[lit()]))
            {
              ShiftIterator shiftIt = a_domain.shiftIterator();
              IntVect shiftMult = a_domain.domainBox().size();
              Box shiftedBox = a_fineBoxes[lit()];

              for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                {
                  IntVect shiftVect(shiftMult*shiftIt());

                  shiftedBox.shift(shiftVect);
                  m_fiinterpIVS -= shiftedBox;
                  shiftedBox.shift(-shiftVect);
                }
            }
        }
    }

  m_fiinterpIVS.compact();

  if (m_fiinterpIVS.isEmpty())
    {
      m_packed = false;
    }
  else
  {
    Box min = m_fiinterpIVS.minBox();

    if (min.numPts() == m_fiinterpIVS.numPts())
    {
      m_packCount++;
      m_packed = true;
      m_packedBox = min;
    }
    else
    {
      m_sparseCount++;
      m_packed = false;
    }
  }
}


void CFIVS::define(const ProblemDomain&     a_domain,
                   const Box&               a_boxIn,
                   const Vector<Box>&       a_periodicfineBoxes,
                   int                      a_direction,
                   Side::LoHiSide           a_hiorlo)
{
  CH_TIME("CFIVS::define(fast)");
  m_isDefined = true;

  CH_assert(a_direction >= 0);
  CH_assert(a_direction < SpaceDim);
  CH_assert(!a_domain.isEmpty());
  CH_assert(a_domain.contains(a_boxIn));
  CH_assert((a_hiorlo == Side::Lo) || (a_hiorlo == Side::Hi));

  // create fine stencil
  Box finebox = a_boxIn;
  Box edgebox;

  if (a_hiorlo == Side::Lo)
    {
      edgebox = adjCellLo(finebox,a_direction,1);
    }
  else
    {
      edgebox = adjCellHi(finebox,a_direction,1);
    }

  edgebox &= a_domain;

  if(edgebox.isEmpty())
    {
      m_packed = false;
      return;
    }
 
  m_fiinterpIVS.define(edgebox);
  int w1 = edgebox.smallEnd()[0];
  int w2 = edgebox.bigEnd()[0];
  for(int i=0; i<a_periodicfineBoxes.size(); ++i)
    {
      const Box& b = a_periodicfineBoxes[i];
      if(b.bigEnd()[0] >= w1)
        {
          m_fiinterpIVS -= b;
          if(b.smallEnd()[0] > w2){
            i=a_periodicfineBoxes.size();
          }
        }
    }
   

  m_fiinterpIVS.compact();

  if (m_fiinterpIVS.isEmpty())
    {
      m_packed = false;
    }
  else
    {
      Box min = m_fiinterpIVS.minBox();
      
      if (min.numPts() == m_fiinterpIVS.numPts())
        {
          m_packCount++;
          m_packed = true;
          m_packedBox = min;
        }
      else
        {
          m_sparseCount++;
          m_packed = false;
        }
    }
}
long long CFIVS::m_packCount = 0;

long long CFIVS::m_sparseCount = 0;
#include "NamespaceFooter.H"
