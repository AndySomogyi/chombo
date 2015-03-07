#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "DenseIntVectSet.H"
#include "BoxIterator.H"
#include "MayDay.H"
#include "ProblemDomain.H"
#include "SPACE.H"
#include "SPMD.H"
#include "NamespaceHeader.H"

DenseIntVectSet::DenseIntVectSet(const Box& a_domain, bool init)
  :
  m_domain(a_domain),
  m_bits(a_domain.numPts(), init),
  m_minBox(a_domain)
{}

bool DenseIntVectSet::operator[](const IntVect& index) const
{
  if(!m_domain.contains(index)) return false;

  return m_bits[m_domain.index(index)];
}

DenseIntVectSet& DenseIntVectSet::operator-=(const Box& b)
{
  if(isEmpty()) return *this;
  if(!m_domain.intersects(b)) return *this;
  BoxIterator it(b & m_domain);
  for(it.begin(); it.ok(); ++it) this->operator-=(it());
  return *this;
}

DenseIntVectSet& DenseIntVectSet::operator|=(const DenseIntVectSet& d)
{
  if(m_domain.contains(d.m_domain))
        {
          BoxIterator bit(d.m_domain);
          int i=0;
          for(bit.begin(); bit.ok(); ++bit, ++i){
                if(d.m_bits[i]) m_bits.setTrue(m_domain.index(bit()));
          }
        }
  else if(d.m_domain.contains(m_domain))
        {
          DenseIntVectSet newSet = d;
          BoxIterator bit(m_domain);
          int i=0;
          for(bit.begin(); bit.ok(); ++bit, ++i){
                if(m_bits[i]) newSet.m_bits.setTrue(newSet.m_domain.index(bit()));
          }
          *this = newSet;
        }
  else
        {
          Box newDomain = minBox(m_domain, d.m_domain);
          DenseIntVectSet newSet(newDomain, false);
          newSet |= *this;
          newSet |= d;
          *this = newSet;
        }
  return *this;
}

DenseIntVectSet& DenseIntVectSet::operator|=(const Box& b)
{
  //if(!m_domain.contains(b)){
  //      MayDay::Error("Box union with DenseIntVectSet outside m_domain");
  //}
  CH_assert(m_domain.contains(b));

  BoxIterator bit(b);
  for(bit.begin(); bit.ok(); ++bit)
        {
          m_bits.setTrue(m_domain.index(bit()));
        }
  return *this;
}

DenseIntVectSet& DenseIntVectSet::operator|=(const IntVect& intvect)
{
  //if(!m_domain.contains(intvect)){
  //      MayDay::Error("union with DenseIntVectSet outside m_domain");
  //}
  CH_assert(m_domain.contains(intvect));

  m_bits.setTrue(m_domain.index(intvect));
  return *this;
}

DenseIntVectSet& DenseIntVectSet::operator&=(const ProblemDomain& domain)
{
  if(domain.domainBox().contains(m_domain)) return *this;
  DenseIntVectSet tmp(*this);
  tmp -= domain.domainBox();
  if(m_domain.intersects(domain.domainBox()))
    *this &= domain.domainBox();
  else
    for(int i=0; i<m_bits.size(); i++)
      m_bits.setFalse(i);

  DenseIntVectSetIterator it(tmp);
  for(it.begin(); it.ok(); ++it)
    {
      IntVect iv = it();
      if(domain.image(iv)) *this |= it();
    }
  return *this;
}

DenseIntVectSet& DenseIntVectSet::operator&=(const Box& b)
{
  if(b.contains(m_domain) || b == m_domain) return *this;
  if(isEmpty()) return *this;
  if(!m_domain.intersects(b))
    {
      *this = DenseIntVectSet();
      return *this;
    }
  if(m_domain.numPts() < b.numPts())
    {
      BoxIterator bit(m_domain);
      int i=0;
      for(bit.begin(); bit.ok();++i, ++bit)
        {
          if(b.contains(bit())) { }// do nothing;
          else m_bits.setFalse(i);
        }
    }
  else
    {
      Box btmp(b);
      btmp&=m_domain;
      DenseIntVectSet tmp = DenseIntVectSet(btmp, true);
      BoxIterator bit(btmp); int i=0;
      for(bit.begin(); bit.ok(); ++i, ++bit)
        {
          if(!this->operator[](bit())) tmp.m_bits.setFalse(i);
        }
      *this = tmp;
    }
  return *this;

}

Vector<Box> DenseIntVectSet::createBoxes() const
{
  Vector<Box> boxes;
  DenseIntVectSetIterator it(*this);
  for(it.begin(); it.ok(); ++it){
        boxes.push_back(Box(it(), it()));
  }
  return boxes;
}

void DenseIntVectSet::recalcMinBox() const
{
  if(isEmpty()){
    (Box&)m_minBox = Box();
    return;
  }
  DenseIntVectSetIterator it(*this);
  it.begin();
  (Box&)m_minBox = Box(it(), it());
  ++it;
  for(; it.ok(); ++it){
    Box m(it(), it());
    ((Box&)m_minBox).minBox(m);
  }
}

DenseIntVectSet DenseIntVectSet::chop(int dir, int chop_pnt)
{
  if(m_minBox.smallEnd(dir) >= chop_pnt)
    {
      DenseIntVectSet rtn(*this);
      *this = DenseIntVectSet();
      return rtn;
    }
  if(m_minBox.bigEnd(dir) < chop_pnt)
    {
      return DenseIntVectSet();
    }
  Box chop = m_domain;
  Box chopThis = chop.chop(dir, chop_pnt);
  DenseIntVectSet left(chop);
  DenseIntVectSet right(chopThis);
  left.intersect(*this);
  right.intersect(*this);
  *this = right;
  return left;

}

bool DenseIntVectSet::contains(const Box& box) const
{
  if(!m_minBox.contains(box)) return false;
  for(BoxIterator bit(box); bit.ok(); ++bit)
    if(!this->operator[](bit())) return false;
  return true;
}

DenseIntVectSet& DenseIntVectSet::operator&=(const DenseIntVectSet& ivs)
{
  if(isEmpty()) return *this;
  if(ivs.isEmpty())
    {
      *this = DenseIntVectSet();
      return *this;
    }
  if(!m_domain.intersectsNotEmpty(ivs.m_domain))
    {
      *this = DenseIntVectSet();
      return *this;
    }
  return intersect(ivs);
}

DenseIntVectSet& DenseIntVectSet::intersect(const DenseIntVectSet& ivs)
{
  BoxIterator bit(m_domain);
  int i=0;
  for(bit.begin(); bit.ok();++i, ++bit)
    {
      if(m_bits[i])
        {
          if(ivs[bit()]) { } // do nothing
          else
            m_bits.setFalse(i);
        }
    }
  return *this;
}

void DenseIntVectSet::coarsen(int iref)
{
  if(iref == 1) return;
  CH_assert(iref >= 1);
  // int refinements = iref/2;
  CH_assert((iref/2)*2 == iref); // check iref for power of 2

  Box newDomain(m_domain);
  newDomain.coarsen(iref);
  DenseIntVectSet newSet(newDomain, false);
  BoxIterator bit(m_domain);
  int count=0;
  for(bit.begin(); bit.ok(); ++bit, ++count)
    {
      if(m_bits[count])
        {
          IntVect iv(bit());
          iv.coarsen(iref);
          long index = newDomain.index(iv);
          newSet.m_bits.setTrue(index);
        }
    }

  *this = newSet;
}

void DenseIntVectSet::refine(int iref)
{
  if(iref == 1) return;
  if(isEmpty()) return;
  CH_assert(iref >= 1);
  //int refinements = iref/2;
  CH_assert((iref/2)*2 == iref); // check iref for power of 2
  Box newDomain(m_domain);
  newDomain.refine(iref);
  DenseIntVectSet newSet(newDomain, false);
  IntVect iv;
  BoxIterator bit(newDomain);
  int count=0;
  for(bit.begin(); bit.ok(); ++bit, ++count)
    {
      iv = bit();
      iv.coarsen(iref);
      if(this->operator[](iv))
        {
          newSet.m_bits.setTrue(count);
        }
    }
  *this = newSet;
}

void DenseIntVectSet::nestingRegion(int radius, const Box& domain)
{
  CH_assert(radius >= 0);
  if (radius == 0) return;

  DenseIntVectSet tmp(*this);

  {
    IntVect lo = m_domain.smallEnd();
    IntVect hi = m_domain.bigEnd();
    for(int i=0; i<SpaceDim; ++i)
      {
        if(lo[i] != domain.smallEnd()[i]) lo[i] += radius;
        if(hi[i] != domain.bigEnd()[i])   hi[i] -= radius;
      }
    Box shrink(lo, hi);
    *this &= shrink;
  }

  Box clobberBox(IntVect::Zero, IntVect::Zero);
  IntVect center(IntVect::Zero);
  clobberBox.grow(radius);
  BoxIterator bit(m_domain);
  int i=0;
  for(bit.begin(); bit.ok();++i, ++bit)
    {
      if(!(tmp.m_bits[i])) // was there a zero at this bit ?
        {
          if(domain.contains(bit()))
            {
              clobberBox.shift(bit()-center);
              *this -= clobberBox;
              center=bit();
            }
        }
    }
}

void DenseIntVectSet::nestingRegion(int radius, const ProblemDomain& a_domain)
{
  CH_assert(radius >= 0);
  if (radius == 0) return;

  DenseIntVectSet tmp;
  Box region = a_domain.domainBox();
  if(!a_domain.isPeriodic())
    {
      tmp = *this;
    }
  else
    {
      D_TERM6(if(a_domain.isPeriodic(0)) region.grow(0, radius);,
              if(a_domain.isPeriodic(1)) region.grow(1, radius);,
              if(a_domain.isPeriodic(2)) region.grow(2, radius);,
              if(a_domain.isPeriodic(3)) region.grow(3, radius);,
              if(a_domain.isPeriodic(4)) region.grow(4, radius);,
              if(a_domain.isPeriodic(5)) region.grow(5, radius);)
        tmp = DenseIntVectSet(region, false);
      tmp |= *this;
      for(int i=0; i<SpaceDim; ++i)
        {
          if(a_domain.isPeriodic(i))
            {
              int size = a_domain.domainBox().size(i);
              Box hiRegion = adjCellHi(a_domain.domainBox(), i, radius);
              for(BoxIterator bit(hiRegion); bit.ok(); ++bit)
                {
                  IntVect image = bit();
                  image[i] -= size;
                  if(this->operator[](image))
                    {
                      tmp |= bit();
                    }
                  image[i] += (size - radius);
                  if(this->operator[](image))
                    {
                      image[i] -= size;
                      tmp |= image;
                    }
                }
            }
        }
    }

  {
    IntVect lo = m_domain.smallEnd();
    IntVect hi = m_domain.bigEnd();
    for(int i=0; i<SpaceDim; ++i)
      {
        if(lo[i] != region.smallEnd()[i]) lo[i] += radius;
        if(hi[i] != region.bigEnd()[i])   hi[i] -= radius;

      }
    Box shrink(lo, hi);
    *this &= shrink;
  }

  Box clobberBox(IntVect::Zero, IntVect::Zero);
  IntVect center(IntVect::Zero);
  clobberBox.grow(radius);
  BoxIterator bit(region);
  int i=0;
  for(bit.begin(); bit.ok();++i, ++bit)
    {
      if(!(tmp.m_bits[i])) // was there a zero at this bit ?
        {
                  clobberBox.shift(bit()-center);
                  *this -= clobberBox;
                  center=bit();
        }
    }
}

void DenseIntVectSet::grow(int igrow)
{
  if(igrow >= 1)
    {
      IntVect range(IntVect::Unit);
      range*=igrow;
      grow(range);
    }
  else if(igrow < 0)
    {
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          grow(idir, igrow);
        }
    }
  else
    {
      //nothing to do if igrow == 0
      return;
    }

}

void DenseIntVectSet::grow(int idir, int igrow)
{
  CH_assert(idir >= 0);
  CH_assert(idir < SpaceDim);
  if(igrow >= 1)
    {
      IntVect range(IntVect::Zero);
      range[idir] = igrow;
      grow(range);
    }
  else if(igrow < 0)
    {
      IntVect shiftvec= igrow*BASISV(idir);
      DenseIntVectSet ivsShiftPlus = *this;
      DenseIntVectSet ivsShiftMinu = *this;
      ivsShiftPlus.shift(shiftvec);
      ivsShiftMinu.shift(-shiftvec);
      *this &= ivsShiftPlus;
      *this &= ivsShiftMinu;
    }
  else
    {
      //nothing to do if igrow == 0
      return;
    }

}

void DenseIntVectSet::grow(const IntVect& range)
{
  Box newDomain(m_domain);
  newDomain.grow(range);
  DenseIntVectSet newSet(newDomain, false);

  int index = 0;
#if CH_SPACEDIM == 1
  for(int i=0; i<m_domain.size(0); ++i, ++index)
    {
      if(m_bits[index])
        {
          int indexNew = i;
          for(int ii = 0; ii<2*range[0] + 1; ++ii, ++indexNew)
            newSet.m_bits.setTrue(indexNew);
        }
    }
#endif

#if CH_SPACEDIM == 2

  for(int j = 0; j<m_domain.size(1); ++j)
    {
      for(int i=0; i<m_domain.size(0); ++i, ++index)
        {
          if(m_bits[index])
            {
              for(int jj=0; jj<2*range[1] +1 ; ++jj)
                {
                  int indexNew = i + newDomain.size(0)*(j+jj);
                  for(int ii = 0; ii<2*range[0] + 1; ++ii, ++indexNew)
                    newSet.m_bits.setTrue(indexNew);
                }
            }
        }
    }
#endif

#if CH_SPACEDIM == 3
  for(int k=0; k<m_domain.size(2); ++k)
    for(int j = 0; j<m_domain.size(1); ++j)
      for(int i = 0; i<m_domain.size(0); ++i, ++index)
        {
          if(m_bits[index])
            {
              for(int kk=0; kk<2*range[2]+1; ++kk)
                for(int jj=0; jj<2*range[1] +1 ; ++jj)
                  {
                    int indexNew = i+
                      newDomain.size(0)*(j+jj) +
                      newDomain.size(1)*newDomain.size(0)*(k+kk);
                    for(int ii = 0; ii<2*range[0] + 1; ++ii, ++indexNew)
                      newSet.m_bits.setTrue(indexNew);
                  }
            }
        }

#endif


#if CH_SPACEDIM == 4
  for(int u=0; u<m_domain.size(3); u++)
    for(int k=0; k<m_domain.size(2); ++k)
      for(int j = 0; j<m_domain.size(1); ++j)
        for(int i = 0; i<m_domain.size(0); ++i, ++index)
          {
            if(m_bits[index])
              {
                for(int uu=0; uu<2*range[3]+1; ++uu)
                  for(int kk=0; kk<2*range[2]+1; ++kk)
                    for(int jj=0; jj<2*range[1] +1 ; ++jj)
                      {
                        int indexNew = i+
                          newDomain.size(0)*(j+jj) +
                          newDomain.size(0)*newDomain.size(1)*(k+kk) +
                          newDomain.size(0)*newDomain.size(1)*newDomain.size(2)*(u+uu);
                        for(int ii = 0; ii<2*range[0] + 1; ++ii, ++indexNew)
                          newSet.m_bits.setTrue(indexNew);
                      }
              }
          }

#endif

#if CH_SPACEDIM == 5
  for(int v=0; v<m_domain.size(4); v++)
    for(int u=0; u<m_domain.size(3); u++)
      for(int k=0; k<m_domain.size(2); ++k)
        for(int j = 0; j<m_domain.size(1); ++j)
          for(int i = 0; i<m_domain.size(0); ++i, ++index)
            {
              if(m_bits[index])
                {
                  for(int vv=0; vv<2*range[4]+1; ++vv)
                    for(int uu=0; uu<2*range[3]+1; ++uu)
                      for(int kk=0; kk<2*range[2]+1; ++kk)
                        for(int jj=0; jj<2*range[1] +1 ; ++jj)
                          {
                            int indexNew = i+
                              newDomain.size(0)*(j+jj) +
                              newDomain.size(0)*newDomain.size(1)*(k+kk) +
                              newDomain.size(0)*newDomain.size(1)*newDomain.size(2)*(u+uu) +
                              newDomain.size(0)*newDomain.size(1)*newDomain.size(2)*newDomain.size(3)*(v+vv);
                            for(int ii = 0; ii<2*range[0] + 1; ++ii, ++indexNew)
                              newSet.m_bits.setTrue(indexNew);
                          }
                }
            }

#endif

  // this is starting to look pretty ugly..
#if CH_SPACEDIM == 6
  for(int w=0; w<m_domain.size(5); w++)
    for(int v=0; v<m_domain.size(4); v++)
      for(int u=0; u<m_domain.size(3); u++)
        for(int k=0; k<m_domain.size(2); ++k)
          for(int j = 0; j<m_domain.size(1); ++j)
            for(int i = 0; i<m_domain.size(0); ++i, ++index)
              {
                if(m_bits[index])
                  {
                    for(int ww=0; ww<2*range[5]+1; ++ww)
                      for(int vv=0; vv<2*range[4]+1; ++vv)
                        for(int uu=0; uu<2*range[3]+1; ++uu)
                          for(int kk=0; kk<2*range[2]+1; ++kk)
                            for(int jj=0; jj<2*range[1] +1 ; ++jj)
                              {
                                int indexNew = i+
                                  newDomain.size(0)*(j+jj) +
                                  newDomain.size(0)*newDomain.size(1)*(k+kk) +
                                  newDomain.size(0)*newDomain.size(1)*newDomain.size(2)*(u+uu) +
                                  newDomain.size(0)*newDomain.size(1)*newDomain.size(2)*newDomain.size(3)*(v+vv) +
                                  newDomain.size(0)*newDomain.size(1)*newDomain.size(2)*newDomain.size(3)*newDomain.size(4)*(v+vv);
                                for(int ii = 0; ii<2*range[0] + 1; ++ii, ++indexNew)
                                  newSet.m_bits.setTrue(indexNew);
                              }
                  }
              }

#endif

#if (CH_SPACEDIM > 6)
  MayDay::Error("DenseIntVectSet::grow undefined for DIM>3");
#endif

  *this = newSet;
}

DenseIntVectSet& DenseIntVectSet::operator-=(const DenseIntVectSet& ivs)
{
  BoxIterator bit(m_domain);
  int i=0;
  for(bit.begin(); bit.ok();++i, ++bit)
    {
      if(m_bits[i])
        {
          if(ivs[bit()])
            m_bits.setFalse(i);
        }
    }
  return *this;
}

bool DenseIntVectSet::isEmpty() const
{
  return m_bits.isEmpty();
}
bool DenseIntVectSet::isFull() const
{
  return m_bits.isFull();
}

int DenseIntVectSet::numPts() const
{
  int n = m_domain.numPts();
  int c=0;
  for(int i=0; i<n; ++i)
    if(m_bits[i]) ++c;
  return c;
}

bool DenseIntVectSet::operator==(const DenseIntVectSet& a_lhs) const
{
  if(numPts() != a_lhs.numPts()) return false;
  DenseIntVectSetIterator it(*this);
  for(; it.ok(); ++it)
    {
      if(!a_lhs[it()]) return false;
    }
  return true;
}


bool DenseIntVectSet::operator<(const DenseIntVectSet& a_ivs) const
{
  // Primary criterion: Box::operator<() applied to m_domain.
  // Secondary criterion: BitSet::operator<() applied to m_bits.

  if( m_domain < a_ivs.m_domain )
  {
    return true;
  } else
  if( a_ivs.m_domain < m_domain )
  {
    return false;
  }

  if( m_bits < a_ivs.m_bits )
  {
    return true;
  }

  return false;
}


int DenseIntVectSet::linearSize() const
{
  if(isEmpty()){
    return CH_XD::linearSize<int>(0);
  }
  return m_bits.linearSize() + 2 * CH_XD::linearSize<Box>(m_minBox);
}

void DenseIntVectSet::linearIn(const void* const inBuf)
{
  static DenseIntVectSet emptySet;
  *this = emptySet;
  int* b = (int*)inBuf;
  if(*b == 0) return;

  char* buf = (char*)inBuf;
  m_bits.linearIn(buf);
  buf+=m_bits.linearSize();
  CH_XD::linearIn<Box>(m_domain, buf);
  buf+=CH_XD::linearSize<Box>(m_domain);
  CH_XD::linearIn<Box>(m_minBox, buf);
}

void DenseIntVectSet::linearOut(void* const a_outBuf) const
{
  if(isEmpty())
    {
      int* b = (int*)a_outBuf;
      *b = 0;
      return;
    }
  char* buf = (char*)a_outBuf;
  m_bits.linearOut(buf);
  buf+=m_bits.linearSize();
  CH_XD::linearOut<Box>(buf, m_domain);
  buf+=CH_XD::linearSize<Box>(m_domain);
  CH_XD::linearOut<Box>(buf, m_minBox);
}

void DenseIntVectSet::compact() const
{
  DenseIntVectSet* nthis = (DenseIntVectSet*)this;
  if(isEmpty()){
    *nthis = DenseIntVectSetIterator::emptyDenseIntVectSet;
    return;
  } else if (isFull()){
    return;
  }
  DenseIntVectSetIterator it(*this);
  Box nDomain(it(), it());
  IntVect& lo = (IntVect&)(nDomain.smallEnd());
  IntVect& hi = (IntVect&)(nDomain.bigEnd());
  for(++it; it.ok(); ++it)
    {
      lo.min(it());
      hi.max(it());
    }
  nDomain.computeBoxLenNotEmpty();
  DenseIntVectSet nset(nDomain);
  nset.intersect(*this);
  *nthis = nset;
}

DenseIntVectSet  DenseIntVectSetIterator::emptyDenseIntVectSet;

void DenseIntVectSetIterator::nextIntVect()
{
  ++(m_current[0]);
#if CH_SPACEDIM == 1
  // I don't think there is anything to do here
#else
  // don't declare these unless they're used
  const IntVect& big = m_ivsPtr->m_domain.bigEnd();
  const IntVect& small = m_ivsPtr->m_domain.smallEnd();

  if(m_current[0] > big[0])
    {
      m_current[0] = small[0];
      ++m_current[1];
    }
#if CH_SPACEDIM == 3
  if(m_current[1] > big[1])
    {
      m_current[1] = small[1];
      ++m_current[2];
    }

#elif (CH_SPACEDIM > 3)
  for (int dir=1; dir<SpaceDim-1; dir++)
    {
      if (m_current[dir] > big[dir])
        {
          m_current[dir] = small[dir];
          ++m_current[dir+1];
        }
    }

#endif
#endif
}

void DenseIntVectSetIterator::nextIntVect(int skip)
{
  if(skip == 1) {nextIntVect(); return;}

#if CH_SPACEDIM == 1
  // 1D is pretty simple
  m_current[0] += skip;
#endif

#if CH_SPACEDIM == 2

  int j_jump = skip / isize;
  m_current[1] += j_jump;
  m_current[0] += skip - j_jump*isize;
  if(m_current[0] > bigi)
    {
      ++(m_current[1]);
      m_current[0] -= isize;
    }
#endif
#if CH_SPACEDIM == 3

  int jump = skip / ijsize;
  m_current[2] += jump;
  skip -=   jump*ijsize;
  jump  =   skip / isize;
  m_current[1] += jump;
  skip -=   jump*isize;
  m_current[0] += skip;

  if(m_current[0] > bigi)
    {
      ++(m_current[1]);
      m_current[0] -= isize;
    }
  if(m_current[1] > bigj)
    {
      ++(m_current[2]);
      m_current[1] -= ijsize/isize;
    }

#endif
#if (CH_SPACEDIM > 3)
  // hack for now to make sure code actually works....
  for(;skip>0; --skip) nextIntVect();
#endif


}

void DenseIntVectSetIterator::begin()
{
  if(m_ivsPtr->m_bits.isEmpty())
    {
      m_iterator.end();
      return;
    }
  isize =  m_ivsPtr->m_domain.size(0);
  bigi  =   m_ivsPtr->m_domain.bigEnd()[0];
#if CH_SPACEDIM > 2
  ijsize = isize * m_ivsPtr->m_domain.size(1);
  bigj  =   m_ivsPtr->m_domain.bigEnd()[1];
#endif
  m_iterator = BitSetIterator(m_ivsPtr->m_bits);
  m_current = m_ivsPtr->m_domain.smallEnd();
  int i=0;
  while(m_iterator.ok() && !m_iterator())
    {
      ++i;
      ++m_iterator;
    }
  nextIntVect(i);
}
#include "NamespaceFooter.H"
