#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "FluxBox.H"
#include "NamespaceHeader.H"

// first do simple access functions
// ---------------------------------------------------------
int
FluxBox::nComp() const
{
  return m_nvar;
}

// ---------------------------------------------------------
const Box&
FluxBox::box() const
{
  return m_bx;
}

// ---------------------------------------------------------
FArrayBox&
FluxBox::getFlux(const int dir)
{
  CH_assert(m_nvar >0);
  CH_assert(dir < SpaceDim);
  CH_assert(m_fluxes[dir] != NULL);

  return *m_fluxes[dir];
}

// ---------------------------------------------------------
const FArrayBox&
FluxBox::getFlux(const int dir) const
{
  CH_assert(m_nvar >0);
  CH_assert(dir < SpaceDim);
  CH_assert(m_fluxes[dir] != NULL);

  return *m_fluxes[dir];
}

// ---------------------------------------------------------
FArrayBox&
FluxBox::operator[] (const int dir)
{
  CH_assert(m_nvar >0);
  CH_assert(dir < SpaceDim);
  CH_assert(m_fluxes[dir] != NULL);

  return *m_fluxes[dir];
}

// ---------------------------------------------------------
const FArrayBox&
FluxBox::operator[] (const int dir)  const
{
  CH_assert(m_nvar >0);
  CH_assert(dir < SpaceDim);
  CH_assert(m_fluxes[dir] != NULL);

  return *m_fluxes[dir];
}

// ---------------------------------------------------------
// constructors and destructors
// ---------------------------------------------------------
FluxBox::FluxBox() : m_fluxes(SpaceDim,NULL)
{
  m_nvar = -1;
}

// ---------------------------------------------------------
FluxBox::FluxBox(const Box& a_bx, int a_nComp)  : m_fluxes(SpaceDim,NULL)
{
  define(a_bx, a_nComp);
}

// ---------------------------------------------------------
FluxBox::~FluxBox()
{
  clear();
}

// ---------------------------------------------------------
void
FluxBox::clear()
{
  // first delete storage
  for (int dir = 0; dir < SpaceDim; dir++)
  {
    if (m_fluxes[dir] != NULL)
    {
      delete m_fluxes[dir];
      m_fluxes[dir] = NULL;
    }
  }

  // now reset all other variables
  m_nvar = -1;

  // set the box to the empty box...
  m_bx = Box();
}

// ---------------------------------------------------------
// define function
void
FluxBox::define(const Box& a_bx, int a_nComp)
{
  CH_assert(a_nComp > 0);

  m_bx = a_bx;
  m_nvar = a_nComp;

  if (m_fluxes.size() == 0)
  {
    m_fluxes.resize(SpaceDim,NULL);
  }

  for (int dir = 0; dir < SpaceDim; dir++)
  {
    Box edgeBox(surroundingNodes(m_bx,dir));

    if (m_fluxes[dir] != NULL)
    {
      delete m_fluxes[dir];
    }

    m_fluxes[dir] = new FArrayBox(edgeBox,m_nvar);
  }
}

// ---------------------------------------------------------
// should resize fluxes in space (could be faster than re-allocating
// storage)
void
FluxBox::resize(const Box& a_bx, int a_nComp)
{
  // if this object has not already been defined, call define fn.
  if (m_nvar < 0)
    {
      define(a_bx, a_nComp);
    }
  else
    {
      CH_assert(a_nComp > 0);
      m_bx = a_bx;
      m_nvar = a_nComp;

      for (int dir = 0; dir < SpaceDim; dir++)
        {
          Box edgeBox(surroundingNodes(m_bx, dir));
          if (m_fluxes[dir] != NULL)
            {
              m_fluxes[dir]->resize(edgeBox,m_nvar);
            }
          else
            {
              FArrayBox* newFabPtr = new FArrayBox(edgeBox, m_nvar);
              m_fluxes[dir] = newFabPtr;
            }
        }
    }
}

// ---------------------------------------------------------
void
FluxBox::setVal(const Real val)
{
  CH_assert(m_nvar > 0);

  for (int dir = 0; dir < SpaceDim; dir++)
    {
      CH_assert(m_fluxes[dir] != NULL);
      m_fluxes[dir]->setVal(val);
    }
}

// ---------------------------------------------------------
void
FluxBox::setVal(const Real val, const int dir)
{
  CH_assert(dir < SpaceDim);

  CH_assert(m_fluxes[dir] != NULL);
  m_fluxes[dir]->setVal(val);
}

// ---------------------------------------------------------
void
FluxBox::setVal(const Real val, const int dir, const int startComp,
                const int nComp)
{
  CH_assert(startComp >-1);
  CH_assert(startComp + nComp <= m_nvar);
  CH_assert(dir < SpaceDim);
  CH_assert(m_fluxes[dir] != NULL);

  for (int comp=startComp; comp < startComp+nComp; comp++)
    {
      m_fluxes[dir]->setVal(val,comp);
    }

}

// ---------------------------------------------------------
void
FluxBox::setVal(const Real val, const Box& bx)
{
  CH_assert(m_bx.contains(bx));

  for (int dir = 0; dir<SpaceDim; dir++) {
    CH_assert(m_fluxes[dir] != NULL);
    // move cell-centered box to appropriate edge
    Box edgeBox(surroundingNodes(bx, dir));
    m_fluxes[dir]->setVal(val,edgeBox,0,m_nvar);
  }

}

// ---------------------------------------------------------
void
FluxBox::setVal(const Real val, const Box& bx, const int dir,
                const int startComp, const int nComp)
{
  CH_assert(m_bx.contains(bx));
  CH_assert(m_fluxes[dir] != NULL);
  CH_assert(startComp > -1);
  CH_assert(startComp +nComp <=m_nvar);

  Box edgeBox(surroundingNodes(bx, dir));
  m_fluxes[dir]->setVal(val, edgeBox, startComp, nComp);
}

// ---------------------------------------------------------
void
FluxBox::copy(const FluxBox& src)
{
  CH_assert(src.box() == m_bx);
  CH_assert(src.nComp() == m_nvar);

  for (int dir=0; dir < SpaceDim; dir++)
    {
      CH_assert(m_fluxes[dir] != NULL);
      m_fluxes[dir]->copy(src[dir]);
    }

}

// ---------------------------------------------------------
void
FluxBox::copy(const FluxBox& src, const int srcComp,
              const int destComp, const int numComp)
{
  // to ensure that neither comp is negative
  CH_assert(srcComp*destComp > -1);
  CH_assert(srcComp+numComp <= src.nComp());
  CH_assert(destComp+numComp <= m_nvar);

  for (int dir=0; dir<SpaceDim; dir++)
    {
      CH_assert(m_fluxes[dir] != 0);
      const FArrayBox& srcFab = src[dir];
      m_fluxes[dir]->copy(srcFab,srcComp, destComp, numComp);
    }
}

// ---------------------------------------------------------
void
FluxBox::copy(const FluxBox& src, const int dir,
              const int srcComp, const int destComp,
              const int numComp)
{
  // to ensure that neither comp is negative
  CH_assert(srcComp*destComp > -1);
  CH_assert(srcComp+numComp <= src.nComp());
  CH_assert(destComp+numComp <= m_nvar);
  CH_assert(dir < SpaceDim);
  CH_assert(dir > -1);
  CH_assert(m_fluxes[dir] != NULL);

  const FArrayBox& srcFab = src[dir];
  m_fluxes[dir]->copy(srcFab,srcComp, destComp, numComp);
}

// ---------------------------------------------------------
void
FluxBox::copy(const FluxBox& a_src,
              const Box&     a_destbox)
{
  CH_assert(a_src.nComp() == m_nvar);
  CH_assert(m_bx.contains(a_destbox));
  CH_assert(a_src.box().contains(a_destbox));

  Interval comps(0, m_nvar-1);
  copy(a_destbox, comps, a_src, comps);
}


// ---------------------------------------------------------
void
FluxBox::copy(const Box& R, const Interval& Cdest, const FluxBox& src,
              const Interval& Csrc)
{

  for (int dir=0; dir<SpaceDim; dir++)
    {
      CH_assert(m_fluxes[dir] != NULL);
      Box Redge(R);
      Redge.surroundingNodes(dir);
      const FArrayBox& srcFab = src[dir];
      // all this intersecting is necessary due to the face-centered
      // nature of things
      Redge &= srcFab.box();
      Redge &= m_fluxes[dir]->box();
      if (!Redge.isEmpty())
        {
          //this is probably wrong in periodic---dtg
          m_fluxes[dir]->copy(Redge,Cdest,Redge,srcFab,Csrc);
        }
    }

}

// ---------------------------------------------------------
void
FluxBox::copy(const Box& srcbox,
              const Interval& destcomps,
              const Box& destbox,
              const FluxBox& src,
              const Interval& srccomps)
{

  for (int dir=0; dir<SpaceDim; dir++)
    {
      CH_assert (m_fluxes[dir] != NULL);
      Box srcEdgeBox(srcbox);
      srcEdgeBox.surroundingNodes(dir);
      Box destEdgeBox(destbox);
      destEdgeBox.surroundingNodes(dir);
      const FArrayBox& srcFab = src[dir];
      // safety check -- due to edge-centered nature of things,
      // destbox may not be contained in m_fluxes[dir]
      destEdgeBox &= (m_fluxes[dir]->box());
      m_fluxes[dir]->copy(srcEdgeBox, destcomps,
                          destEdgeBox, srcFab, srccomps);
    }

}

// ---------------------------------------------------------
FluxBox&
FluxBox::negate()
{
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      CH_assert(m_fluxes[dir] != NULL);
      m_fluxes[dir]->negate();
    }
  return *this;
}

// ---------------------------------------------------------
FluxBox&
FluxBox::negate(int        comp,
                int        numcomp)
{
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      CH_assert(m_fluxes[dir] != NULL);
      m_fluxes[dir]->negate(comp,numcomp);
    }
  return *this;
}

// ---------------------------------------------------------
FluxBox&
FluxBox::negate(const Box& subbox,
                int        comp,
                int        numcomp)
{
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      CH_assert(m_fluxes[dir] != NULL);
      m_fluxes[dir]->negate(subbox,comp,numcomp);
    }
  return *this;
}

// ---------------------------------------------------------
FluxBox&
FluxBox::plus(const FluxBox& a_src,
              const Box&     a_subbox,
              int            a_srccomp,
              int            a_destcomp,
              int            a_numcomp)
{
  CH_assert(m_bx.contains(a_subbox));
  CH_assert(a_src.box().contains(a_subbox));

  for (int dir=0; dir<SpaceDim; dir++)
    {
      Box faceBox(surroundingNodes(a_subbox,dir));
      // now call corresponding FArrayBox function
      m_fluxes[dir]->plus(a_src[dir], faceBox,
                         a_srccomp, a_destcomp, a_numcomp);
    }

  return *this;
}

// ---------------------------------------------------------
FluxBox&
FluxBox::minus(const FluxBox& a_src,
               const Box&     a_subbox,
               int            a_srccomp,
               int            a_destcomp,
               int            a_numcomp)
{
  CH_assert(m_bx.contains(a_subbox));
  CH_assert(a_src.box().contains(a_subbox));

  for (int dir=0; dir<SpaceDim; dir++)
    {
      Box faceBox(surroundingNodes(a_subbox,dir));
      // now call corresponding FArrayBox function
      m_fluxes[dir]->minus(a_src[dir], faceBox,
                           a_srccomp, a_destcomp, a_numcomp);
    }

  return *this;
}


// ---------------------------------------------------------
FluxBox&
FluxBox::mult(const FluxBox& a_src,
              const Box&     a_subbox,
              int            a_srccomp,
              int            a_destcomp,
              int            a_numcomp)
{
  CH_assert(m_bx.contains(a_subbox));
  CH_assert(a_src.box().contains(a_subbox));

  for (int dir=0; dir<SpaceDim; dir++)
    {
      Box faceBox(surroundingNodes(a_subbox,dir));
      // now call corresponding FArrayBox function
      m_fluxes[dir]->mult(a_src[dir], faceBox,
                          a_srccomp, a_destcomp, a_numcomp);
    }

  return *this;
}


// ---------------------------------------------------------
FluxBox&
FluxBox::divide(const FluxBox& a_src,
                const Box&     a_subbox,
                int            a_srccomp,
                int            a_destcomp,
                int            a_numcomp)
{
  CH_assert(m_bx.contains(a_subbox));
  CH_assert(a_src.box().contains(a_subbox));

  for (int dir=0; dir<SpaceDim; dir++)
    {
      Box faceBox(surroundingNodes(a_subbox,dir));
      // now call corresponding FArrayBox function
      m_fluxes[dir]->divide(a_src[dir], faceBox,
                            a_srccomp, a_destcomp, a_numcomp);
    }

  return *this;
}


// ---------------------------------------------------------
FluxBox&
FluxBox::operator+= (Real r)
{
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      CH_assert(m_fluxes[dir] != NULL);
      m_fluxes[dir]->plus(r);
    }
  return *this;
}

// ---------------------------------------------------------
FluxBox&
FluxBox::operator+= (const FluxBox& f)
{
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      CH_assert(m_fluxes[dir] != NULL);
      m_fluxes[dir]->plus(f[dir]);
    }
  return *this;
}

// ---------------------------------------------------------
FluxBox&
FluxBox::operator-= (Real r)
{
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      CH_assert(m_fluxes[dir] != NULL);
      m_fluxes[dir]->plus(-r);
    }
  return *this;
}

// ---------------------------------------------------------
FluxBox&
FluxBox::operator-= (const FluxBox& f)
{
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      CH_assert(m_fluxes[dir] != NULL);
      m_fluxes[dir]->minus(f[dir]);
    }
  return *this;
}

// ---------------------------------------------------------
FluxBox&
FluxBox::operator*= (Real r)
{
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      CH_assert(m_fluxes[dir] != NULL);
      m_fluxes[dir]->mult(r);
    }
  return *this;
}

// ---------------------------------------------------------
FluxBox&
FluxBox::operator*= (const FluxBox& f)
{
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      CH_assert(m_fluxes[dir] != NULL);
      m_fluxes[dir]->mult(f[dir]);
    }
  return *this;
}

// ---------------------------------------------------------
FluxBox&
FluxBox::shift(const IntVect& iv)
{
  m_bx.shift(iv);

  for (int dir=0; dir<SpaceDim; dir++)
    {
      m_fluxes[dir]->shift(iv);
    }

  return *this;

}

// ---------------------------------------------------------
int
FluxBox::size(const Box& bx, const Interval& comps) const
{
  int totalSize = 0;

  FArrayBox tempFab;
  for (int dir =0; dir<SpaceDim; dir++)
    {
      Box edgeBox(surroundingNodes(bx,dir));
      const int dirSize = tempFab.size(edgeBox,comps);
      totalSize += dirSize;
    }

  return totalSize;
}

// ---------------------------------------------------------
void
FluxBox::linearOut(void* buf, const Box& R, const Interval& comps) const
{
  Real* buffer = (Real*) buf;
  for (int dir=0; dir<SpaceDim; dir++)
    {
      CH_assert(m_fluxes[dir] != NULL);
      Box dirBox(surroundingNodes(R,dir));
      int dirSize = m_fluxes[dir]->size(dirBox, comps);
      m_fluxes[dir]->linearOut(buffer, dirBox, comps);
      buffer += dirSize/sizeof(Real);
    }
}

// ---------------------------------------------------------
void
FluxBox::linearIn(void* buf, const Box& R, const Interval& comps)
{
  Real* buffer = (Real*) buf;
  for (int dir=0; dir<SpaceDim; dir++)
    {
      CH_assert(m_fluxes[dir] != NULL);
      Box dirBox(surroundingNodes(R,dir));
      int dirSize = m_fluxes[dir]->size(dirBox, comps);
      m_fluxes[dir]->linearIn(buffer,dirBox,comps);
      buffer += dirSize/sizeof(Real);
    }
}
#include "NamespaceFooter.H"
