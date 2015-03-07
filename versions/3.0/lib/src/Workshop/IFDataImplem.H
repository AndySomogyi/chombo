#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _IFDATAIMPLEM_H_
#define _IFDATAIMPLEM_H_

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include <iostream>

// empty constructor
template <int dim> inline IFData<dim>::IFData()
{
  m_function = NULL;
}

// empty constructor (dim == 1)
inline IFData<1>::IFData()
  :m_dx(0.0)
{
}

template <int dim> inline IFData<dim>::IFData(const IFData<dim>& a_iFData)
  :m_cornerSigns(a_iFData.m_cornerSigns),
   m_intersections(a_iFData.m_intersections),
   m_dx(a_iFData.m_dx),
   m_globalCellCenter(a_iFData.m_globalCellCenter),
   m_normalDir(a_iFData.m_normalDir),
   m_changingDir(a_iFData.m_changingDir),
   m_normal(a_iFData.m_normal),
   m_gradNormal(a_iFData.m_gradNormal),
   m_allVerticesIn(a_iFData.m_allVerticesIn),
   m_allVerticesOut(a_iFData.m_allVerticesOut),
   m_allVerticesOn(a_iFData.m_allVerticesOn)
{
  // Copy the implicit function; m_function is a pointer
  // maps require the ability to copy null constructed objects
  if (a_iFData.m_function != NULL)
  {
    m_function = a_iFData.m_function->newImplicitFunction();
  }
  else
  {
    m_function = NULL;
  }
}

inline IFData<1>::IFData(const IFData<1>& a_iFData)
  :m_cornerSigns(a_iFData.m_cornerSigns),
   m_intersection(a_iFData.m_intersection),
   m_dx(a_iFData.m_dx),
   m_allVerticesIn(a_iFData.m_allVerticesIn),
   m_allVerticesOut(a_iFData.m_allVerticesOut),
   m_allVerticesOn(a_iFData.m_allVerticesOn)
{
}

// Constructor from the implicit function
template <int dim> inline IFData<dim>::IFData(const BaseIF      & a_function,
                                              const Rvect       & a_dx,
                                              const RvgDim      & a_globalCellCenter,
                                              const Vector<int> & a_normalDir)
  :m_dx(a_dx),
   m_globalCellCenter(a_globalCellCenter)
{
  // set the normal directions and the changing directions
  m_normalDir = a_normalDir;
  setChangingDirection();

  // Copy the implicit function; m_function is a pointer
  m_function = a_function.newImplicitFunction();

  // set normal
  setNormal();

  // TODO = reduce Grad Normal for arbitrary directions
  Vector< IndexTM<Real,GLOBALDIM> > globalGradNormal = m_function->gradNormal(m_globalCellCenter);

  if (GLOBALDIM==dim)
  {
    for (int i = 0; i < dim; i++)
    {
      m_gradNormal.push_back(IndexTM<Real,dim>::Zero);
      for (int j = 0; j < dim; j++)
      {
        m_gradNormal[i][j] = globalGradNormal[i][j];
      }
    }
  }
  else if (GLOBALDIM > dim)
  {
    reduceGradNormal(globalGradNormal,m_function->normal(m_globalCellCenter));
  }
  else
  {
    MayDay::Abort("Current dimension is higher than GLOBALDIM");
  }

  // makecornerSigns sets the bools m_allVerticesIn and m_allVerticesOut and
  // m_allVerticesOn
  makeCornerSigns();

  // find an intersection point, if it exists, of the interface with each edge
  findIntersectionPts();
}

template <int dim> inline IFData<dim>::IFData(const IFData<dim+1> & a_hIFData,
                                              const int           & a_idir,
                                              const int           & a_hilo)
{
  CH_TIME("IFData::ConstructorReduceInfo");

  typedef IndexTM<int,dim+1>HDEdgeIndex;
  typedef map<HDEdgeIndex,Real,LexLT<HDEdgeIndex> > HDEdgeIntersections;
  typedef IndexTM<int,dim+1>HDVertex;
  typedef map<HDVertex,int,LexLT<HDVertex> > HDCornerSigns;

  // pass the implicit function
  m_function = a_hIFData.m_function->newImplicitFunction();

  // cellCenter is in physical coordinates (a the global dimension)
  m_globalCellCenter = a_hIFData.m_globalCellCenter;

  // pm = -1 for lo and 1 for hi
  int pm = a_hilo;
  pm *= 2;
  pm -= 1;

  // cell center for this bd cut cell
  m_globalCellCenter[a_idir] += pm*0.5*a_hIFData.m_dx[a_idir];

  // TODO
  m_normalDir.push_back(a_idir);
  setChangingDirection();

  m_allVerticesIn  = true;
  m_allVerticesOut = true;
  m_allVerticesOn  = true;
  for (typename HDCornerSigns::const_iterator it = a_hIFData.m_cornerSigns.begin();
       it != a_hIFData.m_cornerSigns.end(); ++it)
    {
      const HDVertex& hVertex = it->first;
      if (hVertex[a_idir] == a_hilo)
        {
          Vertex vertex;
           for (int j = 0; j < dim; ++j)
             {
               if (j < a_idir)
                 {
                   vertex[j] = hVertex[j];
                 }
               else
                 {
                   vertex[j] = hVertex[j+1];
                 }
             }

           m_cornerSigns[vertex] = it->second;
           if (m_cornerSigns[vertex] == IN)
             {
               m_allVerticesOut = false;
               m_allVerticesOn  = false;
             }
           else if (m_cornerSigns[vertex] == OUT)
             {
               m_allVerticesIn = false;
               m_allVerticesOn = false;
             }
        }
    }

  for (typename HDEdgeIntersections::const_iterator it = a_hIFData.m_intersections.begin();
       it != a_hIFData.m_intersections.end(); ++it)
    {
      const HDEdgeIndex& hEdgeIndex = it->first;
      EdgeIndex edgeIndex;
      int hEdgeDir = hEdgeIndex[0];

      // only interested in edgeDir != a_idir
      // Among these, need to find edges in the higher dim that are a_hilo in the a_idir component.
      // Example: a_idir = 1; a_hilo = 0. Observe that an edge in the x (N.B 0 < 1) direction
      // might look like (edgeDir= 0:yhilo,zhilo ). In other words, y information is in index 1.
      // On the other hand, an edge in the z direction(N.B 2 > 1)
      // might look like (edgeDir= 0:xhilo,yhilo ). In other words, y information is in index 2.
      if (hEdgeDir < a_idir)
        {
          if (hEdgeIndex[a_idir] == a_hilo)
            {
              edgeIndex[0] = hEdgeDir;

              for (int j = 1; j < dim; ++j)
                {
                  if (j < a_idir)
                    {
                      edgeIndex[j] = hEdgeIndex[j];
                    }
                  else
                    {
                      edgeIndex[j] = hEdgeIndex[j+1];
                    }
                }
              m_intersections[edgeIndex] = it->second;
            }
        }
      else if (hEdgeDir > a_idir)
        {
           if (hEdgeIndex[a_idir+1] == a_hilo)
            {
              edgeIndex[0] = hEdgeDir - 1;

              for (int j = 1; j < dim; ++j)
                {
                  if (j < a_idir+1)
                    {
                      edgeIndex[j] = hEdgeIndex[j];
                    }
                  else
                    {
                      edgeIndex[j] = hEdgeIndex[j+1];
                    }
                }
              m_intersections[edgeIndex] = it->second;
            }
        }
    }

  // m_dx reduced for this sub-polytope

  Real mag = 0.0;
  // normal and gradNormal should be taken at the new cell center
  IndexTM<Real,GLOBALDIM> localNormal = m_function->normal(m_globalCellCenter);
  Vector< IndexTM<Real,GLOBALDIM> > localGradNormal = m_function->gradNormal(m_globalCellCenter);

  if (!m_allVerticesOut)
    {
      for (int jdir = 0; jdir < dim+1; ++jdir)
        {
          if (jdir < a_idir)
            {
              m_dx[jdir]              = a_hIFData.m_dx[jdir];
              m_normal[jdir]          = localNormal[jdir];

              mag += m_normal[jdir]*m_normal[jdir];
            }
          else if (jdir > a_idir)
            {
              m_dx[jdir-1]              = a_hIFData.m_dx[jdir];
              m_normal[jdir-1]          = localNormal[jdir];

              mag += m_normal[jdir-1]*m_normal[jdir-1];
            }
        }
      // normal may legitimately be zero vector.
      if (mag > 0.0)
        {
          m_normal /= sqrt(mag);
        }
      reduceGradNormal(localGradNormal,localNormal);
    }
}

inline IFData<1>::IFData(const IFData<2> & a_2DIFData,
                         const int       & a_idir,
                         const int       & a_hilo)
{
  // we want the edge on the a_hilo side of the square with normal in the
  // a_idir direction

  // dx for this edge
  m_dx = a_2DIFData.m_dx[(a_idir + 1)%2];

  IFData<2>localInfo = a_2DIFData;

  // This 2D edgeIndex locates the 1D edge in the edgeIntersection map
  IndexTM<int,2>twoDEdge;
  twoDEdge[0] = (a_idir + 1)%2;
  twoDEdge[1] = a_hilo;

  m_intersection = LARGEREALVAL;
  if (localInfo.m_intersections.find(twoDEdge) != localInfo.m_intersections.end())
  {
    m_intersection = localInfo.m_intersections[twoDEdge];
  }

  // This 2D vertex locates the hi and lo ends of the 1D segment in the
  // cornerSigns map
  IndexTM<int,2>loPt2D;
  loPt2D[(a_idir + 1)%2] = 0;
  loPt2D[a_idir] = a_hilo;

  IndexTM<int,2>hiPt2D;
  hiPt2D[(a_idir+ 1)%2] = 1;
  hiPt2D[a_idir] = a_hilo;

  if (localInfo.m_cornerSigns.find(loPt2D) != localInfo.m_cornerSigns.end())
  {
    m_cornerSigns[0] = localInfo.m_cornerSigns[loPt2D];
  }
  else
  {
    MayDay::Abort("Lo endpoint not in Map");
  }

  if (localInfo.m_cornerSigns.find(hiPt2D) != localInfo.m_cornerSigns.end())
  {
    m_cornerSigns[1] = localInfo.m_cornerSigns[hiPt2D];
  }
  else
  {
    MayDay::Abort("Hi endpoint not in Map");
  }

  // set bools
  m_allVerticesIn  = true;
  m_allVerticesOut = true;
  m_allVerticesOn  = true;

  if (m_cornerSigns[0] != ON || m_cornerSigns[1] != ON)
  {
    m_allVerticesOn  = false;
  }

  if (m_cornerSigns[0] == IN || m_cornerSigns[1] == IN)
  {
    m_allVerticesOut = false;
  }

  if (m_cornerSigns[0] == OUT || m_cornerSigns[1] == OUT)
  {
    m_allVerticesIn = false;
  }
}

// Destructor
template <int dim> inline IFData<dim>::~IFData()
{
  if (m_function != NULL)
  {
    delete m_function;
  }
}

// Destructor (dim == 1)
inline IFData<1>::~IFData()
{
}

template <int dim> inline void IFData<dim>::reduceGradNormal(const Vector< IndexTM<Real,GLOBALDIM> > & a_gradNormal,
                                                             const IndexTM<Real,GLOBALDIM>           & a_normal)
{
  // computes m_gradNormal, given the normal in dim+1 and the gradNormal in
  // dim+1 if we are reducing from the GLOBALDIM to GLOBALDIM-1, the
  // changingDir are m_changingDir

  // generates the denominator
  Real norm = 0.0;
  for (int idir = 0; idir < dim; idir++)
  {
    m_gradNormal.push_back(IndexTM<Real,dim>::Zero);
    norm += a_normal[m_changingDir[idir]] * a_normal[m_changingDir[idir]];
  }

  Real sumNormSquared = norm;
  norm = pow(norm,(Real)1.5);

  for (int j = 0; j < dim; j++)
  {
    for (int k = 0; k < dim; k++)
    {
      if (norm!=0.0)
      {
        for (int i = 0; i < dim; i++)
        {
          m_gradNormal[j][k] += a_normal[m_changingDir[i]]
                              * a_gradNormal[m_changingDir[i]][m_changingDir[k]];
        }

        m_gradNormal[j][k] *= -a_normal[m_changingDir[j]];
        m_gradNormal[j][k] += a_gradNormal[m_changingDir[j]][m_changingDir[k]]
                            * sumNormSquared;

        m_gradNormal[j][k] /= norm;
      }
      else
      {
        m_gradNormal[j][k] = a_gradNormal[m_changingDir[j]][m_changingDir[k]];
      }
    }
  }
}

template <int dim> inline void IFData<dim>::setNormal()
{
  IndexTM<Real,GLOBALDIM> globalNormal = m_function->normal(m_globalCellCenter);

  Real norm = 0.0;
  for (int idir = 0; idir < dim; idir++)
  {

    m_normal[idir] = globalNormal[m_changingDir[idir]];
    norm += m_normal[idir] * m_normal[idir];
  }

  if (norm!=0.0)
  {
    m_normal /= sqrt(norm);
  }
}

template <int dim> inline void IFData<dim>::makeCornerSigns()
{
  m_allVerticesIn  = true;
  m_allVerticesOut = true;
  m_allVerticesOn  = true;

  Vertex vertex;

  // there are 2^globalDim vertices
  int numVertices = 1 << dim;
  for (int i = 0; i < numVertices; ++i)
  {
    int ii = i;
    for (int j = 0; j < dim; ++j)
    {
      vertex[j] = (ii & 1);
      ii = ii >> 1;
    }

    RvgDim offset = relCoord(vertex);
    Real val = m_function->value(m_globalCellCenter + offset);

    // true = negative = in the computational domain = in the fluid
    if (val == 0)
    {
      m_cornerSigns[vertex] = ON;
    }
    else if (val + MACHINEPRECISION < 0)
    {
      m_cornerSigns[vertex] = IN;
      m_allVerticesOut      = false;
      m_allVerticesOn       = false;
    }
    else if (val - MACHINEPRECISION > 0)
    {
      m_cornerSigns[vertex] = OUT;
      m_allVerticesIn       = false;
      m_allVerticesOn       = false;
    }
    else
    {
      m_cornerSigns[vertex] = ON;
    }
  }
}

template <int dim> inline void IFData<dim>::findIntersectionPts()
{
  // double iterate through corner signs map,writing to a edge map if corner1 is connected to corner2
  for (typename CornerSigns::const_iterator it0 = m_cornerSigns.begin();
       it0 != m_cornerSigns.end(); ++it0)
  {
    for (typename CornerSigns::const_iterator it1 = m_cornerSigns.begin();
         it1 != m_cornerSigns.end(); ++it1)
    {
      int edgeDir = LARGEINTVAL;
      if (isConnected(edgeDir,it0->first,it1->first))
      {
        // make edge key:m_intersections[thisEdge] = pt;
        makeEdgeKey(edgeDir,it0->first,it1->first);
      }
    }
  }
}

template <int dim> inline bool IFData<dim>::isConnected(int          & a_edgeDir,
                                                        const Vertex & a_vertex1,
                                                        const Vertex & a_vertex2)
{
  bool connected = true;
  int numDif = 0;

  // connected = true if and only if  (a_vertex1 and a_vertex2 differ in
  // exactly one coordinate.)
  for (int idir = 0; idir < dim; ++idir)
  {
    if (a_vertex1[idir] != a_vertex2[idir])
    {
      // vertices differ in idir direction
      a_edgeDir = idir;
      numDif += 1;
    }
  }

  if (numDif == 1)
  {
    connected = true;
  }
  else
  {
    connected = false;
    a_edgeDir = LARGEINTVAL;
  }

  return connected;
}

template <int dim> inline void IFData<dim>::makeEdgeKey(const int    & a_edgeDir,
                                                        const Vertex & a_vertex0,
                                                        const Vertex & a_vertex1)
{
  EdgeIndex thisEdge;
  thisEdge[0] = a_edgeDir;
  for (int idir = 0; idir < dim; ++idir)
  {
    if (idir < a_edgeDir)
    {
      thisEdge[idir + 1] = a_vertex0[idir];
    }
    else if (idir > a_edgeDir)
    {
      thisEdge[idir] = a_vertex0[idir];
    }
  }

  int lo = LARGEINTVAL;
  int hi = LARGEINTVAL;
  if (m_cornerSigns.find(a_vertex0)!= m_cornerSigns.end())
  {
    lo = m_cornerSigns.find(a_vertex0)->second;
  }
  else
  {
    MayDay::Abort("Vertex not well defined in makeEdgeKey");
  }

  if (m_cornerSigns.find(a_vertex1)!= m_cornerSigns.end())
  {
    hi = m_cornerSigns.find(a_vertex1)->second;
  }
  else
  {
    MayDay::Abort("Vertex not well defined in makeEdgeKey");
  }

  if ((lo==OUT && hi==IN)||(lo==IN && hi==OUT))
  {
    // calculate a value between -1.0 and 1.0
    Real pt = rootFinder(thisEdge);

    // check whether intersection is at the end point
    bool hiOn;
    bool loOn;
    checkIntersection(hiOn,loOn,pt);
    if (hiOn || loOn)
    {
      if (hiOn)
      {
        m_cornerSigns[a_vertex1] = ON;
      }
      else if (loOn)
      {
        m_cornerSigns[a_vertex0] = ON;
      }

      remakeCornerSigns();
    }
    else
    {
      m_intersections[thisEdge] = pt;
    }
  }
}

template <int dim> inline Real IFData<dim>::rootFinder(const EdgeIndex& a_thisEdge)
{
  // this function is only called when dim == globalDim
  Real pt = LARGEREALVAL;

  int edgeDir = a_thisEdge[0];

  Vertex loEnd;
  loEnd[edgeDir] = 0;

  Vertex hiEnd;
  hiEnd[edgeDir] = 1;

  for (int idir = 0; idir < dim; ++idir)
  {
    if (idir < edgeDir)
    {
      loEnd[idir] = a_thisEdge[idir + 1];
      hiEnd[idir] = a_thisEdge[idir + 1];
    }
    else if (idir > edgeDir)
    {
      loEnd[idir] = a_thisEdge[idir];
      hiEnd[idir] = a_thisEdge[idir];
    }
  }

  // make relative to dx and relative to origin at cell center
  RvgDim loOffset = relCoord(loEnd);
  RvgDim hiOffset = relCoord(hiEnd);

  // move to physical cell center (at the globalDim)
  RvgDim loPt = loOffset + m_globalCellCenter;
  RvgDim hiPt = hiOffset + m_globalCellCenter;

  // get the edge direction in GLOBALDIM
  int globalEdgeDir = m_changingDir[edgeDir];

  // returns the distance from loEnd to intersection
  // pt = midPtRootFinder(loPt,hiPt);
  pt = BrentRootFinder(loPt,hiPt,globalEdgeDir);

  // make the intersection independent of dx
  // pt /= m_dx[edgeDir];
  // make the intersection  bewteen (-0.5,0.5)
  // pt -= 0.5;

  return pt;
}

template <int dim> inline Real IFData<dim>::midPtRootFinder(const Rvect & loPt,
                                                            const Rvect & hiPt) const
{
  Real dist = LARGEREALVAL;

  //  Max allowed iterations and floating point precision
  const unsigned int  MAXITER = 100;

#if defined(CH_USE_DOUBLE)
  const Real      EPS   = 3.0e-15;
#elif defined(CH_USE_FLOAT)
  const Real      EPS   = 3.0e-7;
#else
#error Unknown Chombo precision
#endif

  const Real tol = 10.0*EPS;
  Rvect aPt;
  Rvect bPt;
  Rvect xMid;

  Real phi_aPt = m_function->value(loPt);
  Real phi_bPt = m_function->value(hiPt);
  Real phi_xMid = LARGEREALVAL;

  aPt = loPt;
  bPt = hiPt;

  int iter;
  if (phi_aPt*phi_bPt > 0)
  {
    MayDay::Abort("root must be bracketed");
  }
  else
  {
    for (iter = 0; iter < MAXITER; iter++)
    {
      xMid = 0.5*(aPt + bPt);

      phi_xMid = m_function->value(xMid);
      if (Abs(phi_xMid) > tol)
      {
        if (phi_bPt*phi_xMid > 0)
        {
          bPt = xMid;
        }
        else
        {
          aPt = xMid;
        }
      }
      else
      {
        break;
      }
    }

    if (iter >= MAXITER)
    {
      cerr  << "rootFinder: exceeding maximum iterations: "
            << MAXITER << "\n";
    }
    else
    {
      dist = 0.0;
      for (int idir = 0; idir < SpaceDim; ++idir)
      {
        dist += (xMid[idir] - loPt[idir])*(xMid[idir] - loPt[idir]);
      }
      dist = sqrt(dist);
    }
  }

  return dist;
}

template <int dim> inline Real IFData<dim>::BrentRootFinder(const RvgDim & a_loPt,
                                                            const RvgDim & a_hiPt,
                                                            const int    & a_edgeDir) const
{
  // Max allowed iterations and floating point precision
  const unsigned int MAXITER = 100;
  const Real         EPS   = 3.0e-15;

  unsigned int i;
  Real aPt;
  Real bPt;
  Real c, fa, fb, fc;
  Real d, e;
  Real tol1, xm;
  Real p, q, r, s;

  aPt = -1.0;
  bPt =  1.0;

  RvgDim physCoordAPt = a_loPt;
  RvgDim physCoordBPt = a_hiPt;

  fa = -m_function->value(physCoordAPt);
  fb = -m_function->value(physCoordBPt);

  //  Init these to be safe
  c = d = e = 0.0;

  if (fb*fa > 0)
  {
    MayDay::Error("root must be bracketed");
  }

  fc = fb;

  for (i = 0; i < MAXITER; i++)
  {
    if (fb*fc > 0)
    {
      //  Rename a, b, c and adjust bounding interval d
      c = aPt;
      fc  = fa;
      d = bPt - aPt;
      e = d;
    }

    if (Abs(fc) < Abs(fb))
    {
      aPt = bPt;
      bPt = c;
      c = aPt;
      fa  = fb;
      fb  = fc;
      fc  = fa;
    }

    //  Convergence check
    tol1  = 2.0 * EPS * Abs(bPt) + 0.5 * TOLERANCE;
    xm    = 0.5 * (c - bPt);

    if (Abs(xm) <= tol1 || fb == 0.0)
    {
      break;
    }

    if (Abs(e) >= tol1 && Abs(fa) > Abs(fb))
    {
      //  Attempt inverse quadratic interpolation
      s = fb / fa;
      if (aPt == c)
      {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      }
      else
      {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * xm * q * (q-r) - (bPt-aPt) * (r-1.0));
        q = (q-1.0) * (r-1.0) * (s-1.0);
      }

      //  Check whether in bounds
      if (p > 0) q = -q;

      p = Abs(p);

      if (2.0 * p < Min(((Real)3.0)*xm*q-Abs(tol1*q), Abs(e*q)))
      {
        //  Accept interpolation
        e = d;
        d = p / q;
      }
      else
      {
        //  Interpolation failed, use bisection
        d = xm;
        e = d;
      }
    }
    else
    {
      //  Bounds decreasing too slowly, use bisection
      d = xm;
      e = d;
    }

    //  Move last best guess to a
    aPt = bPt;
    fa  = fb;

    //  Evaluate new trial root
    if (Abs(d) > tol1)
    {
      bPt = bPt + d;
    }
    else
    {
      if (xm < 0) bPt = bPt - tol1;
      else        bPt = bPt + tol1;
    }

    physCoordBPt[a_edgeDir] = ((1.0 - bPt)/2.0) * a_loPt[a_edgeDir] + ((1.0 + bPt)/2.0) * a_hiPt[a_edgeDir];
    fb = -m_function->value(physCoordBPt);
  }

  if (i >= MAXITER)
  {
    cerr << "BrentRootFinder: exceeding maximum iterations: "
         << MAXITER << "\n";
  }

  return bPt;
}

template <int dim> inline void IFData<dim>::checkIntersection(bool       & a_hiOn,
                                                              bool       & a_loOn,
                                                              const Real & a_pt) const
{
  a_hiOn = false;
  a_loOn = false;

  if (a_pt >= 1.0 - MACHINEPRECISION)
  {
    a_hiOn = true;
  }
  else if (a_pt <= -1.0 + MACHINEPRECISION)
  {
    a_loOn = true;
  }
}

template <int dim> inline void IFData<dim>::remakeCornerSigns()
{
  m_allVerticesIn  = true;
  m_allVerticesOut = true;
  m_allVerticesOn  = true;

  for (typename CornerSigns::const_iterator it = m_cornerSigns.begin();
       it != m_cornerSigns.end(); ++it)
  {
    if (it->second == IN)
    {
      m_allVerticesOut = false;
    }

    if (it->second == OUT)
    {
      m_allVerticesIn = false;
    }

    if (it->second != ON)
    {
      m_allVerticesOn = false;
    }
  }
}

template <int dim> inline void IFData<dim>::print(ostream& a_out) const
{
  for (typename CornerSigns::const_iterator it = m_cornerSigns.begin();
       it != m_cornerSigns.end(); ++it)
  {
    a_out << "Vertex " << "(" << it->first << ") = " << it->second << "\n";
  }

  for (typename EdgeIntersections::const_iterator it = m_intersections.begin();
       it != m_intersections.end(); ++it)
  {
    a_out << "EdgeIndex " << it->first << " = " << it->second << "\n";
  }

  a_out << "allVerticesIn  = " << m_allVerticesIn  << "\n";
  a_out << "allVerticesOut = " << m_allVerticesOut << "\n";
  a_out << "allVerticesOn  = " << m_allVerticesOn  << "\n";

  if (!m_allVerticesOut)
  {
    a_out << "m_dx = " << m_dx << "\n";
    a_out << "m_normalDir" << m_normalDir << "\n";
    a_out << "m_normal = " << m_normal << "\n";
    a_out << "m_gradNormal = " << m_gradNormal << "\n";
    a_out << "m_globalCellCenter = " << m_globalCellCenter << "\n";
  }
  else
  {
    a_out << "normal and dx not defined : all Vertices Out";
  }
}

inline void IFData<1>::print(ostream& a_out) const
{
  typedef  map<int,int> oneDCornerSigns;
  for (oneDCornerSigns::const_iterator it = m_cornerSigns.begin();
       it != m_cornerSigns.end(); ++it)
  {
    a_out << "Vertex " << "(" << it->first << ") = " << it->second << "\n";
  }

  a_out << "m_allVerticesIn  = " << m_allVerticesIn << "\n";
  a_out << "m_allVerticesOut = " << m_allVerticesOut << "\n";
  a_out << "m_allVerticesOn  = " << m_allVerticesOn  << "\n";

  if (!m_allVerticesOut)
  {
    a_out << "m_dx  = " << m_dx << "\n";
  }

  int lo = LARGEINTVAL;
  int hi = LARGEINTVAL;
  if (m_cornerSigns.find(0)!= m_cornerSigns.end())
  {
    lo = m_cornerSigns.find(0)->second;
  }
  else
  {
    MayDay::Abort("no lo in cornerSigns");
  }

  if (m_cornerSigns.find(1)!= m_cornerSigns.end())
  {
    hi = m_cornerSigns.find(1)->second;
  }
  else
  {
    MayDay::Abort("no hi in cornerSigns");
  }

  if (lo == OUT && hi == OUT)
  {
    a_out << "Edge Covered" << "\n";
  }
  else if (lo==IN  && hi==IN)
  {
    a_out << "Edge Uncovered" << "\n";
  }
  else if (lo==ON  || hi==ON)
  {
    a_out << "Edge has vertex on the interface" << "\n";
  }
  else
  {
    if (m_intersection == LARGEREALVAL)
    {
      MayDay::Warning("no intersection pt");
    }
    else
    {
      a_out << "Intersection pt is " << m_intersection << "\n";
    }
  }
}

template <int dim> inline void IFData<dim>::dump() const
{
  print(pout());
}

template <int dim> inline ostream& operator<<(ostream           & a_out,
                                              const IFData<dim> & a_iFData)
{
  a_iFData.print(a_out);
  return a_out;
}

// equals operator
template <int dim> inline void IFData<dim>::operator=(const IFData & a_iFData)
{
  m_cornerSigns      = a_iFData.m_cornerSigns;
  m_intersections    = a_iFData.m_intersections;

  // maps require the ability to copy null constructed objects
  if (a_iFData.m_function != NULL)
  {
    m_function     = a_iFData.m_function->newImplicitFunction();
  }
  else
  {
    m_function     = NULL;
  }

  m_dx               = a_iFData.m_dx;
  m_globalCellCenter = a_iFData.m_globalCellCenter;

  m_normalDir        = a_iFData.m_normalDir;
  m_changingDir      = a_iFData.m_changingDir;
  m_normal           = a_iFData.m_normal;

  m_gradNormal       = a_iFData.m_gradNormal;

  m_allVerticesIn    = a_iFData.m_allVerticesIn;
  m_allVerticesOut   = a_iFData.m_allVerticesOut;
  m_allVerticesOn    = a_iFData.m_allVerticesOn;
}

// equals operator
inline void IFData<1>::operator=(const IFData & a_iFData)
{
  m_cornerSigns    = a_iFData.m_cornerSigns;
  m_intersection   = a_iFData.m_intersection;

  m_dx             = a_iFData.m_dx;

  m_allVerticesIn  = a_iFData.m_allVerticesIn;
  m_allVerticesOut = a_iFData.m_allVerticesOut;
  m_allVerticesOn  = a_iFData.m_allVerticesOn;
}

// computes coordinates of a_vertex relative to the cell-center = (0,0...0)
template <int dim> inline IndexTM<Real,GLOBALDIM> IFData<dim>::relCoord(const Vertex & a_vertex)
{
RvgDim retval=IndexTM<Real,GLOBALDIM>::Zero;
Rvect hilo;

for (int idir = 0; idir < dim; ++idir)
{
  hilo[idir] = a_vertex[idir] - 0.5;
  retval[m_changingDir[idir]] = hilo[idir]*m_dx[idir];
}

return retval;
}

template <int dim> inline void IFData<dim>::setChangingDirection()
{
if (m_normalDir.size()==0)
{
  for (int i = 0 ; i < GLOBALDIM ; ++i)
  {
    m_changingDir.push_back(i);
  }
}
else
{
  for (int i = 0 ; i < GLOBALDIM ; ++i)
  {
    if (searchNormalDir(i) == false)
    {
      m_changingDir.push_back(i);
    }
  }
}
}

// this looks for the integer i in the normalDir vector, assuming m_normalDir is ordered
template <int dim> inline bool IFData<dim>::searchNormalDir(const int i) const
{
int size = m_normalDir.size();

for (int j = 0; j < size; j++)
{
  if (m_normalDir[j] == i)
  {
    return true;
  }
}

return false;
}

#endif
