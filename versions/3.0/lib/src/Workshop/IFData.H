#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _IFDATA_H_
#define _IFDATA_H_

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif


#include <map>
using std::map;
#include "Vector.H"
#include "REAL.H"
#include "IndexTM.H"

#include "Notation.H"
#include "BaseIF.H"


using std::endl;
using std::cerr;
#include "NamespaceHeader.H"

template <int dim> class IFData
{
public:
  // data types to record vertex signs and edge intersections
  typedef IndexTM<int,dim>                      EdgeIndex;
  typedef map<EdgeIndex,Real,LexLT<EdgeIndex> > EdgeIntersections;
  typedef IndexTM<int,dim>                      Vertex;
  typedef map<Vertex,int,LexLT<Vertex> >        CornerSigns;
  typedef IndexTM<Real,dim>                     Rvect;

  // member data
  CornerSigns                m_cornerSigns;
  EdgeIntersections          m_intersections;
  BaseIF*                    m_function;
  Rvect                      m_dx;
  RvgDim                     m_globalCellCenter;
  Vector<int>                m_normalDir;
  Vector<int>                m_changingDir;
  Rvect                      m_normal;
  Vector<IndexTM<Real,dim> > m_gradNormal;
  bool                       m_allVerticesIn;
  bool                       m_allVerticesOut;
  bool                       m_allVerticesOn;

  // empty constructor
  inline IFData();

  // copy constructor
  inline IFData(const IFData<dim>& a_info);

  // makes an initial IFData
  inline IFData(const BaseIF      & a_function,
                const Rvect       & a_dx,
                const RvgDim      & a_globalCellCenter,
                const Vector<int> & a_normalDir = Vector<int>());

  // makes a reduced IFData
  inline IFData(const IFData<dim+1> & a_hInfo,
                const int           & a_idir,
                const int           & a_hilo);

  // destructor
  inline ~IFData();

  inline void reduceGradNormal(const Vector< IndexTM<Real,GLOBALDIM> > & a_gradNormal,
                               const IndexTM<Real,GLOBALDIM>           & a_normal);

  inline void setNormal();

  inline void makeCornerSigns();
  inline void findIntersectionPts();

  // used by the constructor
  inline bool isConnected(int          & a_edgeDir,
                          const Vertex & a_vertex1,
                          const Vertex & a_vertex2);

  inline void  makeEdgeKey(const int    & a_edgeDir,
                           const Vertex & a_vertex1,
                           const Vertex & a_vertex2);

  //this function calls the midPtRootFinder
  inline Real rootFinder(const EdgeIndex & a_thisEdge);

  inline Real midPtRootFinder(const Rvect & a_loPt,
                              const Rvect & a_hiPt) const;

  inline Real BrentRootFinder(const RvgDim & a_loPt,
                              const RvgDim & a_hiPt,
                              const int    & a_edgeDir) const;

  inline void checkIntersection(bool       & a_hiOn,
                                bool       & a_loOn,
                                const Real & a_pt) const;

  inline void remakeCornerSigns();

  // output functions
  inline void print(ostream& out) const;

  inline void dump() const;

  // operators
  inline void operator=(const IFData & a_ifData);

  // returns (a_vertex[idir]*a_dx[idir])
  inline RvgDim relCoord(const Vertex & a_vertex);

  // this function sets the changing direction given the normal direction
  inline void setChangingDirection();

  inline bool searchNormalDir(const int i) const;
};

template<> class IFData<1>
{
public:
  // true is in the fluid
  typedef Real         EdgeIntersection;
  typedef int          Vertex;
  typedef map<int,int> CornerSigns;
  typedef Real         Rvect;

  // member data
  CornerSigns      m_cornerSigns;
  EdgeIntersection m_intersection;
  Rvect            m_dx;
  Rvect            m_normal; // this is a dummy normal
  bool             m_allVerticesIn;
  bool             m_allVerticesOut;
  bool             m_allVerticesOn;

  inline IFData();

  inline IFData(const IFData<1> & a_ifData);

  inline IFData(const IFData<2> & a_hIFData,
                const int       & a_idir,
                const int       & a_hilo);

  inline ~IFData();

  inline void print(ostream& a_out) const;

  inline void operator=(const IFData<1> & a_ifData);
};

#include "IFDataImplem.H"

#include "NamespaceFooter.H"

#endif
