#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//  ANAG, LBNL

#ifndef _FACEITERATOR_H_
#define _FACEITERATOR_H_

#include "IntVect.H"
#include "FaceIndex.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "NamespaceHeader.H"
class EBGraph;

#define FACEITER_CTOR_ARGS   const IntVectSet&           a_ivs,       \
                             const EBGraph&              a_ebgraph,   \
                             int                         a_direction, \
                             const FaceStop::WhichFaces& a_location

///
/**
   Ennumeration class to distinguish where FaceIterator will stop.
*/
class FaceStop
{
public:
  ///
  /**
     SurroundingWithBoundary  means stop at all faces on the high and
     low sides of IVS cells. \                                          \
     SurroundingNoBoundary  means stop at all faces on the high and
     low sides of IVS cells, excluding faces on the domain boundary.\   \
     LoWithBoundary  means stop at all faces on the
     low side of IVS cells.\\
     LoNoBoundary  means stop at all faces on the
     low side of IVS cells, excluding faces on the domain boundary.\    \
     HiWithBoundary  means stop at all faces on the
     high side of IVS cells.\\
     HiNoBoundary  means stop at all faces on the
     high side of IVS cells, excluding faces on the domain boundary.\   \
  */
  enum WhichFaces{Invalid=-1,
                  SurroundingWithBoundary=0, SurroundingNoBoundary=1,
                  HiWithBoundary, LoWithBoundary,
                  HiNoBoundary  , LoNoBoundary,
                  LoBoundaryOnly, HiBoundaryOnly, AllBoundaryOnly,
                  NUMTYPES};
};


/// Iterator over faces within an IntVectSet and an Ebgraph
/**
   Iterator over within an IntVectSet and an Ebgraph.
   The enumeration tells which faces on which to stop.
 */
class FaceIterator
{
public:

  /// general construction
  FaceIterator( FACEITER_CTOR_ARGS );

  /// specialized construction for most typical case: Irregular cells inside 'region'
  FaceIterator(const Box& region, const EBGraph& a_ebgraph, int a_direction, const FaceStop::WhichFaces& a_location);

  ///
  FaceIterator();

  /// Going with the default...
  // FaceIterator( const FaceIterator & that);

  /// Going with the default...
  // FaceIterator & operator=( const FaceIterator & that);

  ///
  ~FaceIterator();


  ///
  /**
    Looks in cache first, unless optional arg is true.  Delegates to
    this->doDefine() if FaceIterator not found in cache.
  */
  void define( FACEITER_CTOR_ARGS );

  ///
  void reset();

  ///
  void operator++();

  ///
  const FaceIndex& operator() () const;

  ///
  bool ok() const;

  ///
  bool isDefined() const;

  ///
  const Vector<FaceIndex>& getVector() const;

  ///
  //  static void clearCache();

  ///
  /**
    Enable or disable FaceIterator caching.  Caching is enabled by default.
    Seems to break in parallel.  Caveat emptor.
  */
  //  static void setUseCache( bool );

private:
  ///
  /**
    Does the real work.  this->define() looks in the cache first.
  */
  void doDefine( FACEITER_CTOR_ARGS );

  bool              m_isDefined;
  Vector<FaceIndex> m_faces;
  int               m_iface;
  int               m_direction;


};



//#undef FACEITER_CTOR_ARGS

#include "NamespaceFooter.H"
#endif
