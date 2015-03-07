#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _INDEXTM_H_
#define _INDEXTM_H_

//
// This code supports the chombotomization-of-Tempest project.
//
// What we have here (and in IndexTMI.H) specifically is a templatized
// unification of IntVect and RealVect.
//
// Author: Ted
//


#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include "SPACE.H"
#include "REAL.H"
#include "Misc.H"
#include "GenericArithmetic.H"
#include "NamespaceHeader.H"

class HDF5Handle;

template<typename T, int N> class IndexTM;
template<typename T, int N> IndexTM<T,N> min(const IndexTM<T,N>& p1,
                                             const IndexTM<T,N>& p2);
template<typename T, int N> IndexTM<T,N> max(const IndexTM<T,N>& p1,
                                             const IndexTM<T,N>& p2);
template<typename T, int N> IndexTM<T,N> scale(const IndexTM<T,N>& p,
                                               T                   s);
//template<typename T, int N> IndexTM<T,N> scale(T                   s,
//                                               const IndexTM<T,N>& p );
template<typename T, int N> IndexTM<T,N> reflect(const IndexTM<T,N>& a,
                                                 T                   ref_ix,
                                                 int                 idir);
template<typename T, int N> IndexTM<T,N> diagShift(const IndexTM<T,N>& p,
                                                   T                   s);
template<typename T, int N> IndexTM<T,N> coarsen(const IndexTM<T,N>& p,
                                                 T                   s);
template<typename T, int N> IndexTM<T,N> coarsen(const IndexTM<T,N>& p1,
                                                 const IndexTM<T,N>& p2);
template<typename T, int N> std::ostream& operator<<(std::ostream&       os,
                                                     const IndexTM<T,N>& iv);
template<typename T, int N> std::istream& operator>>(std::istream& os,
                                                     IndexTM<T,N>& iv);

///
/**
   Returns a basis vector in the given coordinate direction.<br>
   In 3-D:
   BASISV(0) == (1,0,0); BASISV(1) == (0,1,0); BASISV(2) == (0,0,1).<br>
   In 2-D:
   BASISV(0) == (1,0); BASISV(1) == (0,1).<br>
   Note that the coordinate directions are based at zero.
*/
template<typename T, int N> IndexTM<T,N> BASISV_TM(int dir);

template<typename T, int N=CH_SPACEDIM> class IndexTM
  : public GenericArithmeticable< T,IndexTM<T,N> >
{
public:
  /**
     \name Constructors and Accessors
    See derived class IndexTM for more constructors.
  */
  /*@{*/

  ///
  /**
     Construct an IndexTM whose components are uninitialized.
  */
  IndexTM() : GenericArithmeticable< T, IndexTM<T,N> >( this )
  {}

  ///
  /**
     Destructor.
  */
  ~IndexTM()
  {}

  // Each of these, in its implementation, has a STATIC_ASSERT to check that
  // N is the same as the number of arguments:
  explicit IndexTM( T i );
  IndexTM( T i, T j );
  IndexTM( T i, T j, T k );
  IndexTM( T i, T j, T k, T m );
  IndexTM( T i, T j, T k, T m, T n );
  IndexTM( T i, T j, T k, T m, T n, T p );


  ///
  /**
     Construct an IndexTM setting the coordinates to the corresponding
     values in the integer array <i>a</i>.
  */
  explicit IndexTM (const T* a);

  ///
  /**
     The copy constructor.
  */
  IndexTM (const IndexTM& rhs);

  IndexTM copy() const {return *this;}

//operator IndexTM<Real,N>();

  ///
  /**
     The assignment operator.
  */
  IndexTM& operator= (const IndexTM& rhs);

  ///
  /**
     Returns a modifiable lvalue reference to the <i>i</i>'th coordinate of the
     IndexTM.
  */
  inline
  T& operator[] (int i);

  ///
  /**
     Returns the <i>i</i>'th coordinate of the IndexTM.
  */
  inline
  T operator[] (int i) const;

  ///
  /**
     Set <i>i</i>'th coordinate of IndexTM to <i>val</i>.
  */
  void setVal (int i,
               T val);

  ///
  /**
     set all values to <i>val</i>
  */
  void setAll(T val);

  /*@}*/

  /**
     \name Data pointer functions
  */
  /*@{*/

  ///
  /**
     Returns a const pointer to an array of coordinates of the IndexTM.
     Useful for arguments to FORTRAN calls.
  */
  const T* getVect () const;

  ///
  /**
     Only for sending to Fortran
   */
  const T*  dataPtr() const;

  ///
  /**
     Only for sending to Fortran
   */
  T*  dataPtr();

  /*@}*/

  /**
     \name Comparison Operators
  */
  /*@{*/

  ///
  /**
     Returns true if this IndexTM is equivalent to argument IndexTM.  All
     comparisons between analogous components must be satisfied.
  */
  bool operator== (const IndexTM& p) const;

  ///
  /**
     Returns true if this IndexTM is different from argument IndexTM.
     All comparisons between analogous components must be satisfied.
  */
  bool operator!= (const IndexTM& p) const;

  ///
  /**
     Returns true if this IndexTM is lexically less than the argument.
     An IndexTM MUST BE either lexically less than, lexically greater
     than, or equal to another IndexTM.

     iv1 is lexically less than iv2 if:

     in 2-D:<br>
     (iv1[0] < iv2[0]) || ((iv1[0] == iv2[0]) && (iv1[1] < iv2[1]));

     in 3-D:<br>
     (iv1[0] < iv2[0]) || (iv1[0]==iv2[0] && ((iv1[1] < iv2[1] || ((iv1[1] == iv2[1]) && (iv1[2] < iv2[2])))));
  */
  bool lexLT (const IndexTM& s) const;

  ///
  /**
     Returns true if this IndexTM is lexically greater than the
     argument.  An IndexTM MUST BE either lexically less than,
     lexically greater than, or equal to another IndexTM.

     iv1 is lexically less than iv2 if:

     in 2-D:<br>
     (iv1[0] > iv2[0]) || ((iv1[0] == iv2[0]) && (iv1[1] > iv2[1]));

     in 3-D:<br>
     (iv1[0] > iv2[0]) || (iv1[0]==iv2[0] && ((iv1[1] > iv2[1] || ((iv1[1] == iv2[1]) && (iv1[2] > iv2[2])))));
  */
  bool lexGT (const IndexTM& s) const;

  /*@}*/

  /**
     \name Unary operators
  */
  /*@{*/

  ///
  /**
     Unary plus -- for completeness.
  */
  IndexTM operator+ () const;

  ///
  /**
     Unary minus -- negates all components of this IndexTM.
  */
  IndexTM operator- () const;

  ///
  T dotProduct( const IndexTM& a_rhs ) const;

  ///
  /**
     Sum of all components of this IndexTM.
  */
  T sum () const;

  ///
  /**
     Product of all components of this IndexTM.
  */
  T product () const;

  /*@}*/


  /** Used by GenericArithmeticable to implement all pointwise arithmetic
   *  and comparison functions.
  */
  template<typename OP> bool operatorCompare(const IndexTM<T,N>&, const OP&) const;
  template<typename OP> IndexTM<T,N>& operatorOpEquals(
    const IndexTM<T,N>&, const OP&);
  template<typename OP> IndexTM<T,N>& operatorOpEquals(const T&, const OP&);

  ///
  /**
     Modifies this IndexTM by division of each component into T(1).
  */
  IndexTM& reciprocal();


  ///
  /**
     Component with the minimum value of this Index(returns 0 if they are all the same).
     a_doAbs : if true then take the absolute value before comparing
  */
  int minDir(bool a_doAbs) const;

  ///
  /**
     Component with the maximum value of this Index (returns 0 if they are all the same).
     a_doAbs : if true then take the absolute value before comparing.
  */
  int maxDir(bool a_doAbs) const;


  /*@}*/

  /**
     \name Other arithmetic operators
  */
  /*@{*/

  ///
  /**
     Modifies this IndexTM by taking component-wise min with IndexTM
     argument.
  */
  IndexTM& min (const IndexTM& p);

  ///
  /**
     Modifies this IndexTM by taking component-wise max with IndexTM
     argument.
  */
  IndexTM& max (const IndexTM& p);

  ///
  /**
     Modifies this IndexTM by multiplying each component by a scalar.
  */
  IndexTM& scale (T s);

  ///
  /**
     Modifies IndexTM by reflecting it in the plane defined by the
     index <i>ref_ix</i> and with normal in the direction of <i>idir</i>.
     Directions are based at zero.
  */
  IndexTM& reflect(T ref_ix,
                   int idir);

  ///
  /**
     Modifies this IndexTM by adding <i>s</i> to component in given coordinate
     direction.
  */
  IndexTM& shift (int coord,
                  T s);

  ///
  /**
     Modifies this IndexTM by component-wise addition with IndexTM
     argument.
  */
  IndexTM& shift (const IndexTM& iv);

  ///
  /**
     Modifies this IndexTM by adding a scalar <i>s</i> to each component.

  */
  IndexTM& diagShift (T s);

  ///
  /**
     Modify IndexTM by component-wise integer projection.
  */
  IndexTM& coarsen (const IndexTM& p);

  ///
  /**
     Modify IndexTM by component-wise integer projection.
  */
  IndexTM& coarsen (T p);

  /*@}*/

  /**
     \name I/O Functions
  */
  /*@{*/

  ///
  /**
     Print an IndexTM to the ostream.
  */
  void printOn (std::ostream& os) const;

  ///
  /**
     Print an IndexTM to the pout().
  */
  void p() const;

  ///
  /**
     Print an IndexTM to the ostream a bit more verbosely.
  */
  void dumpOn (std::ostream& os) const;

  ///
  /**
     Print the IndexTM to given output stream in ASCII.
  */
  friend std::ostream& operator<< <>(std::ostream&       os,
                                     const IndexTM& iv);

  ///
  /**
     Read next IndexTM from given input stream.
  */
  friend std::istream& operator>> <> (std::istream& os,
                                      IndexTM& iv);

  /*@}*/

  /**
     \name Constants
  */
  /*@{*/

  /**
     This is an IndexTM all of whose components are equal to zero.
  */
  static const IndexTM Zero;

  /**
     This is an IndexTM all of whose components are equal to one.
  */
  static const IndexTM Unit;

  /**
     Initializes Zero and Unit.
  */
  static int InitStatics();


  /**
    Low-level data copy.
  */
  void linearIn( const void* inBuf );
  void linearOut( void* outBuf ) const;

  /*@}*/

protected:
  friend class HDF5Handle;

  /**
     The individual components of this IndexTM.
  */
  T vect[N];
};


/**
    Useful for preventing division by IntVects, e.g.
        CH_assert( ! IndexTraits<T>::isInt );
    or better yet,
        STATIC_ASSERT( ! IndexTraits<T>::isInt );
*/
template<typename T> struct IndexTraits
{
    static bool const isInt=false;
    static bool const isReal=false;
};
template<> struct IndexTraits<int>
{
    static bool const isInt=true;
    static bool const isReal=false;
};
template<> struct IndexTraits<Real>
{
    static bool const isInt=false;
    static bool const isReal=true;
};

//
// Static initialization.  Gotta make sure there are no static object
// definitions above here (except possibly stuff in headers).  Otherwise,
// the danger is that some static object's constructor calls IndexTM::Zero or
// IndexTM::Unit -- the very things the following definition is supposed to
// initialize.
//
static int s_dummyForIntVectH1( IndexTM<int,1>::InitStatics() );
static int s_dummyForIntVectH2( IndexTM<int,2>::InitStatics() );
static int s_dummyForIntVectH3( IndexTM<int,3>::InitStatics() );
static int s_dummyForIntVectH4( IndexTM<int,4>::InitStatics() );
static int s_dummyForIntVectH5( IndexTM<int,5>::InitStatics() );
static int s_dummyForIntVectH6( IndexTM<int,6>::InitStatics() );
static int s_dummyForRealVectH1( IndexTM<Real,1>::InitStatics() );
static int s_dummyForRealVectH2( IndexTM<Real,2>::InitStatics() );
static int s_dummyForRealVectH3( IndexTM<Real,3>::InitStatics() );
static int s_dummyForRealVectH4( IndexTM<Real,4>::InitStatics() );
static int s_dummyForRealVectH5( IndexTM<Real,5>::InitStatics() );
static int s_dummyForRealVectH6( IndexTM<Real,6>::InitStatics() );

#include "NamespaceFooter.H"
#include "IndexTMI.H"

#endif // include guard
