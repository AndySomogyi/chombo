#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _INDEXTMI_H_
#define _INDEXTMI_H_

#include "SPACE.H"
#include <iostream>
using std::ostream;
using std::istream;
using std::ws;

#include "MayDay.H"
#include "Misc.H"
#include "IndexTM.H"
#include "parstream.H"
#include <cmath>
#include <algorithm>
#include <functional>
#include <numeric>
#include "Metaprograms.H"
#include "NamespaceHeader.H"


template<typename T, int N> ostream&
operator<< (ostream&       os,
            const IndexTM<T,N>& p)
{
    os << '(' << p[0];
    for( int i=1;i<N;++i )
    {
        os << ',' << p[i];
    }
    os << ')';
    if (os.fail())
        MayDay::Error("operator<<(ostream&,Index&) failed");
    return os;
}

//
// Copied from <Utility.H>
//
#define CH_IGNORE_MAX 100000

template<typename T, int N> istream&
operator>> (istream& is,
            IndexTM<T,N>& p)
{
    is >> ws;
    char c;
    is >> c;
    is.putback(c);
    if (c == '(')
    {
        is.ignore(CH_IGNORE_MAX, '(') >> p[0];
        for( int i=1;i<N;++i )
        {
            is.ignore(CH_IGNORE_MAX, ',') >> p[i];
        }
        is.ignore(CH_IGNORE_MAX, ')');
    }
    else if (c == '<')
    {
        is.ignore(CH_IGNORE_MAX, '<') >> p[0];
        for( int i=1;i<N;++i )
        {
            is.ignore(CH_IGNORE_MAX, ',') >> p[1];
        }
        is.ignore(CH_IGNORE_MAX, '>');
    }
    else
        MayDay::Error("operator>>(istream&,Index&): expected \'(\' or \'<\'");

    if (is.fail())
        MayDay::Error("operator>>(istream&,Index&) failed");

    return is;
}


template<typename T, int N> void
IndexTM<T,N>::printOn (ostream& os) const
{
    os << "Index: " << *this << '\n';
}

template<typename T, int N> void
IndexTM<T,N>::p() const
{
    pout() << *this << '\n';
}

template<typename T, int N> void
IndexTM<T,N>::dumpOn (ostream& os) const
{
    os << "Index " << *this << '\n';
}

//
// Static object initialization.
//
template<typename T, int N> int
IndexTM<T,N>::InitStatics()
{
  IndexTM<T,N>& pz = const_cast<IndexTM<T,N>&>( IndexTM<T,N>::Zero );

  pz.setAll(0);

  IndexTM<T,N>& pu = const_cast<IndexTM<T,N>&>( IndexTM<T,N>::Unit );

  pu.setAll(1);

  // No danger of Index::Zero and Unit not having been allocated, as ARM section
  // 3.4 says "The initialization of nonlocal static objects in a translation unit
  // is done before the first use of any function...defined in that translation
  // unit."
  //
  // Had to go through the const_cast stuff because it's nice to be able to declare
  // Index::Zero and Index::Unit as const.

  return 0; // arbitrary
}

//
// Inlines.
//
template<typename T, int N> inline
IndexTM<T,N>::IndexTM(const T *a)
  : GenericArithmeticable< T, IndexTM<T,N> >( this )
{
  memcpy( vect, a, N*sizeof(T) );
}

template<typename T, int N> inline
IndexTM<T,N>::IndexTM(const IndexTM<T,N> &iv)
  : GenericArithmeticable< T, IndexTM<T,N> >( this )
{
  memcpy(vect, iv.vect, N*sizeof(T) );
}

template<typename T, int N> inline
IndexTM<T,N>&
IndexTM<T,N>::operator=(const IndexTM<T,N> &iv)
{
  memcpy( vect, iv.vect, N*sizeof(T) );
  return *this;
}

template<typename T, int N> inline
T&
IndexTM<T,N>::operator[] (int i)
{
  CH_assert(i>=0 && i < N);
  return vect[i];
}

template<typename T, int N> inline
T
IndexTM<T,N>::operator[] (int i) const
{
  CH_assert(i>=0 && i < N);
  return vect[i];
}

template<typename T, int N> inline
void
IndexTM<T,N>::setVal (int i,
                      T val)
{
  CH_assert(i >=0 && i < N);
  vect[i] = val;
}

template<typename T, int N> inline
void
IndexTM<T,N>::setAll (T val)
{
  for(int i=0; i<N; i++)
    vect[i] = val;
}


template<typename T, int N> inline
const T*
IndexTM<T,N>::dataPtr() const
{
  return vect;
}

template<typename T, int N> inline
T*
IndexTM<T,N>::dataPtr()
{
  return vect;
}

template<typename T, int N> inline
const T*
IndexTM<T,N>::getVect () const
{
  return vect;
}

template<typename T, int N> inline
bool
IndexTM<T,N>::operator== (const IndexTM& p) const
{
  return Metaprograms::pointwiseCompare<N,T,std::equal_to<T> >(vect,p.vect);
}

template<typename T, int N> inline
bool
IndexTM<T,N>::operator!= (const IndexTM& p) const
{
  return ! operator==(p);
}

template<typename T, int N> inline
bool
IndexTM<T,N>::lexLT (const IndexTM &s) const
{
  return Metaprograms::LexLT<N,T>()(vect,s.vect);
}

template<typename T, int N> inline
bool
IndexTM<T,N>::lexGT (const IndexTM& s) const
{
  return ! Metaprograms::LexLT<N,T>()(vect,s.vect);
}

template<typename T, int N> inline
IndexTM<T,N>
IndexTM<T,N>::operator+ () const
{
  return IndexTM<T,N>(*this);
}

template<typename T, int N> inline
IndexTM<T,N>
IndexTM<T,N>::operator- () const
{
  IndexTM<T,N> result(*this);
  for( int i=0;i<N;++i )
    {
      result.vect[i] *= -1;
    }
  return result;
}

template<typename T, int N> inline
T
IndexTM<T,N>::dotProduct ( const IndexTM<T,N>& a_rhs ) const
{
  return Metaprograms::InnerProduct<N,T,T,
                                    std::plus<T>,
                                    std::multiplies<T> >()( vect,
                                                            a_rhs.vect );
}

template<typename T, int N> inline
T
IndexTM<T,N>::sum () const
{
  return Metaprograms::Accum<N,T,std::plus<T> >()( vect );
}

template<typename T, int N> inline
T
IndexTM<T,N>::product () const
{
  return Metaprograms::Accum<N,T,std::multiplies<T> >()( vect );
}

template<typename T, int N> template<typename OP> bool
IndexTM<T,N>::operatorCompare(const IndexTM<T,N>& p, const OP& op) const
{
  return Metaprograms::pointwiseCompare<N,T,OP>(vect,p.vect);
}


template<typename T, int N> template<typename OP> inline
IndexTM<T,N>&
IndexTM<T,N>::operatorOpEquals(const T& s, const OP& op)
{
  Metaprograms::Transform<N,T,OP>()(vect,s);
  return *this;
}

template<typename T, int N> template<typename OP> inline
IndexTM<T,N>&
IndexTM<T,N>::operatorOpEquals(const IndexTM<T,N>& p, const OP& op)
{
  Metaprograms::Transform<N,T,OP>()( vect, p.vect );
  return *this;
}

template<typename T, int N> inline
IndexTM<T,N>&
IndexTM<T,N>::reciprocal()
{
  std::transform( vect, vect+N, vect, bind1st(std::divides<T>(),T(1)) );
  return *this;
}


template<typename T> static bool abscompare(const T& a, const T& b)
{
  return ::fabs(a) < ::fabs(b);
}
template<typename T, int N> inline
int
IndexTM<T,N>::minDir (bool a_doAbs ) const
{
  if( a_doAbs )
    {
      return std::min_element( vect, vect+N, std::ptr_fun(abscompare<T>) ) - vect;
    }
  else
    {
      return std::min_element( vect, vect+N ) - vect;
    }
}
template<typename T, int N> inline
int
IndexTM<T,N>::maxDir (bool a_doAbs ) const
{
  if( a_doAbs )
    {
      return std::max_element( vect, vect+N, std::ptr_fun(abscompare<T>) ) - vect;
    }
  else
    {
      return std::max_element( vect, vect+N ) - vect;
    }
}

template<typename T> static T ourmin(T a, T b) { return a<b ? a:b; }
template<typename T> static T ourmax(T a, T b) { return a>b ? a:b; }

template<typename T, int N> inline
IndexTM<T,N>&
IndexTM<T,N>::min (const IndexTM<T,N>& p)
{
  std::transform( vect, vect+N, p.vect, vect,
                  std::ptr_fun(ourmin<T>) );
  return *this;
}

template<typename T, int N> inline
IndexTM<T,N>&
IndexTM<T,N>::max (const IndexTM<T,N>& p)
{
  std::transform( vect, vect+N, p.vect, vect,
                  std::ptr_fun(ourmax<T>) );
  return *this;
}


template<typename T, int N> inline
IndexTM<T,N>&
IndexTM<T,N>::scale (T s)
{
    return (*this) *= s;
}

template<typename T, int N> inline
IndexTM<T,N>&
IndexTM<T,N>::reflect (T ref_ix,
                       int idir)
{
  CH_assert(idir >= 0 && idir < N);
  vect[idir] = -vect[idir] + 2*ref_ix;
  return *this;
}

template<typename T, int N> inline
IndexTM<T,N>&
IndexTM<T,N>::shift (int coord,
                  T s)
{
  CH_assert(coord >= 0 && coord < N);
  vect[coord] += s;
  return *this;
}

template<typename T, int N> inline
IndexTM<T,N>&
IndexTM<T,N>::shift (const IndexTM<T,N>& iv)
{
  return (*this) += iv;
}

template<typename T, int N> inline
IndexTM<T,N>&
IndexTM<T,N>::diagShift (T s)
{
  return (*this) += s;
}


template<typename T, int N> inline
IndexTM<T,N>
scale (const IndexTM<T,N>& p,
       T                   s)
{
  return p * s;
}

template<typename T, int N> inline
IndexTM<T,N>
diagShift(const IndexTM<T,N> &p, T s)
{
  return p + s;
}

template<typename T, int N> inline
IndexTM<T,N>
min(const IndexTM<T,N>& p1,
    const IndexTM<T,N>& p2)
{
  IndexTM<T,N> result( p1 );
  return result.min(p2);
}

template<typename T, int N> inline
IndexTM<T,N>
max(const IndexTM<T,N>& p1,
    const IndexTM<T,N>& p2)
{
  IndexTM<T,N> result( p1 );
  return result.max(p2);
}

template<typename T, int N> inline
IndexTM<T,N>
BASISV_TM(int dir)
{
  CH_assert(dir >= 0 && dir < N);
  IndexTM<T,N> tmp = IndexTM<T,N>::Zero ;
  tmp.dataPtr()[dir] = T(1);
  return tmp;
}

template<typename T, int N> inline
IndexTM<T,N>
reflect(const IndexTM<T,N>& a,
        T            ref_ix,
        int            idir)
{
  IndexTM<T,N> result( a );
  return result.reflect( ref_ix, idir );
}

template<typename T> static T ourcoarsen( T a, T b )
{
    return a<0 ? T(-::fabs(a+1))/b-1 : a/b;
}

template<typename T, int N> inline
IndexTM<T,N>
coarsen(const IndexTM<T,N>& p,
        T                   s)
{
  IndexTM<T,N> result( p );
  return result.coarsen( s );
}

template<typename T, int N> inline
IndexTM<T,N>
coarsen(const IndexTM<T,N>& p1,
        const IndexTM<T,N>& p2)
{
  IndexTM<T,N> result( p1 );
  return result.coarsen( p2 );
}

/*
template<typename T, int N>
IndexTM<T,N>::operator IndexTM<Real,N>()
{
  IndexTM<Real,N> result;
  for( int i=0;i<N;++i ) result.dataPtr()[i] = Real(vect[i]);
  return result;
}
*/

template<typename T, int N> inline
IndexTM<T,N>&
IndexTM<T,N>::coarsen (T s)
{
  CH_assert(s > 0);
  std::transform( vect, vect+N, vect,
                  std::bind2nd(std::ptr_fun(ourcoarsen<T>),s) );
  return *this;
}

template<typename T, int N> inline
IndexTM<T,N>&
IndexTM<T,N>::coarsen (const IndexTM<T,N>& p)
{
  CH_assert(p > (IndexTM<T,N>::Zero));
  std::transform( vect, vect+N, p.vect, vect, std::ptr_fun(ourcoarsen<T>) );
  return *this;
}

template<typename T, int N>
void
IndexTM<T,N>::linearIn( const void* inBuf )
{
  memcpy( vect, (T*)inBuf, N*sizeof(T) );
}
template<typename T, int N>
void
IndexTM<T,N>::linearOut( void* inBuf ) const
{
  memcpy( (T*)inBuf, vect, N*sizeof(T) );
}


template<typename T, int N>
IndexTM<T,N>::IndexTM( T i )
  : GenericArithmeticable< T, IndexTM<T,N> >( this )
{
  STATIC_ASSERT( N==1 );
  vect[0] = i;
}
template<typename T, int N>
IndexTM<T,N>::IndexTM( T i, T j )
  : GenericArithmeticable< T, IndexTM<T,N> >( this )
{
  STATIC_ASSERT( N==2 );
  vect[0] = i;
  vect[1] = j;
}
template<typename T, int N>
IndexTM<T,N>::IndexTM( T i, T j, T k )
  : GenericArithmeticable< T, IndexTM<T,N> >( this )
{
  STATIC_ASSERT( N==3 );
  vect[0] = i;
  vect[1] = j;
  vect[2] = k;
}
template<typename T, int N>
IndexTM<T,N>::IndexTM( T i, T j, T k, T m )
  : GenericArithmeticable< T, IndexTM<T,N> >( this )
{
  STATIC_ASSERT( N==4 );
  vect[0] = i;
  vect[1] = j;
  vect[2] = k;
  vect[3] = m;
}
template<typename T, int N>
IndexTM<T,N>::IndexTM( T i, T j, T k, T m, T n )
  : GenericArithmeticable< T, IndexTM<T,N> >( this )
{
  STATIC_ASSERT( N==5 );
  vect[0] = i;
  vect[1] = j;
  vect[2] = k;
  vect[3] = m;
  vect[4] = n;
}
template<typename T, int N>
IndexTM<T,N>::IndexTM( T i, T j, T k, T m, T n, T p )
  : GenericArithmeticable< T, IndexTM<T,N> >( this )
{
  STATIC_ASSERT( N==6 );
  vect[0] = i;
  vect[1] = j;
  vect[2] = k;
  vect[3] = m;
  vect[4] = n;
  vect[5] = p;
}


template<typename T, int N> const IndexTM<T,N> IndexTM<T,N>::Zero;
template<typename T, int N> const IndexTM<T,N> IndexTM<T,N>::Unit;
static int s_dummyForIntVectCpp1( IndexTM<int,1>::InitStatics() );
static int s_dummyForIntVectCpp2( IndexTM<int,2>::InitStatics() );
static int s_dummyForIntVectCpp3( IndexTM<int,3>::InitStatics() );
static int s_dummyForIntVectCpp4( IndexTM<int,4>::InitStatics() );
static int s_dummyForIntVectCpp5( IndexTM<int,5>::InitStatics() );
static int s_dummyForIntVectCpp6( IndexTM<int,6>::InitStatics() );
static int s_dummyForRealVectCpp1( IndexTM<Real,1>::InitStatics() );
static int s_dummyForRealVectCpp2( IndexTM<Real,2>::InitStatics() );
static int s_dummyForRealVectCpp3( IndexTM<Real,3>::InitStatics() );
static int s_dummyForRealVectCpp4( IndexTM<Real,4>::InitStatics() );
static int s_dummyForRealVectCpp5( IndexTM<Real,5>::InitStatics() );
static int s_dummyForRealVectCpp6( IndexTM<Real,6>::InitStatics() );
// If Index::Zero and Index::Unit were pointers, we wouldn't need this extra
// static int.  But they're objects, and the danger is that the initializations
// right above here ("Index Index::Zero;" and "Index Index::Unit;") are hit
// after the last call to Index::InitStatics, and in that case the
// Index::Index() constructor could redefine Zero and Unit.  In fact, the way
// things stand now, nothing bad would happen, because the Index::Index()
// constructor doesn't assign anything to any of the data members.  But we don't
// want to count on that always being the case.
#include "NamespaceFooter.H"

#endif // include guard
