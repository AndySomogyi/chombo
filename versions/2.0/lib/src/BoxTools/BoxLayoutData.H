#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _BOXLAYOUTDATA_H_
#define _BOXLAYOUTDATA_H_

#include "LayoutData.H"
#include "Interval.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "Copier.H"
#include "SPMD.H"
#include "memtrack.H"
#include "NamespaceHeader.H"

class MotionItem;

extern int LinearizationTest;

/// Factory object to data members of a BoxLayoutData container
template <class T> class DataFactory
{
public:

  virtual ~DataFactory(){;}
  /// factory function.  creates a new 'T' object
  /** creates a new 'T' object and returns a pointer to it.
   */
  virtual T* create(const Box& box, int ncomps, const DataIndex& a_datInd) const=0;

  virtual bool callDelete() const {return true;}
};

/// Factory object to data members of a BoxLayoutData container
/**
   creates new T object by invoking operator new

 \code
 template <class T>
 T* DefaultDataFactory<T>::create(const Box& box,
                                  int ncomps,
                                  const DataIndex& a_datInd) const
 {
   return new T(box, ncomps);
 }
 \endcode
 */

template <class T> class DefaultDataFactory : public DataFactory<T>
{
public:
  /// factory function.  creates a new 'T' object
  /** creates a new 'T' object and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user. */
  virtual T* create(const Box& box, int ncomps, const DataIndex& a_datInd) const;

};

class FABAliasDataFactory : public DataFactory<FArrayBox>
{
public:
  virtual ~FABAliasDataFactory()
  {}

  FABAliasDataFactory(const LayoutData<Real*>& aliases);
  void define(const LayoutData<Real*>& aliases);
  /// factory function.  creates a new 'T' object using an aliased dataPtr for T
  /** creates a new 'T' object and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user. */
  virtual FArrayBox* create(const Box& box, int ncomps, const DataIndex& a_datInd) const;

protected:
  LayoutData<Real*> aliasPtrs;
};

template<class T> class BoxLayoutData;

template <class T>
class AliasDataFactory : public DataFactory<T>
{
public:
  virtual ~AliasDataFactory()
  {}

  AliasDataFactory(BoxLayoutData<T>* a_original, const Interval& interval);
  void define(BoxLayoutData<T>* a_original, const Interval& interval);
    /** creates a new 'T' object and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user. */
  virtual T* create(const Box& box, int ncomps, const DataIndex& a_datInd) const;

protected:
  BoxLayoutData<T>* m_origPointer;
  Interval          m_interval;
};

template<class T> class LevelData;

template <class T>
class LDOperator
{
public:
  virtual ~LDOperator(){}    //to keep the compiler happy

  virtual int size(const T& arg, const Box& b, const Interval& comps) const
  {
    return arg.size(b, comps);
  }
  virtual void linearOut(const T& arg, void* buf, const Box& R,
                         const Interval& comps) const
  {
    arg.linearOut(buf, R, comps);
  }
  virtual void linearIn(T& arg,  void* buf, const Box& R,
                        const Interval& comps)const
  {
    arg.linearIn(buf, R, comps);
  }
  virtual void op(T& dest,
          const Box& RegionFrom,
          const Interval& Cdest,
          const Box& RegionTo,
          const T& src,
          const Interval& Csrc) const
  {
    if(LinearizationTest == 0){
      dest.copy(RegionFrom, Cdest,RegionTo, src, Csrc);
    } else {
      int sizeSource = size(src,  RegionFrom, Csrc);
      int sizeDest   = size(dest, RegionTo,   Cdest);
      
      if(T::preAllocatable() == 0 || T::preAllocatable() == 1){
        if(sizeSource != sizeDest)
          {
            MayDay::Abort("LinearizationTest failure: dest and source have different sizes");
          }
      }
      Vector<char> buffer(sizeSource);
      void* b = (void*)&(buffer[0]);
      linearOut(src, b, RegionFrom, Csrc);
      linearIn(dest, b, RegionTo,   Cdest);
    }
  }
};

/// Data on a BoxLayout
/**

  A collection of Box-oriented objects.  The arrangement
  of Boxes is given by the underlying BoxLayout object.
  BoxLayoutData attempt to prevent users from maniupulating
  the template class T to become out of synch with the
  boxes in the BoxLayout.  Caveat emptor.

  All access to the data in a BoxLayoutData is forced to be
  data-parallel, by use of the DataIterator.

  class T must provide the following methods:
    - <PRE> T() </PRE>
    - <PRE> T(const Box& box, int comps) </PRE>
    - <PRE> define(const Box& box, int comps) </PRE>
    - <PRE> void copy(const Box& R, const Interval& Cd, const T& source, const Interval Cs)
       Cs and Cd must have same length, but can be displaced </PRE>

    - <PRE> static int preAllocatable()
     returns 0 if the size(...) function is strictly a function of Box and
        Interval, and does not depend on the current state of the T object.
     return 1 if size(..) is symmetric, in that sender and receiver T object
        can size their message buffers, but a static object cannot.
     returns 2 if the object is truly dynamic.  the message size is subject to
        unique object data. </PRE>
    - <PRE> int size(const Box& R, const Interval& comps) </PRE>
    - <PRE> void linearOut(void* buf, const Box& R, const Interval& comps) </PRE>
    - <PRE> void linearIn(void* buf, const Box& R, const Interval& comps) </PRE>

  In sequential codes, many of these functions will not be called, but in the
  general parallel case, they all will be.

*/

template<class T>
class BoxLayoutData : public LayoutData<T>
{
public:
  ///
  BoxLayoutData();

  virtual ~BoxLayoutData();
  ///
  BoxLayoutData(const BoxLayout& boxes, int comps,
                const DataFactory<T>& factory = DefaultDataFactory<T>());

  ///
  virtual void define(const BoxLayout& boxes, int comps,
                      const DataFactory<T>& factory = DefaultDataFactory<T>());

  ///
  virtual void define(const BoxLayoutData<T>& da,
                      const DataFactory<T>& factory = DefaultDataFactory<T>());

  ///
  /** Define this BoxLayoutData to be the subset of the
      data in da defined by the Interval comps.
  */
  virtual void define(const BoxLayoutData<T>& da, const Interval& comps,
                      const DataFactory<T>& factory = DefaultDataFactory<T>());

  /// overridden and throws an error.
  virtual void define(const BoxLayout& boxes);

  ///
  int nComp() const { return m_comps;}

  ///
  Interval interval() const
  {
    Interval outint(0, m_comps-1);
    return(outint);
  }

  /// General data copying operation.
  /**
     @param a_destGrids BoxLayout for the destination data holder
     @param a_dest return argument.  upon completion contains a Vector of T objects
       for each Box in a_destGrids that overlaps this->boxLayout()
     @param a_interval range of components to source
     @param a_domain ProblemDomain that this grid and a_destGrids reside on.
     @param factory optional data factory for template classes that have non-trivial construction
  */
  void generalCopyTo(const BoxLayout& a_destGrids,
                     LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                     const Interval& a_interval,
                     const ProblemDomain& a_domain,
                     const DataFactory<T>& factory = DefaultDataFactory<T>()) const ;

  ///
  void generalCopyTo(const BoxLayout& a_destGrids,
                     LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                     const Interval& a_interval,
                     const ProblemDomain& a_domain,
                     const Copier& a_copier,
                     const DataFactory<T>& factory = DefaultDataFactory<T>()) const ;


    ///
  /* User writes a function with the signature:

     <PRE>
     void myfunction(const Box& box, int comps, T& t){ your code here;}
     </PRE>

     They can then hand this off to LayoutData::apply.  This class
     then cycles through all the T objects and invokes this function.  Function
     must not be inline. (I'm still trying to figure out a nice way to send
     in non-static member functions).
     */
  virtual void apply(void (*a_Function)(const Box& box, int comps, T& t));

  ///
  virtual bool isDefined() const;


protected:
  int             m_comps;
  bool            m_isdefined;

  friend class LevelData<T>;

  void setVector(const BoxLayoutData<T>& da,
                 const Interval& srcComps,
                 const Interval& destComps);

  void allocateGhostVector(const DataFactory<T>& factory,
                           const IntVect& ghost = IntVect::Zero);

  void makeItSo(const Interval&     a_srcComps,
                const BoxLayoutData<T>& a_src,
                BoxLayoutData<T>&   a_dest,
                const Interval&     a_destComps,
                const Copier&       a_copier,
                const LDOperator<T>& a_op = LDOperator<T>()) const;
  void makeItSoBegin(
                const Interval&     a_srcComps,
                const BoxLayoutData<T>& a_src,
                BoxLayoutData<T>&   a_dest,
                const Interval&     a_destComps,
                const Copier&       a_copier,
                const LDOperator<T>& a_op = LDOperator<T>()) const;
  void makeItSoEnd(
                BoxLayoutData<T>&   a_dest,
                const Interval&     a_destComps,
                const LDOperator<T>& a_op = LDOperator<T>()) const;

  //========================================================================
  //
  // data structures used by makeItSo when we have some
  // data that needs to be moved (ie. there are entries
  // in the 'FROM' or 'TO' CopyIterators)
  //
  void completePendingSends() const;

  void allocateBuffers(const BoxLayoutData<T>& a_src,
                       const Interval& a_srcComps,
                       const BoxLayoutData<T>& a_dest,
                       const Interval& a_destComps,
                       const Copier&   a_copier,
                       const LDOperator<T>& a_op) const;

  void writeSendDataFromMeIntoBuffers(const BoxLayoutData<T>& a_src,
                                      const Interval& a_srcComps,
                                      const LDOperator<T>& a_op) const;

  void postSendsFromMe() const ;

  void postReceivesToMe() const ;

  void unpackReceivesToMe(BoxLayoutData<T>& a_dest,
                          const Interval&   a_destComps,
                          const LDOperator<T>& a_op) const ;

  void unpackReceivesToMe_append(LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                 const Interval&   a_destComps,
                                 int ncomp,
                                 const DataFactory<T>& factory,
                                 const LDOperator<T>& a_op) const;

  /** \name Parallel messaging members */
  /*@{*/
  mutable void*  m_sendbuffer; // pointer member OK here,
                               // since LevelData<T> has no copy
  mutable size_t m_sendcapacity;
  mutable void*  m_recbuffer;  // pointer member OK here,
                               // since LevelData<T> has no copy
  mutable size_t m_reccapacity;

#ifdef CH_MPI

#ifndef DOXYGEN

  struct bufEntry
  {
    void* bufPtr; // pointer into contiguous m_buffer
    size_t size;
    const MotionItem* item;
    unsigned int procID;
    bool operator < (const bufEntry& rhs) const
      {
        if(procID == rhs.procID)
          {
            const Box& left = item->toRegion;
            const Box& right= rhs.item->toRegion;
            if(left.smallEnd() == right.smallEnd())
              {
                return left.bigEnd().lexLT(right.bigEnd());
              }
            else
              {
                return item->toRegion < rhs.item->toRegion;
              }
          }
        //else
        return procID < rhs.procID;
      }
  };

#endif
  mutable std::vector<bufEntry> m_fromMe;
  mutable std::vector<bufEntry> m_toMe;

  mutable Vector<MPI_Request>  m_sendRequests,  m_receiveRequests;
  mutable Vector<MPI_Status>   m_receiveStatus, m_sendStatus;
  mutable int numSends, numReceives;
#endif

};

///
/** not actually L-p norms, since it doesn't take into
  account the dx of the system.  A user can take that
  into account or not.

  For p != 0, returns pth root of sum of pth powers over all
  points in all fabs and all components in the interval:

  ( sum [ |A[i][pt,var]|^p : FArrayBox A[i], point pt in A[i].box(), var in interval ] )^(1/p)

  To turn into an L-p norm, one needs to multiply this by dx^(SpaceDim/p).

  For p == 0, returns global max over all points in all fabs and all
  components in the interval:

  max [ |A[i][pt,var]| : FArrayBox A[i], point pt in A[i].box(), var in interval ]

  Some people don't like that this norm is not normalized based on
  number of points in A.  Normalization is your problem.
  */
Real norm(const BoxLayoutData<FArrayBox>& A,
          const Interval& interval,
          const int& p = 2);

//======================================================================
template < >
BaseFab<int>* DefaultDataFactory<BaseFab<int> >::create(const Box& box,
                                                        int ncomps,
                                                        const DataIndex& a_datInd) const;

template < >
FArrayBox* DefaultDataFactory<FArrayBox>::create(const Box& box,
                                                 int ncomps,
                                                 const DataIndex& a_datInd) const;

template <class T>
T* DefaultDataFactory<T>::create(const Box& box,
                                 int ncomps,
                                 const DataIndex& a_datInd) const
{
  return new T(box, ncomps);
}

template<class T>
inline bool BoxLayoutData<T>::isDefined() const
{
  return m_isdefined;
}

template <class T>
inline void BoxLayoutData<T>::setVector(const BoxLayoutData<T>& da,
                                        const Interval& srcComps,
                                        const Interval& destComps)
{
  if(&da != this)
    {
      for(DataIterator it(this->dataIterator()); it.ok(); ++it)
        {
          this->m_vector[this->m_boxLayout.index(it())]->copy( this->box(it()), destComps,
                                                   this->box(it()), da[it()], srcComps);
        }
    }
}

template<class T>
inline void BoxLayoutData<T>::define(const BoxLayoutData<T>& da, const Interval& comps,
                                     const DataFactory<T>& factory)
{
  if(this == &da){
    MayDay::Error("BoxLayoutData<T>::define(const LayoutData<T>& da,.....) called with 'this'");
  }
 CH_assert(comps.size()>0);
 CH_assert(comps.end()<=m_comps);
 CH_assert(comps.begin()>=0);
  this->m_boxLayout = da.boxLayout();

  this->m_comps = comps.size();

  Interval dest(0, m_comps-1);
  allocateGhostVector(factory);
  setVector(da, comps, dest);
#ifdef CH_MPI
  this->m_fromMe.resize(0);
  this->m_toMe.resize(0);
#endif

}

template<class T>
inline void BoxLayoutData<T>::define(const BoxLayout& boxes, int comps,
                                     const DataFactory<T>& factory)
{
 CH_assert(boxes.isClosed());
  this->m_boxLayout = boxes;
  m_comps = comps;
  m_isdefined = true;
  allocateGhostVector(factory);
#ifdef CH_MPI
  m_fromMe.resize(0);
  m_toMe.resize(0);
#endif
}

template<class T>
inline void BoxLayoutData<T>::define(const BoxLayout& boxes)
{
  MayDay::Error("BoxLayoutData<T>::define(const BoxLayout& boxes)...needs comps");
}

template <class T>
inline   BoxLayoutData<T>::BoxLayoutData():m_comps(0) , m_sendbuffer(NULL),
    m_sendcapacity(0),  m_recbuffer(NULL),m_reccapacity(0)
{
  m_isdefined = false;
#ifdef CH_MPI
  this->numSends = 0;
  this->numReceives = 0;
#endif
}
template<class T>
inline BoxLayoutData<T>::BoxLayoutData(const BoxLayout& boxes, int comps,
                                       const DataFactory<T>& factory)
  :m_comps(comps), m_sendbuffer(NULL),
    m_sendcapacity(0),  m_recbuffer(NULL),m_reccapacity(0)
{
 CH_assert(boxes.isClosed());
  this->m_boxLayout = boxes;
  m_isdefined = true;
  allocateGhostVector(factory);
#ifdef CH_MPI
  m_fromMe.resize(0);
  m_toMe.resize(0);
  this->numSends = 0;
  this->numReceives = 0;
#endif
}

template<class T>
BoxLayoutData<T>::~BoxLayoutData()
{
  completePendingSends();
  free(m_sendbuffer);
  free(m_recbuffer);
}

template<class T>
inline void BoxLayoutData<T>::define(const BoxLayoutData<T>& da,
                                     const DataFactory<T>& factory)
{
  if(this != &da){
    m_isdefined = da.m_isdefined;
    this->m_boxLayout = da.boxLayout();
    m_comps    = da.nComp();
    Interval srcAnddest(0, m_comps-1);
    allocateGhostVector(factory);
    setVector(da, srcAnddest, srcAnddest);
  }
#ifdef CH_MPI
  m_fromMe.resize(0);
  m_toMe.resize(0);
#endif
}

template<class T>
inline void BoxLayoutData<T>::allocateGhostVector(const DataFactory<T>& factory, const IntVect& ghost)
{
  if(this->m_callDelete == true){
  for(unsigned int i=0; i<this->m_vector.size(); ++i)
    {
      delete this->m_vector[i];
      this->m_vector[i] = NULL;
    }
  }

  this->m_callDelete = factory.callDelete();

  this->m_vector.resize(this->m_boxLayout.size(), NULL);

  for(DataIterator it(this->dataIterator()); it.ok(); ++it)
    {
      unsigned int index = this->m_boxLayout.index(it());
      Box abox = this->box(it());
      abox.grow(ghost);
      this->m_vector[index] = factory.create(abox, m_comps, it());
      if(this->m_vector[index] == NULL)
        {
          MayDay::Error("OutOfMemory in BoxLayoutData::allocateGhostVector");
        }
    }
}

template<class T>
inline void BoxLayoutData<T>::apply(void (*a_func)(const Box& box, int comps, T& t))
{
  for(DataIterator it(this->dataIterator()); it.ok(); ++it)
    {
      unsigned int index = this->m_boxLayout.index(it());
      a_func(this->box(it()), m_comps, *(this->m_vector[index]));
    }
}

//======================================================================
template <class T>
AliasDataFactory<T>::AliasDataFactory(BoxLayoutData<T>* a_original, const Interval& interval)
{
  define(a_original, interval);
}

template <class T>
void AliasDataFactory<T>::define(BoxLayoutData<T>* a_original, const Interval& interval)
{
  m_origPointer = a_original;
  m_interval    = interval;
}

template <class T>
T* AliasDataFactory<T>::create(const Box& a_box, int ncomps, const DataIndex& a_dataInd) const
{
  //CH_assert(this->box(a_dataInd) == a_box);
 CH_assert(ncomps = m_interval.size());
  T* rtn = new T(m_interval, m_origPointer->operator[](a_dataInd));
  return rtn;
}

#include "NamespaceFooter.H"
#include "BoxLayoutDataI.H"

#endif //BOXLAYOUTDATA
