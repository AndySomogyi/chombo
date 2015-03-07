#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif


#ifndef _BOXLAYOUT_H_
#define _BOXLAYOUT_H_

#include "Box.H"
#include "Vector.H"
#include "RefCountedPtr.H"
#include "DataIndex.H"
#include "SPMD.H"
#include "NamespaceHeader.H"

class DataIterator;
class TimedDataIterator;
class LayoutIterator;

#ifndef DOXYGEN

struct Entry
{
  Entry()
    :
    m_procID(procID())
  {}

  Entry(const Box& a_box)
    :
    box(a_box),
    m_procID(procID())
  {}

  Entry(const Box& a_box, const unsigned int a_index)
    :
    box(a_box),
    index(a_index),
    m_procID(procID())
  {}

  bool operator < (const Entry& rhs) const
  {
    return box < rhs.box;
  }

  Box box;
  unsigned int index;
  unsigned int m_procID;// not used in serial code.
};

#endif

///A not-necessarily-disjoint collective of boxes.
/**
  A BoxLayout is a collection of Box objects that are assigned
  to process numbers.  Each box is associated with only one
  process.  Processes are numbered from 0 to n-1 (for a job with
  n processes).

  A BoxLayout can be either open or closed.

  Open BoxLayout:
  - Created by null construction or deepCopy.
  - Boxes may be added to it.
  - Non-const operations may be performed on the boxes in it.

  Closed BoxLayout:
  - Created by constructor with vectors of Boxes and processors given explicitly.
  - Cannot be modified.
  - Represented as sorted boxes.
  - Many uses of BoxLayouts require a closed BoxLayout.

  <b>Ref-counting</b>

  BoxLayout is an explicitly ref-counted object.

  Assignment and copy are compiler-generated.  They increment the refcount
  on the contained data members.  They perform shallow, ref-counted
  operations.

  Refcounting is a process whereby multiple instantiations make
  use of a single implementation of that object and keep a tally
  of how many instantiations are sharing.  Thus:
  <PRE>

  BoxLayout b1(boxes, procIDs);

            b1  ----> refcount = 1
                ---->  m_boxes
                ---->  m_processors

  BoxLayout b2(b1)

            b1  ----> refcount = 2  <---- b2
                ---->  m_boxes      <----
                ---->  m_processors <----

  BoxLayout b3;
  b3 = b2;

            b1  ----> refcount = 3  <---- b2
                ---->  m_boxes      <----
                ---->  m_processors <----
                        ^^^
                        |||
                         b3
  </PRE>
*/
class BoxLayout
{
public:

  /**
     \name Constructors, destructors, assignments, defines
  */
  /*@{*/

  ///
  /**
     Construct BoxLayout with no boxes.
   */
  BoxLayout();

  ///
  /** Construct from a Vector of Boxes and a Vector of
      processor assignments.  On exit, the BoxLayout will be closed.
  */
  BoxLayout(const Vector<Box>& a_boxes,
            const Vector<int>& a_procIDs);


  ///
  /**
    Ref-counted destruction.  Once the last reference to the
    implementation of this class is destroyed, the data members
    are cleaned up
  */
  virtual
  ~BoxLayout();

  ///
  /**
     Ref-counted assignment.
   */
  BoxLayout& operator=(const BoxLayout& a_rhs);

  ///
  /** Define this BoxLayout from a Vector of Boxes and a
      Vector of processor assignments.  Any pre-existing layout will
      be lost (ref-count dropped by one).  The processor assignment Vector
      must be the same length
      as the Box Vector.  On exit, the BoxLayout will be closed.
  */
  virtual void
  define(const Vector<Box>& a_boxes,
         const Vector<int>& a_procIDs);



  /*@}*/

  /**
     \name Accessors
  */
  /*@{*/

  ///
  /** const accessor operator.  See also get(const LayoutIndex&).
   */
  const Box&
  operator[](const LayoutIndex& it) const;

   ///
   /** accessor operator.  See also get(const LayoutIndex&).
       not overloaded with operator[] on purpose.
    */
  Box& ref(const LayoutIndex& it);

  ///
  /** accessor operator.  See also get(const LayoutIndex&).
   */
  const Box&
  operator[](const LayoutIterator& it) const;

  ///
  /** accessor operator.  See also get(const LayoutIndex&).
   */
  const Box&
  operator[](const DataIterator& it) const;

  ///
  /** Get box indexed by <i>it</i>.

      As a consequence of the C++ compiler being free to choose which
      version of operator[] when the object is technically non-const, we very
      often get 'BoxLayout closed' errors.
      This is a non-overloaded get method.
  */
  Box get(const LayoutIndex& it) const;


  ///
  /** Get box indexed by <i>it</i>. equivalent to  get(it()), just does the extra() for you.
   */
  Box get(const DataIterator& it) const;

  ///
  /** Get box indexed by <i>it</i>. equivalent to  get(it()), just does the extra() for you.
   */
  Box get(const LayoutIterator& it) const;

  ///
  /** Returns the processor to which this box has been assigned.
      Not a user function, at least, not a new user function.  It can
      be used safely at anytime, closed or open.  A person needing this
      level of knowledge of the processor assignment should have non-trivial
      needs, like writing your own load balancer or such.
      Most user-level parallel

  */
  unsigned int
  procID(const LayoutIndex& a_index) const ;





  /// Returns the number of boxes assigned to a given procID.
  /** Returns the number of boxes assigned to a given procID.
   */
  int numBoxes(const int procID) const;


  /// Return number of cells in all boxes of entire box layout
  long long numCells() const;

  /// Returns the total number of boxes in the BoxLayout.
  /** Returns the total number of boxes in the BoxLayout.
   */
  unsigned int
  size() const ;

  /** Not a user function.  Used in I/O routine.
   */
  unsigned int index(const LayoutIndex& index) const;

  unsigned int lindex(const DataIndex& index) const;
  /*@}*/

  /**
     \name Checks
  */
  /*@{*/

  ///
  /** Refcounted pointer check.  Return <tt>true</tt> if these two objects
      share the same implementation.
    */
  inline bool
  operator==(const BoxLayout& rhs) const ;

  /**
     RefCounted pointer compare.  Need this to use DisjointBoxLayout as a Key in a std::map
  */
  inline bool operator<(const BoxLayout& rhs) const ;

  /**
      current ref count for this object
  */
  int refCount() const
  {
    return m_boxes.refCount();
  }
  ///
  /** Refcounted pointer check.  Return <tt>true</tt> if these two objects
      share the same implementation.
    */
  bool eq(const BoxLayout& rhs) const { return *this == rhs;}

  /** Return <tt>true</tt> if close() has been called.
      Closed BoxLayout is always sorted.
   */
  bool
  isClosed() const;

  ///
  /** Return <tt>true</tt> if sort() has been called.
   */
  bool
  isSorted() const {return isClosed();}

  /** not a user function
   */
  bool check(const LayoutIndex& index) const
    { return index.m_layoutIntPtr == m_layout;}


  ///
  /**
     returns 'true' if you can use the same LayoutIterator and DataIterator on these
     two layouts and data holders built on top of these layouts
  */
  bool compatible(const BoxLayout& a_rhs) const
  {
    return m_layout == a_rhs.m_layout;
  }

  /*@}*/

  /**
     \name Modification functions
  */
  /*@{*/

  ///
  /** Add a box to this BoxLayout. Takes the processor
      assignment. The input <i>box</i> is copied.  BoxLayout must be "open"
      for this operation, or else a runtime error will occur.
      */
 //  DataIndex
//   addBox(const Box& box, int procID);
  void addBox(const Box& box, int procID);

  ///
  /** Mark this BoxLayout as complete and unchangeable.
   */
  virtual void
  close();

  // not for general consumption.
  virtual void
  closeNoSort();


  ///
  /** Actual deep copy operation.  New object created with copied
      data.  This object disassociates itself with original implementation
      safely.  This object now is considered 'open' and can be non-const
      modified.  There is no assurance that the order in which this BoxLayout
      is indexed corresponds to the indexing of <i>a_source</i>.
<PRE>
       BoxLayout b1(boxes, procIDs);

            b1  ----> refcount = 1
                ---->  m_boxes
                ---->  m_processors

       BoxLayout b2(b1)

            b1  ----> refcount = 2  <---- b2
                ---->  m_boxes      <----
                ---->  m_processors <----

       BoxLayout b3;
       b3.deepCopy(b2);
            b1  ----> refcount = 2  <---- b2  b3 ----> refcount = 1
                ---->  m_boxes      <----        ---->  m_boxes
                ---->  m_processors <----        ---->  m_processors
</PRE>
  */
  virtual void
  deepCopy(const BoxLayout& a_source);

  ///
  /**
     returns true iff:
     - every Box in the BoxLayout can be coarsened by refRatio and return back
     to the original Box when refined by refRatio.
     - refRatio must be a positive non-zero integer.
  */
  bool coarsenable(int refRatio) const;

  ///
  /**
     Coarsen a BoxLayout:
     - <i>output</i> must be open
     - <i>input</i> must be closed
     - <i>refinement</i> must be a positive non-zero integer
     - <i>output</i> and <i>input</i> do not share an
     implementation.

     <i>output</i> is first deepCopy'ed from <i>input</i>,
     then coarsen(refinement) is called on each box of <i>output</i>.

     <i>output</i> returns from this function closed.

     LayoutIterators and DataIterators from <i>input</i> and <i>output</i>
     can be used interchangeably.
  */
  friend void coarsen(BoxLayout& output,
                      const BoxLayout& input,
                      int refinement);

  ///
  /**
     Refine a BoxLayout:
     - <i>output</i> must be open
     - <i>input</i> must be closed
     - <i>refinement</i> must be a positive non-zero integer
     - <i>output</i> and <i>input</i> do not share an
     implementation.

     <i>output</i> is first deepCopy'ed from <i>input</i>,
     then refine(refinement) is called on each box of <i>output</i>.

     <i>output</i> returns from this function closed.

     LayoutIterators and DataIterators from <i>input</i> and <i>output</i>
     can be used interchangeably.
  */
  friend void refine(BoxLayout& output,
                     const BoxLayout& input,
                     int refinement);

  ///
  /** Assign a Box in the BoxLayout to a processor.
      Requires the BoxLayout to be open.
  */
  void
  setProcID(const LayoutIndex& a_index, unsigned int a_procID);

  //const or non-const operation ?.....I can think of usages either way...bvs
  ///
  /**
     Sort the boxes of an open BoxLayout.
     No change if the BoxLayout is closed.
  */
  void
  sort();

  /** Not a user function, used in Chombo interface to AMRData.H
      interoperability tools. */
  void aliasAddBox(const Box& box);

  /** Not a user function, used in Chombo interface to AMRData.H
      interoperability tools. */
  void aliasClose();

  /*@}*/

  /**
     \name Iterators
  */
  /*@{*/

  ///Parallel iterator.
  /** Parallel iterator.

      Returns DataIndex object that corresponds to boxes
      on the local processor only. */
  DataIterator
  dataIterator() const;

  ///
  TimedDataIterator
  timedDataIterator() const;

  ///Iterator that processes through ALL the boxes in a BoxLayout.
  /** Iterator that processes through ALL the boxes in a BoxLayout.

      If BoxLayout is closed, then LayoutIterator will return
      the Boxes in sorted order. */
  LayoutIterator
  layoutIterator() const;

  /*@}*/

  /**
     \name I/O functions
  */
  /*@{*/

#ifndef WRAPPER
  ///
  /**
     Invokes cout<<*this; pretty-print dump of BoxLayout.
   */
  void
  print() const;
#endif
  ///
  /**
     Invokes cout<<*this; pretty-print dump of BoxLayout.
   */
  void p() const { print();}

  /*@}*/

  /* Return all the constituent boxes. */
  Vector<Box> boxArray() const;

  /** Return the processor id numbers corresponding to the boxes as returned by
   *  this->boxArray().
  */
  Vector<int> procIDs() const;

  inline unsigned int indexI(const LayoutIndex&) const;

protected:

  void setIndexVector();
  void buildDataIndex();
  friend class LayoutIterator;
  friend class DataIterator;

  RefCountedPtr<Vector<Entry> >        m_boxes;
  RefCountedPtr<Vector<unsigned int> > m_index;
  RefCountedPtr<Vector<unsigned int> > m_localIndex;
  RefCountedPtr<int>                   m_layout;
  RefCountedPtr<bool>                  m_closed;
  RefCountedPtr<DataIterator>          m_dataIterator;

#ifdef CH_MPI
  RefCountedPtr<Vector<DataIndex> >    m_dataIndex;
#endif

  void checkDefine(const Vector<Box>& a_boxes, const Vector<int>& procIDs);
private:

};

void coarsen_bl(BoxLayout& output, const BoxLayout& input, int refinement);

void refine_bl(BoxLayout& output, const BoxLayout& input, int refinement);

#ifndef WRAPPER

inline
void coarsen_bl(BoxLayout& output, const BoxLayout& input, int refinement)
{ coarsen(output, input, refinement);}

inline
void refine_bl(BoxLayout& output, const BoxLayout& input, int refinement)
{ refine(output, input, refinement);}

void mortonOrdering(Vector<Box>& a_boxes);

//========================================================

// inlined functions
// =================

// constructors
// ============

/*
inline
BoxLayout::BoxLayout(int a_boxes)
  :m_boxes(new std::vector<Entry>(a_boxes)),
   m_index(new std::vector<unsigned int>(a_boxes)),
   m_layout(new int),
   m_closed(new bool(false))
{
  std::vector<Entry>& boxes = *m_boxes;
  std::vector<unsigned int>& index = *m_index;
  for(unsigned int i=0; i<boxes.size(); ++i)
    {
      boxes[i].index = i;
      index[i] = i;
    }
}
*/

// operators
// =========

inline unsigned int BoxLayout::indexI(const LayoutIndex& index) const
{
  return m_index->operator[](index.m_index);
}

inline unsigned int BoxLayout::lindex(const DataIndex& index) const
{
  return m_localIndex->operator[](index.m_index);
}

inline const Box&
BoxLayout::operator[](const LayoutIndex& index) const
{
  CH_assert(check(index));// make sure this LayoutIndex came from my own iterator
  return m_boxes->operator[](m_index->operator[](index.m_index)).box;
}


inline Box&
BoxLayout::ref(const LayoutIndex& index)
{
  if(*m_closed)
    {
      MayDay::Error("attempt to modify closed BoxLayout");
    }
  CH_assert(check(index)); // make sure this LayoutIndex came from my own iterator
  return m_boxes->operator[](m_index->operator[](index.m_index)).box;
}

inline bool
BoxLayout::operator==(const BoxLayout& rhs) const
{
  return m_boxes == rhs.m_boxes;
}

inline bool
BoxLayout::operator<(const BoxLayout& rhs) const
{
  return m_boxes < rhs.m_boxes;
}

// member functions
// ================

inline  Box
BoxLayout::get(const LayoutIndex& index) const
{
  CH_assert(check(index)); // make sure this LayoutIndex came from my own iterator
  return m_boxes->operator[](m_index->operator[](index.m_index)).box;
}



inline unsigned int
BoxLayout::index(const LayoutIndex& a_index) const
{
  return m_index->operator[](a_index.m_index);
}

inline bool
BoxLayout::isClosed() const
{
  return *m_closed;
}

inline unsigned int
BoxLayout::procID(const LayoutIndex& a_index) const
{
  CH_assert(check(a_index));
  return m_boxes->operator[](m_index->operator[](a_index.m_index)).m_procID;
}


inline void
BoxLayout::setProcID(const LayoutIndex& a_index, unsigned int a_procID)
{
  CH_assert(check(a_index));
  m_boxes->operator[](m_index->operator[](a_index.m_index)).m_procID = a_procID;
}

// global functions
// ================

std::ostream& operator<<(std::ostream& os, const BoxLayout& a_layout);

#endif /*WRAPPER*/

#include "NamespaceFooter.H"
#endif
