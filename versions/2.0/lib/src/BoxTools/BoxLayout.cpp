#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <iostream>
#include <algorithm>
#include "BoxLayout.H"
#include "DataIterator.H"
#include "LayoutIterator.H"
#include "NamespaceHeader.H"

using std::ostream;

//need at least one non-inlined function, otherwise
//some compilers don't build a class description in the
//object file.
Box BoxLayout::get(const LayoutIterator& it) const
{
  return get(it());
}
Box BoxLayout::get(const DataIterator& it) const
{
  return get(it());
}

const Box&
BoxLayout::operator[](const LayoutIterator& index) const
{
  return this->operator[](index());
}

const Box&
BoxLayout::operator[](const DataIterator& index) const
{
  return this->operator[](index());
}
BoxLayout::~BoxLayout()
{}

BoxLayout::BoxLayout()
  :m_boxes(new Vector<Entry>()),
   m_index(new Vector<unsigned int>()),

   m_layout(new int),
   m_closed(new bool(false)),
   m_dataIterator(RefCountedPtr<DataIterator>())
{
}

BoxLayout& BoxLayout::operator=(const BoxLayout& a_rhs)
{
  if(this == &a_rhs) return *this;
  m_boxes = a_rhs.m_boxes;
  m_index = a_rhs.m_index;
  m_layout = a_rhs.m_layout;
  m_closed = a_rhs.m_closed;
  m_dataIterator = a_rhs.m_dataIterator;
  return *this;
}

void BoxLayout::sort()
{
  if(!*m_closed)
    {
      //      std::sort(m_boxes->begin(), m_boxes->end());
      m_boxes->sort();
      setIndexVector();
    }
}

void BoxLayout::closeNoSort()
{
  if(!*m_closed){
    //sort();  there, no sort.
    *m_closed = true;
    m_dataIterator = RefCountedPtr<DataIterator>(new DataIterator(*this, m_layout));
  }
}
void BoxLayout::close()
{
  //  CH_assert(m_index->size() != 0);
  if(!*m_closed){
    sort();
    *m_closed = true;
    m_dataIterator = RefCountedPtr<DataIterator>(new DataIterator(*this, m_layout));
  }
}
bool BoxLayout::coarsenable(int refRatio) const
{
  for(int i=0; i<size(); i++)
    {
      Box b =  m_boxes->operator[](m_index->operator[](i)).box;
      b.coarsen(refRatio);
      b.refine(refRatio);
      if(b !=  m_boxes->operator[](m_index->operator[](i)).box)
        return false;
    }
  return true;
}

void BoxLayout::setIndexVector()
{
  Vector<Entry>& boxes = *m_boxes;
  Vector<unsigned int>& index = *m_index;
  index.resize(boxes.size());
  for(unsigned int i=0; i<boxes.size(); ++i)
    {
      index[boxes[i].index] = i;
    }
}

// Constructors and such
// =====================

DataIterator BoxLayout::dataIterator() const
{
  CH_assert(*m_closed);
  //DataIterator rtn(*this, m_layout);
  return *m_dataIterator;
  //return rtn;
}


LayoutIterator BoxLayout::layoutIterator() const
{
  return LayoutIterator(*this, m_layout);
}

BoxLayout::BoxLayout(const Vector<Box>& a_boxes, const Vector<int>& assignments)
  :m_boxes(new Vector<Entry>()),
   m_index(new Vector<unsigned int>()),
   m_layout(new int),
   m_closed(new bool(false))
{
  define(a_boxes, assignments);
}

BoxLayout::BoxLayout(const Vector<Box>& a_boxes, 
                     const Vector<int>& assignments,
                     const Vector<int>& a_blockIDs)
  :m_boxes(new Vector<Entry>()),
   m_index(new Vector<unsigned int>()),
   m_layout(new int),
   m_closed(new bool(false))
{
  define(a_boxes, assignments, a_blockIDs);
}

void BoxLayout::checkDefine(const Vector<Box>& a_boxes, const Vector<int>& a_procIDs)
{

  if(*m_closed)
    {
      MayDay::Error("attempt to define(..) a closed BoxLayout");
    }
  const int num_boxes = a_boxes.size();
  const int num_procs = a_procIDs.size();
  if( num_procs > 0 && num_boxes != num_procs )
    {
      MayDay::Error("BoxLayout::define(): vector of processor assignments is different length from vector of boxes");
    }
  // Check for negative proc ID's and ID's larger than total number of procs.
  for (unsigned int i = 0; i < num_procs; ++i)
    {
      if (a_procIDs[i] < 0)
        {
          MayDay::Error("BoxLayout::define(): Negative processor assignments not allowed");
        }
   //    if (a_procIDs[i] >= numProc())
//         {
//           MayDay::Error("BoxLayout::define(): Attempting to assign data to processor ID larger than total number of processors available");
//         }
    }
}

void
BoxLayout::define(const Vector<Box>& a_boxes, const Vector<int>& a_procIDs)
{
  checkDefine(a_boxes, a_procIDs);
  const int num_boxes = a_boxes.size();
  const int num_procs = a_procIDs.size();
  m_boxes->resize(num_boxes);
  for (unsigned int i = 0; i < num_boxes; ++i)
    {
      m_boxes->operator[](i) = a_boxes[i];
      m_boxes->operator[](i).index = i;
      if( num_procs > 0 ) m_boxes->operator[](i).m_procID = a_procIDs[i];
    }
  m_index->resize(1); //signal for close routine.
  close();
}

void
BoxLayout::define(const Vector<Box>& a_boxes, 
                  const Vector<int>& a_procIDs,
                  const Vector<int>& a_blockIDs)
{
  checkDefine(a_boxes, a_procIDs);
  const int num_boxes = a_boxes.size();
  const int num_procs = a_procIDs.size();
  m_boxes->resize(num_boxes);
  for (unsigned int i = 0; i < num_boxes; ++i)
    {
      m_boxes->operator[](i) = a_boxes[i];
      m_boxes->operator[](i).index = i;
      if( num_procs > 0 ) m_boxes->operator[](i).m_procID = a_procIDs[i];
      m_boxes->operator[](i).m_blockID = a_blockIDs[i];
    }
  m_index->resize(1); //signal for close routine.
  close();
}

// Other member functions
// ======================

//DataIndex
void BoxLayout::addBox(const Box& box, int procID)
{
  if(*m_closed)
    {
      MayDay::Error("attempt to addBox to closed BoxLayout");
    }
  if(m_layout.isNonUnique())
    {
      m_layout = RefCountedPtr<int>(new int);
    }
  Vector<Entry>& boxes = *m_boxes;
  unsigned int i = boxes.size();
  boxes.push_back(box);
  boxes[i].index = i;
  boxes[i].m_procID = procID;
  (*m_index).push_back(i);
  //return DataIndex(i, m_layout);
}

void BoxLayout::aliasAddBox(const Box& box)
{
  if(*m_closed)
    {
      MayDay::Error("attempt to addBox to closed BoxLayout");
    }
  if(m_layout.isNonUnique())
    {
      m_layout = RefCountedPtr<int>(new int);
    }
  m_boxes->push_back(box);
  CH_assert(m_index->size() == 0);
}

void
BoxLayout::deepCopy(const BoxLayout& a_source)
{
  m_boxes =  RefCountedPtr<Vector<Entry> >(
                new Vector<Entry>(*(a_source.m_boxes)));
  m_index =  RefCountedPtr<Vector<unsigned int> >(
                new Vector<unsigned int>(*(a_source.m_index)));
  m_layout = a_source.m_layout;
  *m_closed = false;
}

/*
void
BoxLayout::resize(unsigned int n, const Box& box)
{
  if(*m_closed)
    {
      MayDay::Error("attempt to resize closed BoxLayout");
    }
  if(n == m_boxes->size()) return;
  if(m_layout.isNonUnique())
    {
      m_layout = new int;
    }
  unsigned int size = m_boxes->size();
  m_boxes->resize(n, box);
  m_index->resize(n);
  if(n > size)
    {
      for(unsigned int i=size; i<n; ++i)
        m_boxes->operator[](i).index = i;
    }
  setIndexVector();
}

*/
void BoxLayout::computeNeighbors(){;}

void BoxLayout::aliasClose()
{

  CH_assert(m_index->size()==0); //check to make sure the other addBox function was never called

  Vector<Vector<Box> > allBoxes;
  Vector<Box> vbox(m_boxes->size());

  for(int i=0; i<m_boxes->size(); i++)
    {
      vbox[i] = m_boxes->operator[](i).box;
    }

  gather(allBoxes, vbox, 0);
  broadcast(allBoxes, 0);
  vbox.resize(0);
  Vector<int> procID;
  for(int p=0; p<allBoxes.size(); p++)
    {
      Vector<Box>& b = allBoxes[p];
      for(int i=0; i<b.size(); i++)
        {
          procID.push_back(p);
          vbox.push_back(b[i]);
        }
    }

  m_boxes->resize(vbox.size());
  for (unsigned int i = 0; i < vbox.size(); ++i)
    {
      m_boxes->operator[](i) = vbox[i];
      m_boxes->operator[](i).index = i;
      m_boxes->operator[](i).m_procID = procID[i];
    }
  setIndexVector();
  *m_closed = true;
  m_dataIterator = RefCountedPtr<DataIterator>(
                    new DataIterator(*this, m_layout));

}

// Global functions
// ================

// For now, we can just have the one coarsen funtion.  If a DisjointBoxLayout
// enters this function, is coarsened, and then doesn't remain disjoint, it
// will be caught here at the call to close().  Debugging should not be

void
coarsen(BoxLayout& a_output, const BoxLayout& a_input, int a_refinement)
{
   if(!a_input.isClosed())
    {
      MayDay::Error("input to coarsen must be called with closed BoxLayout");
    }
  if(a_output.isClosed())
    {
      MayDay::Error("output of coarsen must be called on open BoxLayout");
    }
  a_output.deepCopy(a_input);

  for(LayoutIterator it(a_input.layoutIterator()); it.ok(); ++it)
    {
      a_output.ref(it()).coarsen(a_refinement);
    }
  a_output.close();
}

// we have an easier time with refine, since we know that refinement will
// not change the state of a sort, but, we will play it safe for now
// until this function shows up in the profiler.

void refine(BoxLayout& a_output, const BoxLayout& a_input, int a_refinement)
{
  if(!a_input.isClosed())
    {
      MayDay::Error("input to refine must be called with closed BoxLayout");
    }
  if(a_output.isClosed())
    {
      MayDay::Error("output of refine must be called on open BoxLayout");
    }
  a_output.deepCopy(a_input);

  for(LayoutIterator it(a_input.layoutIterator()); it.ok(); ++it)
    {
      a_output.ref(it()).refine(a_refinement);
    }
  a_output.close();
}

ostream& operator<<(ostream& os, const BoxLayout& a_layout)
{
  int i=0;
  for(LayoutIterator it(a_layout.layoutIterator()); it.ok(); ++it)
    {
      os << a_layout.get(it())<<"["<<a_layout.procID(it())<<"]";
      ++i;
      if(i==4){ os <<"\n"; i=0;}
      else {os <<" # ";}
    }

  os <<"\n";
  return os;
}

void BoxLayout::print() const { std::cout << *this;}

int BoxLayout::numBoxes(const int procID) const
{
  int num = 0;
  for(int i=0; i<m_boxes->size(); ++i)
    {
      if(m_boxes->operator[](i).m_procID == procID) ++num;
    }
  return num;
}

Vector<Box> BoxLayout::boxArray() const
{
  Vector<Box> result( m_boxes->size() );

  // FIXME: maybe we want the order to be the one LayoutIterator would take?
  for( int i=0;i<m_boxes->size();++i )
  {
    result[i] = (*m_boxes)[i].box;
  }
  return result;
}

Vector<int> BoxLayout::procIDs() const
{
  Vector<int> result( m_boxes->size() );
  for( int i=0;i<m_boxes->size();++i )
  {
    result[i] = (*m_boxes)[i].m_procID;
  }
  return result;
}

class MortonOrdering
{
public:
  MortonOrdering(int a_maxSize):maxSize(a_maxSize){;}
  MortonOrdering():maxSize(8*sizeof(int)-2){;}
  inline bool operator()(const Box& lhs, const Box& rhs) const;
  int maxSize;
};

inline bool MortonOrdering::operator()(const Box& lhs, const Box& rhs) const
{
  const IntVect l = lhs.smallEnd();
  const IntVect r = rhs.smallEnd();
  for(int i = maxSize; i>0; i--)
    {
      const int N = (1<<i); // march from most significant bit to least.
      for(int dir=CH_SPACEDIM-1; dir>=0; dir--)
        {
          if      ((l[dir]/N) < (r[dir]/N)) return true;
          else if ((l[dir]/N) > (r[dir]/N)) return false;
        }
    }
  return false ;
}

void mortonOrdering(Vector<Box>& a_boxes)
{
  int maxSize=0;
  std::vector<Box>& b = a_boxes.stdVector();
  for(std::vector<Box>::iterator p=b.begin(); p<b.end(); ++p)
    {
      IntVect small = p->smallEnd();
      D_EXPR( maxSize = Max(maxSize, Abs(small[0])),
              maxSize = Max(maxSize, Abs(small[1])),
              maxSize = Max(maxSize, Abs(small[2])));
    }
  int bits;
  for(bits=8*sizeof(int)-2; bits>0; bits--)
    {
      const int N = (1<<bits);
      int rem = maxSize/N;
      if (rem > 0) break;
    }
  bits++;
  std::sort(b.begin(), b.end(), MortonOrdering(bits));

  
}

#include "NamespaceFooter.H"
