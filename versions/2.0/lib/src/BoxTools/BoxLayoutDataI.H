#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _BOXLAYOUTDATAI_H_
#define _BOXLAYOUTDATAI_H_

#include <cstdlib>
#include <algorithm>
#include <limits.h>
#include <list>

#include "parstream.H"
#include "memtrack.H"
#include "Misc.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

using std::sort;


template<class T>
void BoxLayoutData<T>::makeItSo(const Interval&   a_srcComps,
                                const BoxLayoutData<T>& a_src,
                                BoxLayoutData<T>& a_dest,
                                const Interval&   a_destComps,
                                const Copier&     a_copier,
                                const LDOperator<T>& a_op) const
{
  makeItSoBegin(a_srcComps, a_src, a_dest, a_destComps, a_copier, a_op);
  makeItSoEnd(a_dest, a_destComps, a_op);
}

template<class T>
void BoxLayoutData<T>::makeItSoBegin(const Interval&   a_srcComps,
                            const BoxLayoutData<T>& a_src,
                            BoxLayoutData<T>& a_dest,
                            const Interval&   a_destComps,
                            const Copier&     a_copier,
                            const LDOperator<T>& a_op) const
{
  // The following five functions are nullOps in uniprocessor mode

  // Instead of doing this here, do it an end of makeItSo (ndk)
  //completePendingSends(); // wait for sends from possible previous operation

  // new evil logic to determine under what conditions we can just
  // re-use our messaging pattern and buffers from the last call to
  // makeItSo.  pretty elaborate, I know.  bvs
#ifdef CH_MPI
  static Copier* lastCopier=NULL;

  if(T::preAllocatable() == 2 || !a_copier.bufferAllocated() ||
     (m_fromMe.size() + m_toMe.size() == 0) ||lastCopier != &a_copier){
    allocateBuffers(a_src,  a_srcComps,
                    a_dest, a_destComps,
                    a_copier,
                    a_op);  //monkey with buffers, set up 'fromMe' and 'toMe' queues
    a_copier.setBufferAllocated(true);
  }
  //lastCopier = (Copier*)(&a_copier);

#endif

  writeSendDataFromMeIntoBuffers(a_src, a_srcComps, a_op);

  // If there is nothing to recv/send, don't go into these functions
  // and allocate memory that will not be freed later.  (ndk)
  // The #ifdef CH_MPI is for the m_toMe and m_fromMe
#ifdef CH_MPI
  this->numReceives = m_toMe.size();
  if (this->numReceives > 0) {
    postReceivesToMe(); // all non-blocking
  }

  this->numSends = m_fromMe.size();
  if (this->numSends > 0) {
    postSendsFromMe();  // all non-blocking
  }
#endif

  //  computation that could occur during communication should really
  //  go here somehow.  while all the non-blocking sending and receiving is
  //  going on.
  //
  //  my thought is to return from this function at this point an object
  //  that encapsulates the argument list above.
  //  a "ChomboMessaging" object.
  //  The user can keep a reference
  //  to this object and do computations.  When they reach the limit of what
  //  they can compute without this communication completing, they call the
  //  "finalize()" function of their ChomboMessaging object and the rest of this
  //  code below gets executed.
  //  a real question though is: is there really enough computation to do while
  //  messaging is going on to justify the effort, and what machines really have
  //  good asynchronous messaging to make the work worthwhile.
  //
  //  the other approach is to more finely decompose the overlapping of
  //  messaging and computation by using the ChomboMessaging object in the
  //  DataIterator construction.  The DataIterator returns T objects as they
  //  are completed from messaging.  This preserves almost all of the Chombo
  //  code as is but would be mucho tricky to actually implement and might only
  //  gain little.  This would not be a thing to try unitl Chombo is
  //  heavily instrumented for performance measuring.  in this design, unpackRecievesToMe()
  //  would become a complicated process interwoven with a DataIterator.


  // perform local copy
  {
    CH_TIME("local copying");
    for(CopyIterator it(a_copier, CopyIterator::LOCAL); it.ok(); ++it)
      {
        const MotionItem& item = it();
        a_op.op(a_dest[item.toIndex], item.fromRegion,
                a_destComps,
                item.toRegion,
                a_src[item.fromIndex],
                a_srcComps);
        
      }
  }
}

template<class T>
void BoxLayoutData<T>::makeItSoEnd(BoxLayoutData<T>& a_dest,
                                   const Interval&   a_destComps,
                                   const LDOperator<T>& a_op) const
{
  // Uncomment and Move this out of unpackReceivesToMe()  (ndk)
  completePendingSends(); // wait for sends from possible previous operation
  
  unpackReceivesToMe(a_dest, a_destComps, a_op); // nullOp in uniprocessor mode
  
}


#ifndef CH_MPI
// uniprocessor version of all these nullop functions.
template<class T>
void BoxLayoutData<T>::completePendingSends() const
{}

template<class T>
void BoxLayoutData<T>::allocateBuffers(const BoxLayoutData<T>& a_src,
                                   const Interval& a_srcComps,
                                   const BoxLayoutData<T>& a_dest,
                                   const Interval& a_destComps,
                                   const Copier&   a_copier,
                                   const LDOperator<T>& a_op
                                   ) const
{}

template<class T>
void BoxLayoutData<T>::writeSendDataFromMeIntoBuffers(const BoxLayoutData<T>& a_src,
                                                  const Interval&     a_srcComps,
                                                  const LDOperator<T>& a_op) const
{}

template<class T>
void BoxLayoutData<T>::postSendsFromMe() const
{}

template<class T>
void BoxLayoutData<T>::postReceivesToMe() const
{}

template<class T>
void BoxLayoutData<T>::unpackReceivesToMe(BoxLayoutData<T>& a_dest,
                                      const Interval&   a_destComps,
                                      const LDOperator<T>& a_op) const
{}

template<class T>
void BoxLayoutData<T>::unpackReceivesToMe_append(LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                                 const Interval&   a_destComps,
                                                 int ncomp,
                                                 const DataFactory<T>& factory,
                                                 const LDOperator<T>& a_op) const
{}

#else

// MPI versions of the above codes.

template<class T>
void BoxLayoutData<T>::completePendingSends() const
{
  if(this->numSends > 0){
    m_sendStatus.resize(this->numSends);
    int result = MPI_Waitall(this->numSends, &(m_sendRequests[0]), &(m_sendStatus[0]));
    if(result != MPI_SUCCESS)
      {
        //hell if I know what to do about failed messaging here
      }


  }
  this->numSends = 0;
}


template<class T>
void BoxLayoutData<T>::allocateBuffers(const BoxLayoutData<T>& a_src,
                                   const Interval& a_srcComps,
                                   const BoxLayoutData<T>& a_dest,
                                   const Interval& a_destComps,
                                   const Copier&   a_copier,
                                   const LDOperator<T>& a_op) const
{
  CH_TIME("MPI_allocateBuffers");
  m_fromMe.resize(0);
  m_toMe.resize(0);
  size_t sendBufferSize = 0;
  size_t recBufferSize  = 0;
  // two versions of code here.  one for preAllocatable T, one not.

  T dummy;
  for(CopyIterator it(a_copier, CopyIterator::FROM); it.ok(); ++it)
    {
      const MotionItem& item = it();
      bufEntry b;
      b.item = &item;
      b.size = a_op.size(a_src[item.fromIndex], item.fromRegion, a_srcComps);
      sendBufferSize+=b.size;
      b.procID = item.procID;
      m_fromMe.push_back(b);
    }
  sort(m_fromMe.begin(), m_fromMe.end());
  for(CopyIterator it(a_copier, CopyIterator::TO); it.ok(); ++it)
    {
      const MotionItem& item = it();
      bufEntry b;
      b.item = &item;
      if(T::preAllocatable() == 0)
        {
          b.size = a_op.size(dummy, item.toRegion, a_destComps);
          recBufferSize+=b.size;
        }
      else if (T::preAllocatable() == 1)
        {
          b.size = a_op.size(a_dest[item.toIndex], item.toRegion, a_destComps);
          recBufferSize+=b.size;
        }
      b.procID = item.procID;
      m_toMe.push_back(b);
    }
  sort(m_toMe.begin(), m_toMe.end());

  if(T::preAllocatable() == 2) // dynamic allocatable, need two pass
    {
      CH_TIME("MPI_ Phase 1 of 2 Phase: preAllocatable==2");
      // in the non-preallocatable case, I need to message the
      // values for the m_toMe[*].size
      Vector<unsigned long> fdata;
      Vector<unsigned long> tdata;
      int count = 1;
      int scount = 1;
      if(m_toMe.size() > 0)
        {
          tdata.resize(m_toMe.size(), ULONG_MAX);
          m_receiveRequests.resize(numProc()-1);
          m_receiveStatus.resize(numProc()-1);
          MPI_Request* Rptr = &(m_receiveRequests[0]);

          int lastProc = m_toMe[0].procID;
          int messageSize = 1;
          unsigned long * dataPtr = &(tdata[0]);
          int i = 1;
          
          for(;i<m_toMe.size(); ++i)
            {
              bufEntry& b = m_toMe[i];
              if(b.procID == lastProc) 
                messageSize++;
              else
                {                     
                 
                  MPI_Irecv(dataPtr, messageSize, MPI_UNSIGNED_LONG, lastProc,
                            1, Chombo_MPI::comm, Rptr);
                  Rptr++;
              
                  lastProc = b.procID;
                  messageSize = 1;
                  dataPtr = &(tdata[i]);
                  count++;
                }
            }
          
          MPI_Irecv(dataPtr, messageSize, MPI_UNSIGNED_LONG, lastProc,
                    1, Chombo_MPI::comm, Rptr );

         
        }
      if(m_fromMe.size() > 0)
        {
          fdata.resize(m_fromMe.size());
          fdata[0]=m_fromMe[0].size;
          m_sendRequests.resize(numProc()-1);
          m_sendStatus.resize(numProc()-1);
          MPI_Request* Rptr = &(m_sendRequests[0]);

          int lastProc = m_fromMe[0].procID;
          int messageSize = 1;
          unsigned long * dataPtr = &(fdata[0]);
          int i = 1;
          for(;i<m_fromMe.size(); ++i)
            {
              fdata[i]    = m_fromMe[i].size;
              bufEntry& b = m_fromMe[i];
              if(b.procID == lastProc) 
                messageSize++;
              else
                {                     
                
                  MPI_Isend(dataPtr, messageSize, MPI_UNSIGNED_LONG, lastProc,
                            1, Chombo_MPI::comm, Rptr);
                 
                
                  Rptr++;
                  lastProc = b.procID;
                  messageSize = 1;
                  dataPtr = &(fdata[i]);
                  scount++;
                }
            }
   
          MPI_Isend(dataPtr, messageSize, MPI_UNSIGNED_LONG, lastProc,
                    1, Chombo_MPI::comm, Rptr);
     
          
        }

      if(m_toMe.size() > 0)
        {
          
          int result = MPI_Waitall(count, &(m_receiveRequests[0]), &(m_receiveStatus[0]));
          if(result != MPI_SUCCESS)
            {
              MayDay::Error("First pass of two-phase communication failed");
            }
          
          for(int i=0; i<m_toMe.size(); ++i) {
            CH_assert(tdata[i] != ULONG_MAX);
            m_toMe[i].size = tdata[i];
            recBufferSize+= tdata[i];
          }
          
  
        }
      if(m_fromMe.size() > 0)
        {
                   
          int result = MPI_Waitall(scount, &(m_sendRequests[0]), &(m_sendStatus[0]));
          if(result != MPI_SUCCESS)
            {
              MayDay::Error("First pass of two-phase communication failed");
            }
     
        }
    }

  // allocate send and receveive buffer space.

  if(sendBufferSize > m_sendcapacity)
    {
      free(m_sendbuffer);
      m_sendbuffer = malloc(sendBufferSize);
      if(m_sendbuffer == NULL)
        {
          MayDay::Error("Out of memory in BoxLayoutData::allocatebuffers");
        }
      m_sendcapacity = sendBufferSize;
    }

  if(recBufferSize > m_reccapacity)
    {
      free(m_recbuffer);
      m_recbuffer = malloc(recBufferSize);
      if(m_recbuffer == NULL)
        {
          MayDay::Error("Out of memory in BoxLayoutData::allocatebuffers");
        }
      m_reccapacity = recBufferSize;
    }

  /*
    pout()<<"\n";
    for(int i=0; i<m_fromMe.size(); i++)
    pout()<<m_fromMe[i].item->region<<"{"<<m_fromMe[i].procID<<"}"<<" ";
    pout() <<"::::";
    for(int i=0; i<m_toMe.size(); i++)
    pout()<<m_toMe[i].item->region<<"{"<<m_toMe[i].procID<<"}"<<" ";
    pout() << endl;
  */

  char* nextFree = (char*)m_sendbuffer;
  if(m_fromMe.size() > 0)
    {
      for(unsigned int i=0; i<m_fromMe.size(); ++i)
        {
          m_fromMe[i].bufPtr = nextFree;
          nextFree += m_fromMe[i].size;
        }
    }

  nextFree = (char*)m_recbuffer;
  if(m_toMe.size() > 0)
    {
      for(unsigned int i=0; i<m_toMe.size(); ++i)
        {
          m_toMe[i].bufPtr = nextFree;
          nextFree += m_toMe[i].size;
        }
    }

  // since fromMe and toMe are sorted based on procID, messages can now be grouped
  // together on a per-processor basis.

}

template<class T>
void BoxLayoutData<T>::writeSendDataFromMeIntoBuffers(const BoxLayoutData<T>& a_src,
                                                  const Interval&     a_srcComps,
                                                  const LDOperator<T>& a_op) const
{
  CH_TIME("write Data to buffers");
  for(unsigned int i=0; i<m_fromMe.size(); ++i)
    {
      const bufEntry& entry = m_fromMe[i];
      a_op.linearOut(a_src[entry.item->fromIndex], entry.bufPtr,
                     entry.item->fromRegion, a_srcComps);
    }
}

template<class T>
void BoxLayoutData<T>::postSendsFromMe() const
{
  CH_TIME("post Sends");
  // now we get the magic of message coalescence
  // fromMe has already been sorted in the allocateBuffers() step.

  this->numSends = m_fromMe.size();

  if(this->numSends > 1){
    for(unsigned int i=m_fromMe.size()-1; i>0; --i)
      {
        if(m_fromMe[i].procID == m_fromMe[i-1].procID)
          {   
            this->numSends--;
            m_fromMe[i-1].size = m_fromMe[i-1].size + m_fromMe[i].size;
            m_fromMe[i].size = 0;
          }
      }
  }
  m_sendRequests.resize(this->numSends);
  std::list<MPI_Request> extraRequests;

  unsigned int next=0;
  long long maxSize = 0;
  for(int i=0; i<this->numSends; ++i)
    {
      const bufEntry& entry = m_fromMe[next];
      char*  buffer = (char*)entry.bufPtr;
      size_t bsize = entry.size;
      int idtag=0;
      while(bsize > CH_MAX_MPI_MESSAGE_SIZE){
        extraRequests.push_back(MPI_Request());
        {
          CH_TIME("MPI_Isend");
          MPI_Isend(buffer, CH_MAX_MPI_MESSAGE_SIZE, MPI_BYTE, entry.procID,
                    idtag, Chombo_MPI::comm, &(extraRequests.back()));
        }
        maxSize = CH_MAX_MPI_MESSAGE_SIZE;
        bsize -= CH_MAX_MPI_MESSAGE_SIZE;
        buffer+=CH_MAX_MPI_MESSAGE_SIZE;
        idtag++;
      }
      {
        CH_TIME("MPI_Isend");
        MPI_Isend(buffer, bsize, MPI_BYTE, entry.procID,
                idtag, Chombo_MPI::comm, &(m_sendRequests[i]));
      }
      maxSize = Max<long long>(bsize, maxSize);
      ++next;
      while(next < m_fromMe.size() && m_fromMe[next].size == 0) ++next;
    }
  for(std::list<MPI_Request>::iterator it = extraRequests.begin(); it != extraRequests.end(); ++it){
    m_sendRequests.push_back(*it);
  }
  this->numSends = m_sendRequests.size();

  CH_MaxMPISendSize = Max<long long>(CH_MaxMPISendSize, maxSize);
  
}

template<class T>
void BoxLayoutData<T>::postReceivesToMe() const
{
  CH_TIME("post Receives");
  this->numReceives = m_toMe.size();

  if(this->numReceives > 1){
    for(unsigned int i=m_toMe.size()-1; i>0; --i)
      {
        if(m_toMe[i].procID == m_toMe[i-1].procID)
          {
            this->numReceives--;
            m_toMe[i-1].size += m_toMe[i].size;
            m_toMe[i].size = 0;
          }
        
      }
  }
  m_receiveRequests.resize(this->numReceives);
  std::list<MPI_Request> extraRequests;
  unsigned int next=0;
  long long maxSize = 0;
  for(int i=0; i<this->numReceives; ++i)
    {
      const bufEntry& entry = m_toMe[next];
      char*  buffer = (char*)entry.bufPtr;
      size_t bsize = entry.size;
      int idtag=0;
      while(bsize > CH_MAX_MPI_MESSAGE_SIZE){
        extraRequests.push_back(MPI_Request());
        {
          CH_TIME("MPI_Irecv");
          MPI_Irecv(buffer, CH_MAX_MPI_MESSAGE_SIZE, MPI_BYTE, entry.procID,
                    idtag, Chombo_MPI::comm, &(extraRequests.back()));
        }
        maxSize = CH_MAX_MPI_MESSAGE_SIZE;
        bsize -= CH_MAX_MPI_MESSAGE_SIZE;
        buffer+=CH_MAX_MPI_MESSAGE_SIZE;
        idtag++;
      }
      {
        CH_TIME("MPI_Irecv");
        MPI_Irecv(buffer, bsize, MPI_BYTE, entry.procID,
                  idtag, Chombo_MPI::comm, &(m_receiveRequests[i]));
      }
      ++next;
      maxSize = Max<long long>(bsize, maxSize);
      while(next < m_toMe.size() && m_toMe[next].size == 0) ++next;
    }
  for(std::list<MPI_Request>::iterator it = extraRequests.begin(); it != extraRequests.end(); ++it){
    m_receiveRequests.push_back(*it);
  }
  this->numReceives = m_receiveRequests.size();

  CH_MaxMPIRecvSize = Max<long long>(CH_MaxMPIRecvSize, maxSize);
  //pout()<<"maxSize="<<maxSize<<" posted "<<this->numReceives<<" receives\n";

}

template<class T>
void BoxLayoutData<T>::unpackReceivesToMe(BoxLayoutData<T>& a_dest,
                                      const Interval&   a_destComps,
                                      const LDOperator<T>& a_op) const
{
 
  CH_TIME("unpack messages");
  if(this->numReceives > 0){
    m_receiveStatus.resize(this->numReceives);
    int result;
    {
      CH_TIME("MPI_Waitall");
      result = MPI_Waitall(this->numReceives, &(m_receiveRequests[0]), 
                             &(m_receiveStatus[0]));
    }
    if(result != MPI_SUCCESS)
      {
        //hell if I know what to do about failed messaging here
      }

    for(unsigned int i=0; i<m_toMe.size(); ++i)
      {
        const bufEntry& entry = m_toMe[i];
        a_op.linearIn(a_dest[entry.item->toIndex], entry.bufPtr, entry.item->toRegion, a_destComps);
      }

 
  }
  this->numReceives = 0;
}

template<class T>
void BoxLayoutData<T>::unpackReceivesToMe_append(LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                                 const Interval&   a_destComps,
                                                 int ncomp,
                                                 const DataFactory<T>& factory,

                                                 const LDOperator<T>& a_op) const
{

  if(this->numReceives > 0){
   m_receiveStatus.resize(this->numReceives);
   int result;
   {
     CH_TIME("MPI_Waitall");
     result = MPI_Waitall(this->numReceives, &(m_receiveRequests[0]), 
                              &(m_receiveStatus[0]));
   }
   if(result != MPI_SUCCESS)
      {
        //hell if I know what to do about failed messaging here
      }

    for(unsigned int i=0; i<m_toMe.size(); ++i)
      {
        const bufEntry& entry = m_toMe[i];
        const MotionItem& item = *(entry.item);
        RefCountedPtr<T> newT( factory.create(item.toRegion, ncomp, item.toIndex) );;

        a_op.linearIn(*newT, entry.bufPtr, item.toRegion, a_destComps);
        a_dest[item.toIndex].push_back(newT);
      }


  }
  this->numReceives = 0;
}
#endif

template <class T>
void BoxLayoutData<T>::generalCopyTo(const BoxLayout& a_destGrids,
                                     LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                     const Interval& a_srcComps,
                                     const ProblemDomain& a_domain,
                                     const Copier& a_copier,
                                     const DataFactory<T>& factory) const
{

 CH_assert(T::preAllocatable() == 0);
  a_dest.define(a_destGrids);

  LDOperator<T> a_op;

  int ncomp = a_srcComps.size();
  Interval destComps(0, ncomp-1);
  allocateBuffers(*this,  a_srcComps,
                  *this,  destComps,
                  a_copier, a_op);

  writeSendDataFromMeIntoBuffers(*this, a_srcComps, a_op);

  // If there is nothing to recv/send, don't go into these functions
  // and allocate memory that will not be freed later.  (ndk)
  // The #ifdef CH_MPI is for the m_toMe and m_fromMe
#ifdef CH_MPI
  this->numReceives = m_toMe.size();
  if (this->numReceives > 0) {
    postReceivesToMe(); // all non-blocking
  }

  this->numSends = m_fromMe.size();
  if (this->numSends > 0) {
    postSendsFromMe();  // all non-blocking
  }
#endif

    // perform local copy
  for(CopyIterator it(a_copier, CopyIterator::LOCAL); it.ok(); ++it)
    {
      const MotionItem& item = it();
      RefCountedPtr<T> newT( factory.create(item.toRegion, ncomp, item.toIndex) );

      a_op.op(*newT, item.fromRegion,
              destComps,
              item.toRegion,
              this->operator[](item.fromIndex),
              a_srcComps);
      a_dest[item.toIndex].push_back(newT);
    }

  // Uncomment and Move this out of unpackReceivesToMe()  (ndk)
  completePendingSends(); // wait for sends from possible previous operation

  unpackReceivesToMe_append(a_dest, destComps, ncomp, factory, a_op); // nullOp in uniprocessor mode
}


template <class T>
void BoxLayoutData<T>::generalCopyTo(const BoxLayout& a_destGrids,
                                     LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                     const Interval& a_srcComps,
                                     const ProblemDomain& a_domain,
                                     const DataFactory<T>& factory) const
{
  Copier copier;
  copier.define(this->m_boxLayout, a_destGrids, a_domain, IntVect::Zero);

  generalCopyTo(a_destGrids, a_dest, a_srcComps, a_domain, copier, factory);
}

#include "NamespaceFooter.H"
#endif
