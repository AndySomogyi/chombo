#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _BITSET_H_
#define _BITSET_H_

//
// Arch dependant setting for what a BITSETWORD should be
// BITSETWORD is the smallest efficiently accessible memory
// object.  BITSETWORD_BITS is the size, in bits, of this object
//
// Later I can make this object use a Pool, and we can
// avoid all the memory system call overhead for creation and
// use placement new.
//

#include "CH_assert.H"
#include <iostream>
#include "BaseNamespaceHeader.H"

#define BITSETWORD unsigned int
const int  BITSETWORDSIZE = 8*sizeof(unsigned int);
// BITSETWORDSIZE is no the same as sizeof(BITSETWORD).  This is the number of bits, not chars

using std::ostream;

class BitSetIterator;

///
/* stores a contiguous memory chunk and represents true with a 1 bit
   and false with a 0 bit.

   example: 35 bits, set to true, BITSETWORDSIZE=8:

   m_index = 5 (need at least 5 BITSETWORDS to hold that many bits)
   m_size  = 35

   +-------+-------+-------+-------+-------+
   1111111111111111111111111111111111111111

   now, set bit 35 to 0

   index = 35/BITSETWORDSIZE = 4       , hence, we are in the fourth BITSETWORD.
   mask  = 35 - 4*BITSETWORDSIZE =  3    hence, we are in the third bit of the fourth word

   now, we can use the masks:
static BITSETWORD trueMasks[BITSETWORDSIZE]; //10000000, 01000000, 00100000, ....

   word[index] = word[index] & ~trueMasks[mask]  (~trueMasks[3] = 11101111)

   now we have:

   +-------+-------+-------+-------+-------+
   1111111111111111111111111111111111101111
*/
class BitSet
{
public:
  ///
  BitSet();

  ///
  BitSet(int bits, bool init);

  ///
  void define(int bits, bool init);

  ///
  BitSet(const BitSet& rhs);

  ///
  BitSet& operator=(const BitSet& rhs);

  ///
  /**
    Primary criterion: m_length.
    Secondary criterion: BITSETWORD-by-BITSETWORD lexicographic comparison
    of *m_bits.
  */
  bool operator<(const BitSet& rhs) const;

  ///
  ~BitSet();

  // somewhat slow random access. Fast iterating done
  // with BitSetIterator
  bool operator[](int i) const;

  /*
      logic operations
  */

  ///
  void setTrue(int i); // set bit to 1

  ///
  void setFalse(int i); // set bit to 0

  ///
  /**
     returns 'true' if the entire bitset is zero
  */
  bool isEmpty() const;

  ///
  /**
     returns 'true' if entire bitset is 1
  */
  bool isFull() const;

  ///
  int size() const;
  static int initialize();

  int linearSize() const;

  void linearIn(const void* const inBuf);

  void linearOut(void* const a_outBuf) const;

  // not for public, used by memory tracker.
  static long int bytes;
  static long int peak;

private:
  friend class BitSetIterator;

  BITSETWORD* m_bits;
  int   m_size;
  int   m_length;  //length of char array, not bit length

  static BITSETWORD trueMasks[BITSETWORDSIZE]; //10000000, 01000000, 00100000, ....
};

// somewhat slow random access. Fast iterating done
// with BitSetIterator
inline bool BitSet::operator[](int i) const
{
  CH_assert(i>=0);
  CH_assert(i<m_size);
  int index = i/BITSETWORDSIZE;
  CH_assert(index < m_length);
  int remainder = i-BITSETWORDSIZE*index;
  BITSETWORD tmp = m_bits[index] & trueMasks[remainder];
  return tmp > 0;
}

inline int BitSet::size() const
{
  return m_size;
}

///
/* Iterate over bits in a BitSet.  return true if bit is 1

   example: 35 bits in bitset, BITSETWORDSIZE=8:

   currently at the 22nd bit: (bit # =21)

   |-m_index = 2---|
   +-------+-------+-------+-------+-------+
   1111111110001111111100011111111100011111
                   ^    ^
       m_remainder |----| = 6      ^  ^
                                   |--| m_partialBits = 3

  returns: false for operator()
.
*/
class BitSetIterator
{
public:
  ///
  BitSetIterator();

  ///
  BitSetIterator(const BitSet& bitset);

  // copy and assign should be fine

  ///
  bool operator()() const;

  ///
  bool ok() const;

  ///
  void operator++();

  ///
  void begin();

  ///
  void end();

private:
  int m_index;
  int m_remainder;
  int m_length;

  int m_partialBits;
  const BitSet* m_bitset;
  static BitSet emptyBitSet;
};

inline
BitSetIterator::BitSetIterator()
  :
  m_index(0),
  m_remainder(0),
  m_length(0),
  m_partialBits(0),
  m_bitset(&emptyBitSet)
{}

inline
BitSetIterator::BitSetIterator(const BitSet& a_bitset)
  :m_index(0), m_remainder(0), m_length(a_bitset.m_length - 1),
   m_partialBits(a_bitset.m_size - BITSETWORDSIZE*(a_bitset.m_length - 1)),
   m_bitset(&a_bitset)
{
  if(m_partialBits == BITSETWORDSIZE)
    {
      m_partialBits = 0;
      m_length++;
    }
}

inline
bool BitSetIterator::operator()() const
{
  return (m_bitset->m_bits[m_index] & BitSet::trueMasks[m_remainder] ) > 0;
}

inline
bool BitSetIterator::ok() const
{
  if(m_index < m_length) return true;
  if(m_remainder < m_partialBits) return true;
  return false;
}

inline
void BitSetIterator::operator++()
{
  ++m_remainder;
  if(m_remainder == BITSETWORDSIZE)
    {
      m_remainder = 0;
      ++m_index;
    }
}

#include "BaseNamespaceFooter.H"
#endif
