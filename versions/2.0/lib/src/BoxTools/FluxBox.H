#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _FLUXBOX_H_
#define _FLUXBOX_H_

#include "Box.H"
#include "Vector.H"
#include "FArrayBox.H"
#include "NamespaceHeader.H"

/// A FArrayBox-like container for face-centered fluxes
/** This is a class to contain face-centered fluxes on a box.
 */
class FluxBox
{

public:
  /// Default constructor
  FluxBox();

  /// Constructs FluxBox on cell-centered box with n components
  FluxBox(const Box& bx, int n=1);

  /// Destructor
  ~FluxBox();

  /// Resize FluxBox similar to BaseFab::resize()
  void resize(const Box& bx, int n=1);

  /// Define function
  void define(const Box& bx, int n=1);

  /// Returns the FluxBox to the undefined state
  void clear();

  //  access functions

  /// Number of components
  int nComp() const;

  /// Returns cell-centered box which defines fluxBox
  const Box& box() const;

  /// Returns face-centered flux in direction \em dir
  FArrayBox& getFlux(const int dir);

  /// Returns const reference to face-centered flux in direction \em dir
  const FArrayBox& getFlux(const int dir) const;

  /// Returns FArrayBox in direction dir
  FArrayBox& operator[] (const int dir);

  /// Constant version
  const FArrayBox& operator[] (const int dir) const;

  //  data modification functions

  /// Set all fluxes to val
  void setVal(const Real val);

  /// Set fluxes in direction dir to val
  void setVal(const Real val, const int dir);

  /// More specific setVal
  void setVal(const Real val, const int dir, const int startComp,
              const int nComp);

  /// Sets fluxes on faces surrounding cell-centered box bx
  void setVal(const Real val, const Box& bx);

  /// Most specific setVal
  /**
     Sets fluxes on faces surrounding cell-centered box bx
  */
  void setVal(const Real val, const Box& bx, const int dir,
              const int startComp, const int nComp);

  /// Copy from src to this FluxBox -- sizes must be identical
  void copy(const FluxBox& src);

  /// Copy on overlap, for all directions
  void copy(const FluxBox& src, const int srcComp,
            const int destComp, const int numComp);

  /// Copy on overlap of FluxBoxes, in direction dir
  void copy(const FluxBox& src, const int dir, const int srcComp,
            const int destComp, const int numComp);

  /// Copies from a subsection of one box into another.
  /**
     Assumes the boxes are both in the same index space, and that
     box R is completely contained in both the src and destination
     boxes.  This function required by BoxLayoutData
  */
  void copy(const Box& R, const Interval& Cdest,  const FluxBox& src,
            const Interval& Csrc);

  /// Modifies this FluxBox by copying the contents of \em src into it.
  /**
     This, the most general form
     of copy, specifies the contents of any sub-box srcbox
     in `FluxBox' src may be copied into a (possibly
     different) destbox in the destination `FluxBox'.  Note
     that although the srcbox and the destbox may be
     disjoint, they must be the same size and shape.  If the
     sizes differ, the copy is undefined and a runtime error
     results.  This copy function is the only one of the
     copy functions to allow a copy between differing
     boxes. The user also specifies how many components are
     copied, starting at component srccomp in src and stored
     starting at component destcomp.  The results are
     UNDEFINED if the src and dest FluxBoxes are the same and
     the srcbox and destbox overlap.
  */
  void copy (const Box&        srcbox,
             const Interval&   destcomps,
             const Box&        destbox,
             const FluxBox&    src,
             const Interval&   srccomps);

  /// Modifies this FluxBox to its additive inverse.
  /**
     Modifies this FluxBox by replacing each value with its additive
     inverse, for the given range of components and within the given subbox.
     Returns *this.
  */
  FluxBox& negate (const Box& subbox,
                   int        comp=0,
                   int        numcomp=1);

  /// Modifies this FluxBox to its additive inverse.
  /**
     Modifies this FluxBox by replacing each value with its additive
     inverse, for the given range of components over the whole domain of
     the FluxBox.  Returns *this.
  */
  FluxBox& negate (int comp,
                   int numcomp=1);

  /// Modifies this FluxBox to its additive inverse.
  /**
     Modifies this FluxBox by replacing each value with its additive
     inverse for all components over the whole domain of the FluxBox.
     Returns *this.
  */
  FluxBox& negate ();

  /// Modifies this FluxBox by adding the scalar Real r to all values.
  FluxBox& operator+= (Real r);

  /// Modifies this FluxBox by incrementing with the argument FluxBox.
  /**
     Modifies this FluxBox by pointwise addition of the values of the
     argument FluxBox.  You might come to grief if the domains of the
     FArrayBoxes don't match, just as in FArrayBox::plus().
  */
  FluxBox& operator+= (const FluxBox& f);

  /// Modifies this FluxBox by subtracting the scalar Real r to all values.
  FluxBox& operator-= (Real r);

  /// Modifies this FluxBox by decrementing with the argument FluxBox.
  /**
     Modifies this FluxBox by pointwise subtraction of the values of the
     argument FluxBox.  You might come to grief if the domains of the
     FluxBoxes don't match, just as in FArrayBox::minus().
  */
  FluxBox& operator-= (const FluxBox& f);

  /// Modifies this FluxBox by multiplying all values by the scalar Real r.
  FluxBox& operator*= (Real r);

  /// Modifies this FluxBox by multiplying by the argument FluxBox.
  /**
     Modifies this FluxBox by pointwise multiplication of the values by the
     argument FluxBox.  You might come to grief if the domains of the
     FluxBoxes don't match, just as in FArrayBox::mult().
  */
  FluxBox& operator*= (const FluxBox& f);

  /// Modifies this FluxBox by shifting its domain box.
  FluxBox& shift (const IntVect& v);

  //  Linearization Functions

  /// Returns size of linearized data over \em bx.
  /**
     Returns size, in number of bytes, of a flat linear
     representation of data in components comps in faces around
     cell-centered box R
  */
  int size(const Box& bx, const Interval& comps) const;

  /// Writes a linear representation of this FluxBox
  /**
     Write a linear representaion of the internal data for
     the faces surrounding cell-centered box R.  Assumes that sufficient
     memory for the buffer has already been allocated by the caller
  */
  void linearOut(void*buf, const Box& R, const Interval& comps) const;

  /// Read a linear representation of the data over the Box \rm R
  /**
      Reads in the output of linearOut
  */
  void linearIn(void*buf, const Box& R, const Interval& comps);

  /// Helper function for linearization
  static int preAllocatable() {return 0;}


protected:
  /**
     Cell-centered Box over which this FluxBox is defined
  */
  Box m_bx;

  /**
     Number of variables on each face
  */
  int m_nvar;

  /**
     CH_SPACEDIM FArrayBoxes which hold fluxes
  */
  Vector<FArrayBox*> m_fluxes;


private:
  // these are disallowed

  /// Disallowed
  FluxBox (const FluxBox&);
  /// Disallowed
  FluxBox& operator = (const FluxBox&);

};

#include "NamespaceFooter.H"
#endif
