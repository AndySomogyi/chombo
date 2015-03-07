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

#ifndef _BASEEBFACEFAB_H_
#define _BASEEBFACEFAB_H_

#include <cmath>
#include <cstdlib>
#include "SPACE.H"
#include "Vector.H"
#include "BaseIFFAB.H"
#include "IntVectSet.H"
#include "BaseFab.H"
#include "EBISBox.H"
#include "NamespaceHeader.H"

///  Array defined at the Faces of an Box in an EBISBox.
/**
   An BaseEBFaceFAB is a templated holder for face-centered
   data over a region which consists of the intersection
   of a face-centered box and an EBISBox.
   At every uncovered VoF in this intersection,
   the BaseEBFaceFAB contains a specified number of data values.
   At singly valued faces,
   the data is stored internally in a BaseFab.  At multiply-valued
   faces, the data is stored internally in an  INFab.
   A multiply-valued face is defined to be a face which
   abuts at least one multiply-valued Face.
   BaseEBFaceFAB provides indexing by Face.
   This class does not provide a copy constructor or assignment operator
   for all the usual reasons.
 */
template <class T>
class BaseEBFaceFAB
{
public:
  ///
  BaseEBFaceFAB();

  ///
  /**
     Box going into this needs to be cell-centered.
     Data will exist over the surrounding nodes of the box.
  */
  BaseEBFaceFAB(const EBISBox& a_ebisBox,
                const Box& a_region,
                int a_iDir, int a_nVar);

  ///
  /**
     Box going into this needs to be cell-centered.
     Data will exist over the surrounding nodes of the box.
  */
  void
  define(const EBISBox& a_ebisBox,
         const Box& a_region,
         int a_iDir, int a_nVar);

  ///
  virtual ~BaseEBFaceFAB();

  ///
  void  clear();

  ///
  void setVal(const T& value);

  ///
  void setVal(int ivar,const T& value);

  ///
  bool isDefined() const;

  ///
  int nComp() const ;

  ///
  const BaseIFFAB<T>&   getMultiValuedFAB() const;

  ///
  BaseIFFAB<T>&  getMultiValuedFAB() ;

  ///
  const BaseFab<T>& getSingleValuedFAB() const;

  ///
  const IntVectSet& getMultiCells() const;

  ///
  BaseFab<T>&   getSingleValuedFAB() ;

  ///
  int  direction() const;

  ///
  /**
     Returns the FACE-CENTERED region of the fab
  */
  const Box& getRegion() const ;

  ///
  /**
     Returns the CELL-CENTERED region of the fab
  */
  const Box& getCellRegion() const ;

  ///
  T& operator() (const FaceIndex& a_facein, int  a_nVarLoc);

  ///
  const T&
  operator() (const FaceIndex& a_facein, int  a_nVarLoc) const;

  ///
  static bool preAllocatable(){ return true;}

  ///
  /**
     Box going into this should be CELL CENTERED.
     copy done over surrounding nodes.
  */
  void copy(const Box& RegionFrom,
            const Interval& destInt,
            const Box& RegionTo,
            const BaseEBFaceFAB<T>& source,
            const Interval& srcInt);

  BaseEBFaceFAB<T>& copy(const BaseEBFaceFAB<T>& source);

  int size(const Box& R, const Interval& comps) const;

  void linearOut(void* buf, const Box& R, const Interval& comps) const;

  void linearIn(void* buf, const Box& R, const Interval& comps);

  /// Invalid but necessary for LevelData<T> to compile
  BaseEBFaceFAB(const Box& a_region, int a_nVar)
  {
    MayDay::Error("invalid constructor called for BaseEBFaceFAB");
  }

  ///
  const EBISBox& getEBISBox() const;

  ///
  void
  setCoveredFaceVal(const T&    a_val,
                    const int&  a_comp,
                    const bool& a_doMulti=true);

protected:
  ///data at faces which abut multi-valued cells
  BaseIFFAB<T> m_irrFAB;
  ///data at faces between two single-valued cells
  BaseFab<T> m_regFAB;

  ///number of data values at each face in BaseEBFaceFAB
  int m_nComp;

  //direction of the faces in the BaseEBFaceFAB
  int m_iDir;

  ///the face-centered region over which the BaseEBFaceFAB lives
  Box m_regionFace;

  ///the cell-centered region over which the BaseEBFaceFAB lives
  Box m_region;

  ///has the full define function been called for the BaseEBFaceFAB?
  bool m_isDefined;

  EBISBox m_ebisBox;
private:
  void
  setDefaultValues();
  //disallowed for performance reasons
  void operator= (const BaseEBFaceFAB<T>& ebcin)
  {
    MayDay::Error("BaseEBFaceFAB<T>::operator=  not defined");
  }
  BaseEBFaceFAB (const BaseEBFaceFAB<T>& ebcin)
  {
    MayDay::Error("BaseEBFaceFAB<T> copy constructor  not defined");
  }
};

#include "NamespaceFooter.H"

#ifndef CH_EXPLICIT_TEMPLATES
#include "BaseEBFaceFABI.H"
#endif // CH_EXPLICIT_TEMPLATES

#endif
