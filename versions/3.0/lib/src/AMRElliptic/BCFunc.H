#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _BCFUNC_H_
#define _BCFUNC_H_

#include "IntVect.H"
#include "RealVect.H"
#include "FArrayBox.H"
#include "ProblemDomain.H"
#include "RefCountedPtr.H"
#include "NamespaceHeader.H"

///
/**
   function interface for ghost cell boundary conditions
   of AMRPoissonOp.   If you are using Neumann or Dirichlet
   boundary conditions,   it is easiest to use the functions
   provided.
 */
typedef void(*BCFunc)(FArrayBox&           a_state,
                      const Box&           a_valid,
                      const ProblemDomain& a_domain,
                      Real                 a_dx,
                      bool                 a_homogeneous);

class BCFunction
{
public:
  virtual ~BCFunction()
  {
  }

  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous) = 0;
};

class BCHolder
{
public:
  BCHolder():m_funcptr(NULL){;}
  BCHolder(BCFunc funcptr):m_funcptr(funcptr){;}
  BCHolder(RefCountedPtr<BCFunction> refptr):m_funcptr(NULL),m_bc(refptr){;}
  void operator()(FArrayBox&           a_state,
                  const Box&           a_valid,
                  const ProblemDomain& a_domain,
                  Real                 a_dx,
                  bool                 a_homogeneous)
  {
    if (m_funcptr != NULL)
    {
      m_funcptr(a_state, a_valid, a_domain, a_dx, a_homogeneous);
    }
    else
    {
      m_bc->operator()(a_state, a_valid, a_domain, a_dx, a_homogeneous);
    }
  }

protected:
  BCFunc m_funcptr;
  RefCountedPtr<BCFunction> m_bc;
};

///
/**
   given
   pos [x,y,z] position on center of cell edge
   int dir direction, x being 0
   int side -1 for low, +1 = high,
   fill in the a_values array
*/

typedef void(*BCValueFunc)(Real*           a_pos,
                           int*            a_dir,
                           Side::LoHiSide* a_side,
                           Real*           a_value);

class BCValueFunction
{
public:
  virtual ~BCValueFunction()
  {
  }

  virtual void operator()(Real*           a_pos,
                          int*            a_dir,
                          Side::LoHiSide* a_side,
                          Real*           a_value) = 0;
};

class BCValueHolder
{
public:
  BCValueHolder():m_funcptr(NULL){;}
  BCValueHolder(BCValueFunc funcptr):m_funcptr(funcptr){;}
  BCValueHolder(RefCountedPtr<BCValueFunction> refptr):m_funcptr(NULL),m_bc(refptr){;}

  virtual ~BCValueHolder(){;}

  virtual void operator()(Real*           a_pos,
                  int*            a_dir,
                  Side::LoHiSide* a_side,
                  Real*           a_value)
  {
    if (m_funcptr != NULL)
    {
      m_funcptr(a_pos, a_dir, a_side, a_value);
    }
    else
    {
      m_bc->operator()(a_pos, a_dir, a_side, a_value);
    }
  }

protected:
  BCValueFunc m_funcptr;
  RefCountedPtr<BCValueFunction> m_bc;
};

// Class used by ConstDiriNeumBC to define BCs which can be passed anywhere a
// BCFunction/BCHolder is needed
class ConstBCFunction: public BCFunction
{
public:
  ConstBCFunction(const IntVect&  a_loSideType,
                  const RealVect& a_loSideValue,
                  const IntVect&  a_hiSideType,
                  const RealVect& a_hiSideValue);

  ~ConstBCFunction();

  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous);

protected:
  IntVect  m_loSideType;
  RealVect m_loSideValue;

  IntVect  m_hiSideType;
  RealVect m_hiSideValue;
};

///
/**
   A helper function to produce the needed object for constant
   Dirichlet/Neumann on all the faces.  The return value can be passed to
   anything expecting a BCFunction/BCHolder.

   "a_loSideType/a_hiSideType" specify the type of boundary condition in a
   given direction by having the enter corresponding to the direction set to
   0 for Neumann or 1 for Dirichlet (on the low or high side, respectively).
   "a_loSideValue/a_hiSideValue" specify the (constant) value for boundary
   condition specified above.

   For example, in 2D if "a_loSideType" = (1,0), "a_hiSideType" = (1,1),
   "a_loSideValue" = (0.0,1.0) and "a_hiSideValue" = (0.0,0.0) then the
   boundary conditions are:

          Low  side x:  Dirichlet = 0.0
          Low  side y:  Neumann   = 1.0
          High side x:  Dirichlet = 0.0
          High side y:  Dirichlet = 0.0
 */
RefCountedPtr<BCFunction> ConstDiriNeumBC(const IntVect&  a_loSideType,
                                          const RealVect& a_loSideValue,
                                          const IntVect&  a_hiSideType,
                                          const RealVect& a_hiSideValue);

///
/**
   Neumann bc for a particular side.
   For use in AMRPoissonOp.
 */
void NeumBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            const BCValueHolder&   a_value,
            int             a_dir,
            Side::LoHiSide  a_side);

///
/**
   Neumann bcs for all sides
   For use in AMRPoissonOp.
 */
void NeumBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            BCValueHolder   a_value);

///
/**
   Dirichlet boundary conditions for a side.
   For use in AMRPoissonOp.
 */
void DiriBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            BCValueHolder   a_value,
            int             a_dir,
            Side::LoHiSide  a_side,
            int             a_order = 1);

///
/**
   Dirchelet boundary conditions for one side.
   For use in AMRPoissonOp.
 */
void DiriBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            BCValueHolder   a_value,
            int             a_order = 1);

///
/**
   No slip vector bc (zero all comps).
   need a_state.ncomp == spacedim
   For use in ResistivityOp, for example.
 */
void NoSlipVectorBC(FArrayBox&     a_state,
                    const Box&     a_valid,
                    Real           a_dx,
                    int            a_dir,
                    Side::LoHiSide a_side,
                    int            a_order = 2);

///
/**
   0 normal comp, reflective for all other comps
   need a_state.ncomp == spacedim
   For use in ResistivityOp, for example.
 */
void ReflectiveVectorBC(FArrayBox&     a_state,
                        const Box&     a_valid,
                        Real           a_dx,
                        int            a_dir,
                        Side::LoHiSide a_side,
                        int            a_order = 2);

#include "NamespaceFooter.H"
#endif
