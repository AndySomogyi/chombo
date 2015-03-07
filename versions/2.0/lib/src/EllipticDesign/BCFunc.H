#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// BVS, July 10, 2003
// revised by dtg in 2006

#ifndef _BCFUNC_H_
#define _BCFUNC_H_

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
                      const Box& a_valid,
                      const ProblemDomain& a_domain,
                      Real                 a_dx,
                      bool                 a_homogeneous);


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


///
/**
   Neumann bc for a particular side.   
   For use in AMRPoissonOp.
 */
void NeumBC(FArrayBox& a_state,
            const Box& a_valid,
            Real a_dx,
            bool a_homogeneous,
            BCValueFunc a_value,
            int a_dir,
            Side::LoHiSide a_side);

///
/**
   Neumann bcs for all sides
   For use in AMRPoissonOp.
 */
void NeumBC(FArrayBox& a_state,
            const Box& a_valid,
            Real a_dx,
            bool a_homogeneous,
            BCValueFunc a_value);

///
/**
   Dirichlet boundary conditions for a side.
   For use in AMRPoissonOp.
 */
void DiriBC(FArrayBox& a_state,
            const Box& valid,
            Real a_dx,
            bool a_homogeneous,
            BCValueFunc a_value,
            int a_dir,
            Side::LoHiSide a_side,
            int order = 1);

///
/**
   Dirchelet boundary conditions for one side.
   For use in AMRPoissonOp.
 */
void DiriBC(FArrayBox& a_state,
            const Box& valid,
            Real a_dx,
            bool a_homogeneous,
            BCValueFunc a_value,
            int order = 1);


///
/**
   No slip vector bc (zero all comps).
   need a_state.ncomp == spacedim
   For use in ResistivityOp, for example.
 */
void NoSlipVectorBC(FArrayBox& a_state,
                    const Box& a_valid,
                    Real a_dx,
                    int a_dir,
                    Side::LoHiSide a_side,
                    int order = 2);

///
/**
   0 normal comp, reflective for all other comps 
   need a_state.ncomp == spacedim
   For use in ResistivityOp, for example.
 */
void ReflectiveVectorBC(FArrayBox& a_state,
                        const Box& a_valid,
                        Real a_dx,
                        int a_dir,
                        Side::LoHiSide a_side,
                        int            a_order = 2);

class BCFunction
{
public:
  virtual ~BCFunction(){;}
  virtual void operator()(FArrayBox&           a_state,
                          const Box& a_valid,
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
                  const Box& a_valid,
                  const ProblemDomain& a_domain,
                  Real                 a_dx,
                  bool                 a_homogeneous)
  {
    if(m_funcptr!= NULL)
      m_funcptr(a_state, a_valid, a_domain, a_dx, a_homogeneous);
    else
      m_bc->operator()(a_state, a_valid, a_domain, a_dx, a_homogeneous);
  }
protected:
  BCFunc m_funcptr;
  RefCountedPtr<BCFunction> m_bc;
};

#include "NamespaceFooter.H"
#endif
