#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif


#ifndef _SPACE_H_
#define _SPACE_H_

#ifdef HAVE_CONFIG_H
#include <autoconfig.H>
#endif

#ifdef CH_LANG_CC
# include "MayDay.H"
#endif

#include "NamespaceHeader.H"

#if !defined(CH_LANG_FORT) || defined(CH_LANG_CC)

/**
  The SpaceDim C++ integer constant specifies the dimension of the
  Chombo library to be built.  It is initialized to the value of
  the macro CH_SPACEDIM.  Hence it will have one of the values 1, 2, or 3.

  C++ code that depends on the dimensionality of the problem that
  you're trying to solve should be cast in terms of this constant,
  instead of the CH_SPACEDIM macro, if possible, as constants are
  visible to debuggers while macros are not.
*/
   const int SpaceDim = CH_SPACEDIM;
#endif

#ifdef CH_LANG_CC
/**
@ingroup macros
compile-time cpp macro that replaces 'assert()' and aborts
correctly when running in parallel. */
# ifndef NDEBUG
#  define CH_assert_str(s) # s
#  define CH_assert_xstr(s) CH_assert_str(s)
#  define CH_assert(cond) if(!(cond)){ CH_XD::MayDay::Abort( __FILE__ ":" CH_assert_xstr(__LINE__) ": Assertion `" #cond "' failed."); }
# else
#  define CH_assert(cond) (void)0
#endif

/** Compile-time assertion */
template<bool> struct static_assert;
template<> struct static_assert<true> {};
#define STATIC_ASSERT( a ) static_assert<a>()


#endif //CH_LANG_CC

#ifdef  CH_LANG_FORT
#ifndef NDEBUG
#  define CH_assert(cond) if(.NOT. cond) call MAYDAYERROR()
#else
#  define CH_assert(cond)
#endif
#endif //CH_LANG_FORT

#ifdef DOXYGEN
/** \defgroup macros Chombo Macros
 */

/**
@ingroup macros
compile-time cpp macro that controls the dimensionality of
    the Chombo package.  */
#define CH_SPACEDIM
#undef  CH_SPACEDIM
#define CH_SPACEDIM 3
#endif

#if ! defined(CH_SPACEDIM)
#error CH_SPACEDIM must be defined
#endif

#if (CH_SPACEDIM != 1 && CH_SPACEDIM !=2 && CH_SPACEDIM != 3)
#error CH_SPACEDIM must be either 1, 2, or 3
#endif

#ifdef DOXYGEN
#undef CH_SPACEDIM
#endif

#include "NamespaceFooter.H"

#endif
#ifdef CH_LANG_CC
// That matches the include guard.  The rest of this file is all macros
// whose definitions depend on CH_SPACEDIM.  In CH_MULTIDIM mode, we want
// to rescan these macros after every time CH_SPACEDIM has changed.
// Of course, if CH_MULTIDIM is not defined, we want to scan just once.
#endif

#ifndef LAST_SPACE_H_SPACEDIM
#define LAST_SPACE_H_SPACEDIM 0
#endif

#if CH_SPACEDIM != LAST_SPACE_H_SPACEDIM
#undef  LAST_SPACE_H_SPACEDIM


#ifdef D_EXPR
#undef D_EXPR
#undef D_DECL
#undef D_TERM
#endif

#if CH_SPACEDIM==1
#  define D_EXPR(a,b,c) ((void)((a),0))
#  define D_DECL(a,b,c) a
#  define D_TERM(a,b,c) a
#  define LAST_SPACE_H_SPACEDIM 1


#elif CH_SPACEDIM==2
#  define D_EXPR(a,b,c) ((void)((a),(b),0))
#  define D_DECL(a,b,c) a,b
#  define D_TERM(a,b,c) a b
#  define LAST_SPACE_H_SPACEDIM 2

#elif CH_SPACEDIM==3
#ifdef DOXYGEN
/**
@ingroup macros

  The D_EXPR(a,b,c) macro expands to a comma expression that will
  evaluate SpaceDim of its arguments, counting from the left; i.e.
  CH_SPACEDIM==1 implies only `a' is evaluated, CH_SPACEDIM==2
  implies both `a' and `b', and CH_SPACEDIM==3 implies all three
  arguments are evaluated.  This macro can be used to enhance
  the portability of code that works with structures whose size
  depends on CH_SPACEDIM.

  For example:

     D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s);

  More explicitly, this macro evaluates to one of

     ((void)((a),0))

     ((void)((a),(b),0))

     ((void)((a),(b),(c),0))

  corresponding to CH_SPACEDIM values of 1, 2, or 3.
*/
#endif
#  define D_EXPR(a,b,c) ((void)((a),(b),(c),0))

#ifdef DOXYGEN
/**
@ingroup macros

  The D_DECL(a,b,c) macro expands to a comma-separated list of
  1, 2, or all 3 of the arguments of the call, depending on the
  value of CH_SPACEDIM.  This can be used to write portable
  function calls that depend on 1, 2, or 3 arguments, corresponding
  to the value of CH_SPACEDIM.

  For example:

     return IntVect(D_DECL(p[0] + s, p[1] + s, p[2] + s));

  More explicitly, this macro evaluates to one of

     a

     a,b

     a,b,c

  corresponding to CH_SPACEDIM values of 1, 2, or 3.
*/
#endif
#  define D_DECL(a,b,c) a,b,c

#ifdef DOXYGEN
/**
@ingroup macros
   define D_TERM(a,b,c) a b c

  More explicitly, this macro evaluates to one of

     a

     a b

     a b c

    corresponding to CH_SPACEDIM values of 1, 2, or 3.

    handy for situations like:

    Real norm = sqrt(D_TERM(u[0]*u[0],+u[1]*u[1],+u[2]*u[2]));
*/
#endif
#  define D_TERM(a,b,c) a b c
#  define LAST_SPACE_H_SPACEDIM 3

#endif // matches elif CH_SPACEDIM==3

#endif // matches CH_SPACEDIM != LAST_SPACE_H_SPACEDIM
