#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>
using namespace std;

#include "MayDay.H"
#include "CHOMBO_VERSION.H"
#ifdef CH_MPI
#include "mpi.h"
#endif

#include "FORT_PROTO.H"
#include "NamespaceHeader.H"

//
// The definition of our NULL string used as default argument.
//
const char * const MayDay::m_nullString = "";

//
// The definition of our version string.
//
#define bl_str(s) # s
#define bl_xstr(s) bl_str(s)

const char * const MayDay::version = "MayDay version " bl_xstr(CHOMBO_VERSION) " built " __DATE__ " at " __TIME__ ;

#undef bl_str
#undef bl_xstr


//
// This is used by MayDay::Error(), MayDay::Abort() and MayDay::Warning()
// to ensure that when writing the message to stderr, that no additional
// heap-based memory is allocated.
//
static void write_to_stderr_without_buffering (const char * const a_str)
{
  //
  // Flush all buffers.
  //
  fflush(NULL);

  if (a_str)
  {
    //
    // Add some `!'s and a newline to the string.
    //
    const char * const end = " !!!\n";

    fwrite(a_str, strlen(a_str), 1, stderr);
    fwrite(end  , strlen(end  ), 1, stderr);
  }
}

void MayDay::Error(const char * const a_msg, int exit_code)
{
  write_to_stderr_without_buffering(a_msg);
#ifdef CH_MPI
  MPI_Abort(MPI_COMM_WORLD ,exit_code);
  // this shouldn't return, but if it does, exit serially
#endif
  exit(exit_code);
}

void MayDay::Abort(const char * const a_msg)
{
  write_to_stderr_without_buffering(a_msg);
#ifdef CH_MPI
  //MPI_Abort(MPI_COMM_WORLD ,CH_DEFAULT_ERROR_CODE);
  // this shouldn't return, but if it does, abort serially
#endif
  abort();
}

void MayDay::Warning(const char * const a_msg)
{
  write_to_stderr_without_buffering(a_msg);
}


extern "C" {

  void
  FORTRAN_NAME(MAYDAYERROR,maydayerror) (void)
  {
    MayDay::Error( "A ChomboFortran routine called MAYDAYERROR().  Rerun with the debugger to find where." ) ;
  }

  void
  FORTRAN_NAME(MAYDAYABORT,maydayabort) (void)
  {
    MayDay::Abort( "A ChomboFortran routine called MAYDAYABORT().  Rerun with the debugger to find where." ) ;
  }

  void
  FORTRAN_NAME(MAYDAYEXIT,maydayexit) (int* exit_code)
  {
    MayDay::Error( "A ChomboFortran routine called MAYDAYEXIT().  Rerun with the debugger to find where." ,*exit_code ) ;
  }

  void
  FORTRAN_NAME(MAYDAY_ERROR,mayday_error) (void)
  {
    MayDay::Error( "A ChomboFortran routine called MAYDAY_ERROR().  Rerun with the debugger to find where." ) ;
  }

  void
  FORTRAN_NAME(MAYDAY_ABORT,mayday_abort) (void)
  {
    MayDay::Abort( "A ChomboFortran routine called MAYDAY_ABORT().  Rerun with the debugger to find where." ) ;
  }

  void
  FORTRAN_NAME(MAYDAY_EXIT,mayday_exit) (int* exit_code)
  {
    MayDay::Error( "A ChomboFortran routine called MAYDAY_EXIT().  Rerun with the debugger to find where." ,*exit_code ) ;
  }

}
#include "NamespaceFooter.H"
