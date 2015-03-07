#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#ifndef CH_DISABLE_SIGNALS
#include <csignal>
#endif

#include "SPMD.H"
#include "CH_Attach.H"
#include "NamespaceHeader.H"

int pfds[2];

bool alreadyAttached = false;

void DebugCont()
{
// #ifndef CH_DISABLE_SIGNALS
//   signal(SIGSEGV, SIG_DFL);
//   signal(SIGABRT, SIG_DFL);
// #endif
  write(pfds[1], "OK", 3);
}


void AttachDebugger(int a_sig)
{
#ifndef CH_DISABLE_SIGNALS
  if(alreadyAttached) return;

  alreadyAttached = true;
  char buf[300];
  char* display;
  FILE *f;
  int proc = getpid();
  char title[75];

  char binaryName[201];
  memset(binaryName, 0, 201); // readlink doesn't insert NULL at end

  char linkName[30];

#ifdef CH_Linux   // note that a machine can be Linux, but not support /proc (Catamount)
  // voodoo specific to linux /proc system to get the full name of the binary file
  sprintf(linkName,"/proc/%d/exe",proc);

  ssize_t s = readlink(linkName, binaryName, 200);
#else
  // you will need to execute a 'load FILENAME' command once in gdb to let
  // it debug proper symbols.
  
  ssize_t s = 0;
#endif

  int a = 50; // used for title of emacs window

  if (s > 150)
  {
    //filename too darn big
    signal(a_sig, SIG_DFL);
    raise(a_sig);
  }

  if (s < a)
  {
    a = s;
  }

  pipe(pfds);

  //System V vs. POSIX naming differences for child signal handler.
#ifdef SIGCHLD
#ifndef SIGCLD
#define SIGCLD SIGCHLD
#endif
#endif

  signal(SIGCLD, SIG_IGN);
  if(!fork()){
    // try to find a display to send this debugger to.
    display = getenv("DISPLAY");
    printf("%s\n%s\n",linkName, display);
    
    if (display != NULL)
      {
        sprintf(title, "'%s %d'", binaryName+s-a, procID());
        sprintf(buf,"emacs -title %s -display %s --eval '(progn (gdb \" gdb -q %s --pid %d\") )'",
                title, display, binaryName, proc);
        printf("%s\n",buf);
        f = popen(buf,"w");
      }
    else
      {
        sprintf(buf,"gdb -q  %s %d", binaryName, proc);
        f=popen(buf, "w");
        fputs("where\ndetach\nquit\n", f);
      }
    
    pclose(f); // child process stalls
    exit(0);
  } else {
    char rbuf[3];
    read(pfds[0], rbuf, 3);
  }
#endif
}

#ifdef CH_MPI
#include "SPMD.H"
#endif

// when running the code inside an actual debugger, the debugger signal handler
// will come in before this one, so the debugger should function correctly.
int registerDebugger()
{
  int rtn = 0;
#ifndef CH_DISABLE_SIGNALS
  signal(SIGSEGV, AttachDebugger);
  signal(SIGABRT, AttachDebugger);

#ifdef CH_MPI
  int r = setChomboMPIErrorHandler();
  rtn+=r;
#endif
#else
  rtn = 2;
#endif
  return rtn;
}

#ifdef CH_MPI

void MPI_Shutdown()
{
  MPI_Abort(MPI_COMM_WORLD, -4);
}


void mpierrorfunction(MPI_Comm* a_comm, int* a_error, ...)
{
#ifndef CH_DISABLE_SIGNALS
  char string[MPI_MAX_ERROR_STRING];
  int resultlen;

  MPI_Error_string(*a_error, string, &resultlen);
  //cout << string << "\n";
  printf("%s\n",string);

  // dumpmemoryatexit();

  AttachDebugger(4);

  MayDay::Warning("MPI error function intercepting naughty MPI_Abort call");
  abort();
  // it would be nice to have a stack dump at this point also.
  //MPI_Abort(*a_comm, *a_error);
#endif
}

int setChomboMPIErrorHandler()
{
  int rtn=1;
#ifndef CH_DISABLE_SIGNALS
  MPI_Errhandler handler;

  // Some MPI implementations are having trouble with this.  In particular,
  // halem at GSFC.  For now, i'm just going to not do this for machines
  // running OSF (which currently, is only halem for us) (ndk)
#ifndef CH_OSF1
  MPI_Errhandler_create(mpierrorfunction, &handler);
  rtn = MPI_Errhandler_set(Chombo_MPI::comm, handler);
  if(rtn == MPI_SUCCESS) rtn = 0;
#else
  rtn = 2;
#endif
#endif
  return rtn;
}
#endif

#include "NamespaceFooter.H"
