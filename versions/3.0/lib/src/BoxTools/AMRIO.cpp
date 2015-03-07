#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// DTGraves, Fri, Dec 3, 1999

#include <fstream>

#include <string>
using std::fstream;
using std::string;

#include <cstdio>
#ifndef __IBMCPP__
using std::sprintf;
#endif
//using std::tempnam;

#include <cstdlib>
using std::system;

#include <cmath>

#ifdef CH_USE_HDF5
#include "CH_HDF5.H"
#endif

#include "AMRIO.H"
#include "BoxIterator.H"
#include "LayoutIterator.H"
#include "VisItChomboDriver.H"
#include "NamespaceHeader.H"

#ifdef CH_USE_HDF5
/*
\\ write out hierarchy of amr data in HDF5 format
\\ filename,  == file to output to
\\ a_vectData == data at each level
\\ a_vectNames== names of variables
\\ a_domain == domain at coarsest level
\\ a_dx     == grid spacing at coarsest level
\\ a_dt     == time step at coarsest level
\\ a_time     == time
\\ a_vectRatio == refinement ratio at all levels
\\ (ith entry is refinement ratio between levels i and i + 1)
\\ a_numLevels == number of levels to output
*/
void
WriteAMRHierarchyHDF5(const string& filename,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<LevelData<FArrayBox>* > & a_vectData,
                      const Vector<string>& a_vectNames,
                      const Box& a_domain,
                      const Real& a_dx,
                      const Real& a_dt,
                      const Real& a_time,
                      const Vector<int>& a_refRatio,
                      const int& a_numLevels)
{
  HDF5Handle handle(filename.c_str(),  HDF5Handle::CREATE);

  WriteAMRHierarchyHDF5(handle, a_vectGrids, a_vectData, a_vectNames,
                        a_domain, a_dx, a_dt, a_time, a_refRatio, a_numLevels);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  handle.close();
}

void
WriteAMRHierarchyHDF5(HDF5Handle& handle,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<LevelData<FArrayBox>* > & a_vectData,
                      const Vector<string>& a_vectNames,
                      const Box& a_domain,
                      const Real& a_dx,
                      const Real& a_dt,
                      const Real& a_time,
                      const Vector<int>& a_refRatio,
                      const int& a_numLevels)
{
  CH_assert(a_numLevels > 0);
  CH_assert(a_vectData.size()  >= a_numLevels);
  CH_assert(a_refRatio.size() >= a_numLevels-1);

  HDF5HeaderData header;
  int nComp = a_vectNames.size();

  string filedescriptor("VanillaAMRFileType");
  header.m_string ["filetype"]      = filedescriptor;
  header.m_int ["num_levels"]       = a_numLevels;
  header.m_int ["num_components"]    = nComp;

  for (int ivar = 0; ivar < nComp; ivar++)
    {
      char labelChSt[80];
      sprintf(labelChSt, "component_%d", ivar);
      string label(labelChSt);
      header.m_string[label] = a_vectNames[ivar];
    }
  header.writeToFile(handle);

  Box domainLevel = a_domain;
  Real dtLevel = a_dt;
  Real dxLevel = a_dx;
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      int refLevel = 1;
      if (ilev != a_numLevels -1)
        {
          refLevel = a_refRatio[ilev];
        }
      if (ilev != 0)
        {
          domainLevel.refine(a_refRatio[ilev-1]);
          dtLevel /= a_refRatio[ilev-1];
          dxLevel /= a_refRatio[ilev-1];
        }
      CH_assert(a_vectData[ilev] != NULL);
      const LevelData<FArrayBox>& dataLevel = *a_vectData[ilev];
      CH_assert(dataLevel.nComp() == nComp);
      Interval comps(0,nComp-1);
      IntVect ghostVect = a_vectData[0]->ghostVect();
      int eek = writeLevel(handle, ilev, dataLevel,
                           dxLevel, dtLevel, a_time,
                           domainLevel, refLevel, ghostVect, comps);
      if (eek != 0)
        {
          MayDay::Error("WriteAMRHierarchyHDF5: Error in writeLevel");
        }
    }
}

void
WriteAMRHierarchyHDF5(const string& filename,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<LevelData<FArrayBox>* > & a_vectData,
                      const Box& a_domain,
                      const Vector<int>& a_refRatio,
                      const int& a_numLevels)
{

  HDF5Handle handle(filename.c_str(),  HDF5Handle::CREATE);
  WriteAMRHierarchyHDF5(handle, a_vectGrids, a_vectData,
                        a_domain, a_refRatio, a_numLevels);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  handle.close();
}

void
WriteAMRHierarchyHDF5(HDF5Handle& handle,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<LevelData<FArrayBox>* > & a_vectData,
                      const Box& a_domain,
                      const Vector<int>& a_refRatio,
                      const int& a_numLevels)
{
  CH_assert(a_numLevels > 0);
  CH_assert(a_vectData.size()  >= a_numLevels);
  CH_assert(a_refRatio.size() >= a_numLevels-1);
  Real dxin = 1.0;
  Real dtin = 1.0;
  Real time = 1.0;

  HDF5HeaderData header;
  int nComp = a_vectData[0]->nComp();

  string filedescriptor("VanillaAMRFileType");
  header.m_string ["filetype"]      = filedescriptor;
  header.m_int ["num_levels"]       = a_numLevels;
  header.m_int ["num_components"]    = nComp;

  for (int ivar = 0; ivar < nComp; ivar++)
    {
      char labelChSt[80];
      sprintf(labelChSt, "component_%d", ivar);
      string label(labelChSt);
      header.m_string[label] = label;
    }
  header.writeToFile(handle);

  Box domainLevel = a_domain;
  Real dtLevel = dtin;
  Real dxLevel = dxin;
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      int refLevel = 1;
      if (ilev != a_numLevels -1)
        {
          refLevel = a_refRatio[ilev];
        }
      if (ilev != 0)
        {
          domainLevel.refine(a_refRatio[ilev-1]);
          dtLevel /= a_refRatio[ilev-1];
          dxLevel /= a_refRatio[ilev-1];
        }
      CH_assert(a_vectData[ilev] != NULL);
      const LevelData<FArrayBox>& dataLevel = *a_vectData[ilev];
      CH_assert(dataLevel.nComp() == nComp);
      Interval comps(0,nComp-1);
      IntVect ghostVect = a_vectData[0]->ghostVect();
      int eek = writeLevel(handle, ilev, dataLevel,
                           dxLevel, dtLevel, time,
                           domainLevel, refLevel, ghostVect, comps);
      if (eek != 0)
        {
          MayDay::Error("WriteAMRHierarchyHDF5: Error in writeLevel");
        }
    }
}

//
/*
\\ Read in hierarchy of amr data in HDF5 format
\\ filename,  == file to output to
\\ a_vectData == data at each level
\\ a_vectNames== names of variables
\\ a_domain == domain at coarsest level
\\ a_dx     == grid spacing at coarsest level
\\ a_dt     == time step at coarsest level
\\ a_time     == time
\\ a_vectRatio == refinement ratio at all levels
\\ (ith entry is refinement ratio between levels i and i + 1)
\\ a_numLevels == number of levels to output

return values:
0: success
-1: bogus number of levels
-2: bogus number of components
-3: error in readlevel
-4: file open failed
*/
int
ReadAMRHierarchyHDF5(const string& filename,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Vector<string>& a_vectNames,
                     Box& a_domain,
                     Real& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     Vector<int>& a_refRatio,
                     int& a_numLevels)
{
  HDF5Handle handle;
  int err = handle.open(filename.c_str(),  HDF5Handle::OPEN_RDONLY);
  if ( err < 0)
    {
      return -4;
    }
  int eekflag = ReadAMRHierarchyHDF5(handle, a_vectGrids, a_vectData,
                                     a_vectNames, a_domain, a_dx, a_dt,
                                     a_time, a_refRatio, a_numLevels);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  handle.close();

  return (eekflag);
}

int
ReadAMRHierarchyHDF5(HDF5Handle& handle,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Vector<string>& a_vectNames,
                     Box& a_domain,
                     Real& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     Vector<int>& a_refRatio,
                     int& a_numLevels)
{

  HDF5HeaderData header;
  header.readFromFile(handle);

  a_numLevels = header.m_int["num_levels"];
  if (a_numLevels <= 0)
  {
    MayDay::Warning("ReadAMRHierarchyHDF5: Bogus number of levels");
    return (-1);
  }
  a_vectData.resize(a_numLevels);
  a_refRatio.resize(a_numLevels);
  a_vectGrids.resize(a_numLevels);

  int nComp = header.m_int["num_components"];
  if (nComp <= 0)
  {
    MayDay::Warning("ReadAMRHierarchyHDF5: Bogus number of Components");
    return (-2);
  }
  a_vectNames.resize(nComp);

  for (int ivar = 0; ivar < nComp; ivar++)
    {
      char labelChSt[80];
      sprintf(labelChSt, "component_%d", ivar);
      string label(labelChSt);
      a_vectNames[ivar] = header.m_string[label];
    }
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      int refLevel = 0;
      Box domainLevel;
      Real dtLevel;
      Real dxLevel;
      a_vectData[ilev] = new LevelData<FArrayBox>();
      int eek = readLevel(handle, ilev, *(a_vectData[ilev]),
                          dxLevel, dtLevel,  a_time,
                          domainLevel, refLevel, Interval(), true);
      if (eek != 0)
      {
        MayDay::Warning("ReadAMRHierarchyHDF5: readLevel failed");
        return (-3);
      }

      const DisjointBoxLayout& dbl = a_vectData[ilev]->getBoxes();
      a_vectGrids[ilev]= dbl;

      if (ilev == 0)
        {
          a_domain = domainLevel;
          a_dt = dtLevel;
          a_dx = dxLevel;
        }
      a_refRatio[ilev] = refLevel;
    }
  return (0);
}

int
ReadAMRHierarchyHDF5(const string& filename,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Box& a_domain,
                     Vector<int>& a_refRatio,
                     int& a_numLevels)
{
  HDF5Handle handle;
  int err = handle.open(filename.c_str(),  HDF5Handle::OPEN_RDONLY);
  if ( err < 0)
  {
    return -4;
  }

  int eekflag = ReadAMRHierarchyHDF5(handle, a_vectGrids, a_vectData,
                                     a_domain, a_refRatio, a_numLevels);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  handle.close();
  return (eekflag);
}

int
ReadAMRHierarchyHDF5(HDF5Handle& handle,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Box& a_domain,
                     Vector<int>& a_refRatio,
                     int& a_numLevels)
{
  HDF5HeaderData header;
  header.readFromFile(handle);

  a_numLevels = header.m_int["num_levels"];
  if (a_numLevels <= 0)
  {
    MayDay::Warning("ReadAMRHierarchyHDF5: Bogus number of levels");
    return (-1);
  }
  a_vectData.resize(a_numLevels);
  a_refRatio.resize(a_numLevels);
  a_vectGrids.resize(a_numLevels);

  //  int nComp = header.m_int["num_components"];
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      int refLevel = 0;
      Box domainLevel;
      Real dtLevel;
      Real dxLevel;
      Real time;
      a_vectData[ilev] = new LevelData<FArrayBox>();
      int eek = readLevel(handle, ilev, *(a_vectData[ilev]),
                          dxLevel, dtLevel,  time,
                          domainLevel, refLevel, Interval(), true);
      if (eek != 0)
      {
        MayDay::Warning("ReadAMRHierarchyHDF5: readLevel failed");
        return (-3);
      }

      const DisjointBoxLayout& dbl = a_vectData[ilev]->getBoxes();
      a_vectGrids[ilev]= dbl;

      if (ilev == 0)
        {
          a_domain = domainLevel;
        }
      a_refRatio[ilev] = refLevel;
    }

  return (0);
}

// put simple debugging functions here, at the end of the hdf5 stuff

// put simple debugging functions here, at the end of the hdf5 stuff
static void
ChomboVisVisualizeFile(const char *fname)
{
  char command[2000];
  sprintf(command,"$CHOMBOVIS_HOME/bin/chombovis debug_level=0 %s &",fname);
  system(command);
}

static void
ChomboBrowserBrowseFile(const char *fname)
{
  char command[2000];
  sprintf(command,"$CHOMBOVIS_HOME/bin/chombobrowser debug_level=0 %s &",fname);
  system(command);
}

static VisItChomboDriver visit;
static void
VisItVisualizeFile(const char *fname)
{
  visit.VisualizeFile(fname);
}

static void
VisItBrowseFile(const char *fname)
{
  visit.BrowseFile(fname);
}


static void
VisualizeFile(const char *fname)
{
  const char *use_visit = getenv("CHOMBO_USE_VISIT");
  if (use_visit)
      VisItVisualizeFile(fname);
  else
      ChomboVisVisualizeFile(fname);
}

static void
BrowseFile(const char *fname)
{
  const char *use_visit = getenv("CHOMBO_USE_VISIT");
  if (use_visit)
      VisItBrowseFile(fname);
  else
      ChomboBrowserBrowseFile(fname);
}

void
writeFAB(const FArrayBox* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const char* fname = "fab.hdf5";
  writeFABname(a_dataPtr, fname);
}

void
viewFAB(const FArrayBox* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const char* fname = tempnam(NULL,NULL);
  writeFABname(a_dataPtr, fname);
  VisualizeFile(fname);
}

void
viewBFI(const BaseFab<int>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  FArrayBox fab(a_dataPtr->box(), a_dataPtr->nComp());
  BoxIterator bit(fab.box());
  for (bit.begin(); bit.ok(); bit.next())
    {
      const IntVect& iv = bit();
      for (int ivar = 0; ivar < fab.nComp(); ivar++)
        {
          fab(iv, ivar) = a_dataPtr->operator()(iv, ivar);
        }
    }

  const char* fname = tempnam(NULL,NULL);
  writeFABname(&fab, fname);
  VisualizeFile(fname);
}

void
viewBFIV(const BaseFab<IntVect>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  FArrayBox fab(a_dataPtr->box(), SpaceDim * a_dataPtr->nComp());
  BoxIterator bit(fab.box());
  for (bit.begin(); bit.ok(); bit.next())
    {
      const IntVect& iv = bit();
      for (int ivar = 0; ivar < a_dataPtr->nComp(); ivar++)
        {
          const IntVect& ivVal = a_dataPtr->operator()(iv, ivar);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              fab(iv, SpaceDim*ivar + idir) = ivVal[idir];
            }
        }
    }

  const char* fname = tempnam(NULL,NULL);
  writeFABname(&fab, fname);
  VisualizeFile(fname);
}

void
viewBFRV(const BaseFab<RealVect>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  FArrayBox fab(a_dataPtr->box(), SpaceDim * a_dataPtr->nComp());
  BoxIterator bit(fab.box());
  for (bit.begin(); bit.ok(); bit.next())
    {
      const IntVect& iv = bit();
      for (int ivar = 0; ivar < a_dataPtr->nComp(); ivar++)
        {
          const RealVect& rvVal = a_dataPtr->operator()(iv, ivar);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              fab(iv, SpaceDim*ivar + idir) = rvVal[idir];
            }
        }
    }

  const char* fname = tempnam(NULL,NULL);
  writeFABname(&fab, fname);
  VisualizeFile(fname);
}

void
browseFAB(const FArrayBox* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const char* fname = tempnam(NULL,NULL);
  writeFABname(a_dataPtr, fname);
  BrowseFile(fname);

}

void
writeBFR(const BaseFab<Real>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  //copy the BaseFab over to a FArrayBox
  Box fabBox = (*a_dataPtr).box();
  int  ncomp = (*a_dataPtr).nComp();
  FArrayBox fab(fabBox,ncomp);
  fab.copy(*a_dataPtr);

  const char* fname = "fab.hdf5";
  writeFABname(&fab, fname);
}

void
viewBFR(const BaseFab<Real>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  //copy the BaseFab over to a FArrayBox
  Box fabBox = (*a_dataPtr).box();
  int  ncomp = (*a_dataPtr).nComp();
  FArrayBox fab(fabBox,ncomp);
  fab.copy(*a_dataPtr);

  const char* fname = tempnam(NULL,NULL);
  writeFABname(&fab, fname);
  VisualizeFile(fname);

}

void
writeFABname(const FArrayBox      * a_dataPtr,
             const char           * a_filename,
             const Vector<string> & a_compNames)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const FArrayBox& data = *a_dataPtr;

  HDF5Handle handle(a_filename, HDF5Handle::CREATE);
  HDF5HeaderData header;

  int numlevels= 1;
  int nComp = data.nComp();

  string filedescriptor("VanillaAMRFileType");
  header.m_string ["filetype"]      = filedescriptor;
  header.m_int ["num_levels"]       = numlevels;
  header.m_int ["num_components"]    = nComp;

  for (int ivar = 0; ivar < nComp; ivar++)
    {
      char labelChSt[80];
      sprintf(labelChSt, "component_%d", ivar);
      string label(labelChSt);
      if ( a_compNames.size() > ivar )
        {
          header.m_string[label] = a_compNames[ivar] ;
        }
      else
        {
          header.m_string[label] = label;
        }
    }
  header.writeToFile(handle);

  Box domainLevel = data.box();
  // put bogus numbers here
  Real dtLevel = 1.0;
  Real dxLevel = 1.0;
  Real time = 1.0;

  int refLevel = 1;

  // build bogus DisjointBoxLayout here
  Vector<Box> boxes(1,domainLevel);
  unsigned int myprocID= procID();
  Vector<int> procAssign(1,myprocID);
  DisjointBoxLayout grids(boxes, procAssign);
  LevelData<FArrayBox> ldf(grids, nComp);
  // now copy fab into ldf
  DataIterator dit = ldf.dataIterator();
  ldf[dit()].copy(data);

  int eek = writeLevel(handle, 0, ldf, dxLevel, dtLevel, time,
                       domainLevel, refLevel);
  if (eek != 0)
  {
    MayDay::Error("writeFABname: error in writeLEvel");
  }

  handle.close();
}

void
viewLevelNoFine(const LevelData<FArrayBox>* a_dataPtr,
                const LevelData<FArrayBox>* a_dataFinePtr,
                int a_refRatio)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  if (a_dataFinePtr == NULL)
  {
    return;
  }

  int nVar = a_dataPtr->nComp();
  const DisjointBoxLayout& layoutOrig = a_dataPtr->disjointBoxLayout();
  const DisjointBoxLayout& layoutFine = a_dataFinePtr->disjointBoxLayout();
  DisjointBoxLayout layoutFineCoarsened;
  coarsen(layoutFineCoarsened, layoutFine, a_refRatio);
  LevelData<FArrayBox>* dataFineCoarsenedPtr =
    new LevelData<FArrayBox>(layoutFineCoarsened, nVar);
  for (DataIterator dit = a_dataFinePtr->dataIterator(); dit.ok(); ++dit)
    {
      dataFineCoarsenedPtr->operator[](dit).setVal(0.);
    }
  LevelData<FArrayBox>* dataNoFinePtr =
    new LevelData<FArrayBox>(layoutOrig, nVar);
  // *dataNoFinePtr := *a_dataPtr
  a_dataPtr->copyTo(*dataNoFinePtr);
  // Zero out *dataNoFinePtr on coarsened grids of *a_dataFinePtr.
  dataFineCoarsenedPtr->copyTo(*dataNoFinePtr);
  delete dataFineCoarsenedPtr;

  viewLevel(dataNoFinePtr);
  delete dataNoFinePtr;
}

void
writeLevel(const LevelData<FArrayBox>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const char* fname = "LDF.hdf5";
  writeLevelname(a_dataPtr, fname);
}

void
viewLevel(const LevelData<FArrayBox>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const char* fname = tempnam(NULL,NULL);
  writeLevelname(a_dataPtr, fname);
  VisualizeFile(fname);

}

void
browseLevel(const LevelData<FArrayBox>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const char* fname = tempnam(NULL,NULL);
  writeLevelname(a_dataPtr, fname);
  BrowseFile(fname);

}

void
viewLevelNoGhost(const LevelData<FArrayBox>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const DisjointBoxLayout& grids = a_dataPtr->getBoxes();
  int nVar = a_dataPtr->nComp();

  LevelData<FArrayBox> temp(grids, nVar);
  a_dataPtr->copyTo(a_dataPtr->interval(), temp, temp.interval());

  const char* fname = tempnam(NULL,NULL);
  writeLevelname(&temp, fname);

  VisualizeFile(fname);

}

void
writeLevelname(const LevelData<FArrayBox>* a_dataPtr,
               const char*                 a_filename)
{
  Vector<LevelData<FArrayBox>*> data(1);
  data[0]=(LevelData<FArrayBox>*)a_dataPtr;
  Vector<int> refRatios(1,1);
  writeVectorLevelName(&data, &refRatios, a_filename);
}

void
writeVectorLevelName(const Vector<LevelData<FArrayBox>*>* a_dataPtr,
                     const Vector<int>*          a_refRatios,
                     const char*                 a_filename)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  HDF5Handle handle(a_filename, HDF5Handle::CREATE);
  HDF5HeaderData header;

  int numlevels = a_dataPtr->size();
  int nComp =(a_dataPtr->operator[](0))->nComp();

  string filedescriptor("VanillaAMRFileType");
  header.m_string ["filetype"]      = filedescriptor;
  header.m_int ["num_levels"]       = numlevels;
  header.m_int ["num_components"]    = nComp;

  for (int ivar = 0; ivar < nComp; ivar++)
    {
      char labelChSt[80];
      sprintf(labelChSt, "component_%d", ivar);
      string label(labelChSt);
      header.m_string[label] = label;
    }
  header.writeToFile(handle);

  // put bogus numbers here
  Real dtLevel = 1.0;
  Real dxLevel = 1.0;
  Real time = 1.0;

  for (int level=0; level<a_dataPtr->size(); level++)
    {
      const LevelData<FArrayBox>& data = *(a_dataPtr->operator[](level));

      // need to figure out what domain will contain this LevelData
      // This must be LayoutIterator instead of DataIterator because
      // we need domain over boxes in ALL procs.
      const DisjointBoxLayout& levelBoxes = data.getBoxes();
      LayoutIterator lit = levelBoxes.layoutIterator();
      lit.reset();
      Box domain = levelBoxes.get(lit());
      for (lit.reset(); lit.ok(); ++lit)
        {
          const Box thisBox = levelBoxes.get(lit());
          D_TERM6(
                  if (thisBox.smallEnd(0)<domain.smallEnd(0))
                    domain.setSmall(0,thisBox.smallEnd(0)); ,
                  if (thisBox.smallEnd(1)<domain.smallEnd(1))
                    domain.setSmall(1,thisBox.smallEnd(1)); ,
                  if (thisBox.smallEnd(2)<domain.smallEnd(2))
                    domain.setSmall(2, thisBox.smallEnd(2)); ,
                  if (thisBox.smallEnd(3)<domain.smallEnd(3))
                    domain.setSmall(3,thisBox.smallEnd(3)); ,
                  if (thisBox.smallEnd(4)<domain.smallEnd(4))
                    domain.setSmall(4,thisBox.smallEnd(4)); ,
                  if (thisBox.smallEnd(5)<domain.smallEnd(5))
                    domain.setSmall(5, thisBox.smallEnd(5)); );

          D_TERM6(
                  if (thisBox.bigEnd(0)>domain.bigEnd(0))
                    domain.setBig(0,thisBox.bigEnd(0)); ,
                  if (thisBox.bigEnd(1)>domain.bigEnd(1))
                    domain.setBig(1,thisBox.bigEnd(1)); ,
                  if (thisBox.bigEnd(2)>domain.bigEnd(2))
                    domain.setBig(2, thisBox.bigEnd(2)); ,
                  if (thisBox.bigEnd(3)>domain.bigEnd(3))
                    domain.setBig(3,thisBox.bigEnd(3)); ,
                  if (thisBox.bigEnd(4)>domain.bigEnd(4))
                    domain.setBig(4,thisBox.bigEnd(4)); ,
                  if (thisBox.bigEnd(5)>domain.bigEnd(5))
                    domain.setBig(5, thisBox.bigEnd(5)); );

        } // end loop over boxes on level to determine "domain"

      int refLevel = 1;
      if (level < a_dataPtr->size()-1)
        {
          refLevel = a_refRatios->operator[](level);
        }

      Interval comps(0,nComp-1);
      IntVect ghostVect = data.ghostVect();
      int eek = writeLevel(handle, level, data, dxLevel, dtLevel, time,
                           domain, refLevel, ghostVect, comps);

      dtLevel /= refLevel;
      dxLevel /= refLevel;

      if (eek != 0)
      {
        MayDay::Error("writeLDFname: error in writeLEvel");
      }
    }
  handle.close();
}

void
writeDBL(const DisjointBoxLayout* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const char* fname = "DBL.hdf5";
  writeDBLname(a_dataPtr, fname);
}

void
viewDBL(const DisjointBoxLayout* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const char* fname = tempnam(NULL,NULL);
  writeDBLname(a_dataPtr, fname);

  VisualizeFile(fname);

}

void writeCopier(const Copier* a_copier)
{
  Vector<Vector<Box> > boxes;
  Vector<Vector<int> > procs;
  Vector<Box> level;
  Vector<int> p;
  Vector<int> refRatio;
  Box domain;
  int last = -1;
  for (CopyIterator it(*a_copier, CopyIterator::LOCAL); it.ok(); ++it)
    {
      const MotionItem& item = it();
      const DataIndex& to = item.toIndex;
      if (last == -1)
        {
          last = to.intCode();
        }
      if (last != to.intCode())
        {
          if (level.size()==0)
            {
              break;
            }
          boxes.push_back(level);
          procs.push_back(p);
          level.resize(0);
          p.resize(0);
          refRatio.push_back(1);
          last = to.intCode();
        }
      p.push_back(0);
      level.push_back(item.toRegion);
      domain.minBox(item.toRegion);
    }

  int numLevels = boxes.size();
  Vector<DisjointBoxLayout> layouts(numLevels);
  Vector<LevelData<FArrayBox>*> ldf(numLevels);
  IntVect ghostVect(IntVect::Zero);
  for (int i=0; i<numLevels; ++i)
    {
      layouts[i].define(boxes[i], procs[i], ProblemDomain(domain));
      ldf[i] = new LevelData<FArrayBox>(layouts[i], 1, ghostVect);
    }
  WriteAMRHierarchyHDF5("copier.hdf5", layouts, ldf, domain, refRatio, numLevels);
  for (int i=0; i<numLevels; ++i)
    {
      delete ldf[i];
    }
}

void viewCopier(const Copier* a_copier)
{
  if (a_copier == NULL)
  {
    return;
  }

  writeCopier(a_copier);

  const char* fname = "copier.hdf5";
  VisualizeFile(fname);
}

void viewIVS(const IntVectSet* a_dataPtr, const Box* a_domain)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  Vector<Box> boxes = a_dataPtr->boxes();
  Vector<int> procs(boxes.size(), 0);
  DisjointBoxLayout dbl(boxes, procs);
  viewDBL(&dbl);
}

void viewVectorBox(const Vector<Box>* boxes, const Box* a_domain)
{
  if (boxes == NULL)
  {
    return;
  }

  Vector<int> procs(boxes->size(), 0);
  DisjointBoxLayout dbl(*boxes, procs);
  viewDBL(&dbl);
}


void
writeDBLname(const DisjointBoxLayout* a_dataPtr,
             const char*                 a_filename)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const DisjointBoxLayout& grids = *a_dataPtr;

  // now create a LevelData<FArrayBox> based on the grids
  IntVect ghostVect(IntVect::Zero);
  LevelData<FArrayBox> data(grids, 1, ghostVect);

  // initialize data to procID.
  Real dataVal = procID();
  DataIterator dit = data.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      data[dit()].setVal(dataVal);
    }

  writeLevelname(&data, a_filename);
}

void
WritePartialAMRHierarchyHDF5(const string& filename,
                             const Vector<DisjointBoxLayout>& a_vectGrids,
                             const Vector<LevelData<FArrayBox>* > & a_vectData,
                             const Vector<string>& a_vectNames,
                             const Box& a_baseDomain,
                             const Real& a_baseDx,
                             const Real& a_dt,
                             const Real& a_time,
                             const Vector<int>& a_vectRatio,
                             const Interval& a_levels)
{
  int numLevels = a_levels.size();

  // now make new dataholders which only have numLevels entries,
  // and which will move the baseLevel to level 0

  Vector<DisjointBoxLayout> newVectGrids(numLevels);
  Vector<LevelData<FArrayBox> * > newVectData(numLevels);
  Vector<int> newVectRatio(numLevels);

  int leveloffset = a_levels.begin();
  for (int srcLevel = a_levels.begin(); srcLevel <= a_levels.end(); srcLevel++)
    {
      int destLevel = srcLevel-leveloffset;
      newVectGrids[destLevel] = a_vectGrids[srcLevel];
      newVectData[destLevel] = a_vectData[srcLevel];
      newVectRatio[destLevel] = a_vectRatio[srcLevel];
    }

  WriteAMRHierarchyHDF5(filename, newVectGrids, newVectData, a_vectNames,
                        a_baseDomain, a_baseDx, a_dt, a_time,
                        newVectRatio, numLevels);
}

#endif // CH_USE_HDF5

void writeAFabASCII(ostream& os, const FArrayBox& fab)
{
    os << fab.box() << "  " << fab.nComp() << endl;;
    BoxIterator bit(fab.box());
    for (bit.begin(); bit.ok(); bit.next())
    {
        for (int ivar = 0; ivar < fab.nComp(); ivar++)
        {
            os << fab(bit(), ivar) << endl;
        }
    }

}
void readAFabASCII(istream& is, FArrayBox& fab)
{
    Box fabbox;
    int ncomp;
    is >> fabbox >> ncomp;
    fab.resize(fabbox, ncomp);
    BoxIterator bit(fabbox);
    for (bit.begin(); bit.ok(); bit.next())
    {
        for (int ivar = 0; ivar < fab.nComp(); ivar++)
        {
            Real temp;
            is >> temp;
            fab(bit(), ivar) = temp;
        }
    }

}

//this outputs stuff to the amrascii format very slowly
/*
\\ write out hierarchy of amr data in HDF5 format
\\ filename,  == file to output to
\\ a_vectData == data at each level
\\ a_vectNames== names of variables
\\ a_domain == domain at coarsest level
\\ a_dx     == grid spacing at coarsest level
\\ a_dt     == time step at coarsest level
\\ a_time     == time
\\ a_refRatio == refinement ratio at all levels
\\ (ith entry is refinement ratio between levels i and i + 1)
\\ a_numLevel s== number of levels to output
*/
void
WriteAMRHierarchyASCII(const string& filename,
                       const Vector<DisjointBoxLayout>& a_vectGrids,
                       const Vector<LevelData<FArrayBox>* > & a_vectData,
                       const Vector<string>& a_vectNames,
                       const Box& a_domain,
                       const Real& a_dx,
                       const Real& a_dt,
                       const Real& a_time,
                       const Vector<int>& a_refRatio,
                       const int& a_numLevels)
{
  CH_assert(a_numLevels > 0);
  CH_assert(a_vectData.size()  >= a_numLevels);
  CH_assert(a_refRatio.size() >= a_numLevels-1);

  //sanity checks
  int nvar = a_vectData[0]->nComp();
  for (int ivec = 0; ivec < a_vectData.size(); ivec++)
    CH_assert(nvar == a_vectData[ivec]->nComp());
  CH_assert(a_vectNames.size() == nvar);
  for (int ivec = 0; ivec < a_refRatio.size(); ivec++)
    CH_assert(a_refRatio[ivec] > 0);

  fstream os;
  os.open(filename.c_str(), ios::out);

  os <<  nvar           << endl; // Number of states dumped
  for (int iname = 0; iname < a_vectNames.size(); iname++)
    {
      os << a_vectNames[iname]  << endl; // Name of state
    }

  os << SpaceDim     << endl; // Dimension of data
  os << a_numLevels  << endl; // number of levels dumped

  for (int irat = 0; irat < a_numLevels-1; irat++)
    os << a_refRatio[irat] << " " ;          //refinement ratios
  os << endl;
  os << a_domain << endl;                   //coarse level domain
  os << a_dx     << endl;                   //coarse level grid spacing
  os << a_dt     << endl;                   //coarse level time step
  os << a_time   << endl;                   // Simulation time of dump

  for (int ilevel = 0; ilevel < a_numLevels; ilevel++)
    {
      os << ilevel << endl;  //level number
      const LevelData<FArrayBox>& mf_lev = *(a_vectData[ilevel]);
      const DisjointBoxLayout& bxa = a_vectGrids[ilevel];

      // dump actual boxes
      LayoutIterator lit = bxa.layoutIterator();
      os << bxa.size() << endl;              // number of boxes
      for (lit.reset(); lit.ok(); ++lit)
        os << bxa[lit()] << endl;              // For each grid, dump box

      // dump actual data
      DataIterator dit = mf_lev.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
        {
          const FArrayBox& dat = mf_lev[dit()];
          writeAFabASCII(os, dat);
        }
    }
  os.close();
}

//this reads stuff in  amrascii format very slowly
/*
\\ write out hierarchy of amr data in HDF5 format
\\ filename,  == file to output to
\\ a_vectData == data at each level
\\ a_vectNames== names of variables
\\ a_domain == domain at coarsest level
\\ a_dx     == grid spacing at coarsest level
\\ a_dt     == time step at coarsest level
\\ a_time     == time
\\ a_refRatio == refinement ratio at all levels
\\ (ith entry is refinement ratio between levels i and i + 1)
\\ a_numLevels == number of levels to output

return values:
0 = success
-1 = failure
*/
int
ReadAMRHierarchyASCII(const string& filename,
                      Vector<DisjointBoxLayout>& a_vectGrids,
                      Vector<LevelData<FArrayBox>* > & a_vectData,
                      Vector<string>& a_vectNames,
                      Box& a_domain,
                      Real& a_dx,
                      Real& a_dt,
                      Real& a_time,
                      Vector<int>& a_refRatio,
                      int& a_numLevels,
                      const IntVect& ghostVector)
{
  fstream is;
  if (numProc() > 1)
    {
      MayDay::Warning("ReadAMRHierarchyASCII is written for one processor");
    }
  is.open(filename.c_str(), ios::in);

  int nvar;
  is  >>  nvar;             // number of states
  a_vectNames.resize(nvar);
  for (int iname = 0; iname < a_vectNames.size(); iname++)
    {
      is >> a_vectNames[iname];   // Name of state
    }

  int dimdat;
  is >> dimdat; //dimension of problem
  if (dimdat != SpaceDim)
    {
      MayDay::Warning("ReadAMRHierarchyASCII: dimension of data != dimension of reader");
      return (-1);
    }
  is >> a_numLevels; // number of levels dumped

  a_refRatio.resize(a_numLevels);
  a_vectGrids.resize(a_numLevels);
  a_vectData.resize(a_numLevels);
  for (int irat = 0; irat < a_numLevels-1; irat++)
    is >> a_refRatio[irat];          //refinement ratios

  is >> a_domain;                //coarse level domain
  is >> a_dx;                   //coarse level grid spacing
  is >> a_dt;                   //coarse level time step
  is >> a_time;                 // Simulation time of dump

  for (int ilevel = 0; ilevel < a_numLevels; ilevel++)
    {
      int itemp;
      is  >> itemp;  //level number
      if (itemp != ilevel)
        {
          MayDay::Warning("ReadAMRHierarchyASCII: Input file wrong");
          return (-1);
        }

      DisjointBoxLayout& dbl = a_vectGrids[ilevel];
      int nboxes;
      is >> nboxes;
      for (int ibox = 0; ibox < nboxes; ibox++)
        {
          Box btemp;
          is >> btemp;
          dbl.addBox(btemp, 0);
        }
      dbl.close();
      LevelData<FArrayBox> mf_lev(dbl,nvar,ghostVector);
      DataIterator dit = mf_lev.dataIterator();
      // read actual data
      for (dit.reset(); dit.ok(); ++dit)
        {
          FArrayBox datnoGhost;
          readAFabASCII(is, datnoGhost);
          FArrayBox& dat = mf_lev[dit()];
          dat.copy(datnoGhost);
        }
      a_vectData[ilevel] = new LevelData<FArrayBox>();
      a_vectData[ilevel]->define(mf_lev);
    }
  is.close();

  return (0);
}

// Read in hierarchy of amr data in ASCII format
int
ReadEBAMRASCII(const string& filename,
               Vector<DisjointBoxLayout>& a_vectGrids,
               Vector<LevelData<FArrayBox>* > & a_vectData,
               Vector<string>& a_vectNames,
               Box& a_domain,
               Vector<IntVectSet>& a_coveredCells,
               Vector<IntVectSet>& a_multiValuedCells,
               Vector<int>& a_refRatio,
               int& a_numLevels,
               const IntVect& a_ghostVector)
{
  fstream is;
  if (numProc() > 1)
    {
      MayDay::Warning("ReadAMRHierarchyASCII is written for one processor");
    }
  is.open(filename.c_str(), ios::in);

  int nvar;
  is  >>  nvar;             // number of states
  a_vectNames.resize(nvar);
  for (int iname = 0; iname < a_vectNames.size(); iname++)
    {
      is >> a_vectNames[iname];   // Name of state
    }

  int dimdat;
  is >> dimdat; //dimension of problem
  if (dimdat != SpaceDim)
    {
      MayDay::Warning("ReadAMRHierarchyASCII: dimension of data != dimension of reader");
      return (-1);
    }
  is >> a_numLevels; // number of levels dumped

  a_coveredCells.resize(a_numLevels);
  a_multiValuedCells.resize(a_numLevels);
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      int ivecsize;
      is >>  ivecsize; // number of covered boxes
      Box coveredBox;
      IntVectSet& coveredCells = a_coveredCells[ilev];
      coveredCells.makeEmpty();
      for (int ibox = 0; ibox < ivecsize; ibox++)
        {
          is >>  coveredBox; //input each covered box
          CH_assert(!coveredBox.isEmpty());
          coveredCells |= coveredBox;
        }

      IntVectSet& multiIVS = a_multiValuedCells[ilev];
      multiIVS.makeEmpty();
      is >> ivecsize ; //number of multi-valued cells
      for (int ivec = 0; ivec < ivecsize; ivec++)
        {
          IntVect iv;
          is >> iv;
          multiIVS |= iv;
        }
    }

  a_refRatio.resize(a_numLevels);
  a_vectGrids.resize(a_numLevels);
  a_vectData.resize(a_numLevels);
  for (int irat = 0; irat < a_numLevels-1; irat++)
    is >> a_refRatio[irat];          //refinement ratios

  is >> a_domain;                //coarse level domain
  for (int ilevel = 0; ilevel < a_numLevels; ilevel++)
    {
      int itemp;
      is  >> itemp;  //level number
      if (itemp != ilevel)
        {
          MayDay::Warning("ReadAMRHierarchyASCII: Input file wrong");
          return (-1);
        }

      DisjointBoxLayout& dbl = a_vectGrids[ilevel];
      int nboxes;
      is >> nboxes;
      for (int ibox = 0; ibox < nboxes; ibox++)
        {
          Box btemp;
          is >> btemp;
          dbl.addBox(btemp, 0);
        }
      dbl.close();
      LevelData<FArrayBox> mf_lev(dbl,nvar,a_ghostVector);
      DataIterator dit = mf_lev.dataIterator();
      // read actual data
      for (dit.reset(); dit.ok(); ++dit)
        {
          FArrayBox datnoGhost;
          //fab was written out as a farraybox.  no need
          //for special read function
          readAFabASCII(is, datnoGhost);
          FArrayBox& dat = mf_lev[dit()];
          dat.copy(datnoGhost);
        }
      a_vectData[ilevel] = new LevelData<FArrayBox>();
      a_vectData[ilevel]->define(mf_lev);
    }
  is.close();
  return (0);
}
#include "NamespaceFooter.H"
