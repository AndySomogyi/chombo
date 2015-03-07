#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EBDebugOut.H"
#include "EBCellFAB.H"
#include "LayoutIterator.H"
#include "parstream.H"
#include "FaceIndex.H"
#include "Vector.H"
#include "VolIndex.H"
#include "LoHiSide.H"
#include "VoFIterator.H"
#include "DebugOut.H"
#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
#include <iomanip>
#include "NamespaceHeader.H"

void dumpEBFABIVS(const EBCellFAB* a_fab, const IntVectSet* a_ivs)
{
  const EBGraph& ebgraph = a_fab->getEBISBox().getEBGraph();
  const int ncomp = a_fab->nComp();
  for(VoFIterator vofit(*a_ivs, ebgraph); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      pout() << "vof= " << vof.gridIndex();

      pout() << ";      data=";
      for(int ivar = 0; ivar < ncomp; ivar++)
        {
          pout() << " "
                 << setprecision(8)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << (*a_fab)(vof, ivar);
         }
      pout() << endl;
    }
}
void dumpEBFAB(const EBCellFAB* a_fab)
{
  const EBCellFAB& fab = *a_fab;
  Box box = fab.getRegion();
  IntVectSet ivs(box);
  pout() << "valid and ghost data in ebcellfab" << endl;
  dumpEBFABIVS(a_fab, &ivs);
}
void dumpEBLevelIrreg(const LevelData<EBCellFAB>* a_level)
{
  const LevelData<EBCellFAB>& lev = *a_level;
  pout() << "valid irregular data in eblevel" << endl;
  for(DataIterator dit = lev.dataIterator(); dit.ok(); ++dit)
    {
      Box box = lev.disjointBoxLayout().get(dit());
      IntVectSet ivs = lev[dit()].getEBISBox().getIrregIVS(box);
      dumpEBFABIVS(&(lev[dit()]), &ivs); 
    }
}

void dumpEBLevel(const LevelData<EBCellFAB>* a_level)
{
  const LevelData<EBCellFAB>& lev = *a_level;
  pout() << "valid data in eblevel" << endl;
  for(DataIterator dit = lev.dataIterator(); dit.ok(); ++dit)
    {
      Box box = lev.disjointBoxLayout().get(dit());
      IntVectSet ivs(box);
      dumpEBFABIVS(&(lev[dit()]), &ivs);
    }
}

void dumpEBFABIrreg(const EBCellFAB* a_fab)
{
  const EBCellFAB& fab = *a_fab;
  const EBGraph& ebgraph = fab.getEBISBox().getEBGraph();
  Box box = fab.getRegion();
  IntVectSet ivs= ebgraph.getIrregCells(box);

  pout() << "valid and ghost irregular cell data in ebcellfab" << endl;
  dumpEBFABIVS(a_fab, &ivs);
}

void
getMaxEBFAB(const EBCellFAB*  ldptr)
{

  const EBCellFAB& fab = *ldptr;
  for(int icomp = 0; icomp < fab.nComp(); icomp++)
    {
      Real maxval = -1.0e16;
      Real minval =  1.0e16;
      VolIndex vofmax, vofmin;
      pout() << "c = " << icomp << ", ";

      const EBGraph&   ebg = fab.getEBISBox().getEBGraph();
      const Box& grid = fab.box();
      IntVectSet ivs(grid);
      for(VoFIterator vofit(ivs, ebg); vofit.ok(); ++vofit)
        {
          Real datval = fab(vofit(), icomp);
          if(datval > maxval)
            {
              maxval = datval;
              vofmax = vofit();
            }
          if(datval < minval)
            {
              minval = datval;
              vofmin = vofit();
            }
        }
      pout() << "max=" <<  maxval << " at " << vofmax << ", ";
      pout() << "min=" <<  minval << " at " << vofmin << endl;
    }
}

void
getMaxEBLevel(const LevelData<EBCellFAB>*  ldptr)
{

  const LevelData<EBCellFAB>& ld = *ldptr;
  for(int icomp = 0; icomp < ld.nComp(); icomp++)
    {
      Real maxval = -1.0e16;
      Real minval =  1.0e16;
      VolIndex vofmax, vofmin;
      pout() << "c = " << icomp << ", ";
      const DisjointBoxLayout& dbl = ld.disjointBoxLayout();
      for(DataIterator dit = ld.dataIterator(); dit.ok(); ++dit)
        {
          const EBCellFAB& fab = ld[dit()];
          const EBGraph&   ebg = fab.getEBISBox().getEBGraph();
          const Box& grid = dbl.get(dit());
          IntVectSet ivs(grid);
          for(VoFIterator vofit(ivs, ebg); vofit.ok(); ++vofit)
            {
              Real datval = ld[dit()](vofit(), icomp);
              if(datval > maxval)
                {
                  maxval = datval;
                  vofmax = vofit();
                }
              if(datval < minval)
                {
                  minval = datval;
                  vofmin = vofit();
                }
            }
        }
      pout() << "max=" <<  maxval << " at " << vofmax << ", ";
      pout() << "min=" <<  minval << " at " << vofmin << endl;
    }
}
void
dumpLDEBCF(const LevelData<EBCellFAB>*  ldptr)
{
  const LevelData<EBCellFAB>& ld = *ldptr;
  for(DataIterator dit = ld.dataIterator(); dit.ok(); ++dit)
    {
      dumpEBFAB(&(ld[dit()]));
    }
}

void
dumpEBLDDBL(const LevelData<EBCellFAB>*  memLDF_Ptr)
{
  const DisjointBoxLayout& dbl = memLDF_Ptr->disjointBoxLayout();
  dumpDBL(&dbl);
}

void dumpLDBIVF(const LayoutData< BaseIVFAB<Real> >* a_ldptr)
{
  const LayoutData< BaseIVFAB<Real> >& ld = *a_ldptr;
  for(DataIterator dit = ld.dataIterator(); dit.ok(); ++dit)
    {
      dumpIVFAB(&ld[dit()]);
    }
}

void dumpIVFAB(const BaseIVFAB<Real>* a_vectPtr)
{
  const BaseIVFAB<Real>& ivfab = *a_vectPtr;
  const int ncomp = ivfab.nComp();
  const EBGraph& ebgraph = ivfab.getEBGraph();
  const IntVectSet& ivs = ivfab.getIVS();

  pout() << "data in base ivfab" << endl;
  for(VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      pout() << "vof= " << vof.gridIndex();
      //RealVect centroid = ebgraph.centroid(vof);
      //pout() << ";      cent=";
      //for(int idir = 0; idir < SpaceDim; idir++)
      //  {
      //    pout() << " "
      //           << setprecision(8)
      //           << setiosflags(ios::showpoint)
      //           << setiosflags(ios::scientific)
      //           << centroid[idir];
      //  }
      pout() << ";      data=";
      for(int ivar = 0; ivar < ncomp; ivar++)
        {
          pout() << " "
                 << setprecision(8)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << ivfab(vof, ivar);
        }
      pout() << endl;
    }
}
void dumpIFFAB(const BaseIFFAB<Real>* a_vectPtr)
{
  const BaseIFFAB<Real>& iffab = *a_vectPtr;
  const int ncomp = iffab.nComp();
  const int dir = iffab.direction();
  const EBGraph& ebgraph = iffab.getEBGraph();
  const IntVectSet& ivs = iffab.getIVS();
  FaceIterator faceit(ivs, ebgraph, dir,
                      FaceStop::SurroundingWithBoundary);
  pout() << "data in base iffab" << endl;
  for(faceit.reset();faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      pout() << "face= ";
      for(SideIterator sit; sit.ok(); ++sit)
        {
          pout() << " " << face.gridIndex(sit());
        }
      for(int ivar = 0; ivar < ncomp; ivar++)
        {
          pout() << " "
                 << setprecision(8)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << iffab(face, ivar);
        }
      pout() << endl;
    }
}
void dumpVVoF(const Vector<VolIndex>* a_vectPtr)
{
  if(a_vectPtr == NULL) return;
  Vector<VolIndex> vect = *a_vectPtr;
  pout() << "vector of vofs contains" << endl;
  for(int iveco = 0; iveco < vect.size(); iveco++)
    {
      const VolIndex& vof = vect[iveco];
      pout() <<  vof.gridIndex() << "   "  << vof.cellIndex() << "  ";;
    }
  pout() << endl;
}

void dumpVFace(const Vector<FaceIndex>* a_vectPtr)
{
  if(a_vectPtr == NULL) return;
  Vector<FaceIndex> vect = *a_vectPtr;
  pout() << "vector of faces contains" << endl;
  for(int iveco = 0; iveco < vect.size(); iveco++)
    {
      for(SideIterator sit; sit.ok(); ++sit)
        {
          VolIndex vof = vect[iveco].getVoF(sit());
          pout() <<  vof.gridIndex() << "   " << vof.cellIndex() << "  ";;
        }
      if(vect[iveco].isBoundary())
        {
          pout() <<  "on boundary";
        }
      pout() << endl;
    }
}

void dumpFace(const FaceIndex* a_vectPtr)
{
  if(a_vectPtr == NULL) return;
  const FaceIndex& face = *a_vectPtr;
  pout() << "faces contains" << endl;
  for(SideIterator sit; sit.ok(); ++sit)
    {
      VolIndex vof = face.getVoF(sit());
      pout() <<  vof.gridIndex() << "   " << vof.cellIndex() << "  ";
    }
  pout() << endl;
  if(face.isBoundary())
    {
      pout() <<  " face is on on boundary";
      pout() << endl;
    }
}

void dumpFaceSten(const FaceStencil* a_vectPtr)
{
  if(a_vectPtr == NULL) return;
  const FaceStencil& sten = *a_vectPtr;
  pout() << "stencil contains:" << endl;
  for(int isten = 0; isten <  sten.size(); isten++)
    {
      const FaceIndex& face = sten.face(isten);
      const Real& weight = sten.weight(isten);
      pout() << "face: ";
      for(SideIterator sit; sit.ok(); ++sit)
        {
          VolIndex vof = face.getVoF(sit());
          pout() <<  vof.gridIndex() << "   " << vof.cellIndex() << "  ";;
        }
      pout() << "weight: " << weight;
      pout() << endl;
    }
}

void dumpVoFSten(const VoFStencil* a_vectPtr)
{
  if(a_vectPtr == NULL) return;
  const VoFStencil& sten = *a_vectPtr;
  pout() << "stencil contains:" << endl;
  for(int isten = 0; isten <  sten.size(); isten++)
    {
      const VolIndex& vof = sten.vof(isten);
      const Real& weight = sten.weight(isten);
      pout() << "vof: ";
      pout() <<  vof.gridIndex() << "   " << vof.cellIndex() << "  ";;
      pout() << "weight: " << weight;
      pout() << endl;
    }
}
#include "NamespaceFooter.H"
