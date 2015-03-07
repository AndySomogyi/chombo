#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "DebugOut.H"
#include "LayoutIterator.H"
#include "IntVectSet.H"
#include "NamespaceHeader.H"

static const char* indent2 = "      " ;

void dumpNodeLDFPar(const LevelData<NodeFArrayBox>* a_ldfabPtr)
{
  const LevelData<NodeFArrayBox>& ldfab = *a_ldfabPtr;
  Vector<Box> boxes;

  const DisjointBoxLayout& memDBL = ldfab.getBoxes();
  LayoutIterator lit  = memDBL.layoutIterator();
  for (lit.reset(); lit.ok(); ++lit)
  {
    boxes.push_back(memDBL.get(lit()));
  }

  Vector<int> assign(boxes.size(), 0);

  DisjointBoxLayout locDBL(boxes, assign);
  locDBL.close();

  LevelData<NodeFArrayBox> locLDF(locDBL,
                                  ldfab.nComp(),
                                  ldfab.ghostVect());

  const Interval& interv = ldfab.interval();
  ldfab.copyTo(interv,
               locLDF,
               interv);

  DataIterator dit = locLDF.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    const FArrayBox& fab = locLDF[dit()].getFab();
    Box nodeGrid = surroundingNodes(locDBL.get(dit()));
    BoxIterator bit(nodeGrid);
    for (bit.reset(); bit.ok(); ++bit)
    {
      pout() << "\t" << bit();
      for(int ivar = 0; ivar < fab.nComp(); ivar++)
        {
          pout() << "\t" << fab(bit(),ivar);
        }
      pout() << endl;
    }
  }
}

void dumpNodeLDFLoc(const LevelData<NodeFArrayBox>* a_ldfabPtr)
{
  const LevelData<NodeFArrayBox>& locLDF = *a_ldfabPtr;

  DataIterator dit = locLDF.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    const FArrayBox& fab = locLDF[dit()].getFab();

    DisjointBoxLayout dbl =  locLDF.disjointBoxLayout();
    Box nodeGrid = surroundingNodes(dbl.get(dit()));
    BoxIterator bit(nodeGrid);
    for (bit.reset(); bit.ok(); ++bit)
    {
      pout() << "\t" << bit() ;
      for(int ivar = 0; ivar < fab.nComp(); ivar++)
        {
          pout() << "\t" << fab(bit(),ivar);
        }
      pout() << endl;
    }
  }
}

void dumpNodeFAB(const NodeFArrayBox* a_fabPtr)
{
  const FArrayBox& fab = a_fabPtr->getFab();
  BoxIterator bit(fab.box());
  for (bit.reset(); bit.ok(); ++bit)
    {
      pout() << "\t" << bit() ;
      for(int ivar = 0; ivar < fab.nComp(); ivar++)
        {
          pout() << "\t" << fab(bit(),ivar);
        }
      pout() << endl;
    }
}
void dumpFAB(const FArrayBox* a_fabPtr)
{
  const FArrayBox& fab = *a_fabPtr;
  BoxIterator bit(fab.box());
  for (bit.reset(); bit.ok(); ++bit)
    {
      pout() << "\t" << bit() ;
      for(int ivar = 0; ivar < fab.nComp(); ivar++)
        {
          pout() << "\t" << fab(bit(),ivar);
        }
      pout() << endl;
    }
}
void dumpLDFPar(const LevelData<FArrayBox>* a_ldfabPtr)
{
  const LevelData<FArrayBox>& ldfab = *a_ldfabPtr;
  Vector<Box> boxes;

  const DisjointBoxLayout& memDBL = ldfab.getBoxes();
  LayoutIterator lit  = memDBL.layoutIterator();
  for (lit.reset(); lit.ok(); ++lit)
  {
    boxes.push_back(memDBL.get(lit()));
  }

  Vector<int> assign(boxes.size(), 0);

  DisjointBoxLayout locDBL(boxes, assign);
  locDBL.close();

  LevelData<FArrayBox> locLDF(locDBL,
                              ldfab.nComp(),
                              ldfab.ghostVect());

  const Interval& interv = ldfab.interval();
  ldfab.copyTo(interv,
               locLDF,
               interv);

  DataIterator dit = locLDF.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    const FArrayBox& fab = locLDF[dit()];
    BoxIterator bit(locDBL.get(dit()));
    for (bit.reset(); bit.ok(); ++bit)
    {
      pout() << "\t" << bit() ;
      for(int ivar = 0; ivar < fab.nComp(); ivar++)
        {
          pout() << "\t" << fab(bit(),ivar);
        }
      pout() << endl;
    }
  }
}

void dumpLDFLoc(const LevelData<FArrayBox>* a_ldfabPtr)
{
  const LevelData<FArrayBox>& locLDF = *a_ldfabPtr;


  DataIterator dit = locLDF.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    const FArrayBox& fab = locLDF[dit()];
    BoxIterator bit(fab.box());
    for (bit.reset(); bit.ok(); ++bit)
    {
      pout() << "\t" << bit() ;
      for(int ivar = 0; ivar < fab.nComp(); ivar++)
        {
          pout() << "\t" << fab(bit(),ivar);
        }
      pout() << endl;
    }
  }
}

void dumpIVSFAB(const IVSFAB<Real>* a_ivsfabPtr)
{
  const IVSFAB<Real>& ivsfab = *a_ivsfabPtr;
  const IntVectSet& ivs = ivsfab.getIVS();
  for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit)
  {
    const IntVect& iv = ivsit();

    for(int icomp = 0; icomp < ivsfab.nComp(); icomp++)
    {
      pout() << ivsfab(iv, icomp) << "  " << endl;
    }
  }
}

void dumpDBL(const DisjointBoxLayout* a_dblInPtr)
{
  pout() << indent2 <<"DisjointBoxLayout ";

  if (a_dblInPtr->isClosed())
  {
    pout() << "closed\n";
  }
  else
  {
    pout() << "open\n";
  }

  pout() << *a_dblInPtr << endl;
}

void dumpBL(const BoxLayout* a_dblInPtr)
{
  pout() << indent2 <<"BoxLayout ";

  if (a_dblInPtr->isClosed())
  {
    pout() << "closed\n";
  }
  else
  {
    pout() << "open\n";
  }

  pout() << *a_dblInPtr << endl;
}

void dumpIVS(const IntVectSet* a_ivsPtr)
{
  const IntVectSet& ivs = *a_ivsPtr;
  IVSIterator it(ivs);

  pout() << indent2 << ": IntVects in the IVS are:" << endl;
  pout() << indent2;

  for (it.begin(); it.ok(); ++it)
  {
    pout() << it() << "  ";
  }

  pout() << endl;
}

void dumpBox(const Box* a_boxPtr)
{
  pout() << indent2 << *a_boxPtr << endl;
}

void dumpVBox(const Vector<Box>* a_vectPtr)
{
  Vector<Box> vect = *a_vectPtr;
  for (int ivec = 0; ivec < vect.size(); ivec++)
  {
    pout() << indent2 << vect[ivec] << endl;
  }
}

void dumpVVBox(const Vector<Vector<Box> >* a_vectPtr)
{
  Vector<Vector<Box> > vect = *a_vectPtr;
  for (int iveco = 0; iveco < vect.size(); iveco++)
  {
    pout() << indent2;

    for (int iveci = 0; iveci < vect[iveco].size(); iveci++)
    {
      pout() <<  vect[iveco][iveci] << "   ";
    }

    pout() << endl;
  }
}

void dumpLDDBL(const LevelData<FArrayBox>* a_ldfabPtr)
{
  const DisjointBoxLayout& dbl = a_ldfabPtr->disjointBoxLayout();
  dumpDBL(&dbl);
}
void dumpNodeLDDBL(const LevelData<NodeFArrayBox>* a_ldfabPtr)
{
  const DisjointBoxLayout& dbl = a_ldfabPtr->disjointBoxLayout();
  dumpDBL(&dbl);
}
#include "NamespaceFooter.H"
