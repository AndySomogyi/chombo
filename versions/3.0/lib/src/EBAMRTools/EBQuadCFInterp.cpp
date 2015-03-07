#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EBQuadCFInterp.H"
#include "EBInterpolateF_F.H"
#include "VoFIterator.H"
#include "DisjointBoxLayout.H"
#include "EBCellFactory.H"
#include "LayoutIterator.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "EBAlias.H"
#include <algorithm>

#include "NamespaceHeader.H"

/***********************/
bool
EBQuadCFInterp::isDefined() const
{
  return m_isDefined;
}
/************************************/
EBQuadCFInterp::EBQuadCFInterp()
{
  m_isDefined = false;
  m_nComp = -1;
  m_refRat = -1;
}
/************************************/
EBQuadCFInterp::~EBQuadCFInterp()
{
}
/************************************/
EBQuadCFInterp::
EBQuadCFInterp(const DisjointBoxLayout& a_dblFine,
               const DisjointBoxLayout& a_dblCoar,
               const EBISLayout&        a_ebislFine,
               const EBISLayout&        a_ebislCoar,
               const ProblemDomain&     a_domainCoar,
               const int&               a_nref,
               const int&               a_nvar,
               const LayoutData<IntVectSet>&  a_cfivs,
               const EBIndexSpace* const a_ebisPtr,
               bool a_doEBCFCrossing)
  :QuadCFInterp()

{
  define(a_dblFine,   a_dblCoar,  a_ebislFine,
         a_ebislCoar, a_domainCoar, a_nref, a_nvar, a_cfivs, a_ebisPtr, a_doEBCFCrossing);
}
/************************************/
void
EBQuadCFInterp::
define(const DisjointBoxLayout& a_fineBoxes,
       const DisjointBoxLayout& a_coarBoxes,
       const EBISLayout&        a_ebislFine,
       const EBISLayout&        a_ebislCoar,
       const ProblemDomain&     a_domainCoar,
       const int&               a_nref,
       const int&               a_nvar,
       const LayoutData<IntVectSet>&  a_cfivs,
       const EBIndexSpace* const a_ebisPtr,
       bool a_doEBCFCrossing)
{
  QuadCFInterp::clear();
  m_nComp  = a_nvar;
  m_refRat = a_nref;
  m_gridsFine = a_fineBoxes;
  m_gridsCoar = a_coarBoxes;
  m_ebislFine = a_ebislFine;
  m_ebislCoar = a_ebislCoar;

  Real dxfine = 1;
  m_domainCoar = a_domainCoar;
  m_domainFine = refine(m_domainCoar, m_refRat);

  QuadCFInterp::define(a_fineBoxes, &a_coarBoxes, dxfine, a_nref, a_nvar, m_domainFine);

  m_isDefined = true;
  CH_assert(a_nvar > 0);
  CH_assert(a_nref > 0);
  CH_assert(a_fineBoxes.isClosed());
  CH_assert(a_coarBoxes.isClosed());
  CH_assert(!a_domainCoar.isEmpty());

  if(m_ebislFine.getMaxCoarseningRatio() < m_refRat)
    {
      m_ebislFine.setMaxCoarseningRatio(m_refRat);
    }

  //do EB-Crossing-CF interface tricks
  int nghostEBISL = 4;
  m_doEBCFCrossing = (a_doEBCFCrossing && (!m_fineCoversCoarse));
  if(m_doEBCFCrossing)
    {
      m_gridsCoarsenedFine = DisjointBoxLayout();
      coarsen(m_gridsCoarsenedFine, m_gridsFine, m_refRat);
      a_ebisPtr->fillEBISLayout(m_ebislCoarsenedFine,
                                m_gridsCoarsenedFine,
                                m_domainCoar, nghostEBISL);

      //define coarse ebis layouts buffers
      EBCellFactory coarFact(m_ebislCoarsenedFine);
      m_ebBufferCoarsenedFine.define(m_gridsCoarsenedFine, m_nComp, 2*IntVect::Unit, coarFact);

      //figure out where we need to do our EB-crossing-tricks (and if we need them at all)
      m_doEBCFCrossing = getEBCFIVS(a_cfivs);

      if(m_doEBCFCrossing)
        {
          //build EBCF stencils
          buildEBCFStencils();
        }

    }
  //there are cases where no crossing occurs that you will still need the corners
  //build the ebcf corner stencils
  buildEBCFCornerStencils(a_cfivs);

  //make a corner copier
  m_cornerCopier.define(m_gridsFine, m_gridsFine, m_domainFine, IntVect::Unit);

}
/************************************/
void
getEdgeAndCornerIVS(IntVectSet& a_edgeIVS, IntVectSet& a_cornerIVS,
                    const Box& a_grid, const ProblemDomain& a_domain,
                    const IntVectSet& a_cfivsGrid)
{
  Box grownBox = grow(a_grid, 1);
  grownBox &= a_domain;
  IntVectSet ivsAll(grownBox);
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      Box growBoxSide = a_grid;
      growBoxSide.grow(idir, 1);
      ivsAll -= growBoxSide;
    }

  a_cornerIVS = ivsAll;
  if(SpaceDim == 2)
    {
      //in 2d there are no edges
      a_edgeIVS = IntVectSet();
    }
  else
    {
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          for(SideIterator sit1; sit1.ok(); ++sit1)
            {
              Box faceBox = adjCellBox(a_grid, idir, sit1(), 1);
              for(int jdir = 0; jdir < SpaceDim; jdir++)
                {
                  if(jdir != idir)
                    {
                      for(SideIterator sit2; sit2.ok(); ++sit2)
                        {
                          Box edgeBox = adjCellBox(faceBox, jdir, sit2(), 1);
                          a_cornerIVS -= edgeBox;
                        }
                    }
                }
            }
        }
      a_edgeIVS = ivsAll - a_cornerIVS;
    }

  IntVectSet grownCFIVS= grow(a_cfivsGrid, 1);
  a_cornerIVS &= grownCFIVS;
  a_edgeIVS   &= grownCFIVS;
}
/************************************/
bool
EBQuadCFInterp::
getEBCFIVS(const LayoutData<IntVectSet>&  a_cfivs)
{
  bool doEBCrossing = false;

  //now for the side ivs in each direction.   If nothing is found to be crossing,
  // then there is no EB crossing
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      m_ebcfivsLo[idir].define(m_gridsFine);
      m_ebcfivsHi[idir].define(m_gridsFine);
      for(DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
        {

          const IntVectSet& fineIVSLo = m_loQCFS[idir][dit()].getFineIVS();
          const IntVectSet& fineIVSHi = m_hiQCFS[idir][dit()].getFineIVS();
          Box grownBox = grow(m_gridsFine.get(dit()), 1);
          grownBox &= m_domainFine;
          IntVectSet ivsIrreg = m_ebislFine[dit()].getIrregIVS(grownBox);
          IntVectSet ivsLo = fineIVSLo & ivsIrreg;
          IntVectSet ivsHi = fineIVSHi & ivsIrreg;
          //need to grow by the refinement ratio because otherwise cell next to
          //non-eb quadcfinterp can do the wrong thing when it takes derivatvies
          //at coarse cells
          for(int jdir = 0; jdir < SpaceDim; jdir++)
            {
              if(jdir != idir)
                {
                  ivsLo.grow(jdir, 2*m_refRat);
                  ivsLo &= fineIVSLo;
                  ivsHi.grow(jdir, 2*m_refRat);
                  ivsHi &= fineIVSHi;
                }
            }
          m_ebcfivsLo[idir][dit()] = ivsLo;
          m_ebcfivsHi[idir][dit()] = ivsHi;
          //keep track if this  crossing ever happens
          if((!m_ebcfivsLo[idir][dit()].isEmpty()) ||
             (!m_ebcfivsHi[idir][dit()].isEmpty()))
            {
              doEBCrossing = true;
            }
        }
    }
  //in the case of parallel, need to check if ANY of the procs
  //have ebcrossing
#ifdef CH_MPI

  int gatherint = 0;
  if(doEBCrossing) gatherint = 1;     
  int idoebcf;
  MPI_Allreduce(&gatherint, &idoebcf, 1, MPI_INT,
                MPI_MAX, Chombo_MPI::comm);
  doEBCrossing = (idoebcf==1);

#endif

  return doEBCrossing;
}

/*********************/
int
getPhiStarStencil(VoFStencil&     a_stencil,
                  const VolIndex& a_ghostVoFFine,
                  const VolIndex& a_ghostVoFCoar,
                  const EBISBox&  a_ebisBoxCoar,
                  int a_idir, Side::LoHiSide a_sd,
                  int a_refRat)
{
  RealVect dxFine = IntVect::Unit;
  RealVect dxCoar = IntVect::Unit;
  dxCoar *= a_refRat;

  RealVect fineLoc = EBArith::getVofLocation(a_ghostVoFFine, dxFine, IntVect::Zero);
  RealVect coarLoc = EBArith::getVofLocation(a_ghostVoFCoar, dxCoar, IntVect::Zero);

  RealVect dist = fineLoc - coarLoc;
  //don't want any extrapolation in a_idir direction
  //phistar lives in the line connnecting coarse cell centers
  //see algorithm doc
  dist[a_idir] = 0.0;


  //don't want any extrapolation in a_idir direction
  //this gets the stencil to extrapolation from phiCoar to phiStar
  int order  = EBArith::getExtrapolationStencil(a_stencil, dist, dxCoar, a_ghostVoFCoar, a_ebisBoxCoar, a_idir);

  return order;
}

/*********************/
int
getStencils(VoFStencil&     a_fineStencil,
            VoFStencil&     a_coarStencil,
            const VolIndex& a_ghostVoFFine,
            const VolIndex& a_ghostVoFCoar,
            const EBISBox&  a_ebisBoxFine,
            const EBISBox&  a_ebisBoxCoar,
            int a_idir, Side::LoHiSide a_sd,
            int a_refRat)
{
  a_fineStencil.clear();
  a_coarStencil.clear();
  bool hasVoF1 = false;
  bool hasVoF2 = false;
  VolIndex vof1, vof2;
  Side::LoHiSide flipsd = flip(a_sd);
  //coarse stencil is phiStencil*factor (the factor depends upon drop order effects)
  Real phiStarFactor;

  int order;
  Vector<FaceIndex> closeFaces = a_ebisBoxFine.getFaces(a_ghostVoFFine, a_idir, flipsd);
  hasVoF1 = (closeFaces.size() == 1) && (!closeFaces[0].isBoundary());
  if(hasVoF1)
    {
      vof1= closeFaces[0].getVoF(flipsd);
      Vector<FaceIndex> farFaces = a_ebisBoxFine.getFaces(vof1, a_idir, flipsd);
      hasVoF2 = (farFaces.size() == 1) && (!farFaces[0].isBoundary());
      if(hasVoF2)
        {
          vof2 = farFaces[0].getVoF(flipsd);
        }
    }
  if(hasVoF2)
    {
      order = 2;
      Real denom =   Real( 3 + 4*a_refRat +   a_refRat*a_refRat);
      Real phi1Fac = Real(-6 + 4*a_refRat + 2*a_refRat*a_refRat)/denom;
      Real phi2Fac = Real(1 - a_refRat*a_refRat)/denom;

      a_fineStencil.add(vof1, phi1Fac);
      a_fineStencil.add(vof2, phi2Fac);
      phiStarFactor = 8.0/denom;
    }
  else if(hasVoF1)
    {
      order = 1;
      Real denom   = Real(a_refRat + 1);
      Real phi1Fac = Real(a_refRat - 1)/denom;
      a_fineStencil.add(vof1, phi1Fac);
      phiStarFactor = 2.0/denom;
    }
  else
    {
      //no fine values--- use coarse value = phiStar
      order = 0;
      phiStarFactor = 1.0;
    }

  //coarse stencil is phiStar stencil*factor (the factor depends upon drop order effects)
  int phiStarOrder = getPhiStarStencil(a_coarStencil,
                                       a_ghostVoFFine, a_ghostVoFCoar,
                                       a_ebisBoxCoar, a_idir, a_sd, a_refRat);
  order = Min(order, phiStarOrder);

  a_coarStencil *= phiStarFactor;

  return order;

}
/*********************/
void
EBQuadCFInterp::
buildEBCFStencils()
{
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      m_coarStencilLo[idir].define(m_gridsFine);
      m_coarStencilHi[idir].define(m_gridsFine);
      m_fineStencilLo[idir].define(m_gridsFine);
      m_fineStencilHi[idir].define(m_gridsFine);
      m_vofItEBCFLo[idir].define(m_gridsFine);
      m_vofItEBCFHi[idir].define(m_gridsFine);
      for(DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
        {
          //define all the holders
          const IntVectSet& loEBCFIVS = m_ebcfivsLo[idir][dit()];
          const IntVectSet& hiEBCFIVS = m_ebcfivsHi[idir][dit()];

          m_vofItEBCFLo[idir][dit()].define(loEBCFIVS, m_ebislFine[dit()].getEBGraph());
          m_vofItEBCFHi[idir][dit()].define(hiEBCFIVS, m_ebislFine[dit()].getEBGraph());

          m_coarStencilLo[idir][dit()].define(loEBCFIVS, m_ebislFine[dit()].getEBGraph(), 1);
          m_coarStencilHi[idir][dit()].define(hiEBCFIVS, m_ebislFine[dit()].getEBGraph(), 1);
          m_fineStencilLo[idir][dit()].define(loEBCFIVS, m_ebislFine[dit()].getEBGraph(), 1);
          m_fineStencilHi[idir][dit()].define(hiEBCFIVS, m_ebislFine[dit()].getEBGraph(), 1);

          for(m_vofItEBCFLo[idir][dit()].reset(); m_vofItEBCFLo[idir][dit()].ok(); ++m_vofItEBCFLo[idir][dit()])
            {
              const VolIndex& vofFine = (m_vofItEBCFLo[idir][dit()])();
              VolIndex vofCoar = m_ebislFine.coarsen(vofFine, m_refRat, dit());
              //phistar only depends upon coarse vofs.
              //the line interp stencil depends upon fine vofs
              //but they interact depending on drop order considerations

              getStencils(m_fineStencilLo[idir][dit()](vofFine, 0),
                          m_coarStencilLo[idir][dit()](vofFine, 0),
                          vofFine, vofCoar,
                          m_ebislFine[dit()], m_ebislCoarsenedFine[dit()],
                          idir, Side::Lo, m_refRat);
            }
          for(m_vofItEBCFHi[idir][dit()].reset(); m_vofItEBCFHi[idir][dit()].ok(); ++m_vofItEBCFHi[idir][dit()])
            {
              const VolIndex& vofFine = (m_vofItEBCFHi[idir][dit()])();
              VolIndex vofCoar = m_ebislFine.coarsen(vofFine, m_refRat, dit());
              //phistar only depends upon coarse vofs.
              //the line interp stencil depends upon fine vofs
              //but they interact depending on drop order considerations
              getStencils(m_fineStencilHi[idir][dit()](vofFine, 0),
                          m_coarStencilHi[idir][dit()](vofFine, 0),
                          vofFine, vofCoar,
                          m_ebislFine[dit()], m_ebislCoarsenedFine[dit()],
                          idir, Side::Hi, m_refRat);
            }
        }
    }
}
/*********************/
void
getExtrapSigns(IntVect& a_signs, const IntVect& a_corner, const Box& a_grid)
{
  const IntVect& loEnd = a_grid.smallEnd();
  const IntVect& hiEnd = a_grid.bigEnd();
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      if(a_corner[idir] > hiEnd[idir])
        {
          a_signs[idir] = -1; //estrapolate from the low side if on the high end
        }
      else if(a_corner[idir] < loEnd[idir])
        {
          a_signs[idir] =  1; //estrapolate from the high side if on the low end
        }
      else
        {
#if CH_SPACEDIM==2
          //in 2d, one should not get here by construction
          pout() << "getExtrapSigns a_corner = " << a_corner << ", a_grid = " << a_grid << endl;
          MayDay::Error("getExtrapSigns: corner vof apparently not in the corner");

#elif CH_SPACEDIM==3
          //in 3d you can be along an edge.  only extrapolate in the other directions.
          a_signs[idir] = 0;
#else
          bogus_spacedim_compilation_error();
#endif
        }

    }
}
//returns the order of the extrapolation
//if it cannot find any points to extrap, return -1
int
getCornerExtrapStencil(VoFStencil& a_stencil, const VolIndex& a_cornerVoF, const EBISBox& a_ebisBox, const Box& a_grid)
{
  int order = 3;
  IntVect extrapSigns;
  getExtrapSigns(extrapSigns,a_cornerVoF.gridIndex(), a_grid);
  a_stencil.clear();
  int numExtraps = 0;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      bool skipThisDir;
      Side::LoHiSide sd;
      if(extrapSigns[idir] == -1)
        {
          skipThisDir = false;
          sd = Side::Lo;
        }
      else if(extrapSigns[idir] ==  1)
        {
          skipThisDir = false;
          sd = Side::Hi;
        }
      else
        {
          //extrapSigns == 0 means do not extrapolate in that direction
          skipThisDir = true;
          sd = Side::Invalid;
        }
      if(!skipThisDir)
        {
          Vector<FaceIndex> closeFaces = a_ebisBox.getFaces(a_cornerVoF, idir, sd);
          VolIndex closeVoF, farVoF, veryFarVoF;
          bool hasClose= false;
          bool hasFar = false;
          bool hasVeryFar = false;
          hasClose = (closeFaces.size() == 1) && (!closeFaces[0].isBoundary());
          if(hasClose)
            {
              //can extrap since we have at least one point
              numExtraps++;
              closeVoF = closeFaces[0].getVoF(sd);
              Vector<FaceIndex> farFaces = a_ebisBox.getFaces(closeVoF, idir, sd);
              hasFar = (farFaces.size() == 1) && (!farFaces[0].isBoundary());
              if(hasFar)
                {
                  farVoF = farFaces[0].getVoF(sd);
                  Vector<FaceIndex> veryFarFaces = a_ebisBox.getFaces(farVoF, idir, sd);
                  hasVeryFar = (veryFarFaces.size() == 1) && (!veryFarFaces[0].isBoundary());
                  if(hasVeryFar)
                    {
                      veryFarVoF = veryFarFaces[0].getVoF(sd);
                    }
                }
            }
          if(hasVeryFar)
            {
              //has all three points
              a_stencil.add(closeVoF,  3.0);
              a_stencil.add(farVoF,   -3.0);
              a_stencil.add(veryFarVoF,1.0);
            }
          else if(hasFar)
            {
              //only has two points
              a_stencil.add(closeVoF,  2.0);
              a_stencil.add(farVoF,   -1.0);
              order = Min(order, 2);
            }
          else if(hasClose)
            {
              //only one point in this direction
              a_stencil.add(closeVoF, 1.0);
              order = Min(order, 1);
            }
        }
    }
  if(numExtraps > 0)
    {
      a_stencil *= (1.0/Real(numExtraps));
    }
  if(numExtraps == 0)
    {
      //we will use the coarse value.
      order = -1;
    }
  return order;
}

void EBQuadCFInterp::
buildEBCFCornerStencils(const LayoutData<IntVectSet>& a_cfivs)
{
  // let us get the corner IVS
  m_vofItCorners.define(m_gridsFine);
  m_vofItEdges.define(m_gridsFine);
  m_stencilCorners.define(m_gridsFine);
  m_stencilEdges.define(m_gridsFine);

  for(DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      const Box& gridFine = m_gridsFine.get(dit());
      IntVectSet cornerIVS,  edgeIVS;
      getEdgeAndCornerIVS(edgeIVS, cornerIVS, gridFine, m_domainFine, a_cfivs[dit()]);

      m_vofItCorners[dit()].define(cornerIVS, m_ebislFine[dit()].getEBGraph());
      m_vofItEdges[dit()].define(    edgeIVS, m_ebislFine[dit()].getEBGraph());

      m_stencilCorners[dit()].define(cornerIVS, m_ebislFine[dit()].getEBGraph(),1);
      m_stencilEdges[dit()].define(    edgeIVS, m_ebislFine[dit()].getEBGraph(),1);

      for(m_vofItCorners[dit()].reset(); m_vofItCorners[dit()].ok(); ++m_vofItCorners[dit()])
        {
          const VolIndex& vof = m_vofItCorners[dit()]();
          getCornerExtrapStencil(m_stencilCorners[dit()](vof,0), vof, m_ebislFine[dit()], m_gridsFine.get(dit()));
        }
      for(m_vofItEdges[dit()].reset(); m_vofItEdges[dit()].ok(); ++m_vofItEdges[dit()])
        {
          const VolIndex& vof = m_vofItEdges[dit()]();
          getCornerExtrapStencil(m_stencilEdges[dit()](vof,0), vof, m_ebislFine[dit()], m_gridsFine.get(dit()));
        }
    }
}
/***********************/
void
EBQuadCFInterp::interpolate(LevelData<EBCellFAB>&       a_fineData,
                            const LevelData<EBCellFAB>& a_coarData,
                            const Interval&             a_variables,
                            bool a_doOnlyRegularInterp)
{
  Interval totInterv(0, m_nComp-1);
  CH_assert(a_fineData.nComp() == m_nComp);
  CH_assert(a_coarData.nComp() == m_nComp);
  CH_assert(a_variables.begin() >= 0);
  CH_assert(a_variables.begin() < m_nComp);
  if(totInterv.size() != a_variables.size())
    {
      MayDay::Warning("EBQuadCFInterp variables mismatch");
    }
  LevelData<FArrayBox> fineDataLDFAB, coarDataLDFAB;

  aliasEB(fineDataLDFAB, a_fineData);
  aliasEB(coarDataLDFAB, (LevelData<EBCellFAB>&)a_coarData);
  QuadCFInterp::coarseFineInterp(fineDataLDFAB, coarDataLDFAB);

  if(!a_doOnlyRegularInterp)
    {
      if(m_doEBCFCrossing)
        {
          //          pout() << "going into doEBCFCrossing" << endl;
          interpEBCFCrossing(a_fineData,a_coarData,a_variables);
          //          pout() << "coming out doEBCFCrossing" << endl;
        }

      //      pout() << "going into interpEBCFCorners" << endl;
      //fix corner vofs --this includes exchanges
      interpEBCFCorners(a_fineData, a_coarData, a_variables);
      //      pout() << "coming out interpEBCFCorners" << endl;
    }
}
/***********************/
void
EBQuadCFInterp::
interpEBCFCrossing(LevelData<EBCellFAB>&       a_fineData,
                   const LevelData<EBCellFAB>& a_coarData,
                   const Interval&             a_variables)
{
  // coarsen to coarse vofs covered by finer HAS TO BE ALREADY DONE
  //took out call to ebcoarsen here

  //copy the coarse data to our coarsened fine grids data holders at CF boundaries
  // (this part is blocking)
  a_coarData.copyTo(a_variables, m_ebBufferCoarsenedFine, a_variables);
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      for(DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
        {
          for(m_vofItEBCFLo[idir][dit()].reset(); m_vofItEBCFLo[idir][dit()].ok(); ++m_vofItEBCFLo[idir][dit()])
            {
              const   VolIndex& vofGhost =  (m_vofItEBCFLo[idir][dit()])();
              const VoFStencil& fineSten = m_fineStencilLo[idir][dit()](vofGhost, 0);
              const VoFStencil& coarSten = m_coarStencilLo[idir][dit()](vofGhost, 0);
              for(int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++)
                {
                  Real fineContrib = applyVoFStencil(fineSten,              a_fineData[dit()], icomp);
                  Real coarContrib = applyVoFStencil(coarSten, m_ebBufferCoarsenedFine[dit()], icomp);
                  a_fineData[dit()](vofGhost, icomp) = fineContrib + coarContrib;
                }
            }
          for(m_vofItEBCFHi[idir][dit()].reset(); m_vofItEBCFHi[idir][dit()].ok(); ++m_vofItEBCFHi[idir][dit()])
            {
              const   VolIndex& vofGhost =  (m_vofItEBCFHi[idir][dit()])();
              const VoFStencil& fineSten = m_fineStencilHi[idir][dit()](vofGhost, 0);
              const VoFStencil& coarSten = m_coarStencilHi[idir][dit()](vofGhost, 0);
              for(int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++)
                {
                  Real fineContrib = applyVoFStencil(fineSten,              a_fineData[dit()], icomp);
                  Real coarContrib = applyVoFStencil(coarSten, m_ebBufferCoarsenedFine[dit()], icomp);
                  a_fineData[dit()](vofGhost, icomp) = fineContrib + coarContrib;
                }
            }
        }
    }
}
/***********************/
void
EBQuadCFInterp::
interpEBCFCorners(LevelData<EBCellFAB>&       a_fineData,
                  const LevelData<EBCellFAB>& a_coarData,
                  const Interval&             a_variables)
{
  //exchange
  a_fineData.exchange(a_variables);

  //fills the edges with extrapolations from neighbors not along edge
  //fills the corners with extrapolations from neighbors (which are edgea
  for(DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      for(m_vofItEdges[dit()].reset(); m_vofItEdges[dit()].ok(); ++m_vofItEdges[dit()])
        {
          const   VolIndex& vofEdge = (m_vofItEdges[dit()])();
          const VoFStencil& stencil = m_stencilEdges[dit()](vofEdge, 0);
          if(stencil.size() > 0)  //there are fine vofs to extrapolate from
            {
              for(int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++)
                {
                  Real edgeVal = applyVoFStencil(stencil, a_fineData[dit()], icomp);
                  a_fineData[dit()](vofEdge, icomp) = edgeVal;
                }
            }
          else //just use the coarse value
            {
              VolIndex coarVoF = m_ebislFine.coarsen(vofEdge, m_refRat, dit());
              for(int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++)
                {
                  Real edgeVal = m_ebBufferCoarsenedFine[dit()](coarVoF, icomp);
                  a_fineData[dit()](vofEdge, icomp) = edgeVal;
                }

            }
        }

      for(m_vofItCorners[dit()].reset(); m_vofItCorners[dit()].ok(); ++m_vofItCorners[dit()])
        {
          const   VolIndex& vofCorn = (m_vofItCorners[dit()])();
          const VoFStencil& stencil = m_stencilCorners[dit()](vofCorn, 0);
          if(stencil.size() > 0)  //there are fine vofs to extrapolate from
            {
              for(int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++)
                {
                  Real cornerVal = applyVoFStencil(stencil, a_fineData[dit()], icomp);
                  a_fineData[dit()](vofCorn, icomp) = cornerVal;
                }
            }
          else //just use the coarse value
            {
              VolIndex coarVoF = m_ebislFine.coarsen(vofCorn, m_refRat, dit());
              for(int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++)
                {
                  Real cornerVal = m_ebBufferCoarsenedFine[dit()](coarVoF, icomp);
                  a_fineData[dit()](vofCorn, icomp) = cornerVal;
                }

            }
        }
    }


  a_fineData.exchange(a_variables, m_cornerCopier);
}

#include "NamespaceFooter.H"
