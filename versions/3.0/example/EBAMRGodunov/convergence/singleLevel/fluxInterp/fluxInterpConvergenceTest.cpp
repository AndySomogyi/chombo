#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <cmath>
#include <cstdio>
#include <iostream>

#include "PolyGeom.H"
#include "TiltedCylinderIF.H"
#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "EBPatchPolytropicFactory.H"
#include "EBPatchPolytropic.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "GeometryShop.H"
#include "FaceIterator.H"
#include "EBIndexSpace.H"
#include "PlaneIF.H"
#include "IntersectionIF.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#if CH_SPACEDIM==2
IntVect ivdebug(2,5);
#endif
#if CH_SPACEDIM==3
IntVect ivdebug(2,5,0);
#endif

/***************/
/***************/
void
makeGeometry(Box& a_domain,
             Real& a_dx)
{
  //parse input file.  single level
  ParmParse pp;
  RealVect origin = RealVect::Zero;
  Vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);

 CH_assert(n_cell.size() == SpaceDim);
  IntVect lo = IntVect::Zero;
  IntVect hi;
  for(int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if(n_cell[ivec] <= 0)
        {
          pout() << " bogus number of cells input = " << n_cell[ivec];
          MayDay::Error();
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  a_domain = Box(lo, hi);

  Real prob_hi;
  int numOpen = n_cell[0];
  pp.get("domain_length",prob_hi);
  a_dx = prob_hi/numOpen;

  int whichGeom;
  pp.get("which_geom", whichGeom);
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  int biggridsize;
  char ebfilename[80];
  int  fineLength = a_domain.size(0);
  pp.get("max_grid_size", biggridsize);
  Real channelRadius;
  pp.get("channel_radius", channelRadius);

  if(whichGeom == 2)
    {
      RealVect channelNormal;
      vector<Real>  channelNormalVect(SpaceDim);
      pp.getarr("channel_normal",channelNormalVect, 0, SpaceDim);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          channelNormal[idir] = channelNormalVect[idir];
        }

      Real norm = 0.0;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          norm += channelNormal[idir] * channelNormal[idir];
        }
      norm = sqrt(norm);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          channelNormal[idir] /= norm;
        }

      RealVect botPoint = RealVect::Zero;
      botPoint -= channelNormal;
      botPoint *= channelRadius;

      RealVect topPoint = RealVect::Zero;
      topPoint += channelNormal;
      topPoint *= channelRadius;

      PlaneIF bot(channelNormal,botPoint,true);
      PlaneIF top(channelNormal,topPoint,false);
      IntersectionIF channel(top,bot);

      pout() << "using a channel" << endl;

      RealVect vectDx = RealVect::Unit;
      vectDx *= a_dx;

      GeometryShop workshop(channel,0,vectDx);
      ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize, 1);
      sprintf(ebfilename, "channel%d.%dd.hdf5",fineLength, SpaceDim);
    }
  else if(whichGeom == 3)
    {
      pout() << "using a tilted cylinder" << endl;
      vector<Real>  cylinderAxisVect(SpaceDim);
      pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
      RealVect cylinderAxis;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          cylinderAxis[idir] = cylinderAxisVect[idir];
        }
      Real sum;
      PolyGeom::unifyVector(cylinderAxis, sum);

      RealVect corner = RealVect::Zero;
      bool negativeInside = true;
      TiltedCylinderIF tunnel(channelRadius, cylinderAxis, corner, negativeInside);

      RealVect vectDx = RealVect::Unit;
      vectDx *= a_dx;

      GeometryShop workshop(tunnel,0,vectDx);
      ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize);


      sprintf(ebfilename, "cyl%d.%dd.hdf5",fineLength, SpaceDim);
    }
  else
    {
      MayDay::Error("bogus geometry flag");
    }
}
/***************/
/***************/
Real
norm(const BaseIFFAB<Real>& a_error,
     const EBISBox& a_ebisBox,
     const IntVectSet& a_set,
     const int& a_idir,
     const int& a_comp,
     const int& a_normtype)
{
 CH_assert(a_normtype >= 0);
 CH_assert(a_normtype <= 2);
  Real retval = 0.;
  FaceIterator faceit(a_set, a_ebisBox.getEBGraph(), a_idir,
                      FaceStop::SurroundingWithBoundary);
  if(a_normtype == 0)
    {
      //max norm
      for(faceit.reset(); faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          Real valFace = a_error(face, a_comp);
          retval = Max(Abs(valFace), retval);
        }
    }
  else
    {
      //integral norm
      Real areaTot = 0.0;
      Real normTot = 0.0;
      for(faceit.reset(); faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          Real areaFrac = a_ebisBox.areaFrac(face);
          Real valFace = a_error(face, a_comp);
          areaTot += areaFrac;
          if(a_normtype == 1)
            normTot += Abs(valFace)*areaFrac;
          else
            normTot += valFace*valFace*areaFrac;
        }
      if(areaTot > 1.0e-8)
        normTot /= areaTot;
      if(a_normtype == 2)
        normTot = sqrt(normTot);

      retval = normTot;
    }
  return retval;
}

void
compareError(const BaseIFFAB<Real> a_errorFine[SpaceDim],
             const EBISBox& a_ebisBoxFine,
             const IntVectSet& a_setsFine,
             const BaseIFFAB<Real> a_errorCoar[SpaceDim],
             const EBISBox& a_ebisBoxCoar,
             const IntVectSet& a_setsCoar)
{
  int ncomp = a_errorFine[0].nComp();
  pout() << "==============================================" << endl;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      pout() << "Direction = " << idir << endl;
      pout() << "==============================================" << endl;
      for(int comp = 0; comp < ncomp; comp++)
        {
          pout() << "component = " << comp << endl;
          pout() << "==============================================" << endl;
          for(int inorm = 0; inorm < 3; inorm++)
            {

              if(inorm == 0)
                {
                  pout() << endl << "Using max norm." << endl;
                }
              else
                {
                  pout() << endl << "Using L-" << inorm << "norm." << endl;
                }
              Real coarnorm = norm(a_errorCoar[idir],
                                   a_ebisBoxCoar, a_setsCoar,
                                   idir, comp, inorm);
              Real finenorm = norm(a_errorFine[idir],
                                   a_ebisBoxFine, a_setsFine,
                                   idir, comp, inorm);
              pout() << "Coarse Error Norm = " << coarnorm << endl;
              pout() << "Fine   Error Norm = " << finenorm << endl;
              if((Abs(finenorm) > 1.0e-8) && (Abs(coarnorm) > 1.0e-8))
                {
                  Real order = log(coarnorm/finenorm)/log(2.0);
                  pout() << "Order of scheme = " << order << endl;
                }
            }
        }
      pout() << "==============================================" << endl;
    }
}
/***************/
/***************/
Real exactFunc(const RealVect& a_xval)
{

  Real retval;

#if CH_SPACEDIM==2

  Real x = a_xval[0];
  Real y = a_xval[1];
  //this ought to be exact
  retval = x + y;
  retval = x*x + y*y + x + y + x*y + x*x*x + y*y*y;
  //  retval = y;
#elif CH_SPACEDIM==3

  Real x = a_xval[0];
  Real y = a_xval[1];
  Real z = a_xval[2];
  retval = x + y + z + x*y + y*z + x*z + x*x + y*y + z*z + x*x*x + y*y*y + z*z*z;
  //  retval = y;

#else
  bogus_ch_spacedim_compilation_macro()
#endif

  //debug
  //retval = x*x + y*y + x + y + x*y + x*x*x + y*y*y;
  //  retval = x+ y + z;
  return retval;
}
/***************/
/***************/
void
getError(BaseIFFAB<Real> a_error[SpaceDim],
         IntVectSet&     a_ivs,
         const EBISBox&  a_ebisBox,
         const Box&      a_domain,
         const Real&     a_dx)
{
  a_ivs = a_ebisBox.getIrregIVS(a_domain);
  BaseIFFAB<Real> fluxInterpol[SpaceDim];
  BaseIFFAB<Real> centroidFlux[SpaceDim];
  for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      //the interpolant must live on a grown IVS
      IntVectSet grownIVS  = a_ivs;
      for(int jdir = 0; jdir < SpaceDim; jdir++)
        {
          if(faceDir != jdir)
            {
              grownIVS.grow(jdir, 1);
            }
        }
      grownIVS &= a_domain;
      BaseIFFAB<Real>& interpol = fluxInterpol[faceDir];

      //define the various data holders
      interpol.define(grownIVS, a_ebisBox.getEBGraph(), faceDir, 1);
      centroidFlux[faceDir].define(a_ivs, a_ebisBox.getEBGraph(), faceDir, 1);
      a_error[faceDir].define(a_ivs, a_ebisBox.getEBGraph(), faceDir, 1);


      //fill the interpolant with the exact answer
      //evaluated at centers of cell faces
      FaceIterator faceit(grownIVS, a_ebisBox.getEBGraph(), faceDir,
                          FaceStop::SurroundingWithBoundary);
      for(faceit.reset(); faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          IntVect ivhi = face.gridIndex(Side::Hi);
          RealVect xval;
          for(int jdir = 0; jdir < SpaceDim; jdir++)
            {
              if(faceDir != jdir)
                {
                  xval[jdir] = a_dx*(Real(ivhi[jdir]) + 0.5);
                }
              else
                {
                  xval[jdir] = a_dx*Real(ivhi[jdir]);
                }
            }

          Real exactVal = exactFunc(xval);
          interpol(face, 0) = exactVal;
        }
    }

  //now interpolate the whole mess to the input flux
   bool is_periodic[SpaceDim];
   for (int i = 0; i < SpaceDim; ++i)
     {
       is_periodic[i] = false;
     }
  ProblemDomain probdomain(a_domain.smallEnd(), a_domain.bigEnd(), is_periodic);
  EBPatchPolytropic patchInt;
  patchInt.define(probdomain, a_dx);
  //single level implies empty cfivsl
  IntVectSet cfivs;
  //last two args are fake time and dt
  patchInt.setValidBox(a_domain, a_ebisBox, cfivs, 1.0, 1.0);
  const BaseIFFAB<Real>*  interpolantGrid[SpaceDim];
  for(int idir = 0; idir < SpaceDim; idir++)
    interpolantGrid[idir] = &fluxInterpol[idir];

  patchInt.interpolateFluxToCentroids(centroidFlux, interpolantGrid, a_ivs);

  //now subtract what we got from the right answer
  for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {

      //fill the interpolant with the exact answer
      FaceIterator faceit(a_ivs, a_ebisBox.getEBGraph(), faceDir,
                          FaceStop::SurroundingWithBoundary);

      const BaseIFFAB<Real>& centroidVal = centroidFlux[faceDir];
      BaseIFFAB<Real>& error = a_error[faceDir];
      //this time evaluated at centroids of faces
      for(faceit.reset(); faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          RealVect centroid = a_ebisBox.centroid(face);
          IntVect ivhi = face.gridIndex(Side::Hi);
          RealVect xval;
          for(int jdir = 0; jdir < SpaceDim; jdir++)
            {
              if(faceDir != jdir)
                {
                  xval[jdir] = a_dx*(Real(ivhi[jdir]) + 0.5 + centroid[jdir]);
                }
              else
                {
                  xval[jdir] = a_dx*Real(ivhi[jdir]);
                }
            }

          Real exactVal = exactFunc(xval);
          Real calcVal  = centroidVal(face, 0);
          error(face, 0) = calcVal - exactVal;
        }
    }


}
/***************/
/***************/
void
fluxInterpTest()
{
  Box domainFine, domainCoar;
  Real dxFine, dxCoar;
  makeGeometry(domainFine, dxFine);
  dxCoar = 2.*dxFine;
  domainCoar = coarsen(domainFine, 2);

  Vector<int> pmap(1, 0);
  Vector<Box> coarBoxes(1, domainCoar);
  Vector<Box> fineBoxes(1, domainFine);
  DisjointBoxLayout dblFine(fineBoxes, pmap);
  DisjointBoxLayout dblCoar;
  //so they will both with with the same iterator
  coarsen(dblCoar, dblFine, 2);
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  EBISLayout ebislFine, ebislCoar;
  ebisPtr->fillEBISLayout(ebislFine, dblFine, domainFine, 0);
  ebisPtr->fillEBISLayout(ebislCoar, dblCoar, domainCoar, 0);

  //just looping through one box but the syntax demands it.
  for(DataIterator dit = dblFine.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBoxFine = ebislFine[dit()];
      const EBISBox& ebisBoxCoar = ebislCoar[dit()];
      IntVectSet ivsFine;
      IntVectSet ivsCoar;
      BaseIFFAB<Real> errorFine[SpaceDim];
      BaseIFFAB<Real> errorCoar[SpaceDim];

      getError(errorFine, ivsFine, ebisBoxFine, domainFine, dxFine);
      getError(errorCoar, ivsCoar, ebisBoxCoar, domainCoar, dxCoar);
      compareError(errorFine, ebisBoxFine, ivsFine,
                   errorCoar, ebisBoxCoar, ivsCoar);
    }
}
/***************/
/***************/
int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
  {
    // setChomboMPIErrorHandler();
#endif
    char* inFile = NULL;


    if (a_argc > 1)
      {
        inFile = a_argv[1];
      }
    else
      {
        pout() << "Usage: <executable name> <inputfile>" << endl;
        pout() << "No input file specified" << endl;
        return -1;
      }
    ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

    fluxInterpTest();

#ifdef CH_MPI
    dumpmemoryatexit();
  }
  MPI_Finalize();
#endif
}
/***************/
/***************/
