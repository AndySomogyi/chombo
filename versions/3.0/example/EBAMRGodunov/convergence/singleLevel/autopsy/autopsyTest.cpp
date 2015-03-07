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
#include <iostream>

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"

#include "EBExplosionIBCFactory.H"
#include "ModianoIBCFactory.H"
#include "EBPatchPolytropicFactory.H"

#include "EBPatchPolytropic.H"

#include "BaseIVFactory.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "VoFIterator.H"
#include "EBArith.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBLGIntegrator.H"
#include "AllRegularService.H"
#include "TiltedCylinderIF.H"
#include "SquareCylinderIF.H"
#include "TransformIF.H"
#include "PlaneIF.H"
#include "IntersectionIF.H"
#include "EBAMRIO.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

//vofs to print out in detail
void getDebugIVS(IntVectSet& a_ivs,
                 const bool& a_thisCalcCoar)
{
  a_ivs = IntVectSet();
  if(!a_thisCalcCoar)
    {
      IntVect ivdebugloFine(D_DECL(66,10,10));
      IntVect ivdebughiFine(D_DECL(68,11,10));
      a_ivs |= ivdebugloFine;
      a_ivs |= ivdebughiFine;
    }
}

/***************/
/***************/
void
makeFinestDomain(Box& a_domain,
                 Real& a_dx)
{
  //parse input file.  single level
  ParmParse pp;
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
}
void
makeGeometry(const Box& a_domain,
             Real& a_dx)
{
  RealVect origin = RealVect::Zero;
  ParmParse  pp;
  int whichGeom;
  pp.get("which_geom", whichGeom);
  Real channelRadius;
  pp.get("channel_radius", channelRadius);
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  int biggridsize;
  pp.get("max_grid_size", biggridsize);
  int maxcoarsen = 0;
  if(whichGeom == 0)
    {
      pout() << "using all reg geom" << endl;
      AllRegularService allreg;
      ebisPtr->define(a_domain, origin, a_dx, allreg, biggridsize, maxcoarsen);
    }
  else if(whichGeom == 1)
    {
      RealVect channelNormal;
      vector<Real>  channelNormalVect(SpaceDim);
      pp.getarr("channel_normal",channelNormalVect, 0, SpaceDim);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          channelNormal[idir] = channelNormalVect[idir];
        }
      Real alpha;
      pp.get("alpha", alpha);

      RealVect channelPoint = RealVect::Zero;
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (channelNormal[idir] != 0.0)
        {
          channelPoint[idir] = alpha / channelNormal[idir];
          break;
        }
      }

      bool inside = false;

      PlaneIF ramp(channelNormal,channelPoint,inside);

      pout() << "using a ramp" << endl;

      RealVect vectDx = RealVect::Unit;
      vectDx *= a_dx;

      GeometryShop workshop(ramp,0,vectDx);
      ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize, maxcoarsen);
    }
  else if(whichGeom == 2)
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
      ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize, maxcoarsen);
    }
  else if(whichGeom == 3)
    {
      vector<Real>  cylinderAxisVect(SpaceDim);
      pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
      RealVect cylinderAxis;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          cylinderAxis[idir] = cylinderAxisVect[idir];
        }
      Real sum;
      PolyGeom::unifyVector(cylinderAxis, sum);

      int cylinderType;
      pp.get("cylinder_type", cylinderType);
      RealVect corner = RealVect::Zero;
      bool negativeInside = true;
      if(cylinderType == 0)
        {
          pout() << "using a tilted circular cylinder" << endl;

          TiltedCylinderIF tunnel(channelRadius, cylinderAxis, corner, negativeInside);
          RealVect vectDx = RealVect::Unit;
          vectDx *= a_dx;
          GeometryShop workshop(tunnel,0,vectDx);
          ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize, 0);
        }
      else if(cylinderType == 1)
        {
          pout() << "using a tilted square cylinder" << endl;
          //make square cylinder around x axis
          SquareCylinderIF sqcyl(channelRadius, negativeInside);
          TransformIF transformif(sqcyl);
          RealVect xAxis= BASISREALV(0);
          transformif.rotate(xAxis, cylinderAxis);
          transformif.translate(corner);
          RealVect vectDx = RealVect::Unit;
          vectDx *= a_dx;
          GeometryShop workshop(transformif,0,vectDx);
          ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize, 0);
        }
      else
        {
          MayDay::Error("bogus cylinder type");
        }
    }
  else
    {
      MayDay::Error("bogus geometry flag");
    }

}
/************/
void
faceNorm(Real& a_fluxfabNorm,
         Real& a_covMinuNorm,
         Real& a_covPlusNorm,
         FaceIndex& a_faceMaxOpen,
         VolIndex&  a_vofMaxMinu,
         VolIndex&  a_vofMaxPlus,
         const EBFaceFAB& a_fferror,
         const BaseIVFAB<Real>& a_coveredErrorMinu,
         const BaseIVFAB<Real>& a_coveredErrorPlus,
         const Vector<VolIndex>& a_coveredFacesMinu,
         const Vector<VolIndex>& a_coveredFacesPlus,
         const EBISBox& a_ebisBox,
         const int& a_comp,
         const int& a_idir,
         const int& a_normtype)
{
 CH_assert(a_normtype >= 0);
 CH_assert(a_normtype <= 2);
  a_fluxfabNorm = 0.0;
  a_covMinuNorm = 0.0;
  a_covPlusNorm = 0.0;
  IntVectSet ivs(a_fferror.getCellRegion());
  FaceIterator faceit(ivs, a_ebisBox.getEBGraph(), a_idir,
                      FaceStop::SurroundingNoBoundary);
  if(a_normtype == 0)
    {
      //max norm
      for(faceit.reset(); faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          Real valFace = Abs(a_fferror(face, a_comp));
          if(valFace > a_fluxfabNorm)
            {
              a_fluxfabNorm = valFace;
              a_faceMaxOpen = face;
            }
        }
      for(int ivof = 0; ivof < a_coveredFacesMinu.size(); ivof++)
        {
          const VolIndex& vof = a_coveredFacesMinu[ivof];
          Real valVoF = Abs(a_coveredErrorMinu(vof, a_comp));
          if(valVoF > a_covMinuNorm)
            {
              a_covMinuNorm = valVoF;
              a_vofMaxMinu = vof;
            }
        }
      for(int ivof = 0; ivof < a_coveredFacesPlus.size(); ivof++)
        {
          const VolIndex& vof = a_coveredFacesPlus[ivof];
          Real valVoF = Abs(a_coveredErrorPlus(vof, a_comp));
          if(valVoF > a_covPlusNorm)
            {
              a_covPlusNorm = valVoF;
              a_vofMaxPlus  = vof;
            }
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
          Real valFace = a_fferror(face, a_comp);
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

      a_fluxfabNorm = normTot;

      areaTot = 0.0;
      normTot = 0.0;
      for(int ivof = 0; ivof < a_coveredFacesMinu.size(); ivof++)
        {
          const VolIndex& vof = a_coveredFacesMinu[ivof];
          Real valVoF = a_coveredErrorMinu(vof, a_comp);
          //area frac one for covered face i guess
          areaTot += 1.0;
          if(a_normtype == 1)
            normTot += Abs(valVoF);
          else
            normTot += valVoF*valVoF;
        }
      if(areaTot > 1.0e-8)
        normTot /= areaTot;
      if(a_normtype == 2)
        normTot = sqrt(normTot);
      a_covMinuNorm = normTot;

      areaTot = 0.0;
      normTot = 0.0;
      for(int ivof = 0; ivof < a_coveredFacesPlus.size(); ivof++)
        {
          const VolIndex& vof = a_coveredFacesPlus[ivof];
          Real valVoF = a_coveredErrorPlus(vof, a_comp);
          //area frac one for covered face i guess
          areaTot += 1.0;
          if(a_normtype == 1)
            normTot += Abs(valVoF);
          else
            normTot += valVoF*valVoF;
        }
      if(areaTot > 1.0e-8)
        normTot /= areaTot;
      if(a_normtype == 2)
        normTot = sqrt(normTot);
      a_covPlusNorm = normTot;

    }

}
/************/
void
compareError(const EBFluxFAB&       a_ffErrorFine,
             const EBFluxFAB&       a_ffErrorCoar,
             const BaseIVFAB<Real>  a_coveredErrorMinuFine[SpaceDim],
             const BaseIVFAB<Real>  a_coveredErrorPlusFine[SpaceDim],
             const BaseIVFAB<Real>  a_coveredErrorMinuCoar[SpaceDim],
             const BaseIVFAB<Real>  a_coveredErrorPlusCoar[SpaceDim],
             const Vector<VolIndex> a_coveredFacesMinuFine[SpaceDim],
             const Vector<VolIndex> a_coveredFacesPlusFine[SpaceDim],
             const Vector<VolIndex> a_coveredFacesMinuCoar[SpaceDim],
             const Vector<VolIndex> a_coveredFacesPlusCoar[SpaceDim],
             const IntVectSet&      a_ivsIrregFine,
             const IntVectSet&      a_ivsIrregCoar,
             const EBISBox&         a_ebisBoxFine,
             const EBISBox&         a_ebisBoxCoar)
{
  int ncomp = a_ffErrorFine[0].nComp();

  pout() << "==============================================" << endl;
  pout() << "face flux test "  << endl;
  for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      pout() << "==============================================" << endl;
      pout() << "Direction = " << faceDir << endl;
      for(int comp = 0; comp < ncomp; comp++)
        {
          pout() << "==============================================" << endl;
          pout() << "component = " << comp << endl;
          pout() << "==============================================" << endl;
          for(int inorm = 0; inorm <= 2; inorm++)
            {
              if(inorm == 0)
                {
                  pout()  << "Using max norm." << endl;
                }
              else
                {
                  pout()  << "Using L-" << inorm << " norm." << endl;
                }

              Real fluxfabNormFine;
              Real covMinuNormFine;
              Real covPlusNormFine;
              Real fluxfabNormCoar;
              Real covMinuNormCoar;
              Real covPlusNormCoar;

              FaceIndex faceMaxOpenFine;
              VolIndex  vofMaxMinuFine;
              VolIndex  vofMaxPlusFine;
              FaceIndex faceMaxOpenCoar;
              VolIndex  vofMaxMinuCoar;
              VolIndex  vofMaxPlusCoar;
              faceNorm(fluxfabNormFine,
                       covMinuNormFine,
                       covPlusNormFine,
                       faceMaxOpenFine,
                       vofMaxMinuFine,
                       vofMaxPlusFine,
                       a_ffErrorFine[faceDir],
                       a_coveredErrorMinuFine[faceDir],
                       a_coveredErrorPlusFine[faceDir],
                       a_coveredFacesMinuFine[faceDir],
                       a_coveredFacesPlusFine[faceDir],
                       a_ebisBoxFine, comp, faceDir, inorm);

              faceNorm(fluxfabNormCoar,
                       covMinuNormCoar,
                       covPlusNormCoar,
                       faceMaxOpenCoar,
                       vofMaxMinuCoar,
                       vofMaxPlusCoar,
                       a_ffErrorCoar[faceDir],
                       a_coveredErrorMinuCoar[faceDir],
                       a_coveredErrorPlusCoar[faceDir],
                       a_coveredFacesMinuCoar[faceDir],
                       a_coveredFacesPlusCoar[faceDir],
                       a_ebisBoxCoar, comp, faceDir, inorm);

              pout() << "Coarse Open Flux Norm = " << fluxfabNormCoar << endl;
              pout() << "Fine   Open Flux Norm = " << fluxfabNormFine << endl;
              if((Abs(fluxfabNormCoar) > 1.0e-12) && (Abs(fluxfabNormFine) > 1.0e-12))
                {
                  Real order = log(fluxfabNormCoar/fluxfabNormFine)/log(2.0);
                  pout() << "Order of open fluxes = " << order << endl << endl;
                  if(inorm == 0)
                    {
                      pout() << "Fine Face With Max Error  = " << faceMaxOpenFine.gridIndex(Side::Lo) << faceMaxOpenFine.gridIndex(Side::Hi)<< endl;
                      pout() << "Coar Face With Max Error  = " << faceMaxOpenCoar.gridIndex(Side::Lo) << faceMaxOpenCoar.gridIndex(Side::Hi)<< endl;

                    }
                }
              if(a_coveredFacesMinuCoar[faceDir].size() == 0)
                {
                  pout() << "no covered minu faces coar" << endl;
                }
              else
                {
                  pout() << "Coarse Covered Minu Error Norm = " << covMinuNormCoar << endl;
                }
              if(a_coveredFacesMinuFine[faceDir].size() == 0)
                {
                  pout() << "no covered minu faces fine" << endl;
                }
              else
                {
                  pout() << "Fine   Covered Minu Error Norm = " << covMinuNormFine << endl;
                }
              if((Abs(covMinuNormCoar) > 1.0e-12) && (Abs(covMinuNormFine) > 1.0e-12))
                {
                  Real order = log(covMinuNormCoar/covMinuNormFine)/log(2.0);
                  pout() << "Order of covered minu  fluxes = " << order << endl  << endl;;
                  if(inorm == 0)
                    {
                      pout() << "fine minu max error vof  = " << vofMaxMinuFine  << endl;
                      pout() << "coar minu max error vof  = " << vofMaxMinuCoar << endl;

                    }
                }
              if(a_coveredFacesPlusCoar[faceDir].size() == 0)
                {
                  pout() << "no covered plus faces coar" << endl;
                }
              else
                {
                  pout() << "Coarse Covered Plus Error Norm = " << covPlusNormCoar << endl;
                }
              if(a_coveredFacesPlusFine[faceDir].size() == 0)
                {
                  pout() << "no covered plus faces fine" << endl;
                }
              else
                {
                  pout() << "Fine   Covered Plus Error Norm = " << covPlusNormFine << endl;
                }
              if((Abs(covPlusNormCoar) > 1.0e-12) && (Abs(covPlusNormFine) > 1.0e-12))
                {
                  Real order = log(covPlusNormCoar/covPlusNormFine)/log(2.0);
                  pout() << "Order of covered plus fluxes = " << order << endl << endl;
                  if(inorm == 0)
                    {
                      pout() << "fine plus max error vof  = " << vofMaxPlusFine  << endl;
                      pout() << "coar plus max error vof  = " << vofMaxPlusCoar << endl;

                    }
                }
            }//end loop over norms
        }// end loop over components
    } //end loop over face directions
}
/************/
void
printOpenSelect(const EBFluxFAB&        a_ffError,
                const EBISBox&          a_ebisBox,
                const IntVectSet&       a_ivs,
                int ncomp)
{
  IntVectSet ivs = a_ivs;
  ivs &= a_ebisBox.getRegion();

  pout() << setw(12)
         << setprecision(6)
         << setiosflags(ios::showpoint)
         << setiosflags(ios::scientific) ;
  for(int idir= 0; idir < SpaceDim; idir++)
    {
      if(!ivs.isEmpty())
        {
          pout() << "debugging faces for direction " << idir << endl;
          for(FaceIterator faceit(ivs, a_ebisBox.getEBGraph(), idir,
                                  FaceStop::SurroundingWithBoundary);
              faceit.ok(); ++faceit)
            {
              const FaceIndex& face = faceit();
              pout() << face.gridIndex(Side::Lo) << "  " << face.gridIndex(Side::Hi) << "   " ;
              for(int icomp = 0; icomp < ncomp; icomp++)
                {
                  Real error = Abs(a_ffError[idir](face, icomp));
                  pout() << error << "  ";
                }
              pout() << endl;
            }
        }
    }
}
/************/
void
printCoveredSelect(const BaseIVFAB<Real>  a_coveredError[SpaceDim],
                   const Vector<VolIndex> a_coveredFaces[SpaceDim],
                   const EBISBox&         a_ebisBox,
                   const IntVectSet&      a_ivsDebug,
                   int ncomp)
{
  pout() << setw(12)
         << setprecision(6)
         << setiosflags(ios::showpoint)
         << setiosflags(ios::scientific) ;
  for(int idir= 0; idir < SpaceDim; idir++)
    {
      if(!a_ivsDebug.isEmpty())
        {
          pout() << "debugging covered faces for direction " << idir << endl;
          for(int ivof= 0; ivof < a_coveredFaces[idir].size(); ivof++)
            {
              const VolIndex& vof = a_coveredFaces[idir][ivof];
              if(a_ivsDebug.contains(vof.gridIndex()))
                {
                  pout() << vof << "   " ;
                  for(int icomp = 0; icomp < ncomp; icomp++)
                    {
                      Real error = Abs(a_coveredError[idir](vof, icomp));
                      pout() << error << "  ";
                    }
                  pout() << endl;
                }
            }
        }
    }
}
/************/
void
subtract(BaseIVFAB<Real>& a_error,
         const BaseIVFAB<Real>& a_calc,
         const BaseIVFAB<Real>& a_exact,
         const Vector<VolIndex> a_faces)
{
  int nCons = a_error.nComp();
 CH_assert(a_calc.nComp() == nCons);
 CH_assert(a_exact.nComp() == nCons);
  for(int ivof = 0; ivof < a_faces.size(); ivof++)
    {
      const VolIndex& vof = a_faces[ivof];
      for(int ivar = 0; ivar < nCons; ivar++)
        {
          a_error(vof, ivar) =  a_calc(vof, ivar) -a_exact(vof, ivar);
        }
    }
}
/************/
void
fluxAutopsy()
{
  //make layouts == domain
  Box domainBoxFine, domainBoxCoar;
  Real dxFine, dxCoar;
  bool is_periodic[SpaceDim];
  for (int i = 0; i < SpaceDim; ++i)
    {
      is_periodic[i] = false;
    }

  makeFinestDomain(domainBoxFine, dxFine);
  domainBoxCoar = coarsen(domainBoxFine, 2);
  dxCoar = 2.*dxFine;
  ProblemDomain domainFine(domainBoxFine.smallEnd(),
                           domainBoxFine.bigEnd(),
                           is_periodic);
  ProblemDomain domainCoar(domainBoxCoar.smallEnd(),
                           domainBoxCoar.bigEnd(),
                           is_periodic);

  Vector<Box> boxFine(1, domainBoxFine);
  Vector<Box> boxCoar(1, domainBoxCoar);
  Vector<int> proc(1, 0);
  DisjointBoxLayout dblFine(boxFine, proc);
  DisjointBoxLayout dblCoar;
  coarsen(dblCoar, dblFine, 2);
  //source is a placeholder.
  EBCellFAB source;
  EBCellFAB  slopePrimCoar[SpaceDim], slopeNLimCoar[SpaceDim];
  EBCellFAB  slopePrimFine[SpaceDim], slopeNLimFine[SpaceDim];

  RealVect modianoAxis, center;
  Real waveAmp, waveWidth, gamma;

  ParmParse pp;
  pp.get("wave_amplitude", waveAmp);
  pp.get("wave_width", waveWidth);
  pp.get("gamma",gamma);

  int idoNegWave;
  pp.get("use_negative_wave",idoNegWave);
  bool useNegativeWave = (idoNegWave==1);

  if(useNegativeWave)
    {
      pout() << "u - c wave " << endl;
    }
  else
    {
      pout() << "u + c wave " << endl;
    }
  int idoFreeStream;
  pp.get("free_stream_prob",idoFreeStream);
  bool doFreeStream = (idoFreeStream==1);
  if(doFreeStream)
    {
      pout() << "free stream modiano bump " << endl;
    }
  else
    {
      pout() << "simple wave prob " << endl;
    }

  vector<Real> centerpp(SpaceDim,0.5);
  pp.getarr("wave_center",centerpp,0,SpaceDim);

  for(int idir = 0; idir < SpaceDim; idir++)
    {
      center[idir] = centerpp[idir];
    }
  int whichGeom;
  pp.get("which_geom", whichGeom);
  if(whichGeom == 3)
    {
      pout() << "uses cylinder axis for direction" << endl;
      vector<Real>  cylinderAxisVect(SpaceDim);
      pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
      RealVect cylinderAxis;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          cylinderAxis[idir] = cylinderAxisVect[idir];
        }
      Real sum;
      PolyGeom::unifyVector(cylinderAxis, sum);
      modianoAxis = cylinderAxis;
    }
  else
    {
      pout() << "uses channel slope for direction" << endl;
      RealVect channelNormal;
      vector<Real>  channelNormalVect(SpaceDim);
      pp.getarr("channel_normal",channelNormalVect, 0, SpaceDim);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          channelNormal[idir] = channelNormalVect[idir];
        }
      modianoAxis = mibcGetNormalVector(channelNormal);;
    }


  ModianoIBC ibcFine(gamma, waveAmp, waveWidth, center, modianoAxis, doFreeStream, useNegativeWave);
  ModianoIBC ibcCoar(gamma, waveAmp, waveWidth, center, modianoAxis, doFreeStream, useNegativeWave);
  ModianoIBCFactory ibcFactFine(gamma, waveAmp, waveWidth, center, modianoAxis, doFreeStream, useNegativeWave);
  ModianoIBCFactory ibcFactCoar(gamma, waveAmp, waveWidth, center, modianoAxis, doFreeStream, useNegativeWave);

  //make ebislayouts
  //no need for ghost since we are using the domain
  //for the box
  Real cfl;
  pp.get("cfl", cfl);

  EBISLayout ebislFine, ebislCoar;
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();

  makeGeometry(domainBoxFine, dxFine);
  ebisPtr->fillEBISLayout(ebislFine, dblFine, domainBoxFine, 0);
  makeGeometry(domainBoxCoar, dxCoar);
  ebisPtr->fillEBISLayout(ebislCoar, dblCoar, domainBoxCoar, 0);

  int nCons = CNUM;
  int ifourth, iflatten, iartvisc;
  pp.get("use_fourth_order_slopes", ifourth);
  pp.get("use_flattening"         , iflatten);
  pp.get("use_art_visc"           , iartvisc);
  bool useFourthOrderSlopes = (ifourth  ==1);
  bool useFlattening        = false;
  bool useArtificialVisc    = (iartvisc ==1);

  for(DataIterator dit = dblFine.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBoxFine = ebislFine[dit()];
      const EBISBox& ebisBoxCoar = ebislCoar[dit()];

      //these get filled by modianoibc
      EBCellFAB consStateFine(ebisBoxFine, domainBoxFine, nCons);
      EBCellFAB consStateCoar(ebisBoxCoar, domainBoxCoar, nCons);

      ibcFine.define(domainFine, dxFine);
      ibcCoar.define(domainCoar, dxCoar);
      ibcFine.initialize(consStateFine, ebisBoxFine);
      ibcCoar.initialize(consStateCoar, ebisBoxCoar);

      EBPatchPolytropic ebppFine, ebppCoar;
      ebppFine.define(domainFine, dxFine);
      ebppCoar.define(domainCoar, dxCoar);
      //get the max wave speed and set the time step
      //cfivs is empty because this is a single-level calc
      IntVectSet cfivs;
      ebppFine.setValidBox(domainBoxFine, ebisBoxFine, cfivs, 0., 0.);
      ebppCoar.setValidBox(domainBoxCoar, ebisBoxCoar, cfivs, 0., 0.);
      ebppFine.setGamma(gamma);
      ebppCoar.setGamma(gamma);

      EBCellFAB& primStateFine = ebppFine.getPrimState();
      EBCellFAB& primStateCoar = ebppCoar.getPrimState();
      bool useLimiting = false;
      ebppFine.setSlopeParameters(useFourthOrderSlopes, useFlattening, useLimiting);
      ebppCoar.setSlopeParameters(useFourthOrderSlopes, useFlattening, useLimiting);
      ebppFine.doRZCoords(false);
      ebppCoar.doRZCoords(false);
      ebppFine.artificialViscosity(useArtificialVisc);
      ebppCoar.artificialViscosity(useArtificialVisc);
      ebppFine.setEBPhysIBC(ibcFactFine);
      ebppCoar.setEBPhysIBC(ibcFactCoar);
      Real maxWaveSpeedFine=
        ebppFine.getMaxWaveSpeed(consStateFine, domainBoxFine);
      Real maxWaveSpeedCoar=
        ebppCoar.getMaxWaveSpeed(consStateCoar, domainBoxCoar);
      Real dtFine = cfl*dxFine/maxWaveSpeedFine;
      Real dtCoar = cfl*dxCoar/maxWaveSpeedCoar;
      ebppFine.setValidBox(domainBoxFine, ebisBoxFine, cfivs, 0., dtFine);
      ebppCoar.setValidBox(domainBoxCoar, ebisBoxCoar, cfivs, 0., dtCoar);

      //To make the exact stuff work the same as
      //the  computed stuff, define the EBFluxFAB data
      //before going into routines and define BaseIVFAB
      //data inside routines.
      //computed fluxes on grid
      EBFluxFAB  cfluxFine(ebisBoxFine, domainBoxFine, nCons);
      EBFluxFAB  cfluxCoar(ebisBoxCoar, domainBoxCoar, nCons);
      //exact fluxes on grid
      EBFluxFAB  exactFine(ebisBoxFine, domainBoxFine, nCons);
      EBFluxFAB  exactCoar(ebisBoxCoar, domainBoxCoar, nCons);
      cfluxFine.setVal(0.);
      cfluxCoar.setVal(0.);
      exactFine.setVal(0.);
      exactCoar.setVal(0.);
      //computed covered fluxes
      BaseIVFAB<Real>*  coveredCfluxMinuFine= ebppFine.getCoveredFluxMinu();
      BaseIVFAB<Real>*  coveredCfluxPlusFine= ebppFine.getCoveredFluxPlus();
      BaseIVFAB<Real>*  coveredCfluxMinuCoar= ebppCoar.getCoveredFluxMinu();
      BaseIVFAB<Real>*  coveredCfluxPlusCoar= ebppCoar.getCoveredFluxPlus();

      //vofs over which covered vofs live
      Vector<VolIndex>* coveredFacesMinuFine= ebppFine.getCoveredFaceMinu();
      Vector<VolIndex>* coveredFacesPlusFine= ebppFine.getCoveredFacePlus();
      Vector<VolIndex>* coveredFacesMinuCoar= ebppCoar.getCoveredFaceMinu();
      Vector<VolIndex>* coveredFacesPlusCoar= ebppCoar.getCoveredFacePlus();

      //exact covered fluxes
      BaseIVFAB<Real>  coveredExactMinuFine[SpaceDim];
      BaseIVFAB<Real>  coveredExactPlusFine[SpaceDim];
      BaseIVFAB<Real>  coveredExactMinuCoar[SpaceDim];
      BaseIVFAB<Real>  coveredExactPlusCoar[SpaceDim];



      //irregular sets for eb fluxes
      IntVectSet ivsIrregFine = ebisBoxFine.getIrregIVS(domainBoxFine);
      IntVectSet ivsIrregCoar = ebisBoxCoar.getIrregIVS(domainBoxCoar);

      BaseIVFAB<Real> ebPressExactFine(ivsIrregFine, ebisBoxFine.getEBGraph(), 1);
      BaseIVFAB<Real> ebPressExactCoar(ivsIrregCoar, ebisBoxCoar.getEBGraph(), 1);
      ebPressExactFine.setVal(0.);
      ebPressExactCoar.setVal(0.);

      EBCellFAB flat;
      bool verbose = false;

      pout() << "computing fine fluxes " << endl;
      ebppFine.computeFluxes(cfluxFine,
                             coveredCfluxMinuFine, coveredCfluxPlusFine,
                             coveredFacesMinuFine, coveredFacesPlusFine,
                             primStateFine, slopePrimFine, slopeNLimFine,
                             flat, consStateFine, source, domainBoxFine,
                             dit(),verbose);

      pout() << "computing coarse fluxes " << endl;
      ebppCoar.computeFluxes(cfluxCoar,
                             coveredCfluxMinuCoar, coveredCfluxPlusCoar,
                             coveredFacesMinuCoar, coveredFacesPlusCoar,
                             primStateCoar, slopePrimCoar, slopeNLimCoar,
                             flat, consStateCoar, source, domainBoxCoar,
                             dit(),verbose);

      IntVectSet coveredSetsMinuFine[SpaceDim];
      IntVectSet coveredSetsPlusFine[SpaceDim];
      IntVectSet coveredSetsMinuCoar[SpaceDim];
      IntVectSet coveredSetsPlusCoar[SpaceDim];
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          coveredSetsMinuFine[idir]=coveredCfluxMinuFine[idir].getIVS();
          coveredSetsMinuCoar[idir]=coveredCfluxMinuCoar[idir].getIVS();
          coveredSetsPlusFine[idir]=coveredCfluxPlusFine[idir].getIVS();
          coveredSetsPlusCoar[idir]=coveredCfluxPlusCoar[idir].getIVS();
        }
      pout() << "computing exact fluxes fields " << endl;
      //get the exact fluxes
      ibcFine.computeExactFluxes(exactFine,
                                 coveredExactMinuFine, coveredExactPlusFine,
                                 coveredFacesMinuFine, coveredFacesPlusFine,
                                 coveredSetsMinuFine, coveredSetsPlusFine,
                                 ebPressExactFine, ivsIrregFine,
                                 ebisBoxFine, domainBoxFine,
                                 domainBoxFine, 0.5*dtFine);

      ibcCoar.computeExactFluxes(exactCoar,
                                 coveredExactMinuCoar, coveredExactPlusCoar,
                                 coveredFacesMinuCoar, coveredFacesPlusCoar,
                                 coveredSetsMinuCoar, coveredSetsPlusCoar,
                                 ebPressExactCoar, ivsIrregCoar,
                                 ebisBoxCoar, domainBoxCoar,
                                 domainBoxCoar,0.5*dtCoar);

      //subtract computed from exact to get error
      EBFluxFAB fferrorFine(ebisBoxFine, domainBoxFine, nCons);
      EBFluxFAB fferrorCoar(ebisBoxCoar, domainBoxCoar, nCons);
      fferrorFine.setVal(0.);
      fferrorCoar.setVal(0.);
      BaseIVFAB<Real>  coveredErrorMinuFine[SpaceDim];
      BaseIVFAB<Real>  coveredErrorPlusFine[SpaceDim];
      BaseIVFAB<Real>  coveredErrorMinuCoar[SpaceDim];
      BaseIVFAB<Real>  coveredErrorPlusCoar[SpaceDim];
      Interval interv(0, nCons-1);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          //copy from the computed solution
          fferrorFine[idir].copy(domainBoxFine, interv,
                                 domainBoxFine, cfluxFine[idir], interv);
          fferrorCoar[idir].copy(domainBoxCoar, interv,
                                 domainBoxCoar, cfluxCoar[idir], interv);
          //subtract off exact solution
          fferrorFine[idir] -= exactFine[idir];
          fferrorCoar[idir] -= exactCoar[idir];

          //define the sparse data holders
          coveredErrorMinuFine[idir].define(coveredSetsMinuFine[idir], ebisBoxFine.getEBGraph(), nCons);
          coveredErrorPlusFine[idir].define(coveredSetsPlusFine[idir], ebisBoxFine.getEBGraph(), nCons);
          coveredErrorMinuCoar[idir].define(coveredSetsMinuCoar[idir], ebisBoxCoar.getEBGraph(), nCons);
          coveredErrorPlusCoar[idir].define(coveredSetsPlusCoar[idir], ebisBoxCoar.getEBGraph(), nCons);
          coveredErrorMinuFine[idir].setVal(0.);
          coveredErrorPlusFine[idir].setVal(0.);
          coveredErrorMinuCoar[idir].setVal(0.);
          coveredErrorPlusCoar[idir].setVal(0.);

          subtract(coveredErrorMinuFine[idir], coveredCfluxMinuFine[idir], coveredExactMinuFine[idir], coveredFacesMinuFine[idir]);
          subtract(coveredErrorPlusFine[idir], coveredCfluxPlusFine[idir], coveredExactPlusFine[idir], coveredFacesPlusFine[idir]);
          subtract(coveredErrorMinuCoar[idir], coveredCfluxMinuCoar[idir], coveredExactMinuCoar[idir], coveredFacesMinuCoar[idir]);
          subtract(coveredErrorPlusCoar[idir], coveredCfluxPlusCoar[idir], coveredExactPlusCoar[idir], coveredFacesPlusCoar[idir]);

        }


      pout() << "comparing error " << endl;
      //do the whole convergence test thing.
      compareError(fferrorFine, fferrorCoar,
                   coveredErrorMinuFine, coveredErrorPlusFine,
                   coveredErrorMinuCoar, coveredErrorPlusCoar,
                   coveredFacesMinuFine, coveredFacesPlusFine,
                   coveredFacesMinuCoar, coveredFacesPlusCoar,
                   ivsIrregFine, ivsIrregCoar,
                   ebisBoxFine, ebisBoxCoar);

      // debug vofs== vofs to dump open flux errors of
      IntVectSet ivsDebugFine, ivsDebugCoar;
      getDebugIVS(ivsDebugFine, false);
      getDebugIVS(ivsDebugCoar, true);

      pout() << "printing open debug faces fine " << endl;
      printOpenSelect(fferrorFine, ebisBoxFine, ivsDebugFine, nCons);
      pout() << "printing open debug faces coar " << endl;
      printOpenSelect(fferrorCoar, ebisBoxCoar, ivsDebugCoar, nCons);

      pout() << "printing covered debug faces minu fine " << endl;
      printCoveredSelect(coveredErrorMinuFine,
                         coveredFacesMinuFine,
                         ebisBoxFine, ivsDebugFine, nCons);
      pout() << "printing covered debug faces plus fine " << endl;
      printCoveredSelect( coveredErrorPlusFine,
                          coveredFacesPlusFine,
                         ebisBoxFine, ivsDebugFine, nCons);

      pout() << "printing covered debug faces minu coar " << endl;
      printCoveredSelect(coveredErrorMinuCoar,
                         coveredFacesMinuCoar,
                         ebisBoxCoar, ivsDebugCoar, nCons);

      pout() << "printing covered debug faces plus coar " << endl;
      printCoveredSelect(coveredErrorPlusCoar,
                         coveredFacesPlusCoar,
                         ebisBoxCoar, ivsDebugCoar, nCons);

    }

}
/************/
int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
  {
    // setChomboMPIErrorHandler();
#endif
    // Check for an input file
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

    fluxAutopsy();

#ifdef CH_MPI
    dumpmemoryatexit();
  }
  MPI_Finalize();
#endif
}
/************/

