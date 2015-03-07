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
#include "EBLevelDataOps.H"

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
#include "SphereIF.H"
#include "PlaneIF.H"
#include "IntersectionIF.H"
#include "EBAMRIO.H"
#include "covExtrapF_F.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

/***************/
/***************/
void
makeFinestDomain(Box&      a_domain,
                 RealVect& a_dx)
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
  pp.get("domain_length",prob_hi);
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      a_dx[idir] = prob_hi/n_cell[idir];
    }
  pout() << "dx = " << a_dx << "; aspect ratio = " << a_dx[0]/a_dx[SpaceDim-1] << endl;
}
void
makeGeometry(const Box& a_domain,
             RealVect&  a_dx)
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
  const int verbosity = 0;
  if(whichGeom == 0)
    {
      pout() << "using all reg geom" << endl;
      AllRegularService allreg;
      ebisPtr->define(a_domain, origin, a_dx[0], allreg, biggridsize, maxcoarsen);
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
      GeometryShop workshop(ramp,verbosity,a_dx);
      ebisPtr->define(a_domain, origin, a_dx[0], workshop, biggridsize, maxcoarsen);
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
      GeometryShop workshop(channel,verbosity,a_dx);
      ebisPtr->define(a_domain, origin, a_dx[0], workshop, biggridsize, maxcoarsen);
    }
  else if(whichGeom == 3)
    {
      vector<Real>  cylinderAxisVect(SpaceDim);
      vector<Real>  cylinderPoint(SpaceDim);
      pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
      pp.getarr("cylinder_point",cylinderPoint, 0, SpaceDim);
      RealVect cylinderAxis;
      RealVect point;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          cylinderAxis[idir] = cylinderAxisVect[idir];
          point[idir]        = cylinderPoint[idir];
        }
      Real sum;
      PolyGeom::unifyVector(cylinderAxis, sum);

      int cylinderType;
      pp.get("cylinder_type", cylinderType);
      bool negativeInside = false;
      if(cylinderType == 0)
        {
          pout() << "using a tilted circular cylinder" << endl;

          TiltedCylinderIF tunnel(channelRadius, cylinderAxis, point, negativeInside);
          GeometryShop workshop(tunnel,verbosity,a_dx);
          ebisPtr->define(a_domain, origin, a_dx[0], workshop, biggridsize, 0);
        }
      else if(cylinderType == 1)
        {
          pout() << "using a tilted square cylinder" << endl;
          //make square cylinder around x axis
          SquareCylinderIF sqcyl(channelRadius, negativeInside);
          TransformIF transformif(sqcyl);
          RealVect xAxis= BASISREALV(0);
          transformif.rotate(xAxis, cylinderAxis);
          transformif.translate(point);
          GeometryShop workshop(transformif,verbosity,a_dx);
          ebisPtr->define(a_domain, origin, a_dx[0], workshop, biggridsize, 0);
        }
      else
        {
          MayDay::Error("bogus cylinder type");
        }
    }
  else if(whichGeom == 5)
    {
      pout() << "sphere geometry" << endl;
      vector<Real> sphere_center(SpaceDim);
      pp.getarr("sphere_center",sphere_center, 0, SpaceDim);
      RealVect sphereCenter;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          sphereCenter[idir] = sphere_center[idir];
        }
      Real sphereRadius;
      pp.get("sphere_radius", sphereRadius);

      bool     insideRegular = false;

      int query = pp.query("inside",insideRegular);
      if(query==0)
        {
          insideRegular = false;
        }
      SphereIF implicit(sphereRadius,sphereCenter,insideRegular);

      GeometryShop workshop(implicit,verbosity,a_dx);
      //this generates the new EBIS
      ebisPtr->define(a_domain, origin, a_dx[0], workshop, biggridsize, 0);
    }
  else
    {
      MayDay::Error("bogus geometry flag");
    }

}
/************/
/************/
void
faceNorm(Real& a_covNorm,
         const BaseIVFAB<Real>&  a_coveredError,
         const Vector<VolIndex>& a_coveredFaces,
         const EBISBox& a_ebisBox,
         const int& a_comp,
         const int& a_normtype)
{
 CH_assert(a_normtype >= 0);
 CH_assert(a_normtype <= 2);
  a_covNorm = 0.0;
  if(a_normtype == 0)
    {
      for(int ivof = 0; ivof < a_coveredFaces.size(); ivof++)
        {
          const VolIndex& vof = a_coveredFaces[ivof];
          Real valVoF = a_coveredError(vof, a_comp);
          a_covNorm = Max(Abs(valVoF), a_covNorm);
        }
    }
  else
    {
      //integral norm
      Real areaTot = 0.0;
      Real normTot = 0.0;
      for(int ivof = 0; ivof < a_coveredFaces.size(); ivof++)
        {
          const VolIndex& vof = a_coveredFaces[ivof];
          Real valVoF = a_coveredError(vof, a_comp);
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
      a_covNorm = normTot;

    }

}
/************/
void
compareError(const BaseIVFAB<Real>&  a_coveredErrorFine,
             const BaseIVFAB<Real>&  a_coveredErrorCoar,
             const Vector<VolIndex>& a_coveredFacesFine,
             const Vector<VolIndex>& a_coveredFacesCoar,
             const EBISBox&          a_ebisBoxFine,
             const EBISBox&          a_ebisBoxCoar,
             const int&              a_faceDir,
             const Side::LoHiSide&   a_sd)
{
  //int ncomp = a_coveredErrorFine.nComp();
  int ncomp = 1;

  pout() << "==============================================" << endl;
  pout() << "covered face extrapolation test "  ;
  if(a_sd == Side::Hi)
    {
      pout() << "for faces covered on HIGH side of vof " ;
    }
  else
    {
      pout() << "for faces covered on LOW side of vof " ;
    }
  pout() << endl;

  pout() << "==============================================" << endl;
  pout() << "Direction = " << a_faceDir << endl;
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

          Real covNormFine;
          Real covNormCoar;

          faceNorm(covNormFine,
                   a_coveredErrorFine,
                   a_coveredFacesFine,
                   a_ebisBoxFine, comp, inorm);

          faceNorm(covNormCoar,
                   a_coveredErrorCoar,
                   a_coveredFacesCoar,
                   a_ebisBoxCoar, comp, inorm);

          if(a_coveredFacesCoar.size() == 0)
            {
              pout() << "no covered  faces coar" << endl;
            }
          else
            {
              pout() << "Coarse Covered  Error Norm = " << covNormCoar << endl;
            }
          if(a_coveredFacesFine.size() == 0)
            {
              pout() << "no covered  faces fine" << endl;
            }
          else
            {
              pout() << "Fine   Covered  Error Norm = " << covNormFine << endl;
            }
          if((Abs(covNormCoar) > 1.0e-12) && (Abs(covNormFine) > 1.0e-12))
            {
              Real order = log(covNormCoar/covNormFine)/log(2.0);
              pout() << "Order of covered  fluxes = " << order << endl << endl;
            }
        }//end loop over norms
    }
}

/************/
void
getExactPrim(EBCellFAB&      a_primStat,
             EBCellFAB&      a_primMinu,
             EBCellFAB&      a_primPlus,
             const EBISBox&  a_ebisBox,
             const Box&      a_domain,
             const int&      a_faceDir,
             const RealVect& a_dx)
{
  BaseFab<Real>& regPrimStat = a_primStat.getSingleValuedFAB();
  BaseFab<Real>& regPrimMinu = a_primMinu.getSingleValuedFAB();
  BaseFab<Real>& regPrimPlus = a_primPlus.getSingleValuedFAB();
  int ncomp = a_primStat.nComp();

  FORT_POLYNOMIALINIT(CHF_FRA(regPrimStat),
                      CHF_FRA(regPrimMinu),
                      CHF_FRA(regPrimPlus),
                      CHF_CONST_REALVECT(a_dx),
                      CHF_CONST_INT(a_faceDir),
                      CHF_CONST_INT(ncomp),
                      CHF_BOX(a_domain));

  IntVectSet ivsMulti = a_ebisBox.getMultiCells(a_domain);
  for(VoFIterator vofit(ivsMulti, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();

      for(int ivar = 0; ivar < ncomp; ivar++)
        {
          Real valcent, valplus, valminu;
          FORT_POINTPOLYEXACT(CHF_REAL(valcent),
                              CHF_REAL(valminu),
                              CHF_REAL(valplus),
                              CHF_CONST_INTVECT(iv),
                              CHF_CONST_REALVECT(a_dx),
                              CHF_CONST_INT(a_faceDir),
                              CHF_CONST_INT(ivar));
          a_primStat(vof, ivar) = valcent;
          a_primMinu(vof, ivar) = valminu;
          a_primPlus(vof, ivar) = valplus;
        }
    }
}
/************/
void
getError(BaseIVFAB<Real>&        a_error,
         const BaseIVFAB<Real>&  a_extendedPrim,
         const Vector<VolIndex>& a_coveredFace,
         const int&              a_faceDir,
         const Side::LoHiSide&   a_sd,
         const Box&              a_domain,
         const RealVect&         a_dx)
{

  int ncomp = a_extendedPrim.nComp();
  for(int ivof = 0; ivof < a_coveredFace.size(); ivof++)
    {
      const VolIndex& vof = a_coveredFace[ivof];
      const IntVect& iv = vof.gridIndex();
      for(int ivar = 0; ivar < ncomp; ivar++)
        {
          Real error;

          Real exaccent, exacplus, exacminu;
          FORT_POINTPOLYEXACT(CHF_REAL(exaccent),
                              CHF_REAL(exacminu),
                              CHF_REAL(exacplus),
                              CHF_CONST_INTVECT(iv),
                              CHF_CONST_REALVECT(a_dx),
                              CHF_CONST_INT(a_faceDir),
                              CHF_CONST_INT(ivar));

          Real exactVal;
          if(a_sd == Side::Hi)
            {
              exactVal = exacplus;
            }
          else
            {
              exactVal = exacminu;
            }

          Real calcVal  = a_extendedPrim(vof, ivar);
          error =  calcVal - exactVal;
          a_error(vof, ivar) = error;
        }
    }
}
/************/
void
coveredExtrapTest()
{
  //make layouts == domain
  Box domainBoxFine, domainBoxCoar;
  RealVect dxFine, dxCoar;
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

  //make ebislayouts
  //no need for ghost since we are using the domain
  //for the box
  EBISLayout ebislFine, ebislCoar;
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  makeGeometry(domainBoxFine, dxFine);
  ebisPtr->fillEBISLayout(ebislFine, dblFine, domainBoxFine, 0);
  makeGeometry(domainBoxCoar, dxCoar);
  ebisPtr->fillEBISLayout(ebislCoar, dblCoar, domainBoxCoar, 0);
  ParmParse pp;
  int ilimit, ifourth;
  pp.get("use_fourth_order_slopes", ifourth);
  pp.get("use_limiting"         , ilimit);
  bool useFourthOrderSlopes = (ifourth   ==1);
  bool useFlattening        = false;
  bool useLimiting          = (ilimit    ==1);
  Box interiorBoxFine = domainBoxFine;
  Box interiorBoxCoar = domainBoxCoar;
  //  interiorBoxFine.grow(-1);
  //  interiorBoxCoar.grow(-1);

  //make something we can view in ChomboVis (to check the geometry)
  //LevelData<EBCellFAB> geomFine;
  //EBCellFactory ebcellfactFine(ebislFine);
  //geomFine.define(dblFine,1,IntVect::Zero,ebcellfactFine);
  //EBLevelDataOps::setToZero(geomFine);
  //viewEBLevel(&geomFine);
  for(DataIterator dit = dblFine.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBoxFine = ebislFine[dit()];
      const EBISBox& ebisBoxCoar = ebislCoar[dit()];

      EBPatchPolytropic ebppFine, ebppCoar;
      ebppFine.define(domainFine, dxFine);
      ebppCoar.define(domainCoar, dxCoar);
      int nVar = ebppFine.numPrimitives();

      //cfivs is empty because this is a single-level calc
      IntVectSet cfivs;
      ebppFine.setValidBox(domainBoxFine, ebisBoxFine, cfivs, 0., 0.);
      ebppCoar.setValidBox(domainBoxCoar, ebisBoxCoar, cfivs, 0., 0.);
      ebppFine.setSlopeParameters(useFourthOrderSlopes, useFlattening, useLimiting);
      ebppCoar.setSlopeParameters(useFourthOrderSlopes, useFlattening, useLimiting);
      Real dtFine = 1.0;
      Real dtCoar = 2.0;
      ebppFine.setValidBox(domainBoxFine, ebisBoxFine, cfivs, 0., dtFine);
      ebppCoar.setValidBox(domainBoxCoar, ebisBoxCoar, cfivs, 0., dtCoar);

      EBCellFAB   primStatFine(ebisBoxFine, domainBoxFine, nVar);
      EBCellFAB   primStatCoar(ebisBoxCoar, domainBoxCoar, nVar);


      EBCellFAB   primMinuFine[SpaceDim];
      EBCellFAB   primMinuCoar[SpaceDim];
      EBCellFAB   primPlusFine[SpaceDim];
      EBCellFAB   primPlusCoar[SpaceDim];


      for(int faceDir = 0; faceDir < SpaceDim; faceDir++) //loop over face directions
        {
          primMinuFine[faceDir].define(ebisBoxFine, domainBoxFine, nVar);
          primMinuCoar[faceDir].define(ebisBoxCoar, domainBoxCoar, nVar);
          primPlusFine[faceDir].define(ebisBoxFine, domainBoxFine, nVar);
          primPlusCoar[faceDir].define(ebisBoxCoar, domainBoxCoar, nVar);

          primMinuFine[faceDir].setVal(0.);
          primMinuCoar[faceDir].setVal(0.);
          primPlusFine[faceDir].setVal(0.);
          primPlusCoar[faceDir].setVal(0.);

          getExactPrim(primStatFine,  primMinuFine[faceDir], primPlusFine[faceDir],
                       ebisBoxFine, domainBoxFine, faceDir, dxFine);
          getExactPrim(primStatCoar,  primMinuCoar[faceDir], primPlusCoar[faceDir],
                       ebisBoxCoar, domainBoxCoar, faceDir, dxCoar);
        }

      EBCellFAB flattening; //placeholder.  no flattening here

      EBCellFAB     slopesPrimFine[SpaceDim];
      EBCellFAB     slopesSecoFine[SpaceDim];
      EBCellFAB     slopesPrimCoar[SpaceDim];
      EBCellFAB     slopesSecoCoar[SpaceDim];
      for(int faceDir = 0; faceDir < SpaceDim; faceDir++) //loop over face directions
        {
          slopesPrimFine[faceDir].define(ebisBoxFine, domainBoxFine, nVar);
          slopesSecoFine[faceDir].define(ebisBoxFine, domainBoxFine, nVar);
          slopesPrimCoar[faceDir].define(ebisBoxCoar, domainBoxCoar, nVar);
          slopesSecoCoar[faceDir].define(ebisBoxCoar, domainBoxCoar, nVar);

          slopesPrimFine[faceDir].setVal(0.);
          slopesSecoFine[faceDir].setVal(0.);
          slopesPrimCoar[faceDir].setVal(0.);
          slopesSecoCoar[faceDir].setVal(0.);

          ebppFine.slope(slopesPrimFine[faceDir], slopesSecoFine[faceDir], primStatFine,
                         flattening, faceDir, domainBoxFine);
          ebppCoar.slope(slopesPrimCoar[faceDir], slopesSecoCoar[faceDir], primStatCoar,
                         flattening, faceDir, domainBoxCoar);
        }

      for(int faceDir = 0; faceDir < SpaceDim; faceDir++) //loop over face directions
        {
          for(SideIterator sit; sit.ok(); ++sit) //loop over which side is covered
            {
              IntVectSet        fulIrregIVSFine, fulIrregIVSCoar;
              Vector<VolIndex>  coveredFaceFine, coveredFaceCoar;
              IntVectSet        coveredSetsFine, coveredSetsCoar;
              ebppFine.computeCoveredFaces(coveredFaceFine,
                                           coveredSetsFine,
                                           fulIrregIVSFine,
                                           faceDir, sit(), interiorBoxFine);
              ebppCoar.computeCoveredFaces(coveredFaceCoar,
                                           coveredSetsCoar,
                                           fulIrregIVSCoar,
                                           faceDir, sit(), interiorBoxCoar);

              BaseIVFAB<Real>  extendedPrimFine(coveredSetsFine, ebisBoxFine.getEBGraph(), nVar);
              BaseIVFAB<Real>  extendedPrimCoar(coveredSetsCoar, ebisBoxCoar.getEBGraph(), nVar);

              extendedPrimFine.setVal(0.);
              extendedPrimCoar.setVal(0.);

              ebppFine.extrapToCoveredFaces(extendedPrimFine,
                                            primMinuFine[faceDir],
                                            primPlusFine[faceDir],
                                            primStatFine,
                                            coveredFaceFine,
                                            faceDir,
                                            sit(),
                                            interiorBoxFine);

              ebppCoar.extrapToCoveredFaces(extendedPrimCoar,
                                            primMinuCoar[faceDir],
                                            primPlusCoar[faceDir],
                                            primStatCoar,
                                            coveredFaceCoar,
                                            faceDir,
                                            sit(),
                                            interiorBoxCoar);

              BaseIVFAB<Real>  errorFine(coveredSetsFine, ebisBoxFine.getEBGraph(), nVar);
              BaseIVFAB<Real>  errorCoar(coveredSetsCoar, ebisBoxCoar.getEBGraph(), nVar);

              errorFine.setVal(0.);
              errorCoar.setVal(0.);

              getError(errorFine, extendedPrimFine, coveredFaceFine, faceDir, sit(), domainBoxFine,dxFine);
              getError(errorCoar, extendedPrimCoar, coveredFaceCoar, faceDir, sit(), domainBoxCoar,dxCoar);

              compareError(errorFine, errorCoar,
                           coveredFaceFine, coveredFaceCoar,
                           ebisBoxFine, ebisBoxCoar, faceDir, sit());

            }
        }
    }
}
/************/
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

    coveredExtrapTest();

#ifdef CH_MPI
    dumpmemoryatexit();
  }
  MPI_Finalize();
#endif
}

