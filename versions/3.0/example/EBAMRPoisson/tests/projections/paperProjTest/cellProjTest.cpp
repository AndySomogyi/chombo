#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <iostream>
#include <cmath>
using std::cerr;

#include "ParmParse.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"
#include "EBFluxFAB.H"
#include "EBFluxFactory.H"

#include "PoissonUtilities.H"
#include "EBAMRPoissonOp.H"

#include "EBFABView.H"
#include "EBDebugDump.H"
#include "DebugDump.H"
#include "EBSimpleSolver.H"
#include "EBLevelCCProjector.H"
#include "EBCompositeCCProjector.H"
#include "EBGradDivFilter.H"
#include "InflowOutflowIBC.H"
#include "cellProjTestF_F.H"
#include "EBAMRDataOps.H"


#define CHECK_FLOATING_PT 0

#if CHECK_FLOATING_PT==1
#include <fenv.h>
#endif
Real getVelExact(const VolIndex& a_vof, const Real& a_dx, const Real& a_domLen, int a_idir)
{
  Real inflowVel = 1.0;
  ParmParse pp;
  pp.get("inflow_vel", inflowVel);
  Real exactVal = 0;
  if(a_idir==0)
    {
      exactVal = inflowVel;
    }

  return exactVal;
}
/******/
void
getGradMag(Vector<LevelData<EBCellFAB>* >&       a_gmag,
           const Vector<LevelData<EBCellFAB>* >& a_gphi,
           const Vector< DisjointBoxLayout >&    a_grids,
           const Vector< EBISLayout >&           a_ebisl)
{
  for(int ilev = 0; ilev < a_grids.size(); ilev++)
    {
      LevelData<EBCellFAB>& gmag = *a_gmag[ilev];
      const LevelData<EBCellFAB>& gphi = *a_gphi[ilev];
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          BaseFab<Real>& gmagfab        = gmag[dit()].getSingleValuedFAB();
          const  BaseFab<Real>& gphifab = gphi[dit()].getSingleValuedFAB();
          FORT_GRADMAG(CHF_FRA1(gmagfab, 0),
                       CHF_CONST_FRA(gphifab),
                       CHF_BOX(a_grids[ilev][dit()]));
        }
    }
}

void setExactVeloc(LevelData<EBCellFAB>&                 a_veloc,
                   const DisjointBoxLayout&              a_grids,
                   const EBISLayout&                     a_ebisl,
                   const RealVect&                       a_dx,
                   const PoissonParameters&              a_params)
{
  ParmParse pp;
  RealVect domainLength = a_params.domainLength;
  for(DataIterator dit = a_grids.dataIterator(); dit.ok();++dit)
    {
      EBCellFAB& vel = a_veloc[dit()];
      vel.setVal(0.);
      IntVectSet ivsBox = IntVectSet(vel.getRegion());
      const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
      for(VoFIterator vofit(ivsBox, ebgraph); vofit.ok(); ++vofit)
        {
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              Real velexact = getVelExact(vofit(), a_dx[idir], domainLength[idir], idir);
              vel(vofit(), idir) = velexact;
            }
        }
    }
}
void outputGraphValues(Vector<Real>&              a_xval,
                       Vector<Real>&              a_yval,
                       char* a_dataFileName)
{
  char fileOut[1000];
  sprintf(fileOut, "%s", a_dataFileName);
  FILE *file1;
  if ((file1 = fopen(fileOut,"w")) == NULL)
    {
      MayDay::Error("problem opening the graph file");
    }

  for(int ivec = 0; ivec < a_xval.size(); ivec++)
    {
      fprintf(file1, "%e    %e \n", a_xval[ivec], a_yval[ivec]);
    }

  fclose(file1);
}
/******/
void runTest(const Vector< DisjointBoxLayout >&   a_grids,
             const Vector< EBISLayout >&          a_ebisl,
             const PoissonParameters&             a_params)
{
  int nlevels = a_params.numLevels;
  Vector<LevelData<EBCellFAB>* > velo(    nlevels, NULL);//new velo
  Vector<LevelData<EBCellFAB>* > gphi(    nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > divu(    nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > gmag(    nlevels, NULL);

  RealVect dxLevCoarsest = RealVect::Unit;
  dxLevCoarsest *=a_params.coarsestDx;
  ProblemDomain domLevCoarsest(a_params.coarsestDomain);

  RealVect dxLev = dxLevCoarsest;
  ProblemDomain domLev = domLevCoarsest;
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      pout() << "creating data for level " << ilev << endl;
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      divu[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], 1       ,  IntVect::Zero,   ebcellfact);
      gmag[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], 1       ,  IntVect::Zero,   ebcellfact);
      velo[ilev]    = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,  2*IntVect::Unit, ebcellfact);
      gphi[ilev]    = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,  2*IntVect::Unit, ebcellfact);

      //set phi = phiExact, rhs=lphiexact  This makes AMRResidual return lphiexact-Lphi
      EBLevelDataOps::setToZero(*(velo[ilev]));
      EBLevelDataOps::setToZero(*(divu[ilev]));
      EBLevelDataOps::setToZero(*(gphi[ilev]));

      setExactVeloc(*velo[ilev], a_grids[ilev], a_ebisl[ilev], dxLev, a_params);

      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);

    }

  ParmParse ppebamrieuler;

  int flowDir = 0;
  Real inflowVel = 1.0;
  ppebamrieuler.get("inflow_vel", inflowVel);
  bool doslip =true;
  int orderebbc = 2;
  InflowOutflowIBCFactory ibcfact(flowDir, inflowVel, orderebbc, doslip);

  EBIBC* ebibc = ibcfact.create();
  RefCountedPtr<BaseDomainBCFactory> macBCVel = ebibc->getMACVelBC();
  RefCountedPtr<BaseDomainBCFactory> celBCPhi = ebibc->getPressBC();
  RefCountedPtr<BaseEBBCFactory>     ebBCPhi  = ebibc->getPressureEBBC();


  Vector<EBLevelGrid>                      eblg   (a_grids.size());
  Vector<RefCountedPtr<EBQuadCFInterp> >   quadCFI(a_grids.size());
  domLev = domLevCoarsest;
  for(int ilev = 0; ilev < a_grids.size(); ilev++)
    {
      int nvar = 1;
      int nref = a_params.refRatio[ilev];
      eblg[ilev] = EBLevelGrid(a_grids[ilev], a_ebisl[ilev], domLev);
      if(ilev > 0)
        {
          int nrefOld = a_params.refRatio[ilev-1];
          ProblemDomain domLevCoar = coarsen(domLev, nrefOld);
          quadCFI[ilev] = RefCountedPtr<EBQuadCFInterp>(new EBQuadCFInterp(a_grids[ilev  ],
                                                                           a_grids[ilev-1],
                                                                           a_ebisl[ilev  ],
                                                                           a_ebisl[ilev-1],
                                                                           domLevCoar,
                                                                           nrefOld, nvar,
                                                                           *(eblg[ilev].getCFIVS())));

        }
      domLev.refine(nref);
    }


  EBCompositeCCProjector projectinator(eblg,  a_params.refRatio, quadCFI,
                                       a_params.coarsestDx,
                                       RealVect::Zero,
                                       macBCVel,
                                       celBCPhi,
                                       ebBCPhi,
                                       true, -1, 3 ,40,1.e99, 1,0);


  //get starting divergence
  Real divunorm[3];
  Real gmagnorm[3];
  projectinator.kappaDivergence(divu, velo);
  for(int inorm = 0; inorm < 3; inorm++)
    {
      divunorm[inorm] = EBArith::norm(divu, a_grids, a_ebisl,a_params.refRatio, 0, inorm, EBNormType::OverBoth);
      gmagnorm[inorm] = EBArith::norm(divu, a_grids, a_ebisl,a_params.refRatio, 0, inorm, EBNormType::OverBoth);
    }
  pout() << "div(vel) before cell project: " << "L_inf = " << divunorm[0]  << ", L_1 = " << divunorm[1] << ", L_2 = " << divunorm[2] << endl;
  pout() << "mag(grad) before cell project: " << "L_inf = " << gmagnorm[0]  << ", L_1 = " << gmagnorm[1] << ", L_2 = " << gmagnorm[2] << endl;

  ParmParse pp;
  int numProj = 10;
  pp.get("num_proj",numProj);
  int numIters = 10;
  pp.get("num_iters",numIters);

  Vector<Real> xval(numIters);
  Vector< Vector<Real> > divuNorms(3, Vector<Real>(numIters,0.0));
  Vector< Vector<Real> > gmagNorms(3, Vector<Real>(numIters,0.0));
  for(int iter = 0;iter<numIters;iter++)
    {
      xval[iter] = Real(iter);
      for(int iproj = 0;iproj<numProj;iproj++)
        {
          projectinator.project(velo, gphi);
          EBAMRDataOps::setCoveredVal(velo,0.0);
          EBAMRDataOps::setCoveredVal(gphi,0.0);

          {//compute: div(velo)
            projectinator.kappaDivergence(divu, velo);
            getGradMag(gmag, gphi, a_grids, a_ebisl);

            EBAMRDataOps::setCoveredVal(divu, 0.0);
            EBAMRDataOps::setCoveredVal(gmag, 0.0);
            for(int inorm = 0; inorm < 3; inorm++)
              {
                divunorm[inorm] = EBArith::norm(divu, a_grids, a_ebisl,a_params.refRatio, 0, inorm, EBNormType::OverBoth);
                gmagnorm[inorm] = EBArith::norm(gmag, a_grids, a_ebisl,a_params.refRatio, 0, inorm, EBNormType::OverBoth);

                divuNorms[inorm][iter] = divunorm[inorm];
                gmagNorms[inorm][iter] = gmagnorm[inorm];
              }
            pout() << setw(12)
                   << setprecision(6)
                   << setiosflags(ios::showpoint)
                   << setiosflags(ios::scientific) ;
            pout() << "div_proj: div(vel)  after cell project: " << iter  << " =L_inf = " << divunorm[0]  << ", L_1 = " << divunorm[1] << ", L_2 = " << divunorm[2] << endl;
            pout() << "mag_grad: mag(gphi) after cell project: " << iter  << " =L_inf = " << gmagnorm[0]  << ", L_1 = " << gmagnorm[1] << ", L_2 = " << gmagnorm[2] << endl;
          }
        }
    }

  for(int inorm = 0; inorm < 3; inorm++)
    {
      char divuName[100];
      char gmagName[100];

      if(inorm == 0)
        {
          sprintf(divuName,"L_inf(div(velo))");
          sprintf(gmagName,"L_inf(mag(grad))");
        }
      else
        {
          sprintf(divuName,"L_%d(div(velo))",inorm);
          sprintf(gmagName,"L_%d(mag(grad))",inorm);
        }
      outputGraphValues(xval,
                        divuNorms[inorm],
                        divuName);
      outputGraphValues(xval,
                        gmagNorms[inorm],
                        gmagName);
    }

  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete divu[ilev];
      delete velo[ilev];
      delete gphi[ilev];
      delete gmag[ilev];
    }

}
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // Scoping trick
  {
    if (argc < 2)
      {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }

    char* inFile = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,inFile);

    PoissonParameters param;
    Vector<DisjointBoxLayout> grids;
    Vector<EBISLayout>        ebisl;

    //read params from file
    getPoissonParameters(param);

    //define geometry from given params
    definePoissonGeometry(param);

    getAllIrregRefinedLayouts(grids, ebisl, param);

    runTest(grids, ebisl, param);
 }// End scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif
}
