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

#include "parstream.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "PoissonUtilities.H"
#include "functionsF_F.H"
#include "PoissProbF_F.H"
#include "AMRIO.H"
#include "BoxIterator.H"
#include "BiCGStabSolver.H"
#include "computeNorm.H"
#include "BCFunc.H"
#include "NodeBCFunc.H"
#include "AMRNodeOp.H"
#include "ResistivityOp.H"
#include "ViscousTensorOp.H"
std::vector<bool> GlobalBCRS::s_printedThatLo = std::vector<bool>(SpaceDim, false);
std::vector<bool> GlobalBCRS::s_printedThatHi = std::vector<bool>(SpaceDim, false);
std::vector<int>  GlobalBCRS::s_bcLo = std::vector<int>();
std::vector<int>  GlobalBCRS::s_bcHi = std::vector<int>();
RealVect          GlobalBCRS::s_trigvec = RealVect::Zero;
bool              GlobalBCRS::s_areBCsParsed= false;
bool              GlobalBCRS::s_valueParsed= false;
bool              GlobalBCRS::s_trigParsed= false;

/***************/
RealVect& getTrigRV()
{
  Real pi = 4.*atan(1.0);
  if(!GlobalBCRS::s_trigParsed)
    {
      GlobalBCRS::s_trigParsed = true;
      ParmParse pp;
      std::vector<Real> trigvec(SpaceDim);
      pp.getarr("trig",trigvec,0,SpaceDim);
      
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          GlobalBCRS::s_trigvec[idir] = pi*trigvec[idir];
        }
    }
  return GlobalBCRS::s_trigvec;
}
void TrigValueNeum(Real* pos,
                   int* dir,
                   Side::LoHiSide* side,
                   Real* a_values)
{
  RealVect& trig = getTrigRV();
  RealVect xval;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      xval[idir] = pos[idir];
    }
  RealVect gradPhi;
  FORT_GETGRADPHIPOINT(CHF_REALVECT(gradPhi),
                       CHF_CONST_REALVECT(trig),
                       CHF_CONST_REALVECT(xval));
  
  a_values[0] = gradPhi[*dir];
}

void ResistDiri(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  RealVect& trig = getTrigRV();
  RealVect xval;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      xval[idir] = pos[idir];
    }
  ParmParse pp;
  int whichMag;
  pp.get("which_mag", whichMag);
  for(int icomp = 0; icomp < SpaceDim; icomp++)
    {
      Real value;
      FORT_GETMAGPOINTRESIST(CHF_REAL(value),
                             CHF_CONST_REALVECT(trig),
                             CHF_CONST_REALVECT(xval),
                             CHF_INT(icomp),
                             CHF_INT(whichMag));
      a_values[icomp] = value;
    }
}


void TrigValueDiri(Real* pos,
                   int* dir,
                   Side::LoHiSide* side,
                   Real* a_values)
{
  RealVect& trig = getTrigRV();
  RealVect xval;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      xval[idir] = pos[idir];
    }
  Real value;
  FORT_GETPHIPOINT(CHF_REAL(value),
                   CHF_CONST_REALVECT(trig),
                   CHF_CONST_REALVECT(xval));
  a_values[0] = value;
}

void 
nodeOutputData(const Vector<LevelData<NodeFArrayBox>* >& vectPhi,
               const Vector<DisjointBoxLayout>& vectGrids,
               const ProblemDomain&     domain,
               const Vector<int>& vectRatio,
               Real dxCoarsest,
               int numlevels,
               string filename,
               string varname)
{
  pout() << "no output yet" << endl ;
}

void 
outputData(const Vector<LevelData<FArrayBox>* >& vectPhi,
           const Vector<DisjointBoxLayout>& vectGrids,
           const ProblemDomain&     domain,
           const Vector<int>& vectRatio,
           Real dxCoarsest,
           int numlevels,
           string filename,
           string varname)
{
#ifdef CH_USE_HDF5
  Vector<string> vectName(vectPhi[0]->nComp(), varname);
  Real time = 1;  //placeholder
  WriteAMRHierarchyHDF5(filename,                             
                        vectGrids, 
                        vectPhi, 
                        vectName,                 
                        domain.domainBox(),                               
                        dxCoarsest, time, time,                                         
                        vectRatio,                           
                        numlevels);
#endif
}

int setRHS(Vector<LevelData<FArrayBox>* >& vectRhs,
           PoissonParameters&              a_params,
           int  a_numLevels)
           
{
  int numlevels = a_numLevels;
  if(numlevels < 0)
    {
      numlevels = a_params.numLevels;
    }
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  Real rhono, rno;
  int iprob;
  ParmParse pp;
  pp.get("iprob", iprob);

  rhono = 0.75;
  rno = 0.5;

  Real dxlev = a_params.coarsestDx;
  ProblemDomain domlev = a_params.coarsestDomain;
  for(int ilev = 0; ilev < numlevels; ilev++)
    {
      if(iprob == -1)
        {
          //trig prob for convergence tests
          RealVect vectDx = dxlev*RealVect::Unit;
          setTrigKappaLOfPhi(*vectRhs[ilev], vectDx, a_params);
        }
      else
        {
          LevelData<FArrayBox>& rhsLD = *vectRhs[ilev];
          DataIterator dit =  rhsLD.dataIterator();
          for(dit.reset(); dit.ok(); ++dit)
            {
              FArrayBox& rhsFab = rhsLD[dit()];
              rhsFab.setVal(0.);
              bool useEBGrids;
              pp.get("use_eb_grids", useEBGrids);

              if(useEBGrids)
                {
                  rhsFab.setVal(1.0);
                }
              else
                {
                  FORT_GETRHSPOIS(CHF_FRA(rhsFab),
                                  CHF_BOX(rhsFab.box()),
                                  CHF_BOX(domlev.domainBox()),
                                  CHF_CONST_REAL(dxlev),
                                  CHF_CONST_REAL(rhono),
                                  CHF_CONST_REAL(rno),
                                  CHF_CONST_INT(iprob));          
                }
            }
        }
      dxlev /= a_params.refRatio[ilev];
      domlev.refine(a_params.refRatio[ilev]);
      
    }
  return 0;
}

/*
  tag cells for refinement based on magnitude(RHS)
*/
void
tagCells(Vector<LevelData<FArrayBox>* >& vectRHS,
         Vector<IntVectSet>& tagVect,
         Vector<Real>& vectDx,
         Vector<ProblemDomain>& vectDomain,
         const Real refine_thresh, 
         const int tags_grow,
         const int baseLevel, 
         int numLevels) 
{
  for (int lev=baseLevel; lev!= numLevels; lev++) 
    {
      IntVectSet local_tags;
      LevelData<FArrayBox> & levelRhs = *vectRHS[lev];
      DisjointBoxLayout level_domain = levelRhs.getBoxes();
      DataIterator dit = levelRhs.dataIterator();

      Real maxRHS = 0;

      maxRHS = norm(levelRhs, levelRhs.interval(), 0);

      Real tagVal = maxRHS * refine_thresh;
      
      // now loop through grids and tag cells where RHS > tagVal
      for (dit.reset(); dit.ok(); ++dit) 
        {
          const Box thisBox = level_domain.get(dit());
          const FArrayBox& thisRhs = levelRhs[dit()];
          BoxIterator bit(thisBox);
          for (bit.begin(); bit.ok(); ++bit) 
            {
              const IntVect& iv = bit();
              if (abs(thisRhs(iv)) >= tagVal)
                local_tags |= iv;
            }
        } // end loop over grids on this level
      
      local_tags.grow(tags_grow);
      const Box& domainBox = vectDomain[lev].domainBox();
      local_tags &= domainBox;

      tagVect[lev] = local_tags;

    } // end loop over levels
}  

/*
  Set grid hierarchy from input file
*/
void getDomainsAndDxes(  Vector<ProblemDomain>&     vectDomain,
                         Vector<Real>&              vectDx,
                         PoissonParameters&         a_params)
{

  vectDomain.resize(a_params.numLevels);
  vectDx.resize(    a_params.numLevels);
  vectDx[0] = a_params.coarsestDx;
  for(int ilev = 1; ilev < a_params.numLevels; ilev++)
    {
      vectDx[ilev] = vectDx[ilev-1]/a_params.refRatio[ilev-1];
    }


  vectDomain[0] = a_params.coarsestDomain;
  for(int ilev = 1;ilev < a_params.numLevels; ilev++)
    {
      vectDomain[ilev] = refine(vectDomain[ilev-1],a_params.refRatio[ilev-1]);
    }
}

int setGrids(Vector<DisjointBoxLayout>& vectGrids,
             PoissonParameters&         a_params)
{
  Vector<ProblemDomain>     vectDomain;
  Vector<Real>              vectDx;
  getDomainsAndDxes(vectDomain, vectDx, a_params);
  
  int numlevels = a_params.numLevels;

  ParmParse pp;
  bool useEBGrids;
  // grid generation parameters

  vectGrids.resize(numlevels);
  bool readInGrids = false;
  pp.query("use_eb_grids", useEBGrids);
  if(pp.contains("read_in_grids"))
    {
      pp.get("read_in_grids", readInGrids);
    }

  if(readInGrids)
    {

      ProblemDomain levDomain = a_params.coarsestDomain;
      for(int ilev = 0; ilev < a_params.numLevels; ilev++)
        {
          Vector<Box>   boxes;
          char boxCountVar[100];
          int boxCount;
          sprintf(boxCountVar, "level_%d_box_count", ilev);
          pp.get(boxCountVar, boxCount);
          boxes.resize(boxCount);
          for(int ibox = 0; ibox < boxCount; ibox++)
            {
              char boxLoVar[100];
              char boxHiVar[100];
              sprintf(boxLoVar, "level_%d_box_%d_lo", ilev, ibox);
              sprintf(boxHiVar, "level_%d_box_%d_hi", ilev, ibox);
              Vector<int> boxLo, boxHi;
              pp.getarr(boxLoVar, boxLo, 0, SpaceDim);
              pp.getarr(boxHiVar, boxHi, 0, SpaceDim);
              IntVect ivLo(D_DECL(boxLo[0], boxLo[1], boxLo[2]));
              IntVect ivHi(D_DECL(boxHi[0], boxHi[1], boxHi[2]));
              boxes[ibox] = Box(ivLo, ivHi);
              if(!levDomain.contains(boxes[ibox]))
                {
                  MayDay::Error("box outside of domain");
                }
            }
          //check to see if level 0 domain is covered
          if(ilev == 0)
            {
              IntVectSet ivDom(levDomain.domainBox());
              for(int ibox = 0; ibox < boxes.size(); ibox++)
                {
                  ivDom -= boxes[ibox];
                }
              if(!ivDom.isEmpty())
                {
                  MayDay::Error("level 0 boxes must cover the domain");
                }
            }
          Vector<int>  proc(a_params.numLevels);
          LoadBalance(proc,boxes);
          vectGrids[ilev] = DisjointBoxLayout(boxes, proc, levDomain);
          levDomain.refine(a_params.refRatio[ilev]);
        }

    }
  else if(useEBGrids)
    {
      pout() << "all regular geometry" << endl;
      pout() << "ignoring grid parameters and making simple  grids" << endl;
      Vector<int> proc(1, 0);
      Box coarBox = vectDomain[0].domainBox();
      Vector<Box> coarBoxes(1, coarBox);
      vectGrids[0] = DisjointBoxLayout(coarBoxes, proc, vectDomain[0]);

      for(int ilev = 1; ilev < numlevels; ilev++)
        {
          int iboxShrink = coarBox.size(0);
          iboxShrink /= 4;
          if(iboxShrink < 2)
            {
              MayDay::Error("wacky DBL generation technique failed, try making base box bigger");
            }
          coarBox.grow(-iboxShrink);
          coarBox.refine(a_params.refRatio[ilev-1]);
          Vector<Box> refBoxes(1, coarBox);
          vectGrids[ilev] = DisjointBoxLayout(refBoxes, proc, 
                                              vectDomain[ilev]);
        }
    }
  else
    {
      pout() << "tagging on gradient of RHS" << endl;
      int maxLevel = numlevels-1;
      Vector<Vector<Box> > newBoxes(numlevels);
      Vector<Vector<Box> > oldBoxes(numlevels);

      // determine grids dynamically, based on grad(RHS)
      // will need temp storage for RHS
      Vector<LevelData<FArrayBox>* > vectRHS(maxLevel+1,NULL);
      int ncomps = 1;

      // define base level first
      Vector< Vector<int> > procAssign(maxLevel+1);
      domainSplit(vectDomain[0], oldBoxes[0], a_params.maxGridSize, a_params.blockFactor);
      procAssign[0].resize(oldBoxes[0].size());
      LoadBalance(procAssign[0],oldBoxes[0]);
      
      vectGrids[0].define(oldBoxes[0],procAssign[0],vectDomain[0]);
      
      vectRHS[0] = new LevelData<FArrayBox>(vectGrids[0], ncomps, 
                                            IntVect::Zero);

      int topLevel = 0;

      bool moreLevels = (maxLevel > 0);
      
      int nesting_radius = 2;
      // create grid generation object
      BRMeshRefine meshrefine(vectDomain[0], a_params.refRatio, 
                              a_params.fillRatio,
                              a_params.blockFactor, nesting_radius, 
                              a_params.maxGridSize);

      while (moreLevels) 
        {
          // default is moreLevels = false
          // (only repeat loop in the case where a new level 
          // is generated which is still less than maxLevel)
          moreLevels = false;

          int baseLevel = 0;
          int oldTopLevel = topLevel;

          // now initialize RHS for this existing hierarchy
          setRHS(vectRHS, a_params, topLevel+1);

          Vector<IntVectSet> tagVect(topLevel+1);
          int tags_grow = 1;
          tagCells(vectRHS, tagVect, vectDx, vectDomain, 
                   a_params.refineThresh,
                   tags_grow, baseLevel, topLevel+1);

          int new_finest = meshrefine.regrid(newBoxes, tagVect, 
                                             baseLevel, 
                                             topLevel, oldBoxes);
                  
          if (new_finest > topLevel) 
            {
              topLevel++;
            }
                  
          oldBoxes = newBoxes;

          //  no need to do this for the base level (already done)
          for (int lev=1; lev<= topLevel; lev++) 
            {
              // do load balancing
              procAssign[lev].resize(newBoxes[lev].size());
              LoadBalance(procAssign[lev], newBoxes[lev]);
              const DisjointBoxLayout newDBL(newBoxes[lev], procAssign[lev],
                                             vectDomain[lev]);
              vectGrids[lev] = newDBL;
              delete vectRHS[lev];
              vectRHS[lev] = new LevelData<FArrayBox>(vectGrids[lev], ncomps,
                                                      IntVect::Zero);
            } // end loop over levels for initialization
          
          // figure out whether we need another pass through grid generation
          if ((topLevel<maxLevel) && (topLevel > oldTopLevel))
            moreLevels = true;

        } // end while moreLevels loop
      // clean up temp storage
      for (int ilev=0; ilev <vectRHS.size(); ilev++)
        {
          if (vectRHS[ilev] != NULL) 
            {
              delete vectRHS[ilev];
              vectRHS[ilev] = NULL;
            }
        }
    }
  return 0;
}




void ParseValue(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  ParmParse pp;
  Real bcVal;
  pp.get("bc_value",bcVal);
  a_values[0]=bcVal;
}


void ParseBC(FArrayBox& a_state,
             const Box& a_valid,
             const ProblemDomain& a_domain,
             Real a_dx,
             bool a_homogeneous)
{
  if(!a_domain.domainBox().contains(a_state.box()))
    {
 
      if(!GlobalBCRS::s_areBCsParsed)
        {
          ParmParse pp;
          pp.getarr("bc_lo", GlobalBCRS::s_bcLo, 0, SpaceDim);
          pp.getarr("bc_hi", GlobalBCRS::s_bcHi, 0, SpaceDim);
          GlobalBCRS::s_areBCsParsed = true;
        }

      Box valid = a_valid;
      for(int i=0; i<CH_SPACEDIM; ++i)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(i))
            {
              Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
              if(!a_domain.domainBox().contains(ghostBoxLo))
                {
                  if(GlobalBCRS::s_bcLo[i] == 1)
                    {
                      if(!GlobalBCRS::s_printedThatLo[i])
                        {
                          if(a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "const neum bcs lo for direction " << i << endl;
                        }
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Lo);
                    }
                  else if(GlobalBCRS::s_bcLo[i] == 2)
                    {
                      if(!GlobalBCRS::s_printedThatLo[i])
                        {
                          if(a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "trig neum bcs lo for direction " << i << endl;
                        }
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             TrigValueNeum,
                             i,
                             Side::Lo);
                    }
                  else if(GlobalBCRS::s_bcLo[i] == 0)
                    {
                      if(!GlobalBCRS::s_printedThatLo[i])
                        {
                          if(a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "const diri bcs lo for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Lo);
                      
                    }
                  else if(GlobalBCRS::s_bcLo[i] == 3)
                    {
                      if(!GlobalBCRS::s_printedThatLo[i])
                        {
                          if(a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "trig diri bcs lo for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             TrigValueDiri,
                             i,
                             Side::Lo);
                      
                    }
                  else if(GlobalBCRS::s_bcLo[i] == 4)
                    {
                      if(!GlobalBCRS::s_printedThatLo[i])
                        {
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "periodic bcs lo for direction " << i << endl;
                        }
                      
                    }
                  else if(GlobalBCRS::s_bcLo[i] == 5)
                    {
                      if(!GlobalBCRS::s_printedThatLo[i])
                        {
                          if(a_state.nComp() == 1)
                            {
                              MayDay::Error("using vector bc function for scalar");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "reflective slip bcs lo for direction " << i << endl;
                        }
                      ReflectiveVectorBC(a_state,
                                         valid,
                                         a_dx,
                                         i,
                                         Side::Lo);
                    }
                  else if(GlobalBCRS::s_bcLo[i] == 6)
                    {
                      if(!GlobalBCRS::s_printedThatLo[i])
                        {
                          if(a_state.nComp() == 1)
                            {
                              MayDay::Error("using vector bc function for scalar");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "no slip bcs lo for direction " << i << endl;
                        }
                      NoSlipVectorBC(a_state,
                                     valid,
                                     a_dx,
                                     i,
                                     Side::Lo, 1);
                    }
                  else if(GlobalBCRS::s_bcLo[i] == 7)
                    {
                      if(!GlobalBCRS::s_printedThatLo[i])
                        {
                          if(a_state.nComp() != SpaceDim)
                            {
                              MayDay::Error("bc function hardwired for ncomp == spacedim");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "Resistive diri bcs lo for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ResistDiri,
                             i,
                             Side::Lo);
                      
                    }
                  else
                    {
                      MayDay::Error("bogus bc flag lo");
                    }
                }
              
              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
                  if(GlobalBCRS::s_bcHi[i] == 1)
                    {
                      if(!GlobalBCRS::s_printedThatHi[i])
                        {
                          if(a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "const neum bcs hi for direction " << i << endl;
                        }
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Hi);
                    }
                  else if(GlobalBCRS::s_bcHi[i] == 2)
                    {
                      if(!GlobalBCRS::s_printedThatHi[i])
                        {
                          if(a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "trig neum bcs hi for direction " << i << endl;
                        }
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             TrigValueNeum,
                             i,
                             Side::Hi);
                    }
                  else if(GlobalBCRS::s_bcHi[i] == 0)
                    {
                      if(!GlobalBCRS::s_printedThatHi[i])
                        {
                          if(a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "const diri bcs hi for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Hi);
                    }
                  else if(GlobalBCRS::s_bcHi[i] == 3)
                    {
                      if(!GlobalBCRS::s_printedThatHi[i])
                        {
                          if(a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "trig diri bcs hi for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             TrigValueDiri,
                             i,
                             Side::Hi);
                    }
                  else if(GlobalBCRS::s_bcHi[i] == 4)
                    {
                      if(!GlobalBCRS::s_printedThatHi[i])
                        {
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "periodic bcs hi for direction " << i << endl;
                        }
                    }
                  else if(GlobalBCRS::s_bcHi[i] == 5)
                    {
                      if(!GlobalBCRS::s_printedThatHi[i])
                        {
                          if(a_state.nComp() == 1)
                            {
                              MayDay::Error("using vector bc function for scalar");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "reflective slip bcs hi for direction " << i << endl;
                        }
                      ReflectiveVectorBC(a_state,
                                         valid,
                                         a_dx,
                                         i,
                                         Side::Hi);
                    }
                  else if(GlobalBCRS::s_bcHi[i] == 6)
                    {
                      if(!GlobalBCRS::s_printedThatHi[i])
                        {
                          if(a_state.nComp() == 1)
                            {
                              MayDay::Error("using vector bc function for scalar");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "no slip bcs hi for direction " << i << endl;
                        }
                      NoSlipVectorBC(a_state,
                                     valid,
                                     a_dx,
                                     i,
                                     Side::Hi, 1);
                    }
                  else if(GlobalBCRS::s_bcHi[i] == 7)
                    {
                      if(!GlobalBCRS::s_printedThatHi[i])
                        {
                          if(a_state.nComp() != SpaceDim)
                            {
                              MayDay::Error("this bc hardwired for spacedim ncomp");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "resisitive mhd diri bcs hi for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ResistDiri,
                             i,
                             Side::Hi);
                    }
                  else
                    {
                      MayDay::Error("bogus bc flag hi");
                    }
                }
            } // end if not periodic in ith direction

        } 
    }
  
}

void NodeParseBC(NodeFArrayBox& a_state,
                 const Box& a_valid,
                 const ProblemDomain& a_domain,
                 Real a_dx,
                 bool a_homogeneous)
{
  if(!a_domain.domainBox().contains(a_state.box()))
    {
 
      if(!GlobalBCRS::s_areBCsParsed)
        {
          ParmParse pp;
          pp.getarr("bc_lo", GlobalBCRS::s_bcLo, 0, SpaceDim);
          pp.getarr("bc_hi", GlobalBCRS::s_bcHi, 0, SpaceDim);
          GlobalBCRS::s_areBCsParsed = true;
        }

      Box valid = a_valid;
      for(int i=0; i<CH_SPACEDIM; ++i)
        {
          if (!a_domain.isPeriodic(i))
            {
              Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
              if(!a_domain.domainBox().contains(ghostBoxLo))
                {
                  if(GlobalBCRS::s_bcLo[i] == 1)
                    {
                      if(!GlobalBCRS::s_printedThatLo[i])
                        {
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "const neum bcs lo for direction " << i << endl;
                        }
                      NodeNeumBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 ParseValue,
                                 i,
                                 Side::Lo);
                    }
                  else if(GlobalBCRS::s_bcLo[i] == 2)
                    {
                      if(!GlobalBCRS::s_printedThatLo[i])
                        {
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "trig neum bcs lo for direction " << i << endl;
                        }
                      NodeNeumBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 TrigValueNeum,
                                 i,
                                 Side::Lo);
                    }
                  else if(GlobalBCRS::s_bcLo[i] == 0)
                    {
                      if(!GlobalBCRS::s_printedThatLo[i])
                        {
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "const diri bcs lo for direction " << i << endl;
                        }
                      NodeDiriBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 ParseValue,
                                 i,
                                 Side::Lo);
                      
                    }
                  else if(GlobalBCRS::s_bcLo[i] == 3)
                    {
                      if(!GlobalBCRS::s_printedThatLo[i])
                        {
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "trig diri bcs lo for direction " << i << endl;
                        }
                      NodeDiriBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 TrigValueDiri,
                                 i,
                                 Side::Lo);
                      
                    }
                  else
                    {
                      MayDay::Error("bogus bc flag lo");
                    }
                }
              
              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
                  if(GlobalBCRS::s_bcHi[i] == 1)
                    {
                      if(!GlobalBCRS::s_printedThatHi[i])
                        {
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "const neum bcs hi for direction " << i << endl;
                        }
                      NodeNeumBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 ParseValue,
                                 i,
                                 Side::Hi);
                    }
                  if(GlobalBCRS::s_bcHi[i] == 2)
                    {
                      if(!GlobalBCRS::s_printedThatHi[i])
                        {
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "trig neum bcs hi for direction " << i << endl;
                        }
                      NodeNeumBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 TrigValueNeum,
                                 i,
                                 Side::Hi);
                    }
                  else if(GlobalBCRS::s_bcHi[i] == 0)
                    {
                      if(!GlobalBCRS::s_printedThatHi[i])
                        {
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "const diri bcs hi for direction " << i << endl;
                        }
                      NodeDiriBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 ParseValue,
                                 i,
                                 Side::Hi);
                    }
                  else if(GlobalBCRS::s_bcHi[i] == 3)
                    {
                      if(!GlobalBCRS::s_printedThatHi[i])
                        {
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "trig diri bcs hi for direction " << i << endl; // 
                        }
                      NodeDiriBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 TrigValueDiri,
                                 i,
                                 Side::Hi);
                    }
                  else
                    {
                      MayDay::Error("bogus bc flag hi");
                    }
                }
            } // end if periodic in ith direction
        }
      
    }
}

void
defineSolver(AMRMultiGrid<LevelData<FArrayBox> >&         a_solver,
             const Vector<DisjointBoxLayout>&             a_grids,
             LinearSolver<LevelData<FArrayBox> >&         a_bottomSolver,
             const PoissonParameters&                     a_params)
{
  ParmParse pp2;
  AMRPoissonOpFactory opFactory;
  opFactory.define(a_params.coarsestDomain,
                   a_grids,
                   a_params.refRatio,
                   a_params.coarsestDx, 
                   &ParseBC, a_params.alpha, a_params.beta);

  ProblemDomain coarsestDomain(a_params.coarsestDomain);
  a_solver.define(coarsestDomain, opFactory,  &a_bottomSolver, a_params.numLevels);

  int numSmooth, numMG, maxIter;
  Real eps, hang;
  pp2.get("num_smooth", numSmooth);
  pp2.get("num_mg",     numMG);
  pp2.get("max_iterations", maxIter);
  pp2.get("tolerance", eps);
  pp2.get("hang",      hang);
  Real normThresh = 1.0e-30;
  a_solver.setSolverParameters(numSmooth, numSmooth, numSmooth, 
                               numMG, maxIter, eps, hang, normThresh);
  a_solver.m_verbosity = a_params.verbosity;

}
/****/
/****/
void
defineResistivityCoef(Vector<RefCountedPtr<LevelData<FluxBox> > >&  a_eta,
                      const Vector<DisjointBoxLayout>&              a_grids,
                      const PoissonParameters&                      a_params)
{
  a_eta.resize(a_params.numLevels);
  RealVect dxLev = RealVect::Unit;
  dxLev *=a_params.coarsestDx;
  for(int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      a_eta[ilev] = RefCountedPtr<LevelData<FluxBox> >(
        new LevelData<FluxBox>(a_grids[ilev], 1, IntVect::Zero));
      setEtaResistive(*a_eta[ilev],   dxLev, a_params);
      dxLev /=      a_params.refRatio[ilev];
  }
}
void
defineViscousTensorCoef(Vector<RefCountedPtr<LevelData<FluxBox> > >&    a_eta,
                        Vector<RefCountedPtr<LevelData<FluxBox> > >&    a_lambda,
                        Vector<RefCountedPtr<LevelData<FArrayBox> > >&  a_beta,
                        const Vector<DisjointBoxLayout>&                a_grids,
                        const PoissonParameters&                        a_params)
{
  a_eta.resize(   a_params.numLevels);
  a_lambda.resize(a_params.numLevels);
  a_beta.resize(  a_params.numLevels);
  RealVect dxLev = RealVect::Unit;
  dxLev *=a_params.coarsestDx;
  for(int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      a_eta[ilev] = RefCountedPtr<LevelData<FluxBox> >(
        new LevelData<FluxBox>(a_grids[ilev], 1, IntVect::Zero));
      a_lambda[ilev] = RefCountedPtr<LevelData<FluxBox> >(
        new LevelData<FluxBox>(a_grids[ilev], 1, IntVect::Zero));
      a_beta[ilev] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(a_grids[ilev], 1, IntVect::Zero));

      setViscousCoefs(*a_eta[ilev], *a_lambda[ilev], *a_beta[ilev], dxLev, a_params);
      dxLev /=      a_params.refRatio[ilev];
  }
}
/****/
/****/
void
defineResistivitySolver(AMRMultiGrid<LevelData<FArrayBox> >&         a_solver,
                        const Vector<DisjointBoxLayout>&             a_grids,
                        LinearSolver<LevelData<FArrayBox> >&         a_bottomSolver,
                        const PoissonParameters&                     a_params)
{
  ParmParse pp2;
  Vector<RefCountedPtr<LevelData<FluxBox> > >  eta;
  defineResistivityCoef(eta, a_grids, a_params);
  
  ResistivityOpFactory opFactory(a_grids, eta, a_params.alpha, a_params.beta,
                                 a_params.refRatio, a_params.coarsestDomain,
                                 a_params.coarsestDx, &ParseBC);

  ProblemDomain coarsestDomain(a_params.coarsestDomain);
  a_solver.define(coarsestDomain, opFactory,  &a_bottomSolver, a_params.numLevels);

  int numSmooth, numMG, maxIter;
  Real eps, hang;
  pp2.get("num_smooth", numSmooth);
  pp2.get("num_mg",     numMG);
  pp2.get("max_iterations", maxIter);
  pp2.get("tolerance", eps);
  pp2.get("hang",      hang);
  Real normThresh = 1.0e-30;
  a_solver.setSolverParameters(numSmooth, numSmooth, numSmooth, 
                               numMG, maxIter, eps, hang, normThresh);
  a_solver.m_verbosity = 3;

}
void
defineViscousTensorSolver(AMRMultiGrid<LevelData<FArrayBox> >&         a_solver,
                          const Vector<DisjointBoxLayout>&             a_grids,
                          LinearSolver<LevelData<FArrayBox> >&         a_bottomSolver,
                          const PoissonParameters&                     a_params)
{
  ParmParse pp2;
  Vector<RefCountedPtr<LevelData<FluxBox> > >    eta;
  Vector<RefCountedPtr<LevelData<FluxBox> > >    lambda;
  Vector<RefCountedPtr<LevelData<FArrayBox> > >  beta;
  defineViscousTensorCoef(eta, lambda, beta, a_grids, a_params);
  
  ViscousTensorOpFactory opFactory(a_grids, eta, lambda, a_params.alpha, beta,
                                   a_params.refRatio, a_params.coarsestDomain,
                                   a_params.coarsestDx, &ParseBC);

  ProblemDomain coarsestDomain(a_params.coarsestDomain);
  a_solver.define(coarsestDomain, opFactory,  &a_bottomSolver, a_params.numLevels);

  int numSmooth, numMG, maxIter;
  Real eps, hang;
  pp2.get("num_smooth", numSmooth);
  pp2.get("num_mg",     numMG);
  pp2.get("max_iterations", maxIter);
  pp2.get("tolerance", eps);
  pp2.get("hang",      hang);
  Real normThresh = 1.0e-30;
  a_solver.setSolverParameters(numSmooth, numSmooth, numSmooth, 
                               numMG, maxIter, eps, hang, normThresh);
  a_solver.m_verbosity = 3;

}

void
nodeDefineSolver(AMRMultiGrid<LevelData<NodeFArrayBox> >&         a_solver,
                 const Vector<DisjointBoxLayout>&                 a_grids,
                 LinearSolver<LevelData<NodeFArrayBox> >&         a_bottomSolver,
                 const PoissonParameters&                         a_params)
{
  ParmParse pp2;
  AMRNodeOpFactory opFactory;
  opFactory.define(a_params.coarsestDomain,
                   a_grids,
                   a_params.refRatio,
                   a_params.coarsestDx, 
                   &NodeParseBC, a_params.alpha, a_params.beta);

  ProblemDomain coarsestDomain(a_params.coarsestDomain);
  a_solver.define(coarsestDomain, opFactory,  &a_bottomSolver, a_params.numLevels);

  int numSmooth, numMG, maxIter;
  Real eps, hang;
  pp2.get("num_smooth", numSmooth);
  pp2.get("num_mg",     numMG);
  pp2.get("max_iterations", maxIter);
  pp2.get("tolerance", eps);
  pp2.get("hang",      hang);
  Real normThresh = 1.0e-30;
  a_solver.setSolverParameters(numSmooth, numSmooth, numSmooth, 
                               numMG, maxIter, eps, hang, normThresh);
  a_solver.m_verbosity = 3;

}

void setVelViscous(LevelData<FArrayBox>&    a_mag,
                   const RealVect&          a_dx,
                   const PoissonParameters& a_params)
{
  int whichVel;
  ParmParse pp;
  pp.get("which_vel", whichVel);
  //for now just pass through.   might want to do something more than this tho
  setMagResistive(a_mag, a_dx, a_params, whichVel);
}

void setKLVViscous(LevelData<FArrayBox>&    a_klv,
                   const RealVect&          a_dx,
                   const PoissonParameters& a_params)
{
  int whichVel, whichLambda, whichEta,  whichBeta;
  Real eps;
  ParmParse pp;
  pp.get("which_vel", whichVel);
  pp.get("which_eta", whichEta);
  pp.get("which_lambda", whichLambda);
  pp.get("which_beta", whichBeta);
  pp.get("eta_eps",   eps);

  //which_lambda == 2 forces lambda = -factor*eta
  Real factor = 1;
  if(whichLambda == 2)
    {
      pp.get("lambda_factor", factor);
    }
  for (DataIterator dit = a_klv.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& curFAB = a_klv[dit()];
      Box curBox = curFAB.box();
      
      const RealVect&     trig = getTrigRV();

      for(int comp = 0; comp < a_klv.nComp(); comp++)
        {
          FORT_GETKLVVISCOUS(CHF_FRA1(curFAB,comp),
                             CHF_CONST_REALVECT(trig),
                             CHF_CONST_REALVECT(a_dx),
                             CHF_CONST_REALVECT(a_params.probLo),
                             CHF_CONST_REAL(a_params.alpha),
                             CHF_BOX(curBox),
                             CHF_INT(comp),
                             CHF_REAL(eps),
                             CHF_INT(whichVel),
                             CHF_INT(whichEta),
                             CHF_INT(whichLambda),
                             CHF_INT(whichBeta),
                             CHF_REAL(factor)
                             );
        }
    }
  
}

void setViscousCoefs(LevelData<FluxBox>&       a_eta,
                     LevelData<FluxBox>&       a_lambda,
                     LevelData<FArrayBox>&     a_beta,
                     const  RealVect&          a_dx,
                     const  PoissonParameters& a_params)
{
  int whichVel, whichLambda, whichEta, whichBeta;
  Real eps;
  ParmParse pp;
  pp.get("which_vel", whichVel);
  pp.get("which_eta", whichEta);
  pp.get("which_lambda", whichLambda);
  pp.get("which_beta", whichBeta);
  pp.get("eta_eps",   eps);
  
  setEtaResistive(a_eta,    a_dx, a_params, whichEta);
  if(whichLambda != 2)
    {
      setEtaResistive(a_lambda, a_dx, a_params, whichLambda);
    }
  else
    {
      //whichlambda == 2 forces lambda = -(2/3)*eta
      Real factor;
      pp.get("lambda_factor", factor);
      for (DataIterator dit = a_eta.dataIterator(); dit.ok(); ++dit)
        {
          a_lambda[dit()].copy(a_eta[dit()]);
          a_lambda[dit()] *= -factor;
        }
    }
  
  for (DataIterator dit = a_beta.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& curFAB = a_beta[dit()];
      Box curBox = curFAB.box();
      
      const RealVect&     trig = getTrigRV();

      FORT_GETBETAVISCOUS(CHF_FRA1(curFAB,0),
                          CHF_CONST_REALVECT(trig),
                          CHF_CONST_REALVECT(a_dx),
                          CHF_CONST_REALVECT(a_params.probLo),
                          CHF_REAL(eps),
                          CHF_BOX(curBox),
                          CHF_INT(whichBeta));
    }
}
/********/
void setMagResistive(LevelData<FArrayBox>&    a_mag,
                     const RealVect&          a_dx,
                     const PoissonParameters& a_params,
                     int a_whichMag)
{
  CH_assert(a_mag.nComp() ==SpaceDim);

  int  whichMag;
  if(a_whichMag < 0)
    {
      ParmParse pp;
      pp.get("which_mag", whichMag);
    }
  else
    {
      whichMag = a_whichMag;
    }
  for (DataIterator dit = a_mag.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& curFAB = a_mag[dit()];
      Box curBox = curFAB.box();
      
      const RealVect&     trig = getTrigRV();

      for(int comp = 0; comp < a_mag.nComp(); comp++)
        {
          FORT_GETMAGRESIST(CHF_FRA1(curFAB,comp),
                            CHF_CONST_REALVECT(trig),
                            CHF_CONST_REALVECT(a_dx),
                            CHF_CONST_REALVECT(a_params.probLo),
                            CHF_BOX(curBox),
                            CHF_INT(comp),
                            CHF_INT(whichMag));
        }
    }
}

/********/
void setKLBResistive(LevelData<FArrayBox>&    a_klb,
                     const RealVect&          a_dx,
                     const PoissonParameters& a_params)
{
  CH_assert(a_klb.nComp() == SpaceDim);

  int whichEta, whichMag;
  Real eps;
  ParmParse pp;
  pp.get("which_eta", whichEta);
  pp.get("which_mag", whichMag);
  pp.get("eta_eps",   eps);

  for (DataIterator dit = a_klb.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& curFAB = a_klb[dit()];
      Box curBox = curFAB.box();
      
      const RealVect&     trig = getTrigRV();

      for(int comp = 0; comp < a_klb.nComp(); comp++)
        {
          FORT_GETKLBRESIST(CHF_FRA1(curFAB,comp),
                            CHF_CONST_REALVECT(trig),
                            CHF_CONST_REALVECT(a_dx),
                            CHF_CONST_REALVECT(a_params.probLo),
                            CHF_CONST_REAL(a_params.alpha),
                            CHF_CONST_REAL(a_params.beta),
                            CHF_BOX(curBox),
                            CHF_INT(comp),
                            CHF_REAL(eps),
                            CHF_INT(whichMag),
                            CHF_INT(whichEta));
        }
    }
}
/********/
void setEtaResistive(LevelData<FluxBox>&      a_eta,
                     const RealVect&          a_dx,
                     const PoissonParameters& a_params,
                     int a_whichEta)
{
  CH_assert(a_eta.nComp() == 1);

  int whichEta;
  ParmParse pp;
  if(a_whichEta < 0)
    {
      pp.get("which_eta", whichEta);
    }
  else
    {
      whichEta = a_whichEta;
    }
  Real eps;
  pp.get("eta_eps", eps);

  for (DataIterator dit = a_eta.dataIterator(); dit.ok(); ++dit)
    {
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          FArrayBox& curFAB = a_eta[dit()][idir];
      
          if(whichEta == 0)
            {
              curFAB.setVal(1.0);
            }
          else
            {

              Box curBox = curFAB.box();
              const RealVect&     trig = getTrigRV();

              FORT_GETETARESIST(CHF_FRA1(curFAB, 0),
                                CHF_CONST_REALVECT(trig),
                                CHF_CONST_REALVECT(a_dx),
                                CHF_CONST_REALVECT(a_params.probLo),
                                CHF_BOX(curBox),
                                CHF_INT(idir),
                                CHF_REAL(eps),
                                CHF_INT(whichEta));
            }
        }
    }
}

/********/
void setTrigPhi(LevelData<FArrayBox>&    a_phi,
                const RealVect&          a_dx,
                const PoissonParameters& a_params)
{
  CH_assert(a_phi.nComp() == 1);

  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& curPhiFAB = a_phi[dit()];
      Box curPhiBox = curPhiFAB.box();
      
      const RealVect&     trig = getTrigRV();

      for(int comp = 0; comp < a_phi.nComp(); comp++)
        {
          FORT_GETPHI(CHF_FRA1(curPhiFAB,comp),
                      CHF_CONST_REALVECT(trig),
                      CHF_CONST_REALVECT(a_dx),
                      CHF_CONST_REALVECT(a_params.probLo),
                      CHF_CONST_REALVECT(a_params.probHi),
                      CHF_BOX(curPhiBox));
        }
    }
}

/********/
void setTrigKappaLOfPhi(LevelData<FArrayBox>&    a_kappaLOfPhi,
                        const RealVect&          a_dx,
                        const PoissonParameters& a_params)
{
  CH_assert(a_kappaLOfPhi.nComp() == 1);

  const RealVect&  trig = getTrigRV();

  for (DataIterator dit = a_kappaLOfPhi.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& curKappaLOfPhiFAB = a_kappaLOfPhi[dit()];
      Box curKappaLOfPhiBox = curKappaLOfPhiFAB.box();
      for(int comp = 0; comp < a_kappaLOfPhi.nComp(); comp++)
        {
          FORT_GETLOFPHI(CHF_FRA1(curKappaLOfPhiFAB,comp),
                         CHF_CONST_REALVECT(trig),
                         CHF_CONST_REALVECT(a_dx),
                         CHF_CONST_REALVECT(a_params.probLo),
                         CHF_CONST_REALVECT(a_params.probHi),
                         CHF_CONST_REAL(a_params.alpha),
                         CHF_CONST_REAL(a_params.beta),
                         CHF_BOX(curKappaLOfPhiBox));
        }

    }
}
/******/
int iscript(int icomp, int inorm, int ncomp)
{
  return icomp + inorm*ncomp;
}
/******/
void compareError(const Vector< LevelData<FArrayBox>* >&   a_errorFine,
                  const Vector< LevelData<FArrayBox>* >&   a_errorCoar,
                  const Vector< DisjointBoxLayout >&       a_gridsFine, 
                  const Vector< DisjointBoxLayout >&       a_gridsCoar, 
                  const PoissonParameters&                 a_paramsFine,
                  const PoissonParameters&                 a_paramsCoar,
                  const string& a_testName)
{
  const Vector<int> refRat = a_paramsCoar.refRatio;
  const int ncomp = a_errorFine[0]->nComp();
  const int nnorm = 3;
  Real* normsCoar = new Real[ncomp*nnorm];
  Real* normsFine = new Real[ncomp*nnorm];
  Real* orders    = new Real[ncomp*nnorm];
  for(int icomp = 0; icomp < ncomp; icomp++)
    {
      orders[icomp] = 0.0;
      for(int inorm = 0; inorm < nnorm; inorm++)
        {
          normsCoar[iscript(icomp, inorm, ncomp)] = 0;
          normsFine[iscript(icomp, inorm, ncomp)] = 0;
        }
    }
  ParmParse pp;
  pout() << "==============================================" << endl;
  for(int comp = 0; comp < ncomp; comp++)
    {
      pout() << "Comparing error in variable  " << comp << endl;
      pout() << "==============================================" << endl;
      for(int inorm = 0; inorm <= 2; inorm++)
        {

          if(inorm == 0)
            {
              pout() << endl << "Using max norm." << endl;
            }
          else
            {
              pout() << endl << "Using L-" << inorm << "norm." << endl;
            }
          Real dxCoar = a_paramsCoar.coarsestDx;
          Real dxFine = a_paramsFine.coarsestDx;
          Interval comps(comp,comp);
          int lbase = 0;
          Real coarnorm = computeNorm(a_errorCoar, refRat, dxCoar, comps, inorm, lbase);
                                      
          Real finenorm = computeNorm(a_errorFine, refRat, dxFine, comps, inorm, lbase);
                                      
          pout() << "Coarse Error Norm = " << coarnorm << endl;
          pout() << "Fine   Error Norm = " << finenorm << endl;

          normsCoar[iscript(comp,inorm,ncomp)] = coarnorm;
          normsFine[iscript(comp,inorm,ncomp)] = finenorm;

          if((Abs(finenorm) > 1.0e-10) && (Abs(coarnorm) > 1.0e-10)) 
            {
              Real order = log(Abs(coarnorm/finenorm))/log(2.0);
              //pout() << "Order of scheme = " << order << endl;
              orders[iscript(comp,inorm,ncomp)] = order;
            }
        }
      pout() << "==============================================" << endl ;;
    }
  

  //output in latex format to be safe
  int nfine = a_paramsFine.coarsestDomain.size(0);
  pout() << setw(12)
         << setprecision(6)
         << setiosflags(ios::showpoint)
         << setiosflags(ios::scientific) ;

  for (int inorm = 0; inorm <= 2; inorm++)
    {
      pout() << "\\begin{table}[p]" << endl;
      pout() << "\\begin{center}" << endl;
      pout() << "\\begin{tabular}{|c|c|c|c|} \\hline" << endl;
      pout() << "Variable & Coarse Error & Fine Error & Order\\\\" << endl;;
      pout() << "\\hline \\hline " << endl;
      for(int icomp = 0; icomp < ncomp; icomp++)
        {
          int iindex = iscript(icomp,inorm,ncomp);
          pout() << "var" << icomp << " &    \t "
                 << setw(12)
                 << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific) 
                 << normsCoar[iindex]  << " & "
                 << setw(12)
                 << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific) 
                 << normsFine[iindex] << " & "
                 << setw(12)
                 << setprecision(2)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific) 
                 << orders[iindex];
          pout() << " \\\\ " << endl <<   "\\hline"  <<  endl;
        }
      pout() << "\\end{tabular}" << endl;
      pout() << "\\end{center}" << endl;
      pout() << "\\caption{";
      pout() << a_testName ;
      pout() << " convergence rates using L-" << inorm << " norm. " << endl;
      pout() << "$h_f = \\frac{1}{" << nfine << "}$ and $h_c = 2 h_f$, $D = " << SpaceDim << "$ }"  << endl;
      pout() << "\\end{table}" << endl;
      pout() << endl << endl;
    }
  pout() << "data for paper:" << endl;
  int nx = a_paramsCoar.coarsestDomain.domainBox().size(0);
  pout() << nx << " ";
  for (int inorm = 0; inorm <= 2; inorm++)
    {
      int icomp = 0;
      int iindex = iscript(icomp,inorm,ncomp);
      pout()     << setw(12)
                 << setprecision(6)
                 << normsCoar[iindex]  << "  ";
    }
  pout () << endl;
  delete[] normsCoar;
  delete[] normsFine;
  delete[] orders   ;
}
/********/
void getCoarseLayoutsFromFine(Vector<DisjointBoxLayout>&       a_gridsCoar,
                              const Vector<DisjointBoxLayout>& a_gridsFine,
                              const PoissonParameters&         a_paramsCoar)
{
  int nlevels = a_paramsCoar.numLevels;
  a_gridsCoar.resize(nlevels);
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      CH_assert(a_gridsFine[ilev].coarsenable(2));
      coarsen(a_gridsCoar[ilev], a_gridsFine[ilev], 2);
    }
  
}
/********/
void PoissonParameters::coarsen(int a_factor)
{
  coarsestDx *= a_factor;
  coarsestDomain.coarsen(a_factor);
}
/********/
void PoissonParameters::refine(int a_factor)
{
  coarsestDx /= a_factor;
  coarsestDomain.refine(a_factor);
}
/********/
void getPoissonParameters(PoissonParameters&  a_params)
{
  ParmParse pp;

  std::vector<int> nCellsArray(SpaceDim);
  pp.getarr("n_cells",nCellsArray,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.nCells[idir] = nCellsArray[idir];
    }

  Vector<int> is_periodic(SpaceDim, false);
  pp.queryarr("periodic", is_periodic, 0, SpaceDim);

  pp.get("refine_threshold",a_params.refineThresh);
  pp.get("block_factor",a_params.blockFactor);
  pp.get("fill_ratio",a_params.fillRatio);
  pp.get("buffer_size",a_params.bufferSize);
  pp.get("alpha",a_params.alpha);
  pp.get("beta", a_params.beta);

  pp.get("verbose", a_params.verbosity);

  pp.get("max_level", a_params.maxLevel);
  a_params.numLevels = a_params.maxLevel + 1;
  pp.getarr("ref_ratio",a_params.refRatio,0,a_params.numLevels);

  IntVect lo = IntVect::Zero;
  IntVect hi = a_params.nCells;
  hi -= IntVect::Unit;

  Box crseDomBox(lo,hi);
  ProblemDomain crseDom(crseDomBox);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      crseDom.setPeriodic(dir, is_periodic[dir]);
    }
  a_params.coarsestDomain = crseDom;

  std::vector<Real> dLArray(SpaceDim);
  pp.getarr("domain_length",dLArray,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.domainLength[idir] = dLArray[idir];
    }

  pp.get("max_grid_size",a_params.maxGridSize);

  //derived stuff
  a_params.coarsestDx = a_params.domainLength[0]/a_params.nCells[0];

  a_params.probLo = RealVect::Zero;
  a_params.probHi = RealVect::Zero;
  a_params.probHi += a_params.domainLength;


}

