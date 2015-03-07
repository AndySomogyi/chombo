#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EBIndexSpace.H"
#include "parstream.H"
#include "LoadBalance.H"
#include "EBGraphFactory.H"
#include "EBDataFactory.H"
#include "BoxIterator.H"
#include "EBISLayout.H"
#include "LayoutIterator.H"
#include "BRMeshRefine.H"
#include "VoFIterator.H"
#include "IrregNode.H"
#include "memtrack.H"
#include "AllRegularService.H"
#include "AMRIO.H"
#include "PolyGeom.H"
#include "CH_Attach.H"
#include "EBLevelDataOps.H"
#include "NamespaceHeader.H"

Real EBIndexSpace::s_tolerance = 1.0e-12;
Real    EBISLevel::s_tolerance = 1.0e-12;
bool EBIndexSpace::s_verbose   = false;
bool    EBISLevel::s_verbose   = false;
bool EBIndexSpace::s_MFSingleBox=false;


int countIrreg(const GeometryService& a_geoserver, const Box& a_region,
               const ProblemDomain&   a_domain,
               const RealVect&        a_origin,
               const Real&            a_dx)
{
  int count = 0;
  int longdir;
  int length = a_region.longside(longdir);
  if(length <= 2)
    {
      return 1;
    }
  else
    {
      GeometryService::InOut  inout = a_geoserver.InsideOutside(a_region, a_domain, a_origin, a_dx);
      if(inout == GeometryService::Irregular)
        {
          int n = length/2;
          //CH_assert(n*2==length);
          Box low(a_region), high;
          high = low.chop(longdir, a_region.smallEnd(longdir)+n);
          count += countIrreg(a_geoserver, low,  a_domain, a_origin, a_dx);
          count += countIrreg(a_geoserver, high, a_domain, a_origin, a_dx);
        }
      else
        {
          return 0;
        }

    }
  return count;
}



/******************/
int
EBIndexSpace::numVoFs(const ProblemDomain& a_domain) const
{
  CH_TIME("EBIndexSpace::numVoFs");
  int ilev = getLevel(a_domain);
  const  EBISLevel&  ebisLevel = *m_ebisLevel[ilev];
  int numPtsLocal = ebisLevel.numVoFsOnProc();
  int  numPtsTot = EBLevelDataOps::parallelSum(numPtsLocal);
  return numPtsTot;
}
/******************/
int
EBISLevel::numVoFsOnProc() const
{
  CH_TIME("EBISLevel::numVoFsOnProc");
  int retval = 0;
  for(DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph& ebgraph = m_graph[dit()];
      int numVoFsBox = ebgraph.numVoFs(m_grids.get(dit()));
      retval += numVoFsBox;
    }

  return retval;
}
/******************/
void
EBISLevel::makeBoxes(Vector<Box>& a_boxes,
                     Vector<long>& a_loads,
                     const Box& a_region,
                     const ProblemDomain& a_domain,
                     const GeometryService& a_geoserver,
                     const RealVect&        a_origin,
                     const Real&            a_dx,
                     const int              a_ncellmax,
                     const EBIndexSpace* const ebisPtr)
{
  if(EBIndexSpace::s_MFSingleBox)
    {
      a_boxes.resize(1);
      a_loads.resize(1);
      a_boxes[0] = a_region;
      a_loads[0] = 1;
      return;
    }
  std::list<Box> boxes;
  makeBoxes(boxes, a_region, a_domain, a_geoserver, a_origin, a_dx, a_ncellmax);
  a_boxes.resize(boxes.size());
  a_loads.resize(boxes.size());
  std::list<Box>::iterator it = boxes.begin();
  for(int i=0; i<a_boxes.size(); ++i, ++it)
    {
      a_boxes[i]=*it;
    }
  mortonOrdering(a_boxes);
  for(int i=0; i<a_boxes.size(); i++)
    {
      if(a_boxes[i].ixType() ==  IndexType::TheNodeType())
        {
          a_boxes[i].convert(IndexType::TheCellType());
          a_loads[i] = 8;
        }
      else
        {
          a_loads[i] = 1;
        }
    }
}
void
EBISLevel::makeBoxes(std::list<Box>& a_boxes,
                     const Box& a_region,
                     const ProblemDomain& a_domain,
                     const GeometryService& a_geoserver,
                     const RealVect&        a_origin,
                     const Real&            a_dx,
                     const int              a_ncellmax)
{

  int longdir;
  int length = a_region.longside(longdir);

  if(length > a_ncellmax)
    {
      int n = length/2;
      //CH_assert(n*2==length);
      Box low(a_region), high;
      high = low.chop(longdir, a_region.smallEnd(longdir)+n);
      makeBoxes(a_boxes, low,  a_domain,
                a_geoserver, a_origin, a_dx, a_ncellmax);
      makeBoxes(a_boxes, high, a_domain,
                a_geoserver, a_origin, a_dx, a_ncellmax);
    }
  else
    {
      if(a_geoserver.InsideOutside(a_region, a_domain, a_origin, a_dx) == GeometryService::Irregular  )
        {
          if(length < 16)
            {
              Box n = a_region;
              n.convert(IndexType::TheNodeType());
              a_boxes.push_back(n);
            }
          else
            {
              int n = length/2;
              //CH_assert(n*2==length);
              Box low(a_region), high;
              high = low.chop(longdir, a_region.smallEnd(longdir)+n);
              makeBoxes(a_boxes, low,  a_domain,
                        a_geoserver, a_origin, a_dx, a_ncellmax);
              makeBoxes(a_boxes, high, a_domain,
                        a_geoserver, a_origin, a_dx, a_ncellmax);
            }
        }
      else
        {
          a_boxes.push_back(a_region);
        }
    }
}

/******************/
#ifdef CH_USE_HDF5
EBISLevel::EBISLevel(HDF5Handle& a_handle)
{
  CH_TIME("EBISLevel::EBISLevel_hdf5");
  m_cacheMisses=0; m_cacheHits=0; m_cacheStale=0;
  HDF5HeaderData header;
  header.readFromFile(a_handle);
  m_origin = header.m_realvect["EBIS_origin"];
  m_domain = header.m_box     ["EBIS_domain"];
  m_dx =     header.m_real    ["EBIS_dx"]    ;
  m_tolerance = m_dx*1E-4;

  //read in the grids
  Vector<Box> boxes;
  read(a_handle,boxes);
  Vector<int> procAssign;
  LoadBalance(procAssign, boxes);
  m_grids.define(boxes, procAssign);//this should use m_domain for periodic...
  EBGraphFactory graphfact(m_domain);
  m_graph.define(m_grids, 1, IntVect::Zero, graphfact);

  //read the graph  in from the file
  std::string graphName("EBIS_graph");
  int eekflag = read(a_handle, m_graph, graphName, m_grids, Interval(), false);
  if(eekflag != 0)
    {
      MayDay::Error("error in writing graph");
    }

  //need a ghosted layout so that the data can be defined properly
  LevelData<EBGraph> ghostGraph(m_grids, 1, IntVect::Unit, graphfact);
  Interval interv(0,0);
  m_graph.copyTo(interv, ghostGraph, interv);

  //now the data for the graph
  EBDataFactory dataFact;
  m_data.define(m_grids, 1, IntVect::Zero, dataFact);
  for(DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      m_data[dit()].defineVoFData( ghostGraph[dit()], m_grids.get(dit()));
      m_data[dit()].defineFaceData(ghostGraph[dit()], m_grids.get(dit()));
    }
  //read the data  in from the file
  std::string  dataName("EBIS_data");
  eekflag = read(a_handle, m_data ,  dataName, m_grids, Interval(), false);
  if(eekflag != 0)
    {
      MayDay::Error("error in writing data");
    }
}
/******************/
void
EBISLevel::write(HDF5Handle& a_handle)
{
  CH_TIME("EBISLevel::write");
  HDF5HeaderData header;
  //this naming stuff kinda depends on the fact
  //that we are only outputting the finest level.
  //we could be slick and incorporate the
  //level number in there if we wanted.
  header.m_int["num_levels"] = 1;
  header.m_int["num_components"] = 1;
  header.m_string["component_0"] = "phi0";
  header.m_realvect["EBIS_origin"] = m_origin;
  header.m_box     ["EBIS_domain"] = m_domain.domainBox();
  header.m_real    ["EBIS_dx"]     = m_dx;
  header.writeToFile(a_handle);
  //write the grids to the file
  CH_XD::write(a_handle, m_grids);

  std::string graphName("EBIS_graph");
  std::string  dataName("EBIS_data");
  int eekflag = CH_XD::write(a_handle, m_graph, graphName);
  if(eekflag != 0)
    {
      MayDay::Error("error in writing graph");
    }
  eekflag = CH_XD::write(a_handle, m_data ,  dataName);
  if(eekflag != 0)
    {
      MayDay::Error("error in writing data");
    }
}
#endif
/******************/
EBISLevel::EBISLevel(const ProblemDomain& a_domain,
                     const RealVect& a_origin,
                     const Real& a_dx,
                     const GeometryService& a_geoserver,
                     const EBIndexSpace* const ebisPtr)
{

  CH_TIME("EBISLevel::EBISLevel_geoserver_domain");
  m_cacheMisses=0; m_cacheHits=0; m_cacheStale=0;
  m_domain = a_domain;
  m_dx = a_dx;
  m_tolerance = a_dx*1E-4;
  m_origin = a_origin;
  //divide up the domain into a layout
  Vector<Box> vbox;
  Vector<long> irregCount;
  {
    CH_TIME("EBISLevel::EBISLevel_makeboxes");
    makeBoxes(vbox, irregCount, a_domain.domainBox(), a_domain, a_geoserver,
              a_origin, a_dx, ebisPtr->getNCellMax(), ebisPtr);

  }

  // pout()<<vbox<<"\n\n";
  //load balance the boxes
  Vector<int> procAssign;
  LoadBalance(procAssign, irregCount, vbox);
  //   pout()<<irregCount<<std::endl;
  //   pout()<<procAssign<<std::endl;
  m_grids.define(vbox, procAssign);//this should use a_domain for periodic

  {
    //both domain and box set by the factory when leveldata construction is used
    EBGraphFactory graphfact(m_domain);
    m_graph.define(m_grids, 1, IntVect::Unit, graphfact);
    LayoutData<Vector<IrregNode> > allNodes(m_grids);

    //define the graph stuff
    for(DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
      {
        Box region = m_grids.get(dit());
        region.grow(1);
        Box ghostRegion = grow(region,1);
        ghostRegion &= m_domain;
        region &= m_domain;

        EBGraph& ebgraph = m_graph[dit()];
        GeometryService::InOut inout = a_geoserver.InsideOutside(region, m_domain, m_origin, m_dx);
        if(inout == GeometryService::Regular)
          {
            ebgraph.setToAllRegular();
          }
        else if(inout == GeometryService::Covered)
          {
            ebgraph.setToAllCovered();
          }
        else
          {
            BaseFab<int>       regIrregCovered(ghostRegion, 1);
            Vector<IrregNode>&  nodes = allNodes[dit()];
 
            a_geoserver.fillGraph(regIrregCovered, nodes, region,
                                  ghostRegion, m_domain,
                                  m_origin, m_dx);
            //pout()<<nodes<<"\n";
#ifndef NDEBUG
            EBGraphImplem::checkGraph(regIrregCovered, nodes, region, m_domain);
#endif
            ebgraph.buildGraph(regIrregCovered, nodes, region, m_domain);

          }
      }

    //overallMemoryUsage();
    //need a ghosted graph to make data so we can define baseiffabs
    //LevelData<EBGraph> ghostGraph(m_grids, 1, IntVect::Unit, graphfact);
    //overallMemoryUsage();
    //Interval interv(0,0);
    //m_graph.copyTo(interv, ghostGraph, interv);
    //overallMemoryUsage();
    //now the data for the graph
    EBDataFactory dataFact;
    m_data.define(m_grids, 1, IntVect::Zero, dataFact);
    for(DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
      {
        //m_data[dit()].define(ghostGraph[dit()], allNodes[dit()], m_grids.get(dit()));
        m_data[dit()].define(m_graph[dit()], allNodes[dit()], m_grids.get(dit()));

      }
    //overallMemoryUsage();
  }


  //
  if(a_geoserver.canGenerateMultiCells())
    {
      fixRegularNextToMultiValued();
    }



}
/************/
//now fix the multivalued next to regular thing for the graph and the data
//the oldgraph/newgraph thing is necessary because the graphs are
//reference counted and they have to be kept consistent with the data
/***********/
void EBISLevel::
fixRegularNextToMultiValued()
{
  CH_TIME("EBISLevel::fixRegularNextToMultiValued");
  EBGraphFactory graphfact(m_domain);
  LayoutData<IntVectSet> vofsToChange(m_grids);
  LevelData<EBGraph> oldGhostGraph(m_grids, 1, 2*IntVect::Unit, graphfact);
  Interval interv(0,0);
  m_graph.copyTo(interv, oldGhostGraph, interv);
  for(DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      m_graph[dit()].getRegNextToMultiValued(vofsToChange[dit()],
                                             oldGhostGraph[dit()]);

      m_graph[dit()].addFullIrregularVoFs(vofsToChange[ dit()],
                                          oldGhostGraph[dit()]);
    }

  EBDataFactory datafact;
  LevelData<EBGraph> newGhostGraph(m_grids, 1, 2*IntVect::Unit, graphfact);
  LevelData<EBData>  newGhostData(m_grids, 1, IntVect::Unit, datafact);
  m_graph.copyTo(interv, newGhostGraph, interv);
  for(DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      Box localBox = m_grids.get(dit());
      localBox.grow(1);
      localBox &= m_domain;
      newGhostData[dit()].defineVoFData(oldGhostGraph[dit()],  localBox);
      newGhostData[dit()].defineFaceData(oldGhostGraph[dit()], localBox);
    }
  m_data.copyTo(interv,  newGhostData,  interv);
  for(DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      m_data[dit() ].addFullIrregularVoFs(vofsToChange[dit()],
                                          newGhostGraph[dit()],
                                          newGhostData[dit()].getVolData(),
                                          oldGhostGraph[dit()]);
    }
}
/******************/
//checks to see the vofs are in the correct cells.
//checks to see that the faces are over the correct cells
//checks that volume fractions, area fractions are positive
//bail out with MayDay::Error if anything fails
/******************/
void
EBISLevel::sanityCheck(  const EBIndexSpace* const ebisPtr)
{
#if 0
  CH_TIME("EBISLevel::sanityCheck");
  EBISLayout ghostLayout;

  ebisPtr->fillEBISLayout(ghostLayout, m_grids, m_domain, 1);

  Real maxcentval = 0.5+ s_tolerance;
  for(DataIterator dit = m_grids.dataIterator();  dit.ok(); ++dit)
    {
      const Box& thisBox = m_grids.get(dit());
      const EBISBox& ebisBox = ghostLayout[dit()];
      for(BoxIterator bit(thisBox); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          Vector<VolIndex> vofs = ebisBox.getVoFs(iv);
          for(int ivof = 0; ivof < vofs.size(); ivof++)
            {
              const VolIndex& vof = vofs[ivof];
              if(vof.gridIndex() != iv)
                {
                  pout() << "EBISLevel::sanityCheck: Error" << endl;
                  pout() << "VoF at Intvect = " << iv
                         << "has grid index = " << vof.gridIndex() << endl;
                  MayDay::Error("EBISLevel::sanityCheck Error 1");
                }
              if(vof.cellIndex() < 0)
                {
                  pout() << "EBISLevel::sanityCheck: Error" << endl;
                  pout() << "VoF at Intvect = " << iv
                         << "has negative cell index = " << vof.cellIndex() << endl;
                  MayDay::Error("EBISLevel::sanityCheck Error 2");
                }
              Real volFrac = ebisBox.volFrac(vof);
              if(volFrac < 0.0)
                {
                  pout() << "EBISLevel::sanityCheck: Error" << endl;
                  pout() << "VoF at Intvect = " << iv
                         << "has invalid volume fraction = " << volFrac << endl;
                  MayDay::Error("EBISLevel::sanityCheck Error 5");
                }
              RealVect volCentroid = ebisBox.centroid(vof);
              for(int idir = 0; idir < SpaceDim; idir++)
                {
                  Real volcentdir = volCentroid[idir];
                  if(volFrac > s_tolerance)
                    {
                      if(volcentdir > maxcentval || volcentdir < -maxcentval)
                        {
                          pout() << "EBISLevel::sanityCheck: Error" << endl;
                          pout() << "VoF at Intvect = " << iv
                                 << " has invalid vol centroid = " << volcentdir
                                 << " at direction "<< idir << endl;
                          MayDay::Error("EBISLevel::sanityCheck Error 51");
                        }
                    }
                }
              for(int idir = 0; idir < SpaceDim; idir++)
                {
                  for(SideIterator sit; sit.ok(); ++sit)
                    {
                      Vector<FaceIndex> faces = ebisBox.getFaces(vof, idir, sit());
                      IntVect iv2 = iv + sign(sit())*BASISV(idir);
                      ////check for regular next to covered and multivalued next to regular
                      if(m_domain.contains(iv2))
                        {
                          if(ebisBox.isRegular(iv))
                            {
                              if(ebisBox.isCovered(iv2))
                                {
                                  pout() << iv << " is regular and " <<  iv2 << " is covered" << endl;
                                  MayDay::Error("EBISLevel::sanityCheck error 420 " );
                                }
                              else
                                {
                                  Vector<VolIndex> otherVoFs = ebisBox.getVoFs(iv2);
                                  if(otherVoFs.size() > 1)
                                    {
                                      pout() << iv << " is regular and " <<  iv2 << " is multivalued" << endl;
                                      MayDay::Error("EBISLevel::sanityCheck error 420.2 " );
                                    }
                                }
                            }
                        }
                      IntVect ivlo, ivhi;
                      if(sit() == Side::Lo)
                        {
                          ivlo = iv2;
                          ivhi = iv;
                        }
                      else
                        {
                          ivlo = iv;
                          ivhi = iv2;
                        }
                      for(int iface = 0; iface < faces.size(); iface++)
                        {
                          const FaceIndex& face = faces[iface];
                          if(face.gridIndex(Side::Lo) != ivlo)
                            {
                              pout() << "EBISLevel::sanityCheck: Error" << endl;
                              pout() << "face at IntVects = " << ivlo << "  " << ivhi
                                     << "has low IntVect  = " << face.gridIndex(Side::Lo)
                                     << endl;
                              MayDay::Error("EBISLevel::sanityCheck Error 3");
                            }
                          if(face.gridIndex(Side::Hi) != ivhi)
                            {
                              pout() << "EBISLevel::sanityCheck: Error" << endl;
                              pout() << "face at IntVects = " << ivlo << "  " << ivhi
                                     << "has high IntVect = " << face.gridIndex(Side::Hi)
                                     << endl;
                              MayDay::Error("EBISLevel::sanityCheck Error 4");
                            }
                          Real areaFrac = ebisBox.areaFrac(face);
                          if(areaFrac  < 0.0)
                            {
                              pout() << "EBISLevel::sanityCheck: Error" << endl;
                              pout() << "VoF at Intvect = " << iv
                                     << "has invalid area fraction = " << areaFrac << endl;
                              MayDay::Error("EBISLevel::sanityCheck Error 51");
                            }
                          if(areaFrac  >  s_tolerance)
                            {
                              RealVect faceCentroid = ebisBox.centroid(face);
                              for(int idir = 0; idir < SpaceDim; idir++)
                                {
                                  Real facecentdir = faceCentroid[idir];
                                  if(facecentdir > maxcentval || facecentdir < -maxcentval)
                                    {
                                      pout() << "EBISLevel::sanityCheck: Error" << endl;
                                      pout() << "VoF at Intvect = " << iv
                                             << " has invalid face centroid = " << facecentdir
                                             << " at direction "<< idir << endl;
                                      MayDay::Error("EBISLevel::sanityCheck Error 51");
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#endif
}

EBISLevel::EBISLevel()
{
  m_cacheMisses=0; m_cacheHits=0; m_cacheStale=0;
}
/******************/
//steps to coarsen an ebislevel:
//1. coarsen vofs
//1a. make a layout over refine(mydbl,2) and copy
//    fine layout into it
//1b.do connectivity bizbaz to make my vofs, volfrac, vof->fineVofs
//2. make faces doing connectivity jive
//2.23 make coarse geometric stuff (centroids and all that) from fine
//3. make coarse layout from coarsen(finelayout). and copy my data into it
//   to make fine->coarserVoF
// (so finer ebislevel does change in this function)
/******************/
void
EBISLevel::
coarsenVoFs(EBISLevel& a_fineEBIS)
{
  CH_TIME("EBISLevel::coarsenVoFs");
  //so that i can do the vofs and faces of this level.
  DisjointBoxLayout fineFromCoarDBL;
  refine(fineFromCoarDBL, m_grids, 2);
  fineFromCoarDBL.close();

  //no need for ghost cells here except to define the face data
  //you need the graph to be one bigger
  EBGraphFactory ebgraphfactFine(a_fineEBIS.m_domain);
  LevelData<EBGraph> fineFromCoarEBGraph(fineFromCoarDBL,1, IntVect::Unit, ebgraphfactFine);

  Interval interv(0,0);
  a_fineEBIS.m_graph.copyTo(interv, fineFromCoarEBGraph, interv);

  for(DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph& fineEBGraph = fineFromCoarEBGraph[dit()];
      const Box& coarRegion      = m_grids.get(dit());
      EBGraph& coarEBGraph = m_graph[dit()];
      coarEBGraph.coarsenVoFs(fineEBGraph, coarRegion);
    }
  EBGraphFactory ebgraphfactCoar(m_domain);
  LevelData<EBGraph> coarGhostEBGraph(m_grids,1, IntVect::Unit, ebgraphfactCoar);
  m_graph.copyTo(interv, coarGhostEBGraph, interv);

  EBDataFactory ebdatafact;
  LevelData<EBData> fineFromCoarEBData(fineFromCoarDBL,1, IntVect::Zero, ebdatafact);
  for(DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& localBox = fineFromCoarDBL.get(dit());
      fineFromCoarEBData[dit()].defineVoFData(fineFromCoarEBGraph[dit()],  localBox);
      fineFromCoarEBData[dit()].defineFaceData(fineFromCoarEBGraph[dit()], localBox);
    }
  a_fineEBIS.m_data.copyTo(interv, fineFromCoarEBData, interv);

  for(DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph& fineEBGraph =  fineFromCoarEBGraph[dit()];
      const EBData& fineEBData = fineFromCoarEBData[dit()];
      const EBGraph& coarEBGraph = coarGhostEBGraph[dit()];

      m_data[dit()].coarsenVoFs(fineEBData, fineEBGraph, coarEBGraph, m_grids.get(dit()));
    }
}
/******************/
void
EBISLevel::
fixFineToCoarse(EBISLevel& a_fineEBIS)
{
  CH_TIME("EBISLevel::fixFineToCoarse");
  // make a coarse layout from the fine layout so that we
  //can fix the fine->coarseVoF thing.
  DisjointBoxLayout coarFromFineDBL;

  coarsen(coarFromFineDBL, a_fineEBIS.m_graph.getBoxes(), 2);
  coarFromFineDBL.close();
  EBGraphFactory ebgraphfact(m_domain);
  LevelData<EBGraph> coarFromFineEBGraph(coarFromFineDBL,1, IntVect::Zero, ebgraphfact);
  Interval interv(0,0);
  m_graph.copyTo(interv, coarFromFineEBGraph, interv);

  for(DataIterator dit = a_fineEBIS.m_grids.dataIterator(); dit.ok(); ++dit)
    {
      EBGraph& fineEBGraph       =  a_fineEBIS.m_graph[dit()];
      const EBGraph& coarEBGraph = coarFromFineEBGraph[dit()];

      coarEBGraph.fixFineToCoarse(fineEBGraph);
    }

}
/******************/
void
EBISLevel::
coarsenFaces(EBISLevel& a_fineEBIS)
{
  CH_TIME("EBISLevel::coarsenFaces");
  //now make a fine ebislayout with two ghost cell
  //on the same mapping as m_dbl
  //so that i can do the vofs and faces of this level.
  //this one will have the fine from coarse stuff fixed
  DisjointBoxLayout fineFromCoarDBL;
  refine(fineFromCoarDBL, m_grids, 2);
  fineFromCoarDBL.close();

  //no need for ghost cells here
  EBGraphFactory ebgraphfactfine(a_fineEBIS.m_domain);
  EBGraphFactory ebgraphfactcoar(m_domain);
  LevelData<EBGraph> fineEBGraphGhostLD(fineFromCoarDBL,1,3*IntVect::Unit, ebgraphfactfine);
  Interval interv(0,0);
  a_fineEBIS.m_graph.copyTo(interv, fineEBGraphGhostLD, interv);
  LevelData<EBGraph> coarEBGraphGhostLD(m_grids,        1,  IntVect::Unit, ebgraphfactcoar);
  m_graph.copyTo(           interv, coarEBGraphGhostLD, interv);

  for(DataIterator dit =  m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph& fineEBGraphGhost = fineEBGraphGhostLD[dit()];
      const EBGraph& coarEBGraphGhost = coarEBGraphGhostLD[dit()];
      EBGraph& coarEBGraph = m_graph[dit()];
      coarEBGraph.coarsenFaces(coarEBGraphGhost, fineEBGraphGhost);
    }
  //redefine coarebghostgraphld so i can use the faces for the ebdata
  coarEBGraphGhostLD.define(m_grids, 1,  IntVect::Unit, ebgraphfactcoar);
  m_graph.copyTo( interv, coarEBGraphGhostLD, interv);

  EBDataFactory ebdatafact;
  LevelData<EBData> fineEBDataGhostLD(fineFromCoarDBL,1, 2*IntVect::Unit, ebdatafact);
  for(DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      Box localBox = grow(fineFromCoarDBL.get(dit()), 2);
      localBox &= a_fineEBIS.m_domain;
      fineEBDataGhostLD[dit()].defineVoFData( fineEBGraphGhostLD[dit()], localBox);;
      fineEBDataGhostLD[dit()].defineFaceData(fineEBGraphGhostLD[dit()], localBox);
    }
  a_fineEBIS.m_data.copyTo(interv, fineEBDataGhostLD, interv);

  for(DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBData&   fineEBData      = fineEBDataGhostLD[dit()];
      const EBGraph& fineEBGraphGhost = fineEBGraphGhostLD[dit()];
      const EBGraph& coarEBGraphGhost = coarEBGraphGhostLD[dit()];

      EBData& coarEBData   = m_data[dit()];
      coarEBData.coarsenFaces(fineEBData,  fineEBGraphGhost, coarEBGraphGhost, m_grids.get(dit()));
    }

}
/******************/
EBISLevel::EBISLevel(EBISLevel& a_fineEBIS, const GeometryService& a_geoserver,
                     const EBIndexSpace* const ebisPtr)
{
  CH_TIME("EBISLevel::EBISLevel_fineEBIS");
  m_cacheMisses=0; m_cacheHits=0; m_cacheStale=0;
  m_domain = coarsen(a_fineEBIS.m_domain,2);
  m_dx = 2.*a_fineEBIS.m_dx;
  m_tolerance = 2.*a_fineEBIS.m_tolerance;
  m_origin = a_fineEBIS.m_origin;
  Vector<Box> vbox;
  Vector<long> irregCount;
  {
    CH_TIME("EBISLevel::EBISLevel_fineEBIS_makeboxes 2");
    makeBoxes(vbox, irregCount, m_domain.domainBox(), m_domain, a_geoserver,
              m_origin, m_dx, ebisPtr->getNCellMax(), ebisPtr);
  }

  //pout()<<vbox<<"\n\n";
  //load balance the boxes
  Vector<int> procAssign;
  LoadBalance(procAssign, irregCount, vbox);
 
  //pout()<<procAssign<<std::endl;
  //define the layout.  this includes the domain and box stuff
  m_grids.define(vbox, procAssign);//this should use m_domain for periodic

  EBGraphFactory ebgraphfact(m_domain);
  m_graph.define(m_grids, 1, IntVect::Zero, ebgraphfact);

  EBDataFactory ebdatafact;
  m_data.define(m_grids, 1, IntVect::Zero, ebdatafact);

  //create coarsened vofs from fine.
  coarsenVoFs(a_fineEBIS);

  //overallMemoryUsage();
  //create coarse faces from fine
  coarsenFaces(a_fineEBIS);
  //overallMemoryUsage();
  //fix the regular next to the multivalued cells
  //to be full irregular cells
  fixRegularNextToMultiValued();
  //overallMemoryUsage();
  // fix the fine->coarseVoF thing.
  fixFineToCoarse(a_fineEBIS);
}
/******************/
EBISLevel::~EBISLevel()
{
}
/****************/
void
EBISLevel::fillEBISLayout(EBISLayout& a_ebisLayout,
                          const DisjointBoxLayout& a_grids,
                          const int& a_nghost) const
{
  CH_assert(a_nghost >= 0);
a_ebisLayout.define(m_domain, a_grids, a_nghost, m_graph, m_data);
return;
  dmap& dcache = cache[a_nghost];
  EBISLayout& l = dcache[a_grids];
  if(!l.isDefined())
    {
      l.define(m_domain, a_grids, a_nghost, m_graph, m_data);
      m_cacheMisses++;
      m_cacheStale++;
    }
  else
    {
      m_cacheHits++;
    }
  a_ebisLayout = l;
  if(m_cacheStale == 1)
    {
      refreshCache();
      m_cacheStale = 0;
    }
}

void EBISLevel::refreshCache() const
{
  for(std::map<int, dmap>::iterator p = cache.begin(); p!= cache.end(); ++p)
    {
      for(dmap::iterator d = p->second.begin(); d!= p->second.end(); ++d)
        {
          if(d->first.refCount() == 1)
            {
              p->second.erase(d);
            }
        }
    }
}

void EBISLevel::clearMultiBoundaries()
{
  CH_TIME("EBISLevel::clearMultiBoundaries");
  DataIterator dit = m_grids.dataIterator();
  for(dit.begin();dit.ok(); ++dit)
    {
      EBData& data = m_data[dit];
      data.clearMultiBoundaries();
    }
}

void EBISLevel::setBoundaryPhase(int phase)
{
  CH_TIME("EBISLevel::setBoundaryPhase");
  DataIterator dit = m_grids.dataIterator();
  for(dit.begin();dit.ok(); ++dit)
    {
      EBData& data = m_data[dit];
      data.setBoundaryPhase(phase);
    }
}

//throughout this routine A refers to objects belonging to this
//ebindexspace,  B refers to the other one.
VolIndex
getOtherVoFSingleValued(const VolData& a_vol, const VolIndex& a_sourceVoF, bool& a_skipCell)
{
  VolIndex retval = a_sourceVoF;
  Real eps = 1.0e-10;
  Real bitLess = 1.0 - eps;
  a_skipCell = false;
  //if a cell is full and has a unit boundary area, it means that the irregular
  //boundary lives on a cell boundary.  We have to link to a cell over.
  // In this case, in the presence of multi valued cells, to quote Yeats:
  //Mere anarchy is loosed upon the world,
  //The blood-dimmed tide is loosed, and everywhere
  //The ceremony of innocence is drowned.
  if((a_vol.m_volFrac  > bitLess) && (a_vol.m_averageFace.m_bndryArea > bitLess))
    {
      a_skipCell = true;
      //figure out which cell we are moving to.   Then we use cellIndex = 0 because
      //this is just for the special case where the cells involved are single valued.
      int whichWayToGo = -23;
      int sign = 1;
      bool found = false;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          Real normDir = a_vol.m_averageFace.m_normal[idir];
          if(Abs(normDir) > bitLess)
            {
              found = true;
              whichWayToGo = idir;
              sign = 1;
              if(normDir < 0)
                {
                  sign = -1;
                }
            }
        }
      if(!found)
        {
          MayDay::Error("Unit normal direction not found");
        }
      IntVect ivOther = sign*BASISV(whichWayToGo);
      ivOther += a_sourceVoF.gridIndex();
      retval = VolIndex(ivOther, 0);
    }
  return retval;
}

void
EBISLevel::levelStitch(EBISLevel&       a_otherPhase,
                       const EBISLevel* a_finePtrA,
                       const EBISLevel* a_finePtrB)
{
  CH_TIME("EBISLevel::levelStitch");
  //hook to keep the tests working.
  //pout() << "levelstitch?  we don't need no stinking levelstitch!" << endl;
  //return;

  //I have no idea what to do if only one of the inputs is null.
  //either both or neither makes sense
  CH_assert(((a_finePtrA != NULL) && (a_finePtrB != NULL)) ||
            ((a_finePtrA == NULL) && (a_finePtrB == NULL)) );
  EBISLayout ebislFineA, ebislFineB;
  DisjointBoxLayout dblFineA, dblFineB;
  if(a_finePtrA != NULL)
    {
      int nghost = 0; //should not need any as this is all about connectivity within a coarse cell
      refine(dblFineA,              m_grids, 2);
      refine(dblFineB, a_otherPhase.m_grids, 2);
      a_finePtrA->fillEBISLayout(ebislFineA, dblFineA, nghost);
      a_finePtrB->fillEBISLayout(ebislFineB, dblFineB, nghost);
    }

  //both layouts are made independently so will not work
  //with the same data iterator but they ought to because they
  //are both breaking up the domain the same way.   This is yet
  //another subtle thing thing that will break when the
  //number of fluids is greater than two.
  DataIterator dita = m_grids.dataIterator();
  DataIterator ditb= a_otherPhase.m_grids.dataIterator();
  dita.begin(); ditb.begin();
  for(;dita.ok(); ++dita, ++ditb)
    {
      EBGraph& ebgrapCoarA =              m_graph[dita];
      EBGraph& ebgrapCoarB = a_otherPhase.m_graph[ditb];
      EBData&  ebdataCoarA =              m_data[dita];
      EBData&  ebdataCoarB = a_otherPhase.m_data[ditb];
      Box region = ebgrapCoarA.getRegion();

      CH_assert(ebgrapCoarB.getRegion() == region);
      IntVectSet seta = ebgrapCoarA.getIrregCells(region);
      IntVectSet setb = ebgrapCoarB.getIrregCells(region);

      //different EBIS's can (correctly) disagree about which cells are multivalued
      IntVectSet setMultiA = ebgrapCoarA.getMultiCells(region);
      IntVectSet setMultiB = ebgrapCoarB.getMultiCells(region);
      IntVectSet setMulti = setMultiA | setMultiB;
      int aphase =              m_phase;
      int bphase = a_otherPhase.m_phase;

      IntVectSet setANoMulti = seta - setMulti;
      IntVectSet setBNoMulti = setb - setMulti;
      {
        //EBIndexSpace does the right thing with all the geometric information
        //in the case of two single valued vofs.   All that needs to be set is
        //the phase on the other side of the irregular face.
        for(IVSIterator ita(setANoMulti); ita.ok(); ++ita)
          {
            VolIndex v(ita(), 0); //0 because we know this whole set is single valued

            VolData& volDatA = ebdataCoarA.getVolData()(v,0);
            volDatA.m_phaseFaces.resize(1);
            //this sets boundary area and normal and centroid to
            //whatever the ebindexspace set it to.
            volDatA.m_phaseFaces[0]=volDatA.m_averageFace;
            bool skip = false;
            VolIndex vother = getOtherVoFSingleValued(volDatA, v, skip);
            if(skip && setMultiB.contains(vother.gridIndex()))
              {
                MayDay::Error("coordinate face boundary also a multi valued one");
              }
            volDatA.m_phaseFaces[0].m_volIndex   = vother;
            volDatA.m_phaseFaces[0].m_bndryPhase = bphase;
          }//end loop over cells  where both phases are singlevalued

        for(IVSIterator itb(setBNoMulti); itb.ok(); ++itb)
          {
            VolIndex v(itb(), 0); //0 because we know this whole set is single valued

            VolData& volDatB = ebdataCoarB.getVolData()(v,0);
            volDatB.m_phaseFaces.resize(1);
            //this sets boundary area and normal and centroid to
            //whatever the ebindexspace set it to.
            volDatB.m_phaseFaces[0]=volDatB.m_averageFace;
            bool skip = false;
            VolIndex vother = getOtherVoFSingleValued(volDatB, v, skip);
            if(skip && setMultiA.contains(vother.gridIndex()))
              {
                MayDay::Error("coordinate face boundary also a multi valued one");
              }
            volDatB.m_phaseFaces[0].m_volIndex   = vother;
            volDatB.m_phaseFaces[0].m_bndryPhase = aphase;
          }//end loop over cells where both phases are singlevalued
      }
      if(!setMulti.isEmpty())
        {
          //now have to do the union of the each set
          //of multivalued cells.   either phase being multivalued can
          //make either ebindex space get the wrong answer for
          //the geometric information
          const EBISBox& ebisbxFineA = ebislFineA[dita()];
          const EBISBox& ebisbxFineB = ebislFineB[ditb()];

          IVSIterator it(setMulti); //see above derivation.
          for(it.begin(); it.ok(); ++it)
            {
              cellStitch(ebdataCoarA, ebdataCoarB,
                         ebgrapCoarA, ebgrapCoarB,
                         ebisbxFineA, ebisbxFineB,
                         it(), aphase, bphase);
            }//end loop over multivalued cells
        }
    } //end loop over grids
}
/*****/
void
getFinerBoundaryData(Vector<Vector<BoundaryData> > & a_bndryDataFineA,
                     Vector<VolIndex>&               a_otherCellIndex,
                     const VolIndex&                 a_coarVoFA,
                     const EBGraph&                  a_ebgrapCoarA,
                     const EBGraph&                  a_ebgrapCoarB,
                     const EBISBox&                  a_ebisbxFineA)
{
  //divdes fine boundary datas a by which coarse vof b they are connected to.
  // to which they are connected?

  const IntVect ivCoar = a_coarVoFA.gridIndex();

  Vector<VolIndex>          coarVoFsB = a_ebgrapCoarB.getVoFs(ivCoar);
  Vector<Vector<VolIndex> > fineVoFsB(coarVoFsB.size());
  a_bndryDataFineA.resize(coarVoFsB.size());
  a_otherCellIndex.resize(coarVoFsB.size());
  for(int ib = 0; ib < coarVoFsB.size(); ib++)
    {
      fineVoFsB[ib] = a_ebgrapCoarB.refine(coarVoFsB[ib]);
      a_otherCellIndex[ib] = coarVoFsB[ib];
    }

  Vector<VolIndex> allFinerVoFsA = a_ebgrapCoarA.refine(a_coarVoFA);
  for(int ib = 0; ib < coarVoFsB.size(); ib++)
    {
      a_bndryDataFineA[ib].resize(0);
      for(int ifineB = 0; ifineB < fineVoFsB[ib].size(); ifineB++)
        {
          const VolIndex& vofFineB = fineVoFsB[ib][ifineB];
          for(int ifineA = 0; ifineA < allFinerVoFsA.size(); ifineA++)
            {
              const VolIndex& vofFineA = allFinerVoFsA[ifineA];
              if(vofFineA.gridIndex() == vofFineB.gridIndex())
                {
                  if(a_ebisbxFineA.isIrregular(vofFineA.gridIndex()))
                    {
                      const VolData& voldat = a_ebisbxFineA.getEBData().getVolData()(vofFineA, 0);
                      const Vector<BoundaryData>& boundaryDat = voldat.m_phaseFaces;
                      for(int iface = 0; iface < boundaryDat.size(); iface++)
                        {
                          if(boundaryDat[iface].m_volIndex.cellIndex() == vofFineB.cellIndex())
                            {
                              a_bndryDataFineA[ib].push_back(boundaryDat[iface]);
                            }
                        }
                    }
                }
            }
        }
    }
} //
/*****/
void
coarsenBoundaryInfo(EBData&                                a_ebdataCoar,
                    const EBISBox&                         a_ebisbxFine,
                    const Vector< Vector<BoundaryData > >& a_finerBoundaryData,
                    const VolIndex&                        a_vof,
                    const int&                             a_otherPhase,
                    const Vector<VolIndex>&                a_otherCellIndex)
{
  //this coarsening factor for boundary area fractions can be thought of this way
  //(1) take fine area and multiply by dxf^{SD-1}
  //(2) add up areas (now we have the real coarse area).
  //(3) divide by dxc^{SD-1} = dxf{SD-1}*2{SD-1}
  Real faceCoarsenFactor = D_TERM(1.0, * 2.0,* 2.0);

  //the point of this routine is to fix the boundary area and normal
  //stuff in this object
  VolData& coarData = a_ebdataCoar.getVolData()(a_vof, 0);
  coarData.m_phaseFaces.resize(a_finerBoundaryData.size());
  for(int ib = 0; ib < a_finerBoundaryData.size(); ib++)
    {
      //initialization is important  in both of these things.
      //since i am taking the average
      RealVect aveNormal  = RealVect::Zero;
      RealVect aveCentro  = RealVect::Zero;
      Real sumArea=0;
      Real aveArea=0;
      int numfaces=0;
      for(int jb = 0; jb < a_finerBoundaryData[ib].size(); jb++)
        {
          const BoundaryData& boundaryDat = a_finerBoundaryData[ib][jb];
          const Real    & fineArea = boundaryDat.m_bndryArea;
          const RealVect& fineNorm = boundaryDat.m_normal;
          const RealVect& fineCent = boundaryDat.m_bndryCentroid;
          sumArea += fineArea;
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              aveNormal[idir] += fineArea*fineNorm[idir];
              aveCentro[idir] += fineArea*fineCent[idir];
            }
          numfaces++;
        }

      //we are taking an average-weighted average
      //of each normal and centroid components so we need to divide by the area
      if(sumArea > 1.0e-12)
        {
          aveNormal /= sumArea;
          aveCentro /= sumArea;
        }
      Real sumSq;
      PolyGeom::unifyVector(aveNormal, sumSq);
      //this takes into account the fact that boundary areas are normalized
      //by dx^(SD-1).  See note above.
      aveArea  = sumArea/faceCoarsenFactor;

      coarData.m_phaseFaces[ib].m_normal            = aveNormal;
      coarData.m_phaseFaces[ib].m_bndryCentroid     = aveCentro;
      coarData.m_phaseFaces[ib].m_bndryArea         = aveArea;
      coarData.m_phaseFaces[ib].m_bndryPhase        = a_otherPhase;
      coarData.m_phaseFaces[ib].m_volIndex          = a_otherCellIndex[ib];
    }

}

/***/
void
EBISLevel::
cellStitch(EBData&           a_ebdataCoarA,
           EBData&           a_ebdataCoarB,
           const EBGraph&    a_ebgrapCoarA,
           const EBGraph&    a_ebgrapCoarB,
           const EBISBox&    a_ebisbxFineA,
           const EBISBox&    a_ebisbxFineB,
           const IntVect&    a_iv,
           const int&        a_phaseA,
           const int&        a_phaseB)
{
  //pout() << "entering cellStitch" << endl;
  Vector<VolIndex> coarVoFsA = a_ebgrapCoarA.getVoFs(a_iv);
  Vector<VolIndex> coarVoFsB = a_ebgrapCoarB.getVoFs(a_iv);
  for(int ivof = 0; ivof < coarVoFsA.size(); ivof++)
    {
      const VolIndex& vof = coarVoFsA[ivof];
      // each element in finerVoFs is a list of finer vofs whose
      // opposing number on the other side of the irregular boundary
      // are connected to each other.   This will correspond with
      // the Vector<VolIndex> of the coarse vofs in the *other* phase.
      Vector<Vector< BoundaryData> > finerBndryData;
      Vector<VolIndex>  cellIndexB;
      getFinerBoundaryData(finerBndryData, cellIndexB, vof, a_ebgrapCoarA, a_ebgrapCoarB, a_ebisbxFineA);


      //once the above list is made, we can correctly coarsen the boundary
      //information using only the information in this phase.
      coarsenBoundaryInfo(a_ebdataCoarA, a_ebisbxFineA, finerBndryData, vof, a_phaseB, cellIndexB);
    }

  for(int ivof = 0; ivof < coarVoFsB.size(); ivof++)
    {
      const VolIndex& vof = coarVoFsB[ivof];
      // each element in finerVoFs is a list of finer vofs whose
      // opposing number on the other side of the irregular boundary
      // are connected to each other.   This will correspond with
      // the Vector<VolIndex> of the coarse vofs in the *other* phase.
      Vector<Vector< BoundaryData> > finerBndryData;
      Vector<VolIndex>  cellIndexA;
      getFinerBoundaryData(finerBndryData, cellIndexA, vof, a_ebgrapCoarB, a_ebgrapCoarA, a_ebisbxFineB);


      //once the above list is made, we can correctly coarsen the boundary
      //information using only the information in this phase.
      coarsenBoundaryInfo(a_ebdataCoarB, a_ebisbxFineB, finerBndryData, vof, a_phaseA, cellIndexA);
    }

}

/******************/
bool
EBIndexSpace::isDefined() const
{
  return m_isDefined;
}
/******************/
int
EBIndexSpace::getNCellMax() const
{
  return m_nCellMax;
}
/******************/
#ifdef CH_USE_HDF5
void
EBIndexSpace::define(HDF5Handle& a_handle, int a_maxCoarsenings)
{
  CH_TIME("EBIndexSpace::define_hdf5handle");
  clear();

  //read in ncellmax
  HDF5HeaderData header;
  header.readFromFile(a_handle);
  int nCellMax = header.m_int["EBIS_ncellMax"];

  //read in finest level
  //coarser levels are derived from graph coarsening
  EBISLevel* level0 = new EBISLevel(a_handle);

  define(level0, nCellMax, a_maxCoarsenings);
}
/******************/
void
EBIndexSpace::write(HDF5Handle& a_handle,
                    ProblemDomain a_outputLevel)
{
  CH_TIME("EBIndexSpace::write");
  //read in ncellmax
  HDF5HeaderData header;
  std::string filedescriptor("EBIndexSpace");
  header.m_string ["filetype"]      = filedescriptor;
  header.m_int["EBIS_ncellMax"] = m_nCellMax;
  header.writeToFile(a_handle);

  int ilev = 0;
  if(!a_outputLevel.isEmpty())
    {
      ilev = getLevel(a_outputLevel);
    }
  //write out finest level
  //coarser levels are derived from graph coarsening
  m_ebisLevel[ilev]->write(a_handle);
}

#endif

void
EBIndexSpace::define(EBISLevel* level0, int nCellMax, int a_maxCoarsenings)
{
  CH_TIME("EBIndexSpace::define_ebislevel0");
  m_nCellMax = nCellMax;
  m_isDefined = true;
  const ProblemDomain&        fineDomain = level0->getDomain();

  //figure out how deep we can go
  m_nlevels = 1;
  bool canref = (fineDomain == refine(coarsen(fineDomain,2), 2));
  CH_assert(!fineDomain.isEmpty());
  ProblemDomain refbox = fineDomain;
  while(canref)
    {
      ProblemDomain refcoarbox = coarsen(refbox,2);
      refcoarbox.refine(2);
      if(refcoarbox != refbox)
        {
          canref = false;
        }
      else
        {
          m_nlevels++;
          refbox.coarsen(2);
        }
    }
  if(a_maxCoarsenings != -1)
    {
      CH_assert(a_maxCoarsenings >= 0);
      if(m_nlevels > a_maxCoarsenings+1) m_nlevels = a_maxCoarsenings + 1;
    }

  m_ebisLevel.resize(m_nlevels, NULL);
  m_domainLevel.resize(m_nlevels);

  ProblemDomain  domLevel = fineDomain;
  level0->clearMultiBoundaries();
  m_ebisLevel[0] = level0;
  m_domainLevel[0] = domLevel;
  AllRegularService dummy;
  for(int ilev = 1; ilev < m_nlevels; ilev++)
    {
      domLevel.coarsen(2);
      m_domainLevel[ilev] = domLevel;
      m_ebisLevel[ilev] = new EBISLevel(*m_ebisLevel[ilev-1], dummy, this);
      m_ebisLevel[ilev]->clearMultiBoundaries();
    }

#ifndef NDEBUG
  for(int ilev = 0; ilev < m_nlevels; ilev++)
    {
      m_ebisLevel[ilev]->sanityCheck(this);
    }
#endif
}

/******************/
void
EBIndexSpace::define(const ProblemDomain& a_domain,
                     const RealVect& a_origin,
                     const Real& a_dx,
                     const GeometryService& a_geoserver,
                     int a_nCellMax,
                     int a_maxCoarsenings)
{
  CH_TIME("EBIndexSpace::define_geoserver_domain0");
  buildFirstLevel(a_domain, a_origin, a_dx, a_geoserver, a_nCellMax, a_maxCoarsenings);
  m_ebisLevel[0]->clearMultiBoundaries();
  EBISLevel* n =  m_ebisLevel[0];
  while(n)
    {
      n->clearMultiBoundaries();
      n=buildNextLevel(a_geoserver);
    }


}

EBISLevel* EBIndexSpace::buildFirstLevel(const ProblemDomain& a_domain,
                                         const RealVect& a_origin,
                                         const Real& a_dx,
                                         const GeometryService& a_geoserver,
                                         int a_nCellMax,
                                         int a_maxCoarsenings)
{
  CH_TIME("EBIndexSpace::buildFirstLevel");
  clear();
  m_isDefined = true;

  if(a_nCellMax > 0)
    {
      m_nCellMax = a_nCellMax;
    }
  else
    {
      if(SpaceDim == 2)
        {
          m_nCellMax = 32;
        }
      else
        {
          m_nCellMax = 32;
        }
    }
  m_nlevels = 1;
  bool canref = (a_domain == refine(coarsen(a_domain,2), 2));

  CH_assert(!a_domain.isEmpty());
  ProblemDomain refbox = a_domain;

  while(canref)
    {
      ProblemDomain refcoarbox = coarsen(refbox,2);
      refcoarbox.refine(2);
      if(refcoarbox != refbox)
        {
          canref = false;
        }
      else
        {
          m_nlevels++;
          refbox.coarsen(2);
        }
    }
  if(a_maxCoarsenings != -1)
    {
      CH_assert(a_maxCoarsenings >= 0);
      if(m_nlevels > a_maxCoarsenings+1) m_nlevels = a_maxCoarsenings + 1;
    }

  m_ebisLevel.resize(m_nlevels, NULL);
  m_domainLevel.resize(m_nlevels);

  ProblemDomain  domLevel = a_domain;
  m_ebisLevel[0] = new EBISLevel(domLevel, a_origin,
                                 a_dx, a_geoserver, this);
  m_domainLevel[0] = domLevel;
  return m_ebisLevel[0];

}

void EBIndexSpace::resetLevels(int nLevel)
{
  CH_TIME("EBIndexSpace::resetLevels");
  for(int i=nLevel; i<m_nlevels; i++)
    {
      if(m_ebisLevel[i] != NULL) delete m_ebisLevel[i];
      m_ebisLevel[i] = NULL;
    }

  m_nlevels = nLevel;

}

EBISLevel* EBIndexSpace::buildNextLevel(const GeometryService& a_geoserver)
{
  CH_TIME("EBIndexSpace::buildNextLevel");
  int ilev=0;
  for(; ilev <m_ebisLevel.size(); ++ilev)
    {
      if(m_ebisLevel[ilev] == NULL) break;
    }
  if(ilev == m_ebisLevel.size()) return NULL;

  m_domainLevel[ilev] = m_domainLevel[ilev-1];
  m_domainLevel[ilev].coarsen(2);
  m_ebisLevel[ilev] = new EBISLevel(*m_ebisLevel[ilev-1], a_geoserver, this);
  return m_ebisLevel[ilev];
}

/******************/
const ProblemDomain&
EBISLevel::getDomain() const
{
  return m_domain;
}
/******************/
const Real&
EBISLevel::getDX() const
{
  return m_dx;
}
/******************/
const RealVect&
EBISLevel::getOrigin() const
{
  return m_origin;
}
/******************/
const DisjointBoxLayout&
EBISLevel::getGrids() const
{
  return m_grids;
}

DisjointBoxLayout
EBISLevel::getIrregGrids() const
{
  DisjointBoxLayout irregGrids;
  DataIterator dit =   m_graph.dataIterator();
  Vector<Box> localBoxes;

  for(dit.begin(); dit.ok(); ++dit)
    {
      const EBGraph& graph = m_graph[dit()];
      const Box& b         = m_graph.box(dit());
      if(graph.hasIrregular())
        {
          localBoxes.push_back(b);

        }
    }
  Vector<Vector<Box> > allBoxes;

  gather(allBoxes, localBoxes, 0);

  broadcast(allBoxes, 0);

  for(int i=0; i<allBoxes.size(); i++)
    {
      int num=allBoxes[i].size();
      Vector<Box>& boxes = allBoxes[i];
      for(int j=0; j<num; j++) irregGrids.addBox(boxes[j], i);
    }

  irregGrids.close();
  return irregGrids;
}



DisjointBoxLayout
EBISLevel::getFlowGrids() const
{
  DisjointBoxLayout flowGrids;
  DataIterator dit =   m_graph.dataIterator();
  Vector<Box> localBoxes;

  for(dit.begin(); dit.ok(); ++dit)
    {
      const EBGraph& graph = m_graph[dit()];
      const Box& b         = m_graph.box(dit());
      if(graph.hasIrregular() || graph.isAllRegular())
        {
          localBoxes.push_back(b);

        }
    }
  Vector<Vector<Box> > allBoxes;

  gather(allBoxes, localBoxes, 0);

  broadcast(allBoxes, 0);

  for(int i=0; i<allBoxes.size(); i++)
    {
      int num=allBoxes[i].size();
      Vector<Box>& boxes = allBoxes[i];
      for(int j=0; j<num; j++) flowGrids.addBox(boxes[j], i);
    }

  flowGrids.close();
  return flowGrids;
}

IntVectSet
EBISLevel::irregCells() const
{
  DataIterator dit =   m_graph.dataIterator();
  IntVectSet rtn;

  for(dit.begin(); dit.ok(); ++dit)
    {
      const EBGraph& graph = m_graph[dit()];
      const Box& b         = m_graph.box(dit());
      if(graph.hasIrregular())
        {
          rtn |= graph.getIrregCells(b);
        }
    }
  return rtn;
}

#ifdef CH_USE_HDF5
void EBIndexSpace::writeInfo(HDF5Handle& handle) const
{
  Vector<DisjointBoxLayout> vectGrids(m_nlevels);
  Vector<LevelData<FArrayBox>* >  vectData(m_nlevels);
  Vector<string> vectNames(1,"unknown");
  ProblemDomain domain = m_domainLevel[m_nlevels-1];
  const EBISLevel& level = *(m_ebisLevel[m_nlevels-1]);
  Real dx = level.getDX();
  Vector<int> vectRatio(m_nlevels, 2);

  for(int i=0; i<m_nlevels; ++i)
    {
      const EBISLevel& level = *(m_ebisLevel[m_nlevels-1-i]);
      vectGrids[i] = level.getGrids();
      vectData[i] = new LevelData<FArrayBox>(level.getGrids(), 1);
      LevelData<FArrayBox>& ld = *(vectData[i]);
      for(DataIterator dit = ld.dataIterator(); dit.ok(); ++dit)
        {
          ld[dit].setVal(0.0);
        }
    }


  WriteAMRHierarchyHDF5(handle, vectGrids, vectData, vectNames, domain.domainBox(),
                        dx, 1.0, 0.0 , vectRatio, m_nlevels);

  for(int i=0; i<m_nlevels; ++i)
    {

      delete vectData[i];

    }
}
#endif

/******************/
int
EBIndexSpace::numLevels() const
{
  return m_nlevels;
}

void
EBIndexSpace::clear()
{
  for(int ilev = 0; ilev < m_ebisLevel.size(); ilev++)
    {
      delete m_ebisLevel[ilev];
      m_ebisLevel[ilev] = NULL;
    }
  m_ebisLevel.resize(0);
  m_domainLevel.resize(0);
  m_nlevels = 0;
  m_isDefined = false;
}
/******************/
EBIndexSpace::EBIndexSpace()
{
}
/******************/
EBIndexSpace::~EBIndexSpace()
{
  clear();
}
/******************/
int
EBIndexSpace::getLevel(const ProblemDomain& a_domain) const
{
  bool found = false;
  int whichlev = -1;
  for(int ilev = 0; ilev < m_domainLevel.size() && !found; ilev++)
    {
      if(m_domainLevel[ilev] == a_domain)
        {
          found = true;
          whichlev = ilev;
        }
    }
  return whichlev;
}
/****************/
void
EBIndexSpace::fillEBISLayout(EBISLayout& a_ebisLayout,
                             const DisjointBoxLayout& a_grids,
                             const ProblemDomain& a_domain,
                             const int& a_nghost) const
{
  CH_TIME("EBIndexSpace::fillEBISLayout");
  //figure out which level we are on
  int whichlev = getLevel(a_domain);
  if(whichlev < 0)
    {
      pout() << "a_domain = " << a_domain
             << " does not correspond to any refinement of EBIS" << endl;
      MayDay::Error("Bad argument to EBIndexSpace::fillEBISLayout");
    }
  m_ebisLevel[whichlev]->fillEBISLayout(a_ebisLayout, a_grids, a_nghost);
  a_ebisLayout.setEBIS(this); //need for mf
}
/****************/
EBIndexSpace* Chombo_EBIS::s_instance = NULL;
/****************/
EBIndexSpace*
Chombo_EBIS::instance()
{
  if(s_instance == NULL)
    s_instance = new EBIndexSpace();

  return  s_instance;
}
/****************/
#include "NamespaceFooter.H"
