//
//  main.cpp
//  graphClustering
//
//  Created by 曹 巍 on 12/06/05.
//  Copyright (c) 2012年 生物情報工学・東京大学. All rights reserved.
//

#include "bilab/common.h"
#include "bilab/cmdline_opt.h"
#include "bilab/gCluster.h"
#include "bilab/utility.h"
#include "bilab/cGraph.h"

using namespace std;
using namespace boost;


int main(int argc, const char * argv[])
{
  int numOfVertices;
  boost::timer timer;

  //Parse command line arguments
  CmdArg::ParseCmdLine(argc,argv);
  numOfVertices = CmdArg::numOfvertex;

  // Define an object of the graph
  cGraph g(numOfVertices);

  //Ready to read the file
  readfile(g,CmdArg::filename.c_str());
  if (CmdArg::verbose){
    std::cout<<"Constructed graph info :"<<std::endl;
    //readfile(g,CmdArg::filename.c_str());
    //std::cout<<" Over."<<std::endl;
    print_graph_info(g);
  }

/*
// Clustering coefficients (not work for the graph including disconnented component)
  typedef exterior_vertex_property<cGraph, float> ClusteringProperty;	
  typedef ClusteringProperty::container_type ClusteringContainer;
  typedef ClusteringProperty::map_type ClusteringMap;
  ClusteringContainer coefs(num_vertices(g));
  ClusteringMap cm(coefs, g);
  //ClusteringMap cm = calculate_clustering_coefficient<cGraph,float>(g);
  float cc =  all_clustering_coefficients(g,cm);

// Print the clustering coefficient of each vertex.
  graph_traits<cGraph>::vertex_iterator i, end;
  for(boost::tie(i, end) = vertices(g); i != end; ++i) {
    std::cout << std::setw(12) << setiosflags(std::ios::left)
              << g[*i].pdbcode << get(cm, *i) << std::endl;
  }
  std::cout << "mean clustering coefficient: " << cc << std::endl;
*/

  // Access vertices by visitor
  //Vertex s = boost::vertex(g);
  //vector<typename vertex_bundle_type<cGraph>::type> s = boost::get(&VertexProperties::name,g);

  //std::cout<"The second vertex: " << g[s].name <<std::endl;
  //breadth_first_search(
  //    g, s,
  //  visitor(make_bfs_visitor(boost::make_list(VertexProperties(names))))
  //);
  if (!CmdArg::ographname.empty()){
    // Create dynamic_properties
    boost::dynamic_properties dp;
    dp.property("id",boost::get(&VertexProperties::name,g));
    dp.property("weight",boost::get(&EdgeProperties::rmsd,g));
    dp.property("label",boost::get(&EdgeProperties::rmsd,g));
    //output the graph into dot format of Graphviz
    write2file<cGraph>(g,CmdArg::ographname, dp, std::string("id"));
    /* Not work yet.
    write2file(g,
      CmdArg::ofilename,
      PropertiesWrapper<VertexProperties>(),
      PropertiesWrapper<EdgeProperties>()
    );
    */
  }
  //print_dp(g, dp);

  //std::vector<Vertex> found =
  //    find_vertex_in_graph(g,&VertexProperties::chainId,"B");
  /*find_vertex_in_graph<
        cGraph
        ,std::string
        //,VertexProperties
    >(g,&VertexProperties::name,"13pkG");
  */
  /*
  if ( found.size() != 0 ){
    for (auto &x : found){
      std::cout<<"Found:"<< g[x].name << std::endl;
    }
  }else{
    std::cout<<"Not any chainId is B found in the graph\n";
  }
  */
  //kruskal minimum spanning tree(not work for disconnected components)
  /*std::vector < Edge > spanning_tree;
   boost::kruskal_minimum_spanning_tree(
        g,
        std::back_inserter(spanning_tree),
        weight_map(boost::get(&EdgeProperties::rmsd,g)));

  for (std::vector < Edge >::iterator ei = spanning_tree.begin();
       ei != spanning_tree.end(); ++ei) {
    Vertex s = boost::source(*ei, g);
    Vertex t = boost::target(*ei, g);
    std::cout << g[s].name << " <--> " << g[t].name << std::endl;
  }
  //prim minimum spanning tree
  std::vector < Vertex > prim_tree(num_vertices(g));
  boost::prim_minimum_spanning_tree(
        g,
        &prim_tree[0],
        weight_map(boost::get(&EdgeProperties::rmsd,g)));
  for (std::size_t i = 0; i != prim_tree.size(); ++i)
    if (prim_tree[i] != i)
      //std::cout << "parent[" << g[prim_tree[i]].name << "] = " << prim_tree[i] << std::endl;
      std::cout << "parent[" << g[i].name << "] = " << g[prim_tree[i]].name << std::endl;
    else
      std::cout << "parent[" << g[i].name << "] = no parent" << std::endl;
  */
  if ( !CmdArg::ofilename.empty() ) {
    auto o = clustering_by_edge(g,&VertexProperties::name);
    save2file(o,CmdArg::ofilename);
  }

  if(CmdArg::verbose){
    printf("Elasped time [CPU]: %f[sec]\n", timer.elapsed());
    std::cout<<"Program finished."<<std::endl;
  }
  return 0;
}

