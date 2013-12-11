//
//  cGraph.h
//  graphClustering
//
//  Created by 曹 巍 on 12/06/18.
//  Copyright (c) 2012年 生物情報工学研究室・東京大学. All rights reserved.
//

#ifndef graphClustering_cGraph_h
#define graphClustering_cGraph_h

#include "common.h"

//Define a graph using bundle properties
struct VertexProperties : public base_visitor<VertexProperties>
{
  std::string pdbcode;
  std::string chainId;
  std::string name;
  int length;

  VertexProperties() = default;
  VertexProperties (std::string pcode,std::string chId, std::string n) :
          pdbcode(pcode), chainId(chId),name(n)
  {};

  void operator()(){
    std::cout<< pdbcode + chainId <<"[shape=circle color=red]"<<std::endl;
  }

  void operator()(std::ostream& out) const {
    out<< pdbcode + chainId <<"[shape=circle color=red]"<<std::endl;
  }

  /*  template <class Vertex, class Graph>
  inline void operator()(Vertex u, Graph& g){
    std::cout << std::endl << "arriving at " << name[u] << std::endl
              << "  neighbore are: ";
  }*/
};


struct EdgeProperties
//: public base_visitor<EdgeProperties>
{
  std::string method;
  int optLength;
  double optChainRmsd;
  double identities_struct_align;
  double pval;
  double rmsd;
  double similarity_s;
  double similiarity_t;

  void operator()() {
    std::cout<<"[weight="<<rmsd<<"]"<<std::endl;
  }

  void operator()(std::ostream& out) const {
    out<<"[weight="<<rmsd<<"]"<<std::endl;
  }

  template <class Edge>
  void operator()(std::ostream& out, const Edge& v) const {
    out << "[label=\"" << rmsd << "\"]";
  }
};

template <typename Properties>
struct PropertiesWrapper
{
  typedef decltype(Properties()()) result_type;

  void operator()(){
    Properties()();
  }
};

template <typename Graph,typename ReturnType = float>
inline typename exterior_vertex_property<Graph, ReturnType>::map_type
calculate_clustering_coefficient(Graph &g)
{
  typedef exterior_vertex_property<Graph, ReturnType> ClusteringProperty;
  typedef typename ClusteringProperty::container_type ClusteringContainer;
  typedef typename ClusteringProperty::map_type ClusteringMap;
  ClusteringContainer coefs(num_vertices(g));
  ClusteringMap cm(coefs, g);
  //float cc = all_clustering_coefficients(g, cm);
  return cm;
}


template <typename DistanceMap>
class bacon_number_recorder : public default_bfs_visitor
{
  public:
  bacon_number_recorder(DistanceMap dist) : d(dist) { }

  template <typename Edge, typename Graph>
  void tree_edge(Edge e, const Graph& g) const {
    typename graph_traits<Graph>::vertex_descriptor u, v;
    u= source(e, g);
    v = target(e, g);
    d[v] = d[u] + 1;
  }
  private:
  DistanceMap d;
};

// Convenience function
template <typename DistanceMap>
bacon_number_recorder<DistanceMap> record_bacon_number(DistanceMap d)
{
  return bacon_number_recorder<DistanceMap>(d);
}

// Iterators over vertices
//boost::graph_traits < cGraph >::vertex_iterator vi, vend;
//for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi) {
//std::cout << "name " << g[*vi].name << std::endl;
//  std::cout << g[*vi].pdbcode <<":"
//            << num_paths_through_vertex(g, *vi) << ":"
//            << num_triangles_on_vertex(g, *vi)   << ":"
//            << clustering_coefficient(g, *vi)<<std::endl;
//num_triangles_on_vertex(g,vi);
//clustering_coefficient(g,vi);
// Store into cm
// cm[vi] =
//}
//Iterators over edges
//boost::graph_traits < cGraph >::edge_iterator ei, eend;
//for (boost::tie(ei, eend) = edges(g); ei != eend; ++ei) {
//    std::cout << g[source(*ei,g)].name << ":" <<g[target(*ei,g)].name << " weight:"<< g[*ei].weight<<std::endl;
//}
template<
  typename Graph,
  typename VertexAttrType,
  typename T,
  typename Bundle
//typename pVertex::* //= boost::vertex_bundle_type<Graph>::type::*
>
typename graph_traits<Graph>::vertex_descriptor
find_first_vertex_in_graph(
    Graph& g,
    //VertexAttrType VertexProperties::* pToVertexProp,
    //VertexAttrType graph_traits<Graph>::vertex_descriptor::* pToVertexProp,
    //typename vertex_bundle_type<Graph>::type::* pToVertexProp,
    T Bundle::* pToVertexProp,
    VertexAttrType vprop){

  //  typename vertex_bundle_type<Graph>::type vBundler;
  //  typename vertex_bundle_type<Graph>::type::* callee;
  //  typename vertex_bundle_type<Graph>::type vbunlder;

  typename boost::graph_traits < Graph >::vertex_iterator vi, vend;
  for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi) {
    if (g[*vi].*pToVertexProp == vprop ){
      return *vi;
    }
  }
  return -1;
}

template<
  typename Graph,
  typename VertexAttrType,
  typename T,
  typename Bundle>
std::vector<typename graph_traits<Graph>::vertex_descriptor>
  find_vertex_in_graph(
                     Graph& g,
                     T Bundle::* pToVertexProp,
                     VertexAttrType vprop){

  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  std::vector<Vertex> r;

  typename boost::graph_traits < Graph >::vertex_iterator vi, vend;

  for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi) {

    if (g[*vi].*pToVertexProp == vprop ){
      if (CmdArg::verbose) {
        std::cout<< g[*vi].pdbcode <<":"<<g[*vi].chainId <<":"<<g[*vi].name<<":";
        std::cout<< g[*vi].length <<"-->"<<g[*vi].*pToVertexProp<<std::endl;
      }
      r.push_back(*vi);
    }
  }
  return r;
}

template<typename Graph>
void print_dp(Graph& g, boost::dynamic_properties& dp)
{
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::edge_descriptor Edge;

  typedef typename graph_traits<Graph>::vertex_iterator v_iter;
  typedef typename graph_traits<Graph>::edge_iterator e_iter;

  v_iter v, v_end;
  e_iter e, e_end;

  //typename graph_traits<Graph>::edges_size_type edge_count = 0;

  //typedef property_map<Graph, >::type
  //std::vector<std::string> s = boost::get(&VertexProperties::name,g);
  for (tie(v, v_end) = vertices(g); v != v_end; ++v) {
    // Output data
    for (boost::dynamic_properties::const_iterator i = dp.begin(); 
        i != dp.end(); ++i) {
      //dp.property("id",boost::get(&VertexProperties::name,g));
      // i->first: "id"

      if (i->second->key() == typeid(Vertex)) {
        std::cout << "\""<<i->second->get_string(*v) << "\";\n";
      }
    }

  }
  Vertex s,t;
  for (tie(e, e_end) = edges(g); e != e_end; ++e) {
    s = boost::source(*e,g);
    t = boost::target(*e,g);
    std::cout <<"\""<<g[s].name <<"\""<< "---" <<"\""<<g[t].name<<"\""<< " [";
    for (boost::dynamic_properties::const_iterator i = dp.begin(); 
          i != dp.end(); ++i) {
      if (i->second->key() == typeid(Edge)) {
        std::cout <<i->first<<"="<<i->second->get_string(*e);
        if( boost::next(i) != dp.end() ){
          std::cout<<",";
        }
      }
    }
    std::cout<<"];\n";
  }
}

template<
  typename Graph,
  typename T,
  typename Bundle
>
std::multimap<int,T>
clustering_by_edge(Graph& g,T Bundle::* pToVertexProp)
{
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::edge_descriptor Edge;

  typedef typename graph_traits<Graph>::vertex_iterator v_iter;
  typedef typename graph_traits<Graph>::edge_iterator e_iter;

  typedef typename graph_traits<Graph>::adjacency_iterator adj_iter;

  typedef typename std::multimap<int,T> Clusters;
  //typename Clusters::iterator pos;

  std::map<Vertex,int> visited;
  typename std::map<Vertex,int>::iterator visited_iter;
  std::queue<Vertex> vertexQueue;

  Clusters c;

  bool inserted;
  int index=1;

  adj_iter neighbourIt, neighbourEnd;
  v_iter v, v_end;
  Vertex curr;

  for (boost::tie(v, v_end) = vertices(g); v != v_end; ++v) {
    int degree_ = boost::degree(*v,g);
    //std::cout<<g[*v].*pToVertexProp <<": "<< degree_ <<std::endl;
    if (degree_ != 0){
      boost::tie(visited_iter,inserted) = visited.insert(std::make_pair(*v,1));
      // if the current vertex is not visited
      if (inserted) {
        //std::cout<<"start from : "<<g[*v].*pToVertexProp <<"("<< degree_ <<")"<<std::endl;
        //std::cout<<"--";
        c.insert(std::make_pair(index,g[*v].*pToVertexProp));
        vertexQueue.push(*v);
        while (!vertexQueue.empty()){
          // bread first search from the current vertex
          curr = vertexQueue.front();
          //std::cout<<g[curr].*pToVertexProp<<"-->(";
          vertexQueue.pop();
          if (boost::degree(curr,g) >= 1){// has connected components.
            // get adjacent vertices
            boost::tie(neighbourIt, neighbourEnd) = boost::adjacent_vertices(curr, g);
            for(; neighbourIt != neighbourEnd; ++neighbourIt ){
              boost::tie(visited_iter,inserted) = visited.insert(std::make_pair(*neighbourIt,1));
              if (inserted){
                vertexQueue.push(*neighbourIt);
                c.insert(std::make_pair(index,g[*neighbourIt].*pToVertexProp));
                //std::cout<<g[*neighbourIt].*pToVertexProp <<",";
              }
            }
            //std::cout<<")"<<std::endl;
          }else{
            //std::cout<<"END)"<<std::endl;
          }
        }
        index++;
        //std::cout<<std::endl;
      } // End of if inserted
    }else{
      //a disconnected vertex
      //std::cout<<"disconnected:"<<g[*v].*pToVertexProp <<std::endl;
      c.insert(std::make_pair(index,g[*v].*pToVertexProp));
      index++;
    }
  }
  /*
  //multimap getKeys
  std::set<typename Clusters::key_type> uKeys;
  getUniqkeys(c,uKeys);
  for (auto pos = uKeys.begin(); pos != uKeys.end(); ++pos) {
  //  for (auto pos = c.begin(); pos != c.end(); ++pos){
        std::cout << *pos <<"\t";
        std::pair<typename Clusters::iterator,typename Clusters::iterator> range = c.equal_range(*pos);
        for(auto c_iter = range.first; c_iter!=range.second; ++c_iter){
          std::cout<< c_iter->second <<" ";
        }
      std::cout<<std::endl;
      }
    //print_graph(g);
  */
  return c;
}

template <typename multimap_t>
void save2file(multimap_t &c, std::string& out)
{
  std::ofstream file(out.c_str());
  if (!file) {
    std::cout <<"Could not open the file "<<out <<"."<<std::endl;
    exit(0);
  }

  std::cout<< "Ready to save clusters to the file "<<out;
  CmdArg::printConditions(file);
  //multimap getKeys
  std::set<typename multimap_t::key_type> uKeys;
  getUniqkeys(c,uKeys);
  for (auto pos = uKeys.begin(); pos != uKeys.end(); ++pos) {
    //for (auto pos = c.begin(); pos != c.end(); ++pos){
    file << *pos <<"\t";
    std::pair<typename multimap_t::iterator,typename multimap_t::iterator> range = c.equal_range(*pos);
    for(auto c_iter = range.first; c_iter!=range.second; ++c_iter){
      file << c_iter->second <<" ";
    }
    file<<std::endl;
  }
  file.close();
  std::cout<< " ... finished."<<std::endl;
}

// Define a graph
typedef boost::adjacency_list<
  vecS,  // adjacency list container
  vecS,  // adjacency list container
  undirectedS, // graph type
  VertexProperties, // Vertex property
  EdgeProperties,   // Edge property
  boost::no_property // Graph property: default is no_property
> cGraph;

// Define vertices
typedef boost::graph_traits<cGraph> Traits;
typedef Traits::vertex_descriptor Vertex;
typedef Traits::edge_descriptor Edge;


#endif
