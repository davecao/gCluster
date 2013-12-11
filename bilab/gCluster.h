//
//  gCluster.h
//  graphClustering
//
//  Created by 曹 巍 on 12/06/05.
//  Copyright (c) 2012年 生物情報工学研究室・東京大学. All rights reserved.
//

#ifndef graphClustering_gCluster_h
#define graphClustering_gCluster_h

#include "common.h"

using namespace boost;

namespace bilabGraph {

/*************************************************************************
 * Class templates for graph:
 * // Adjacency Matrix Class
 *  template<
 *     typename Directed = directedS,
 *     typename VertexProperty = no_property,
 *     typename EdgeProperty = no_property,
 *     typename GraphProperty = no_property,
 *     typename Allocator = std::allocator<bool> 
 *   >
 *   class adjacency_matrix
 *
 * // Adjacency list Class
 *  template<
 *     class OutEdgeListS = vecS,// a Sequence or an AssociativeContainer
 *     class VertexListS = vecS, // a Sequence or a RandomAccessContainer
 *     class DirectedS = directedS,
 *     class VertexProperty = no_property,
 *     class EdgeProperty = no_property,
 *     class GraphProperty = no_property,
 *     class EdgeListS = listS
 *   >
 *   class adjacency_list
 ******************************************************************************/
//    template <typename Graph>
//    struct vertex_vector
//    {
//        typedef graph_traits<Graph> traits;
//        typedef std::vector<typename traits::vertex_descriptor> type;
//    };
//
//    template <typename Graph>
//    void build_graph(Graph& g, typename vertex_vector<Graph>::type& v, int N)
//    {
//        typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
//
//        // add vertices
//        for(size_t i = 0; i < N; ++i) {
//            v[i] = add_vertex(g);
//        }
//
//        // add edges
//        add_edge(v[0], v[1], g);
//        add_edge(v[1], v[2], g);
//        add_edge(v[2], v[0], g);
//        add_edge(v[3], v[4], g);
//        add_edge(v[4], v[0], g);
//    }

  template<
    typename VertexProperty,
    typename EdgeProperty,
    typename Graph,
    typename ClusteringProperty = boost::exterior_vertex_property<Graph,double>
  >
  class gCluster
  {
    // Define Graph vertex
    typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename graph_traits<Graph>::edge_descriptor   Edge;
    typedef typename vertex_bundle_type<Graph>::type vBundler;
    typedef typename edge_bundle_type<Graph>::type   eBundler;
    //typedef std::map<std::string, Vertex> NameVertexMap;

    typedef typename property_map<Graph, VertexProperty>::type VertexProp;
    typedef typename property_map<Graph, EdgeProperty>::type EdgeProp;

    typedef typename graph_traits<Graph>::vertex_iterator vertex_iter;

    //typedef exterior_vertex_property<Graph, double> ClusteringProperty;
    typedef typename ClusteringProperty::container_type ClusteringContainer;
    typedef typename ClusteringProperty::map_type ClusteringMap;

    private:

      Graph g;
      std::vector<Vertex> v;
      unsigned            numOfvertices;
      VertexProp          vertex_prop;
      EdgeProp            edge_prop;

      ClusteringProperty  cp;
      ClusteringMap       cmp;
      ClusteringContainer ccontainer;

    public:

      gCluster() = default;
      ~gCluster() = default;
      /*Vertex createVertex(){
            return add_vertex(g);
        }

        Edge createEdge(Vertex u, Vertex v){
            return add_edge(u,v,g);
        }
        */
        /*
        template<class T>
        void addVertex(T node)
        {
            Vertex u = add_vertex(g);
            boost::put(node_name,u,node);
            //boost::put(node_name,v);
            //boost::tie(e) = add_edge(u,v,g);
        }

        template<class V>
        void addVertex(V &node) {
            Vertex u;

            //add vertices
            u = add_vertex(g);
            boost::put(vertex_prop,u, node);
        }

        template<class S, class T, class E>
        void addEdge (S &s, T &t, E &e) {
            bool is_added;
            Edge edge;

            add_edge(s,t,g);
        }
         */
        //template<class U = typename VertexProperties, class E = EdgeProperties >
        void add (VertexProperty node1, VertexProperty node2, EdgeProperty edge){

            //addVertex(node1);
            //addVertex(node2);
            //addEdge(node1,node2,edge);

            //bool inserted;
            Vertex u,v;
            Edge e;

            //add vertices
            u = boost::add_vertex(g);
            //boost::put(vertex_prop,u, node1);

            v = boost::add_vertex(g);
            //boost::put(vertex_prop,v,node2);

            //add edge
            //edge_property ep(e);
            //ep.m_value = weight;

            add_edge(node1,node2,edge,g);
            //boost::put(edge_prop,e,edge);

        }

        void ToGraphviz(std::string outfile)
        {
            //std::string name;
            // Access vertices

            //
            // method 1: using pair
            //
            //std::pair<VertexIter,VertexIter> vp;

            //obtain node names
            //for ( vp=vertices(g)  ; vp.first != vp.second; ++vp.first )
            //{
            //    std::cout<< boost::get(vertex_name_t(),g, *vp.first) << std::endl;
            //}

            //
            // method 2: using iterators
            /*VertexIter vi, vi_end,next;
             boost::tie(vi,vi_end) = vertices(g);
             for ( next = vi; vi != vi_end; vi = next ){
             ++next;
             std::cout<< boost::get(vertex_name_t(),g,*vi) << std::endl;
             }
             */
            boost::dynamic_properties dp;
            dp.property("id",boost::get(vertex_name_t(),g));
            dp.property("weight",boost::get(edge_name_t(),g));

            // output ti graphviz format
            std::ofstream file(outfile.c_str());
            if (!file) {
                std::cout <<"Could not open the file "<< outfile <<"."<<std::endl;
                exit(0);
            }
            //boost::write_graphviz(file,g,boost::make_label_writer(name.c_str());
            boost::write_graphviz_dp(file,g,dp,std::string("id"));
            file.close();
        }
    };

}


#endif
