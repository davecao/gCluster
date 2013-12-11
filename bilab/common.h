//
//  common.h
//  graphClustering
//
//  Created by 曹 巍 on 12/06/05.
//  Copyright (c) 2012年 生物情報工学研究室・東京大学. All rights reserved.
//

#ifndef graphClustering_common_h
#define graphClustering_common_h

// all the system #include's we'll ever need
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <map>
#include <stack>
#include <queue>
#include <algorithm>
#include <iterator>
#include <limits>   // numeric_limits
#include <typeinfo> // typeid 
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string>
#include <ctype.h>	
#include <assert.h>
//#include <regex> // Not supported by clang++ with llvm 4.2.1

//Boost
#include <boost/config.hpp> // put this first to suppress some VC++ warnings
#include <boost/program_options.hpp>
#include <boost/thread.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>
#include <boost/format.hpp>
#include <boost/bind.hpp>
#include <boost/pending/stringtok.hpp>
#include <boost/regex.hpp>
#include <boost/any.hpp>
//#include <boost/xpressive/xpressive.hpp>

//Boost Graph Libraries
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/clustering_coefficient.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/bc_clustering.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

inline std::string byteConverter(long long byte)
{
  std::stringstream ss;
  if (byte < 1024){
    ss << byte <<" Bytes";
  }else if (byte <1048576 && byte >= 1024){
    ss << static_cast<float>(byte/1024) <<" KB";
  }else if (byte < 1073741824 && byte >= 1048576) {
    ss << static_cast<float>(byte/1048576) <<" MB";
  }else{
    ss << static_cast<float>(byte/1073741824)<<" GB";
  }
  return ss.str();
}

inline std::string byteConverter_s(long long byte)
{
  std::stringstream ss;
  if (byte < 1000){
    ss << byte <<" Bytes";
  }else if (byte <1000000 && byte >= 1000){
    ss << static_cast<float>(byte)/1000 <<" KB";
  }else if (byte < 1000000000 && byte >= 1000000) {
    ss << static_cast<float>(byte)/1000000 <<" MB";
  }else{
    ss << static_cast<float>(byte)/1000000000<<" GB";
  }
  return ss.str();
}

#endif
