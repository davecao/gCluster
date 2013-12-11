//
//  utility.h
//  graphClustering
//
//  Created by 曹 巍 on 12/06/16.
//  Copyright (c) 2012年 生物情報工学研究室・東京大学. All rights reserved.
//

#ifndef graphClustering_utility_h
#define graphClustering_utility_h

#include "common.h"

bool is_comment(std::string line)
{
  //std::regex txt_gex("//.*|#.*");
  //return std::regex_match(line,txt_gex) ? true :false;
  boost::regex txt_gex("//.*|#.*");
  return boost::regex_match(line,txt_gex) ? true :false;
}

void trim(std::string& str)
{
  std::string::size_type pos = str.find_last_not_of(' ');
  if (pos != std::string::npos){
    str.erase(pos+1);
    pos = str.find_first_not_of(' ');
    if (pos!=std::string::npos)	{
      str.erase(0,pos);
    }
  }else{
    str.erase(str.begin(),str.end());
  }
}

template<typename multimap_t>
void getUniqkeys(const multimap_t &m, std::set<typename multimap_t::key_type> &k)
{
  //typename std::set<typename multimap_t::key_type>::iterator keys_end=k.end();
  auto keys_end = k.end();
  for (auto b(m.begin()),e(m.end()); b!=e; ++b){
    k.insert(keys_end,b->first);
  }
}


template <typename Graph, typename InputStream>
//inline std::map<std:string, typename boost::graph_traits<Graph>::vertex_descriptor>
void read_graph2(Graph &g, InputStream& is)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
  typedef typename std::map<std::string, Vertex> NameVertexMap;
  NameVertexMap vName;
  typename NameVertexMap::iterator pos;

  std::string line;

  while (std::getline(is, line)) {
    // Ignore empty line
    if ( line.empty() ) continue;
      // Ignore comment lines started with // or #
      if ( is_comment(line) ) continue;
      // trim the leading and trailing space
      trim(line);
      // Convert NaN to lower case
      boost::replace_all(line,"NaN","nan");
      boost::replace_all(line,".pdb","");

      //VertexProperties u,v;
      //EdgeProperties e;
      Vertex U,V;
      Edge E;
      bool inserted;

      // Set a line token
      std::vector<std::string> line_toks;
      boost::stringtok(line_toks,line,",");

      //process the current line
      //std::list<std::string>::iterator i = line_toks.begin();
      //std::vector<std::string>::iterator i = line_toks.begin();
      //u.name = *i++;
      //v.name = *i++;
      //e.name = *i;
      //e.weight = (*i).length();
      //for (auto i : line_toks) {
      //  std::cout<< i <<" ";
      //}
      //std::cout<<std::endl;
      std::string name=line_toks[0].append(line_toks[1]);
      boost::tie(pos, inserted) = vName.insert(std::make_pair(name,Vertex()));
      if ( inserted ) {
        U = add_vertex(g);
        g[U].pdbcode = line_toks[0];
        g[U].chainId = line_toks[1];
        g[U].name = name;//line_toks[0].append(line_toks[1]);
        g[U].length  = boost::lexical_cast<int>(line_toks[11]);
        //std::cout<<line_toks[0] <<" "<<line_toks[2]<<std::endl;
        pos->second = U;
      }else{
        U = pos->second;
      }
      name=line_toks[2].append(line_toks[3]);
      boost::tie(pos, inserted) = vName.insert(std::make_pair(name,Vertex()));
      if (inserted) {
        V = add_vertex(g);
        g[V].pdbcode = line_toks[2];
        g[V].chainId = line_toks[3];
        g[V].name = name;//line_toks[2].append(line_toks[3]);
        g[V].length  = boost::lexical_cast<int>(line_toks[12]);
        pos->second = V;
      }else{
        V = pos->second;
      }
      // Check the conditions of creating an edge
      double rmsd,sim1,sim2,identity;
      try {
        rmsd =  boost::lexical_cast<double>(line_toks[10]);
        sim1 =  boost::lexical_cast<double>(line_toks[13]);
        sim2 =  boost::lexical_cast<double>(line_toks[14]);
        identity = boost::lexical_cast<double>(line_toks[8]);
      } catch( boost::bad_lexical_cast const& ) {
        std::cout << "Error: input string was not valid" << std::endl;
      }
      // No afp found
      if ( std::isnan(identity) ) continue;
      if ( ((rmsd >= CmdArg::rminCutoff) && (rmsd <= CmdArg::rmaxCutoff))
            &&  (sim1 >= CmdArg::sim1Cutoff)
            &&	(sim2 >= CmdArg::sim2Cutoff) ) {

        boost::tie(E,inserted) = add_edge(U, V, g);
        if (inserted) {
          g[E].method = line_toks[4];
          g[E].optLength = boost::lexical_cast<int>(line_toks[6]);
          g[E].optChainRmsd = boost::lexical_cast<double>(line_toks[7]);
          //boost::algorithm::to_lower(line_toks[8]);
          g[E].identities_struct_align = identity;//boost::lexical_cast<double>(line_toks[8]);
          //printf("%s %s %s %s %s\n",line_toks[0].c_str(),
          //      line_toks[1].c_str(),line_toks[2].c_str(),line_toks[3].c_str(),
          //      line_toks[8].c_str() );
          g[E].pval = boost::lexical_cast<double>(line_toks[9]);
          g[E].rmsd = rmsd;//boost::lexical_cast<double>(line_toks[10]);
          g[E].similarity_s = sim1;//boost::lexical_cast<double>(line_toks[13]);
          g[E].similiarity_t = sim2; //boost::lexical_cast<double>(line_toks[14]);
          //std::cout <<"Identities :"<<g[E].identities_struct_align<<std::endl;
        }
        //std::cout <<u<<":"<<v<<"->"<<e<<std::endl;
      }
  }// End of While
}

template <typename Graph, typename InputStream>
//inline std::map<std:string, typename boost::graph_traits<Graph>::vertex_descriptor>
void read_graph(Graph &g, InputStream& is)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
  typedef typename std::map<std::string, Vertex> NameVertexMap;
  NameVertexMap vName;
  typename NameVertexMap::iterator pos;

  std::string line;

  while (std::getline(is, line)) {
    // Ignore empty line
    if ( line.empty() ) continue;
    // Ignore comment lines started with // or #
    if ( is_comment(line) ) continue;
    // trim the leading and trailing space
    trim(line);
    // Convert NaN to lower case
    boost::replace_all(line,"NaN","nan");
    boost::replace_all(line,".pdb","");

    //VertexProperties u,v;
    //EdgeProperties e;
    Vertex U,V;
    Edge E;
    bool inserted;

    // Set a line token
    std::vector<std::string> line_toks;
    boost::stringtok(line_toks,line," ");

    //process the current line
    //std::list<std::string>::iterator i = line_toks.begin();
    //std::vector<std::string>::iterator i = line_toks.begin();
    //u.name = *i++;
    //v.name = *i++;
    //e.name = *i;
    //e.weight = (*i).length();
    //for (auto i : line_toks) {
    //  std::cout<< i <<" ";
    //}
    //std::cout<<std::endl;
    std::string name=line_toks[0];
    boost::tie(pos, inserted) = vName.insert(std::make_pair(name,Vertex()));
    if ( inserted ) {
      U = add_vertex(g);
      g[U].pdbcode = name;
      //g[U].chainId = line_toks[1];
      g[U].name = name;//line_toks[0].append(line_toks[1]);
      g[U].length  = boost::lexical_cast<int>(line_toks[5]);
      //std::cout<<line_toks[0] <<" "<<line_toks[2]<<std::endl;
      pos->second = U;
    }else{
      U = pos->second;
    }
    name=line_toks[1];
    boost::tie(pos, inserted) = vName.insert(std::make_pair(name,Vertex()));
    if (inserted) {
      V = add_vertex(g);
      g[V].pdbcode = name;
      //g[V].chainId = line_toks[3];
      g[V].name = name;//line_toks[2].append(line_toks[3]);
      g[V].length  = boost::lexical_cast<int>(line_toks[6]);
      pos->second = V;
    }else{
      V = pos->second;
    }
    // Check the conditions of creating an edge
    double rmsd,sim1,sim2,identity;
    int optLen,alnLen;
    //int optLen,alnLen;
    try {
      rmsd =  boost::lexical_cast<double>(line_toks[2]);
      optLen = boost::lexical_cast<double>(line_toks[3]);
      alnLen = boost::lexical_cast<double>(line_toks[4]);
      sim1 =  boost::lexical_cast<double>(line_toks[7]);
      sim2 =  boost::lexical_cast<double>(line_toks[8]);
      identity = boost::lexical_cast<double>(line_toks[9])*100;

      //optLen/alnLen*100; //boost::lexical_cast<double>(line_toks[4])/boost::lexical_cast<double>(line_toks[5])*100;
    } catch( boost::bad_lexical_cast const& ) {
      std::cout << "Error: input string was not valid" << std::endl;
    }
    // No afp found
    if ( std::isnan(identity) ) continue;
    if ( ((rmsd >= CmdArg::rminCutoff) && (rmsd <= CmdArg::rmaxCutoff))
      && (optLen >=CmdArg::optLenCutoff)
      && (sim1 >= CmdArg::sim1Cutoff)
      && (sim2 >= CmdArg::sim2Cutoff)) {

      boost::tie(E,inserted) = add_edge(U, V, g);
      if (inserted) {
        //g[E].method = line_toks[4];
        g[E].optLength = optLen;//boost::lexical_cast<int>(line_toks[4]);
        //g[E].optChainRmsd = boost::lexical_cast<double>(line_toks[7]);
        //boost::algorithm::to_lower(line_toks[8]);
        g[E].identities_struct_align = identity;//boost::lexical_cast<double>(line_toks[8]);
        //printf("%s %s %s %s %s\n",line_toks[0].c_str(),
        //    line_toks[1].c_str(),line_toks[2].c_str(),line_toks[3].c_str(),
        //    line_toks[8].c_str() );
        //g[E].pval = boost::lexical_cast<double>(line_toks[9]);
        g[E].rmsd = rmsd;//boost::lexical_cast<double>(line_toks[10]);
        g[E].similarity_s = sim1;//boost::lexical_cast<double>(line_toks[13]);
        g[E].similiarity_t = sim2; //boost::lexical_cast<double>(line_toks[14]);
        //std::cout <<"Identities :"<<g[E].identities_struct_align<<std::endl;
      }
      //std::cout <<u<<":"<<v<<"->"<<e<<std::endl;
    }
  }
}


template <typename Graph>
void readfile(Graph& g,std::string data)
{
  // Create the stream for reading the data file
  std::ifstream datafile(data.c_str());
  if (!datafile) {
    std::cerr << "Could not find the file: "<< CmdArg::filename << std::endl;
    std::cout << CmdArg::general << std::endl;
    exit(0);
  }
  read_graph(g,datafile);

  // Close file
  datafile.close();
  if (!CmdArg::inputfile_s.empty()){
    std::cout<<" Fininshed."<<std::endl;
  }
}

template <typename Graph>
void write2file(Graph& g,std::string outfile){
  // output ti graphviz format
  
  std::ofstream file(outfile.c_str());
  if (!file) {
    std::cout <<"Could not open the file "<< outfile <<"."<<std::endl;
     exit(0);
  }
  if (CmdArg::verbose)
    std::cout<<"Simply create a dot file of graphviz\n"; 
  boost::write_graphviz(file,g);
  file.close();
}

template <
  typename Graph,
  typename VertexWriter, //= PropertiesWrapper<VertexProperties>,
  typename EdgeWriter // = PropertiesWrapper<EdgeProperties> 
>
void write2file(Graph& g,std::string outfile, 
                VertexWriter vprop, 
                EdgeWriter eprop){
  std::ofstream file(outfile.c_str());
  if (!file) {
    std::cout <<"Could not open the file "<< outfile <<"."<<std::endl;
    exit(0);
  }
  if (CmdArg::verbose)
    std::cout<<"Using overloaded operator of a struct to create a dot file of graphviz\n"; 
  boost::write_graphviz(file,g,vprop,eprop);
  file.close();
}

template <typename Graph>
void write2file(Graph& g,std::string outfile, 
                const boost::dynamic_properties& dp, 
                std::string const& node_id){
  // output ti graphviz format
  std::ofstream file(outfile.c_str());
  if (!file) {
    std::cout <<"Could not open the file "<< outfile <<"."<<std::endl;
    exit(0);
  }
  if (CmdArg::verbose)
    std::cout<<"Using boost::dynamic_properties to create a dot file of graphviz\n";
  boost::write_graphviz_dp(file,g, dp, node_id);
  file.close();
}

template <typename Graph>
void print_graph_info (Graph &g)
{
  //print_graph(g);
  std::cout<<"  Num of vertices : " <<num_vertices(g) << std::endl;
  std::cout<<"  Num of edges : "    <<num_edges(g) << std::endl;
}

#endif
