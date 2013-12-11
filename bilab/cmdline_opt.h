//
//  cmdline_opt.h
//  graphClustering
//
//  Created by 曹 巍 on 12/06/05.
//  Copyright (c) 2012年 生物情報工学研究室・東京大学. All rights reserved.
//

//1. build settings/Linking/other Linker flags: -lboost_program_options
//2. build settings/search path/Library search path: /Users/davecao/Downloads/boost/stage/lib
// on Linux:
//  all:
//     $(CC) $(INC) $(SOURCE) -o $(OUTPUT) $(LIB) -lboost_program_options -dynamic
//namespace po = boost::program_options;


#ifndef graphClustering_cmdline_opt_h
#define graphClustering_cmdline_opt_h

#include "common.h"

#include <sys/stat.h>


namespace po = boost::program_options;

namespace VersionInfo{
  std::string version("graphClustering Version 1.0 developped by Wei Cao.");
}

namespace CmdArg {

  int numOfvertex = 0;
  int graphType = 1;
  std::string filename;
  std::string ographname;// = "out.dot";
  std::string ofilename = "out.clust";
  double rminCutoff = 0.0;
  double rmaxCutoff = std::numeric_limits<double>::infinity();
  double sim1Cutoff = 0.0;
  double sim2Cutoff = 0.0;
  int optLenCutoff = 35;
  bool verbose = false;
  std::string inputfile_s; // store the size of input file.

  // all the options
  po::options_description general("General options");
  po::options_description opt("Graph options");

  // Function: Initialize opts objects
  //
  void init(){
    general.add_options()
      ("version,V", "Show the version number")
      ("num,n", po::value<int>(),"Specify the total number of vertices. Default is 0.")
      ("type,t",po::value<int>(),"Graph types: 0->directed(not implemented yet) or 1->undirected. Default is 1.")
      ("file,f",po::value<std::string>(),"The input data file.")
      ("out,o", po::value<std::string>(),"The output file of clusters in text format.Default name is 'out.clust'.")
      ("graph,g", po::value<std::string>(),"The output file in the dot format of Graphviz.")
      ("help,h", "print help info.")
      ;
    opt.add_options()
      ("verbose,v","The extra verbose.")
      ("rmin",po::value<double>(),"The minimum rmsd for the cutoff between two pdb structures. For rmsd values of two structures, i.e., rmin<= r <= rmax, the edge will be created for the graph. rmin=0.0.")
      ("rmax",po::value<double>(),"The maximium rmsd for the cutoff between two pdb structures. rmax = +inf.")
      ("sim1,a",po::value<double>(),"The similarity cutoff for the structure A in the pairwised structure alignment. Together with rmsd cutoff. Default is zero.")
      ("sim2,b",po::value<double>(),"The similarity cutoff for the structure B in the pairwised structure alignment. Used with rmsd cutoff together. Default is zero.")
      ("optlen,l",po::value<int>(),"The optimal length cutoff for structure comparison.Default is 35 aa.")
      ;
    //opt.add(general);
    general.add(opt);
  }

  template<typename OutputStream>
  void printConditions(OutputStream& o) {
    o<<"############################################"<<std::endl;
    o<<"# Conditions of clustering data            "<<std::endl;
    o<<"# 1. RMSD: ["<<rminCutoff<<","<<rmaxCutoff<<"]         "<<std::endl;
    o<<"# 2. Similarity1: >="<<sim1Cutoff<<std::endl;
    o<<"# 3. Similarity2: >="<<sim2Cutoff<<std::endl;
    o<<"# 4. optLen: >="<<optLenCutoff<<std::endl;
    o<<"# Row specification:"<<std::endl;
    o<<"#Cluster Number[:tab:]pdbid+chainId[:space:]..."<<std::endl;
    o<<"###########################################"<<std::endl;
  }

  // Function: Check the existance of the file
  //
  bool FileExists(std::string fname) {
    struct stat FileInfo;
    int val;
    val = stat(fname.c_str(), &FileInfo);
    if(FileInfo.st_size >1000000){
      inputfile_s = byteConverter_s(FileInfo.st_size);
      std::cout<<fname<<" is so large ("<<inputfile_s<<")"<<". Take a while to read...";
    }
    return (val == 0) ? true : false;
  }

  // Function: Parse cmdline arguments
  //
  void ParseCmdLine(int argc, const char * argv[]) {
    // initialize options
    init();
    //parse command lines
    po::variables_map vm;
    try{
      po::store(parse_command_line(argc, argv, general),vm);
      po::notify(vm);
    }catch(po::ambiguous_option& e){
      //std::cout<<"Program is terminated for the following reason(s):\n\n";
      std::cout<<e.get_option_name()<<" is an ambiguous option name.\n";
      std::cout << general << std::endl;
      exit(0);
      //BOOST_CHECK_EQUAL(std::string(e.what()), "Unknown option");
    }catch(po::unknown_option& e){
      std::cout<<e.get_option_name()<<" is an unknown option name.\n";
      std::cout << general << std::endl;
      exit(0);
      //BOOST_CHECK_EQUAL(e.get_option_name(), "option name");
      //BOOST_CHECK_EQUAL(std::string(e.what()), "Unknown option");
    }

    if( vm.count("help") ){
      // print help info
      std::cout << general << std::endl;
      exit(0);
    }

    if( vm.count("version") ){
      //print version info
      std::cout << VersionInfo::version << std::endl;
      exit(0);
    }

    if ( vm.count("verbose") ){
      verbose = true;
    }

    if ( vm.count("num") ) {
      numOfvertex = vm["num"].as<int>();
    }

    if ( vm.count("type") ) {
      graphType = vm["type"].as<int>();
    }

    if ( vm.count("file") ) {
      filename = vm["file"].as<std::string>();
      if( !FileExists(filename) ) {
        std::cout << "Could not find the file: "<< filename <<std::endl;
        std::cout << "Program terminated." << std::endl;
        exit(0);
      }
    }

    if ( vm.count("out") ) {
      ofilename = vm["out"].as<std::string>();
    }

    if ( vm.count("graph") ) {
      ographname = vm["graph"].as<std::string>();
    }

    if ( vm.count("rmin") ) {
      rminCutoff = vm["rmin"].as<double>();
    }

    if ( vm.count("rmax") ) {
      rmaxCutoff = vm["rmax"].as<double>();
    }

    if ( vm.count("sim1") ) {
      sim1Cutoff = vm["sim1"].as<double>();
    }

    if ( vm.count("sim2") ) {
      sim2Cutoff = vm["sim2"].as<double>();
    }
    
    if ( vm.count("optlen") ) {
      optLenCutoff = vm["optlen"].as<int>();
    }
    /*{
      std::cout<<"Unknown command line option. See help.\n";
      // print help info
      std::cout << general << std::endl;
      exit(0);
    }
    */
  }
} // end of namespace CmdArg

#endif
