#include "Tree.hpp"
#include "GeneTree.hpp"
#include "SpeciesTree.hpp"
#include <pugixml.hpp>
#include "ReconMap.hpp"
#include <iostream>
#include <fstream>
#include <tclap/CmdLine.h>

struct Dummy {};
int main(int argc, const char *argv[]) {
	std::string input_path, output_path;
	bool write_dot;
	try {
		using namespace TCLAP;
		CmdLine cmd("Command description", ' ', "0.9");
		ValueArg<std::string> inputArg("i", "input", "File input path", true, "", "string", cmd);
		ValueArg<std::string> outputArg("o", "output", "File output path", false, "", "string", cmd);
		SwitchArg dotArg("d", "dot", "write dot file", cmd, false);
		cmd.parse(argc, argv);

		input_path = inputArg.getValue();
		output_path = outputArg.getValue();
		write_dot = dotArg.getValue();

	} catch (TCLAP::ArgException &e) {  // catch any exceptions  
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
		return 0;
	}
	/* tc::GeneTree gt(tc::Event::Speciation); */
	tc::SpeciesTree st;
	tc::GeneTree gt;
	/* std::cout << "Reading Gene Tree" << std::endl; */
	/* pugi::xml_document gdoc; */
	/* auto result = gdoc.load_file("../data/gt2.xml"); */
	/* if (result) { */
	/* 	try { */
	/* 		gt.read_xml(gdoc); */
	/* 	} catch(std::invalid_argument e) { */
	/* 		std::cout << e.what() << std::endl; */
	/* 		return 0; */
	/* 	} */
	/* 	/1* gt.print(); *1/ */
	/* } else { */
	/* 	std::cout << "Error: " << result.description() << std::endl; */
	/* 	std::cout << "Offset: " << result.offset << std::endl; */
	/* 	return 0; */
	/* } */
	std::ifstream ft(input_path);
	std::string str_gt, str_st;
	if ((std::getline(ft, str_st))) {
		try {
			st.read_newick(str_st);
		} catch (std::invalid_argument e) {
			std::cerr << "error: " << e.what() << std::endl;
			return 0;
		}
	} else {
		std::cerr << "missing species tree\n";
		return 0;
	}
	if ((std::getline(ft, str_gt))) {
		try {
			gt.read_newick(str_gt);
		} catch (std::invalid_argument e) {
			std::cerr << "error: " << e.what() << std::endl;
			return 0;
		}
	} else {
		std::cerr << "missing gene tree\n";
		return 0;
	}
	/* std::getline(fst, str_st); */

	/* gt.read_newick(str_gt); */
	/* gt.print(); */
	/* st.read_newick(str_st); */
	/* st.print(); */

	/* std::cout << "Reading Species Tree" << std::endl; */
	/* pugi::xml_document sdoc; */
	/* auto result = sdoc.load_file("../data/st2.xml"); */
	/* if (result) { */
	/* 	try { */
	/* 		st.read_xml(sdoc); */
	/* 	} catch(std::invalid_argument e) { */
	/* 		std::cout << e.what() << std::endl; */
	/* 		return 0; */
	/* 	} */
	/* 	/1* st.print(); *1/ */
	/* } else { */
	/* 	std::cout << "Error: " << result.description() << std::endl; */
	/* 	std::cout << "Offset: " << result.offset << std::endl; */
	/* 	return 0; */
	/* } */
	tc::ReconMap rcm(gt, st);
	/* rcm.print(); */
	auto res = rcm.make_time_consistent();
	if (res.first) {
		if (output_path != "") {
			std::ofstream ft(output_path);
			if (write_dot) {
				rcm.write_dot(ft);
			} else {
				rcm.write_nexus(ft);
			}
		} else {
			if (write_dot) {
				rcm.write_dot(std::cout);
			} else {
				rcm.write_nexus(std::cout);
			}
		}
	} else {
		std::cout << "No Time Consistent Reconciliation map exists.\n";
	}
	
	return 0;
}
