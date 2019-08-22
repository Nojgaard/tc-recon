#ifndef TC_RECON_MAP_HPP
#define TC_RECON_MAP_HPP

#include "GeneTree.hpp"
#include "SpeciesTree.hpp"
#include "TimeMap.hpp"
#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include <vector>

namespace tc {

class ReconMap {
public:
	ReconMap(const GeneTree& gt, const SpeciesTree& st);
	std::pair<bool, TimeMap> make_time_consistent();
	void print();
	void write_dot(std::ostream& os = std::cout);
	void write_nexus(std::ostream& os = std::cout);
	void write_aux_graph(std::ostream& os);
private:
	using AuxGraph = boost::adjacency_list<boost::vecS,boost::vecS,boost::directedS,size_t>;
	std::pair<std::vector<size_t>, std::vector<size_t>> build_time_graph(AuxGraph& g);

	SpeciesTree::Node& g2s(GeneTree::Node gn);
	SpeciesTree::Node& lca_g2s(GeneTree::Node gn);
	void build_leaf_map();
	void build_lca_map(GeneTree::Node n);


	void write_nexus_taxa(std::ostream& os);
	void write_nexus_trees(std::ostream& os);
	void write_nexus_recon(std::ostream& os);


	const GeneTree& _gt;
	const SpeciesTree& _st;
	std::vector<SpeciesTree::Node> _m, _lca_g2s;
};
	
} /* tc */ 

#endif /* ifndef TC_RECON_MAP_HPP */
