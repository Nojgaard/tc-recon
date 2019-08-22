#include "ReconMap.hpp"
#include "Graph.hpp"
#include <unordered_map>
#include <unordered_set>
#include <iostream>

namespace tc {

ReconMap::ReconMap(const GeneTree& gt, const SpeciesTree& st): _gt(gt), _st(st), _m(gt.num_nodes()) 
	, _lca_g2s(gt.num_nodes()) {
	build_leaf_map();
	build_lca_map(_gt.root());
	_m = _lca_g2s;
}

void ReconMap::print() {
	for (size_t i = 0; i < _m.size(); ++i) {
		std::cout << i << " -> " << g2s(i) << " (";
		auto ls = _st.leaf_set(g2s(i));
		for (auto l : ls) {
			std::cout << _st.species(l) << " ";
		}
		std::cout << ")" << std::endl;
	}
}

std::pair<std::vector<size_t>, std::vector<size_t>> ReconMap::build_time_graph(AuxGraph& g) {
	std::vector<size_t> st2a(_st.num_nodes()), gt2a(_gt.num_nodes());
	for (auto x : _st.nodes()) { st2a[x] = add_vertex(g); }
	for (auto u : _gt.nodes()) { 
		if (!_gt.is_leaf(u) && _gt.event(u) != Event::Speciation) {
			gt2a[u] = add_vertex(g);
		} else {
			gt2a[u] = st2a[g2s(u)];
		}
	}

	for (auto x : _st.nodes()) {
		for (auto y : _st.children(x)) {
			add_edge(st2a[x],st2a[y],g);
		}
	}
	for (auto u : _gt.nodes()) {
		for (auto v : _gt.children(u)) {
			add_edge(gt2a[u],gt2a[v],g);
		}
		if (!_gt.is_leaf(u) && _gt.event(u) != Event::Speciation) {
			add_edge(gt2a[u], st2a[lca_g2s(u)], g);
		}
		if (!_gt.is_leaf(u) && _gt.event(u) == Event::HGT) {
			for (auto v : _gt.transfer_nodes(u)) {
				auto x = lca_g2s(u), y = lca_g2s(v);
				add_edge(st2a[_st.lca(x, y)], gt2a[u], g);
			}
			/* auto v = (_gt.is_transfer(_gt.left(u))) ? _gt.left(u) : _gt.right(u); */
			/* auto x = lca_g2s(u), y = lca_g2s(v); */
			/* add_edge(st2a[_st.lca(x, y)], gt2a[u], g); */
		}
	}
	return std::make_pair(gt2a, st2a);
}

std::pair<bool, TimeMap> ReconMap::make_time_consistent() {
	TimeMap tm(_st, _gt);
	std::vector<size_t> gt2a, st2a;
	AuxGraph time_graph;
	std::tie(gt2a, st2a) = build_time_graph(time_graph);
	if (has_cycle(time_graph)) { return std::make_pair(false, tm); }

	auto topo_vec = topo_sort(time_graph);
	size_t time_point = 0;
	for (auto ii = topo_vec.rbegin(); ii != topo_vec.rend(); ++ii) { time_graph[*ii] = time_point++; }
	
	for (auto x : _st.nodes()) {
		tm(x, _st) = time_graph[st2a[x]];
	}
	for (auto u : _gt.nodes()) {
		tm(u, _gt) = time_graph[gt2a[u]];
	}
	for (auto u : _gt.nodes()) {
		if (_gt.is_leaf(u) || _gt.event(u) == Event::Speciation) { continue; }
		int tp = tm(_st.parent(g2s(u)), _st);
		while (tm(u,_gt) < tp) {
			g2s(u) = _st.parent(g2s(u));
			tp = tm(_st.parent(g2s(u)), _st);
		}
	}

	return std::make_pair(true, tm);
}

SpeciesTree::Node& ReconMap::g2s(GeneTree::Node gn) { return _m[gn]; }
SpeciesTree::Node& ReconMap::lca_g2s(GeneTree::Node gn) { return _lca_g2s[gn]; }

void ReconMap::build_leaf_map() {
	std::unordered_map<std::string, SpeciesTree::Node> leaf_map;
	for (auto n : _st.nodes()) {
		if (_st.is_leaf(n)) { leaf_map[_st.species(n)] = n; }
	}

	for (auto n : _gt.nodes()) {
		if (_gt.is_leaf(n)) { g2s(n) = leaf_map[_gt.species(n)]; }
	}
}

void ReconMap::build_lca_map(GeneTree::Node n) {
	if (_gt.is_leaf(n)) {
		lca_g2s(n) = g2s(n);
	} else {
		/* GeneTree::Node v = _gt.left(n), w = _gt.right(n); */
		/* build_lca_map(v); */
		/* build_lca_map(w); */
		std::vector<GeneTree::Node> tn;
		for (auto v : _gt.children(n)) {
			build_lca_map(v);
			if (!_gt.is_transfer(v)) { 
				tn.push_back(lca_g2s(v)); 
			}
		}
		lca_g2s(n) = _st.lca(tn);

		/* if (_gt.event(n) == Event::HGT) { */
		/* 	lca_g2s(n) = (_gt.is_transfer(v)) ? lca_g2s(w) : lca_g2s(v); */
		/* } else { */
		/* 	lca_g2s(n) = _st.lca(lca_g2s(v), lca_g2s(w)); */
		/* } */
	}
}

void ReconMap::write_dot(std::ostream& os) {
	os << "graph g {splines=line;" << std::endl;
	/* os << "subgraph cluster_gt {" << std::endl; */
	_gt.write_dot(os);
	/* os << "}" << std::endl; */
	/* os << "subgraph cluster_st {" << std::endl; */
	_st.write_dot(os);
	/* os << "}" << std::endl; */
	os << "edge[constraint=false,splines=curved, dir=forward, style=dashed];" << std::endl;

	for (auto u : _gt.nodes()) {
		if (!_gt.is_leaf(u) && _gt.event(u) == Event::Speciation) {
			os << "g" << u << " -- s" << g2s(u) << ";" << std::endl;
		} else if (!_gt.is_leaf(u)) {
			os << "g" << u << " -- sb" << _st.parent(g2s(u)) << g2s(u) << ";" << std::endl;
		}
	}

	os << "{ rank=same; s1; g0; }" << std::endl;
	os << "}" << std::endl;
}	

void ReconMap::write_aux_graph(std::ostream& os) {
	std::vector<size_t> gt2a, st2a;
	AuxGraph time_graph;
	std::tie(gt2a, st2a) = build_time_graph(time_graph);
	std::vector<bool> used_nodes(num_vertices(time_graph), false);
	os << "digraph g {" << std::endl;
	for (auto x : _st.nodes()) {
		/* os << st2a[x] << std::endl; */
		used_nodes[st2a[x]] = true;
		std::string shape = "circle";
		std::string label = "\"\"";
		os << st2a[x] << "[shape=" << shape << ",height=.5,width=.5,fixedsize=true,label="<<label<<"];" << std::endl;
	}
	for (auto u : _gt.nodes()) {
		if (used_nodes[gt2a[u]]) { continue; }
		used_nodes[gt2a[u]] = true;
		std::string shape = "circle";
		std::string label = "\"\"";
		if (_gt.is_leaf(u)) {
			label = "\"" + _gt.species(u) + "\"";
		}
		if (!_gt.is_leaf(u) && _gt.event(u) == Event::Duplication) {
			shape = "square";
		} else if (!_gt.is_leaf(u) && _gt.event(u) == Event::HGT) {
			shape = "triangle";
		}
		os << gt2a[u] << "[shape=" << shape << ",height=.5,width=.5,fixedsize=true,label="<<label<<"];" << std::endl;
	}
	boost::graph_traits<AuxGraph>::edge_iterator eit, eit_end;
	for (std::tie(eit, eit_end) = boost::edges(time_graph);
			eit != eit_end; ++eit) {
		auto e = *eit;
		size_t src = source(e, time_graph);
		size_t tar = target(e, time_graph);
		os << src << "->" << tar << std::endl;
	}
	os << "}";
}

void ReconMap::write_nexus_taxa(std::ostream& os) {
	std::unordered_set<std::string> leafs;
	for (auto n : _gt.nodes()) {
		if (_gt.is_leaf(n)) { leafs.insert(_gt.gene(n)); }
	}
	for (auto n : _st.nodes()) {
		if (_st.is_leaf(n)) { leafs.insert(_st.species(n)); }
	}
	os << "BEGIN TAXA;" << std::endl;
	os << "\tDIMENSIONS NTAX = " << leafs.size() << std::endl;
	os << "\tTAXLABELS" << std::endl;
	for (const std::string& s : leafs) {
		os << "\t\t'" << s << "'" << std::endl;
	}
	os << "\t\t;" << std::endl;

	os << "END;" << std::endl;
}

void ReconMap::write_nexus_trees(std::ostream& os) {
	os << "BEGIN TREES;" << std::endl;
	os << "\tTRANSLATE" << std::endl;
	for (auto n : _st.nodes()) {
		os << "\t\t" << _st.id(n) << "\t'";
		_st.write_leafset(os, n);
		os << "'," << std::endl;
	}
	for (auto n : _gt.nodes()) {
		os << "\t\t" << _gt.id(n) << "\t'";
		_gt.write_leafset(os, n);
		os << "'" << std::endl;
	}
	os << "\tTREE SPECIESTREE = ";
	_st.write_nexus(os);
	os << std::endl;
	os << "\tTREE GENETREE = ";
	_gt.write_nexus(os);
	os << std::endl;
	os << "END;" << std::endl;
}

void ReconMap::write_nexus_recon(std::ostream& os) {
	os << "BEGIN RECONCILIATION;" << std::endl;
	os << "[SIGMA represents the species a gene resides in]" << std::endl;
	os << "[Syntax is: gene_leaf_name species_leaf_name]" << std::endl;
	os << "\tSIGMA" << std::endl;
	for (auto n : _gt.nodes()) {
		if (_gt.is_leaf(n)) {
			os << "\t\t'" << _gt.gene(n) << "'\t'" << _gt.species(n) << "'," << std::endl;
		}
	}
	os << "\t\t;" << std::endl << std::endl;

	os << "[RECONCILIATIONMAP represents the reconciliation map]" << std::endl;
	os << "[SIGMA represents the species a gene resides in]" << std::endl;
	os << "[in case of a gene vertex to species vertex mapping then host_edge_from_name and host_edge_to_name are the same]" << std::endl;
	os << "\tRECONCILIATIONMAP" << std::endl;

	for (auto n : _gt.nodes()) {
		if (!_gt.is_leaf(n)) {
			os << "\t\t'" << _gt.id(n) << "' ";
			auto x = g2s(n);
			if (_gt.event(n) == Event::Speciation) {
				os << "'" << _st.id(x) << "' '" << _st.id(x) << "'";
			} else {
				os << "'" << _st.id(_st.parent(x)) << "' '" << _st.id(x) << "'";
			}
			os << "," << std::endl;
		}
	}
	os << "\t\t;" << std::endl;
	os << "END;" << std::endl;
}

void ReconMap::write_nexus(std::ostream& os) {
	os << "#NEXUS" << std::endl;
	write_nexus_taxa(os);
	os << std::endl;
	write_nexus_trees(os);
	os << std::endl;
	write_nexus_recon(os);
	os << std::endl;
}


} /* tc */ 
