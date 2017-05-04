#ifndef TC_GRAPH
#define TC_GRAPH 
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/range/iterator_range.hpp>

#include <iostream>


struct CycleDetector : public boost::dfs_visitor<>
{
 CycleDetector( bool& has_cycle) 
	: _has_cycle(has_cycle) { }

 template <class Edge, class Graph>
 void back_edge(Edge, Graph&) {
	_has_cycle = true;
 }
protected:
 bool& _has_cycle;
};

template <typename Graph>
bool has_cycle(const Graph& g) {
	bool has_cycle = false;
	CycleDetector cycle_detector(has_cycle);
	boost::depth_first_search(g, boost::visitor(cycle_detector));
	return has_cycle;
}

template <typename Graph>
std::vector<size_t> topo_sort(const Graph& g) {
	std::vector<size_t> vs;
	boost::topological_sort(g, std::back_inserter(vs));
	return vs;
}

template <typename Graph>
void print_graph(const Graph& g) {
	for (auto u : boost::make_iterator_range(vertices(g))) {
		std::cout << u << ": ";
		for (auto v : make_iterator_range(adjacent_vertices(u, g))) {
			std::cout << v << " ";
		}
		std::cout << std::endl;
	}
}

#endif /* ifndef TC_GRAPH */
