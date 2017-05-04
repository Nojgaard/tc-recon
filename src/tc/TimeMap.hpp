#ifndef TC_TIME_MAP_HPP
#define TC_TIME_MAP_HPP
#include <vector>
#include "SpeciesTree.hpp"
#include "GeneTree.hpp"

namespace tc {
	class TimeMap {
	public:
		TimeMap(const SpeciesTree& st, const GeneTree& gt): stp(st.num_nodes()), gtp(gt.num_nodes()) {}
		int& operator () (SpeciesTree::Node u, const SpeciesTree&) { return stp[u]; }
		int& operator () (GeneTree::Node u, const GeneTree&) { return gtp[u]; }
	private:
		std::vector<int> stp, gtp;
	};
} /* tc */ 

#endif /* ifndef TC_TIME_MAP_HPP */
