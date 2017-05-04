#ifndef TC_TREE_HPP
#define TC_TREE_HPP
#include <boost/graph/adjacency_list.hpp>
/* #include <pugixml.hpp> */
#include <vector>

namespace tc {

struct EmptyProp {};

template<typename NodeProp = EmptyProp, typename EdgeProp = EmptyProp, typename TreeProp = EmptyProp>
class Tree {
private:
	struct MetaNodeProp {
		NodeProp data;
		int depth;
		size_t parent;
	};
	using Graph = typename boost::adjacency_list < boost::vecS
	                                             , boost::vecS
												 				, boost::directedS
																, MetaNodeProp
																, EdgeProp
																, TreeProp
																>;
public:
	using Node = typename boost::graph_traits<Graph>::vertex_descriptor;
	using Edge = typename boost::graph_traits<Graph>::edge_descriptor;
	using NodeIter = typename Graph::vertex_iterator;
	using ChildNodeIter = typename boost::graph_traits<Graph>::adjacency_iterator;

	Tree(): _graph(), _root(add_vertex(MetaNodeProp{NodeProp(), 0, 0 }, _graph)) { }
	Tree(NodeProp&& np): _graph(), _root(add_vertex(MetaNodeProp{np, 0, 0 }, _graph)) { }
	Tree(const NodeProp& np): _graph(), _root(add_vertex(MetaNodeProp{np, 0, 0 }, _graph)) { }

	Node root() const { return _root; }

	Node parent(Node v) const { return _graph[v].parent; }
	Node lca(Node u, Node v) const {
		while (u != v) {
			if (depth(u) > depth(v)) {
				u = parent(u);
			} else {
				v = parent(v);
			}
		}
	}

	Node lca(const std::vector<Node>& ns) const {
		assert(ns.size() > 0);
		if (ns.size() == 1) { return ns[0]; }
		Node lca_node = ns[0];
		for (size_t i = 1; i < ns.size(); ++i) {
			lca_node = lca(lca_node, ns[i]);
		}
		return lca_node;
	}

	auto children(Node n) const { return make_iterator_range(adjacent_vertices(n, _graph)); }
	auto nodes() const { return boost::make_iterator_range(vertices(_graph)); }
	size_t num_children(Node n) const { return out_degree(n, _graph); }
	size_t num_nodes() const { return num_vertices(_graph); }

	bool is_leaf(Node n) const { return num_children(n) == 0; }


protected:
	Node add_node_base(Node p, NodeProp&& np) {
		Node n = add_vertex(MetaNodeProp {np, (depth(p) + 1), p}, _graph);
		add_edge(p, n, _graph);
		return n;
	}

	const NodeProp& data(Node n) const { return _graph[n].data; }
	NodeProp& data(Node n) { return _graph[n].data; }

	/* template <typename AttrParser> */
	/* void read_xml_node_base(const pugi::xml_node& xn, Node p, AttrParser attr_parser) { */
	/* 	NodeProp np = attr_parser(xn); */
	/* 	Node n = add_node_base(p, std::move(np)); */
	/* 	for (auto xc : xn.children("Node")) { */
	/* 		read_xml_node_base(xc, n, attr_parser); */
	/* 	} */
	/* } */

private:
	int depth(Node v) const { return _graph[v].depth; }

private:
	Graph _graph;
	Node _root;
};
	
} /* namespace tc */
#endif /* ifndef TC_TREE_HPP */
