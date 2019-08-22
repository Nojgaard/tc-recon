#ifndef TC_GENE_TREE_HPP
#define TC_GENE_TREE_HPP
#include "Tree.hpp"
#include <boost/variant.hpp>
/* #include <pugixml.hpp> */
#include <iostream>

namespace tc {
enum Event {Speciation, Duplication, HGT};

inline std::ostream& operator << (std::ostream& os, const Event& e) {
	switch(e) {
		case Speciation: os << "S"; break;
		case Duplication: os << "D"; break;
		case HGT: os << "H"; break;
	}
	return os;
}

namespace impl {

struct Leaf {
	std::string gene, species;
};

struct NodeProp {
	boost::variant<Event, Leaf> label;
	bool is_transfer;
};

using GeneTreeBase = Tree<NodeProp>;
	
} /* impl */ 


class GeneTree: public impl::GeneTreeBase {
public:
	GeneTree ();
	const std::string& gene(Node n) const;
	const std::string& species(Node n) const;
	Event event(Node n) const;
	bool is_transfer(Node n) const;
	Node left(Node n) const;
	Node right(Node n) const;
	std::string id(Node n) const;

	Node add_node(Node p, Event e, bool is_transfer = false);
	Node add_node(Node p, std::string gene, std::string species, bool is_transfer = false);
	std::vector<Node> transfer_nodes(Node n) const;

	/* void read_xml(const pugi::xml_document& doc); */
	void read_newick(const std::string& str);
	void write_nexus(std::ostream& os = std::cout) const;
	void write_dot(std::ostream& os = std::cout, bool make_directed = false) const;
	void write_leafset(std::ostream& os, Node n) const;
	void print();

private:
	void read_newick_node(const std::string& str, size_t& cur_pos, Node parent_node);
	void write_nexus(std::ostream& os, Node n) const;
	void print_node(Node n);
	/* data */
};

	
} /* tc  */ 
#endif /* ifndef TC_GENE_TREE_HPP */
