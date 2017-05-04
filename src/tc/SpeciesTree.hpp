#ifndef TC_SPECIES_TREE
#define TC_SPECIES_TREE
#include <boost/optional.hpp>
#include "Tree.hpp"
#include <vector>
/* #include <pugixml.hpp> */
#include <iostream>

namespace tc {

namespace impl {

struct SNodeProp {
	std::string label;
};

using SpeciesTreeBase = Tree<SNodeProp>;
	
} /* impl */ 

class SpeciesTree : public impl::SpeciesTreeBase {
public:
	SpeciesTree();

	const std::string& species(Node n) const;

	Node add_node(Node p, std::string species);
	Node add_node(Node p);
	std::string id(Node n) const;
	
	std::vector<Node> leaf_set(Node n) const;

	/* void read_xml(const pugi::xml_document& doc); */
	void write_dot(std::ostream& os) const;
	void write_nexus(std::ostream& os = std::cout) const;
	void print();
	void read_newick(const std::string& str);
	void write_leafset(std::ostream& os, Node n) const;
private:
	void print_node(Node n);
	void write_nexus(std::ostream& os, Node n) const;
};
	
} /* tc */ 

#endif /* ifndef TC_SPECIES_TREE */
