#include "SpeciesTree.hpp"
#include <cassert>
#include <pugixml.hpp>
#include <stdexcept>
#include <iostream>

namespace tc {

SpeciesTree::SpeciesTree() : impl::SpeciesTreeBase() {}

const std::string& SpeciesTree::species(Node n) const {
	assert(data(n).label != "");
	return data(n).label;
}

SpeciesTree::Node SpeciesTree::add_node(Node p) {
	return add_node_base(p, impl::SNodeProp{});
}

SpeciesTree::Node SpeciesTree::add_node(Node p, std::string species) {
	return add_node_base(p, impl::SNodeProp{species});
}

void SpeciesTree::read_xml(const pugi::xml_document& doc) {
	auto xst = doc.child("SpeciesTree");
	if (!xst) {
		throw std::invalid_argument("Root node must be 'SpeciesTree'");
	}

	auto xn = xst.child("Node");
	if (!xn) {
		throw std::invalid_argument("No root found in SpeciesTree");
	}

	auto attr_parser = [] (const pugi::xml_node& xn) -> impl::SNodeProp {
		return impl::SNodeProp{ xn.attribute("Species").value()};
	};

	read_xml_node_base(xn, root(), attr_parser);
}

void SpeciesTree::print() {
	print_node(root());
	std::cout << std::endl;
}

void SpeciesTree::print_node(Node n) {
	if (is_leaf(n)) {
		std::cout << species(n);
	} else {
		std::cout << "(";
		for (Node c : children(n)) { 
			print_node(c); 
			if (*(children(n).end() - 1) != c) {
				std::cout << ", ";
			}
		}
		std::cout << ")";
	}
}

std::vector<SpeciesTree::Node> SpeciesTree::leaf_set(Node n) const {
	std::vector<Node> q, ls;
	q.push_back(n);
	while (!q.empty()) {
		Node c = q.back();
		q.pop_back();
		if (is_leaf(c)) { ls.push_back(c); }
		for (auto cc : children(c)) { q.push_back(cc); }
	}
	return ls;
}


void SpeciesTree::write_dot(std::ostream& os) const {
	for (auto x : nodes()) {
		std::string label = "\"\"";
		if (is_leaf(x)) {
			label = "\"" + species(x) + "\"";
		}
		os << "s" << x << "[shape=circle,width=.5,fixedsize=true,label="<<label<<"];" << std::endl;
	}
	/* os << "node[shape=none, width=0, height=0, label=\"\",fixedsize=true];" << std::endl; */
	os << "node[constraint=false,shape=point, width=0.01, width=0.01, label=\"\"];" << std::endl;
	/* os << "edge[arrowhead=none];" << std::endl; */
	for (auto x : nodes()) {
		for (auto y : children(x)) {
			os << "s" << x << " -- sb" << x << y << "-- s" << y << ";" << std::endl;
		}
	}
}
std::string species_read_label(size_t& i, const std::string& s) {
	std::string label;
	while (i < s.size() && s[i] != ')' && s[i] != ',') {
		if (s[i] == '(' || s[i] == ';' || s[i] == '[' || s[i] == ']') { 
			throw std::invalid_argument("parse error: pos " + std::to_string(i));
		}
		label += s[i]; ++i;
	}
	if (i == s.size()) { 
		throw std::invalid_argument("parse error: pos " + std::to_string(i)); 
	}
	return label;
}

void SpeciesTree::read_newick(const std::string& str) {
	assert(str.size() > 0 && str[0] == '(');
	size_t i = 1;
	Node cur_node = add_node(root());
	while (i < str.size() && str[i] != ';') {
		/* std::cout << i << "('" << str[i] << "')" << " "; */
		if (str[i] == '(') {
			++i;
			cur_node = add_node(cur_node);
		} else if (str[i] == ')') {
			++i;
			cur_node = parent(cur_node);
		} else if (str[i] == ',') {
			++i;
		} else {
			auto lbl = species_read_label(i, str);
			Node n = add_node(cur_node);
			data(n).label = lbl;
		}
	}
}

std::string SpeciesTree::id(Node n) const {
	if (is_leaf(n)) { return species(n); }

	return "S" + std::to_string(n);
}

void SpeciesTree::write_nexus(std::ostream& os, Node n) const {
	if (is_leaf(n)) {
		os << species(n);
	} else {
		os << "(";
		Node last = *(children(n).end() - 1);
		for (Node c : children(n)) {
			write_nexus(os, c);
			if (c != last) {
				os << ",";
			}
		}
		os << ")" << id(n);
	}
}

void SpeciesTree::write_nexus(std::ostream& os) const {
	write_nexus(os, root());
	os << ";";
}

void SpeciesTree::write_leafset(std::ostream& os, Node n) const {
	if (is_leaf(n)) {
		os << species(n);
	} else {
		for (Node c : children(n)) {
			write_leafset(os, c);
		}
	}
}
	
} /* tc */ 
