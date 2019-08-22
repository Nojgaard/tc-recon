#include "GeneTree.hpp"
/* #include <pugixml.hpp> */
#include <iostream>

namespace tc {

GeneTree::GeneTree(): impl::GeneTreeBase(impl::NodeProp{Event::Speciation, false}) { }

const std::string& GeneTree::gene(Node n) const {
	return get<impl::Leaf>(data(n).label).gene;
}

const std::string& GeneTree::species(Node n) const {
	return get<impl::Leaf>(data(n).label).species;
}

Event GeneTree::event(Node n) const {
	return get<Event>(data(n).label);
}

GeneTree::Node GeneTree::left(Node n) const {
	assert(!is_leaf(n));
	assert(num_children(n) == 2);
	auto cs = children(n);
	return *cs.begin();
}

GeneTree::Node GeneTree::right(Node n) const {
	assert(!is_leaf(n));
	assert(num_children(n) == 2);
	auto cs = children(n);
	return *(cs.begin() + 1);
}

std::vector<GeneTree::Node> GeneTree::transfer_nodes(Node n) const {
	std::vector<Node> tn;
	for (Node c : children(n)) {
		if (is_transfer(c)) { tn.push_back(c); }
	}
	return tn;
}

bool GeneTree::is_transfer(Node n) const {
	return data(n).is_transfer;
}

GeneTree::Node GeneTree::add_node(Node p, Event e, bool is_transfer) {
	return add_node_base(p, impl::NodeProp{e, is_transfer});
}

GeneTree::Node GeneTree::add_node(Node p, std::string gene, std::string species, bool is_transfer) {
	return add_node_base(p, impl::NodeProp{impl::Leaf{gene,species}, is_transfer});
}

Event parse_event(const std::string& se) {
	if (se == "S") {
		return Event::Speciation;
	} else if (se == "D") {
		return Event::Duplication;
	} else if (se == "H") {
		return Event::HGT;
	} else {
		throw std::invalid_argument("Unknown event symbol: " + se);
	}
}

/* void GeneTree::read_xml(const pugi::xml_document& doc) { */
/* 	auto xgt = doc.child("GeneTree"); */
/* 	if (!xgt) { throw std::invalid_argument("No GeneTree Found"); } */

/* 	auto xn = xgt.child("Node"); */
/* 	if (!xn) { throw std::invalid_argument("No root node found"); } */

/* 	if (!xn.attribute("Event")) { throw std::invalid_argument("Root must be an event node"); } */

/* 	auto attr_parser = [] (const pugi::xml_node& xn) -> impl::NodeProp { */
/* 		bool is_transfer = false; */
/* 		if (xn.attribute("Transfer")) { is_transfer = xn.attribute("Transfer").as_bool(); } */
/* 		if (xn.attribute("Event")) { */
/* 			return impl::NodeProp { parse_event(xn.attribute("Event").value()), is_transfer }; */
/* 		} else { */
/* 			return impl::NodeProp { impl::Leaf { xn.attribute("Gene").value() */
/* 				                                , xn.attribute("Species").value()}, is_transfer }; */
/* 		} */
/* 	}; */

	/* auto rp = attr_parser(xn); */
	/* if (rp.is_transfer) { throw std::invalid_argument("Root cannot be transfer node"); } */
	/* get<Event>(data(root()).label) = get<Event>(rp.label); */

	/* for (auto xc : xn.children("Node")) { read_xml_node_base(xc, root(), attr_parser); } */
/* } */


void GeneTree::print_node(Node n) {
	if (is_leaf(n)) {
		/* std::cout << gene(n) << "[" << species(n) << ", " << is_transfer(n) << "]"; */
		std::cout << gene(n) << "[" << species(n) << "]";
	} else {
		std::cout << "(";
		for (Node c : children(n)) { print_node(c); }
		std::string se;
		switch(event(n)) {
			case Event::Speciation: se = "S"; break;
			case Event::Duplication: se = "D"; break;
			case Event::HGT: se = "H"; break;
		};
		/* std::cout << ")[E=" << se << " T= " << is_transfer(n) << "]"; */
		std::cout << ")" << se;
	}
	if (is_transfer(n)) {
		std::cout << "[1]";
	}
}

void GeneTree::print() {
	print_node(root());
	std::cout << std::endl;
}

void GeneTree::write_dot(std::ostream& os, bool directed) const {
	for (auto u : nodes()) {
		std::string shape = "circle";
		std::string label = "\""+std::to_string(u)+"\"";
		if (is_leaf(u)) {
			label = "\"" + species(u) + "\"";
		}
		if (!is_leaf(u) && event(u) == Event::Duplication) {
			shape = "square";
		} else if (!is_leaf(u) && event(u) == Event::HGT) {
			shape = "triangle";
		}
		os << "g" << u << "[shape=" << shape << ",height=.5,width=.5,fixedsize=true,label="<<label<<"];" << std::endl;
	}
	std::string etype = (directed) ? " -> " : " -- ";
	for (auto u : nodes()) {
		for (auto v : children(u)) {
			std::string style = "solid";
			if (is_transfer(v)) { style = "dashed"; }
			os << "g" << u << etype << "g" << v << "[weight=10,splines=curved, style=" << style << "];" << std::endl;
		}
	}
}

void GeneTree::read_newick_node(const std::string& str, size_t& cur_pos, Node parent_node) {

}

std::string read_label(size_t& i, const std::string& s) {
	std::string label;
	while (i < s.size() && s[i] != ')' && s[i] != ',' && s[i] != '[' && s[i] != ']') {
		if (s[i] == '(' || s[i] == ';') { 
			throw std::invalid_argument("parse error: pos " + std::to_string(i));
		}
		label += s[i]; ++i;
	}
	if (i == s.size()) { 
		throw std::invalid_argument("parse error: pos " + std::to_string(i)); 
	}
	return label;
}

Event read_event(char c) {
	switch(c) {
		case 'S': return Event::Speciation;
		case 'D': return Event::Duplication;
		case 'H': return Event::HGT;
		default: 
		throw std::invalid_argument("Unknown event symbol: " + c);
	}
}

void read_pattern(size_t& i, const std::string& str, std::string pattern) {
	for (char c : pattern) {
		if (str[i] != c) {
			throw std::invalid_argument("Expected: '" + std::to_string(c) + "' at position" + std::to_string(c));
		}
		++i;
	}
}

void GeneTree::read_newick(const std::string& str) {
	assert(str.size() > 0 && str[0] == '(');
	size_t i = 1;
	Node cur_node = root();
	bool should_transfer = false;
	while (i < str.size() && str[i] != ';') {
		/* std::cout << i << "('" << str[i] << "')" << " "; */
		if (str[i] == '(') {
			++i;
			cur_node = add_node(cur_node, Event::Speciation, should_transfer);
			should_transfer = false;
		} else if (str[i] == ')') {
			++i;
			get<Event>(data(cur_node).label) = read_event(str[i]);
			++i;
			cur_node = parent(cur_node);
		} else if (str[i] == ',') {
			++i;
		} else if (str[i] == '[') {
			read_pattern(i, str, "[t]");
			should_transfer = true;
		
		} else {
			auto lbl = read_label(i, str);
			std::string species_lbl = lbl;
			if (str[i] == '[') {
				++i;
				species_lbl = read_label(i,str);
				++i;
			}
			Node n = add_node(cur_node, lbl, species_lbl, should_transfer);
			should_transfer = false;
		}
	}
}

std::string GeneTree::id(Node n) const {
	if (is_leaf(n)) {
		return gene(n);
	}
	return "G" + std::to_string(n);
}

void GeneTree::write_nexus(std::ostream& os, Node n) const {
	if (is_transfer(n)) {
		os << "[t]";
	}
	if (is_leaf(n)) {
		os << gene(n);
		return;
	}
	
	os << "(";
	Node last_child = *(children(n).end() - 1);
	for (Node c : children(n)) {
		write_nexus(os, c);
		if (c != last_child) {
			os << ",";
		}
	}
	os << ")";
	os << id(n) << "[" << event(n) << "]";
}

void GeneTree::write_nexus(std::ostream& os) const {
	write_nexus(os, root());
	os << ";";
}

void GeneTree::write_leafset(std::ostream& os, Node n) const {
	if (is_leaf(n)) {
		os << gene(n);
	} else {
		for (Node c : children(n)) {
			write_leafset(os, c);
		}
	}
}

} /* namespace tc */
