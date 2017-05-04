# tc-recon

# Installation
```
git clone https://github.com/Nojgaard/tc-recon.git
```

# Compilation
Assuming you are in the root folder of the project, run the following command:
```
mkdir build && cd build && cmake .. && make
```
This should compile the program as the executable `tc`.

# Input File Format
The input file should consist of exactly two lines. Each line is a newick
string specifying the species tree and gene tree repspectively.

The species tree is written in classic newick format:
```
((A,B),(C,D))
```

In contrast, specifying the genetree requires a bit of custom syntax. For any
given leaf, it is required to specify the gene name and the species it resides
in as follows:
```
gene-name[species-name]
```

Additionaly every internal node must be given an event, by adding
S(Speciation), D(Duplication), or H(HGT), just after the closing parenthesis of
an internal node:
```
((a1[A],b1[B])S,a2[A])D
```

Note, if some internal node is specified as a HGT event, it is required that at
least one of its children are specified as a transfer-edge, by adding the tag
[t] as follows:
```
((a1[A],[t]b1[B])H,a2[A])D
```

For a larger example see `data/test1_valid.newick`

# Running the program
To run the program: (assuming the executable is in your current folder build)
```
./tc -i ../data/test1_valid.newick -o test.nexus
```

This will output a nexus file, test.nexus, specifying a time consistent
reconsiliation map. If no such map exists, for the given gene and species tree,
no file is written and the user is notified of such in the standard output.

To test such a case, run:
```
./tc -i ../data/test2_invalid.newick
```

If no output argument is given, everything is written to standard output.

To write a dot file, depcting the gene tree, species tree, and reconciliation
map, simply use the argument -d:
```
./tc -i ../data/test1_valid.newick -o test.dot -d
```

Note that an extra node is always appended as the root of the species tree.
