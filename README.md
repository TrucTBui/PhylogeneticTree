# PhylogeneticTree: An R Script for Phylogenetic Trees

## Overview
The *make_phylogenie_ggtree.R* script is designed to generate a phylogenetic tree from a list of species names. The script outputs two key items:

#### A phylogenetic tree plot, which can be displayed in either a rectangular or circular layout.
#### A Newick format file representing the phylogenetic tree structure. 
This file can be easily used in the **iTOL** (Interactive Tree Of Life) website for visualizing and analyzing phylogenetic trees. The Newick format provides a compact way to represent tree structures using nested parentheses and branch lengths, which is widely supported by various tree visualization tools. You can upload this file to iTOL to customize and explore your tree with features such as different layouts, color schemes, and annotations. 

## Requirements
Before running the script, ensure you have the following R packages installed:

**taxize**: For taxonomic data handling and phylogenetic tree construction.

**metacoder**: For querying taxonomic databases. 

**ape**: For working with phylogenetic trees and Newick format files.

**ggtree**: For visualizing phylogenetic trees.

**argparser**: For parsing command-line arguments.


You can install the required packages using:

install.packages(c("ggtree", "taxize", "metacoder", "argparser", "ape"))

## Arguments
**--input**: Path to the input file containing the list of species names (one species per line).

**--layout**: The layout for the tree plot. Options are "rectangular" or "circular".

**--output**: Path to the folder where the output files (tree plot and Newick file) will be saved.

### Example:
Rscript make_phylogeny_ggtree.R --input species_list.txt --layout circular --output path_to_your_folder

## License
This script is provided under the MIT License. See the LICENSE file for more details.

## Contact
For questions or feedback, please contact thanhtruc.bui@tum.de
