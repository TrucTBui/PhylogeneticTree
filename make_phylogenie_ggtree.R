
# Specify your R libraries path here if needed. The R library of the remote server does not contain 
# taxize and metacoder so I have to use my additional library
# .libPaths(c("/mnt/raidbio/biosoft/software/R/R-4.3.0.mkl/library", 
#            "/mnt/cip/home/b/buit/R/x86_64-pc-linux-gnu-library/4.3"))

library(ggtree)
library(taxize)
library(metacoder)
library(argparser)
library(ape)

p <- arg_parser("Parser")

p <- add_argument(p, "--input", help="Input file with species names")
p <- add_argument(p, "--layout", help="type of graph: rectangle or circular")
p <- add_argument(p, "--output", help="Ouput folder where the phylogenetic tree plot and tree text in Newick format should be saved")

argv <- parse_args(p)
input <- argv$input
graph_layout <- argv$layout
output <- argv$output

# From a list of species names to a map with IDs, using metacoder package
species_list <- read.csv(input,sep="\n")
# species_list <- read.csv("./HiWi/species_list.txt",sep="\n")

# Note: This method is usable, but it gives HTTP error 400 sometimes due to the fragility of NCBI, so we won't use it
# taxa_map <- metacoder::lookup_tax_data(species_list,type="taxon_name",database = "ncbi")
# taxize_class <- taxize::classification(taxa_map$data$query_data$taxon_id, db = "ncbi")

# Instead, we can use try-error:

unknown_species <- c()
id_list <- c()
for (sp in species_list$Species) {
  classes_i <- try(metacoder::lookup_tax_data(sp,type="taxon_name",database = "ncbi"))
  while(class(classes_i)[1]=="try-error") {
    Sys.sleep(1)
    classes_i <- try(metacoder::lookup_tax_data(sp,type="taxon_name",database = "ncbi"))
  }
  
  if (names(classes_i$data$query_data) == "unknown") { #Search for synonyms
    synonym <- taxize::synonyms(sp, db = "itis")
    synonyms_list <- synonym[[sp]]$syn_name
    if (!is.null(synonyms_list)) {
      # Print the list of synonyms
      cat(sp, "not found. Searching for potential synonyms...","Synonyms for", sp, ":\n")
      for (i in seq_along(synonyms_list)) {
        cat(i, ": ", synonyms_list[i], "\n", sep = "")
      }
      
      # Prompt user for choice using command-line arguments
      cat("Choose a synonym for", sp, ". Enter the number of your choice or 0 to discard the name:\n")
      choice <- as.integer(readLines(con = "stdin", n = 1))
   
      # Check if a valid choice was made
      if (choice == 0 || choice > length(synonyms_list)) {
        message("The name will be be discared.")
        unknown_species <- c(unknown_species, sp)  #discard the name
        next
      } else {
        # Return the chosen synonym
        cat("The chosen synonym is:", synonyms_list[choice], "\n")
        Sys.sleep(3)
        classes_i <- try(metacoder::lookup_tax_data((synonyms_list[choice]),type="taxon_name",database = "ncbi"))
      }
    }
    
    else {
      unknown_species <- c(unknown_species, sp)  #discard the name
      next
    }
  }
 
  id_list <- c(id_list,names(classes_i$data$query_data)) #add the taxa ID to the list
}

# Print out the unknown names. These names are not included in the phylogenetic tree
cat("Unknown species: ", knitr::combine_words(unknown_species, and = ""))

# Fetch classification and generate phylo tree, using taxize packacge
taxize_class <- try(taxize::classification(id_list, db = "ncbi"))
while(class(taxize_class)[1]=="try-error") {
  Sys.sleep(3)
  taxize_class <- try(taxize::classification(id_list, db = "ncbi"))
}

taxize_tree <- taxize::class2tree(taxize_class, check = TRUE)

# Extract the phylo object
phylo_tree <- taxize_tree$phylo


# Create the ggtree plot from the phylo tree
ggtree_plot <- ggtree(phylo_tree, branch.length = "branch.length", layout=graph_layout) + 
  geom_tiplab(aes(label = label), align = TRUE, linesize = 0.5, size = 6, hjust = -0.1,
              fontface = "italic", color = "blue") +  # adjust justification
  theme_tree2() + geom_tree() +
  geom_nodepoint(color = "red", size = 2, shape = 21, fill = "yellow") +  # customize node points
  scale_x_continuous(expand = c(0, 7))  + # adjust the x-axis scale to control edge length
  geom_nodelab(vjust=-0.5,size=5)
  # geom_hilight(node = 10, fill = "lightblue") +  
  # geom_hilight(node = 15, fill = "lightgreen") 


if (graph_layout == "circular"){
  ggsave(paste0(output,"/phylogenetic_tree_circular.png"), plot = ggtree_plot , width = 25, height = 20)
}

if (graph_layout == "rectangular"){
  ggsave(paste0(output,"/phylogenetic_tree_rectangular.png"), plot = ggtree_plot, width = 35, height = 15)
}

# Save the tree as Newick format. The txt file can be used for visualization in iTOL
ape::write.tree(phylo_tree,file=paste0(output,"/phylogenetic_tree_newick.txt"))


#########################################################################################

