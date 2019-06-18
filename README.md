# ConsensusModules
# Comparing consensus modules using S2B and MODifieR

Overview
Modules are inferred from stand-alone methods (DIAMOnD, MCODE, WGCNA) and serve as inputs for S2B, which is used to create additional consensus modules for comparison. The R package for S2B is available from frpinto/S2B, MODifieR is available from ddeweerd/MODifieRDev, and there are testable datasets available from OMIM (omim.org).

Validation
Validation is inherently difficult in this work because disease modules are largely incomplete today, so there is no universally accepted standard to adopt. This is particularly problematic because measures of sensitivity and specificity cannot be derived. In the absence of an ideal standard, benchmark validation for this study was carried out using Pascal. Pascal is available from https://www2.unil.ch/cbg/index.php?title=Pascal

Summarized here is the method’s general order:
1.	Create a consensus module using three methods: clique-based (MCODE), seed-based (DIAMOnD), and hierarchical (WGCNA).
2.	Create a second consensus module in S2B using combinations of disease modules from MCODE, DIAMOnD, and WGCNA as inputs.
3.	Evaluate the results from 1 and 2 for functional enrichment using Pascal.

Initially, MCODE and DIAMOnD are used to establish a data-processing pipeline. Consensus modules are created by identifying subgraphs of overlapping nodes. This process is then repeated using MCODE and DIAMOnD modules as inputs for S2B. Garcia et al. compared pairs of modules, but in this study combinations of both two- and three-inference modules are derived.

An RDS object (e.g., GSE76262_asthma_sputum_SA.rds) is loaded into the R environment. The object contains a data frame with the following:

•	diff_genes: a data frame of lists, including genes (Entrez IDs), p-values, an annotated expression matrix, and a list of DIAMOnD genes

•	annotation_table: a data frame of lists, including probe IDs, probe symbols, and Entrez IDs

•	group_indicii: an index of controls and patients

Both MCODE and DIAMOnD require loading of a protein-protein interaction (PPI) network in addition to an RDS object. The PPI is a comma-separated (.csv) file that contains the following:

•	entrez1: a list of Entrez identifiers for the network (connections to entrez2)

•	entrez2: a list of Entrez identifiers for the network (connections to entrez1)

•	SCORE: a weighted score of each connection (value range: 700–1,000)

The output from each method is a disease module that contains a data frame of lists, including “module_genes”, which are the significant genes identified by each method. Lists of genes from each module are then written to a text file, which is later evaluated for enrichment of GWAS SNPs using Pascal. The module gene lists serve as inputs for S2B in the next step.

The format of the initial PPI network (.csv) needs to be converted into an igraph data frame using the function igraph::graph.data.frame. The S2B function simpmain is then executed on the igraph data frame to remove loops and redundant edges.

Lists of module genes (from the derived outputs of MCODE, DIAMOnD, and WGCNA) are then assigned to temporary lists arbitrarily named “geneset1” and “geneset2”. The simpmain function is then executed on these temporary lists to produce the module inputs index1 and index2. The number of randomizations are defined using the parameters nrep and nrep2, which returns specificity scores — the default for both parameters is 100, but significant results can be found just by leaving these set to 1, which reduces the overall computation time with little discernible effect on the output. S2B is then called in RStudio accordingly:

S2B::S2B(seed_graph, index1, index2, nrep, nrep2)

The output is a new S2B consensus module. Subsequent combinations of modules are then tested:
1.	DIAMOnD/MCODE
2.	MCODE/WGCNA
3.	WGCNA/DIAMOnD

Each module is then validated using Pascal.

Additionally, consensus modules based on the stand-alone modules are created. MODifieR includes a function for combining stand-alone modules just by taking an overlap of WGCNA, MCODE, and DIAMOnD. The MODifieR function is called create_consensus_list, and it provides several useful ways of assessing similarities and differences among disease modules. For example, it can provide a simple side-by-side comparison of modules and return a comprehensive list of overlapping genes.

Modules outputted to S2B are recombined with stand-alone modules using three permutations:
1.	DIAMOnD/MCODE + WGCNA
2.	MCODE/WGCNA + DIAMOnD
3.	WGCNA/DIAMOnD + MCODE

Upon validation with Pascal, the program delivers a single “meta p-value” for each module.
