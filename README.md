# gene_ordering

This R code simulates random gene loss as postulated for the yeast genome.

### Description

This code simulates two species that derive from the yeast whole genome duplication (WGD) event.  The ancestral yeast genome is duplicated, then genes are lost (deleted) randomly until only the number of genes observed in the modern yeast genome remain.  The model determines all gene pairs in the two genomes (excluding genes that remain as duplicates) and calculates the number of gene pairs that are shared between the two species.

The simulation code is written in R, where it requires only the base installation, and is all contained in a single code file, [gene_ordering.R](https://github.com/mpcox/gene_ordering/blob/master/gene_ordering.R).

Note that the code is not parallelized.  However, the simulations are independent and therefore well suited to being run in an ‘embarrassingly parallel’ setting (i.e., running multiple jobs and manually concatenating the output files).

### Running

There are only five input parameters:

-    The number of iterations (i.e., the number of simulation tests you want to perform) [default: 100]
-    The number of genes in the ancestral yeast genome [default: 4943]
-    The number of genes in the modern yeast genome (post-duplication and after gene loss) [default: 5500]
-    The number of chromosomes in the yeast pre-WGD genome [default: 8]
-    The output filename [default: outfile.txt]


### Output

Each simulation output is given on a single line, with multiple lines listing the results from different simulations.  

The output has 8 columns:

-    n.shared – the number of gene pairs shared between species 1 and species 2
-    shared.gene.pairs – the actual shared gene pairs
-    sim1.duplicates – the genes present as duplicates in species 1
-    sim2.duplicates – the genes present as duplicates in species 2
-    sim1.gen1.pairs – all gene pairs present on the first homeologous set of chromosomes in species 1
-    sim1.gen2.pairs – all gene pairs present on the second homeologous set of chromosomes in species 1
-    sim2.gen1.pairs – all gene pairs present on the first homeologous set of chromosomes in species 2
-    sim2.gen2.pairs – all gene pairs present on the second homeologous set of chromosomes in species 2

The program also produces a histogram of the number of gene pairs shared between species 1 and species 2 (i.e., n.shared).


### Worked example

A worked example is shown in the tab-delimited file, [example_outfile.tab](https://github.com/mpcox/gene_ordering/blob/master/example_outfile.tab).

Murray Cox

2 April 2020
