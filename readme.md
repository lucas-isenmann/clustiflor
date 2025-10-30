# Clustiflor

Clustering of graphs and biclustering of bipartite graphs algorithms with overlaps.



## Clustering

## n m File format

Graph file format used is a line "n m" (where n is the number of vertices and m is the number of edges)followed by the list of the edges with or without weights.
Without weights, the default weight is 1.

    n m
    i j
    i j w

For example this graph has 4 vertices and 5 edges and the weight of (2,3) is 0.5

    4 5
    0 1
    1 2
    2 3 0.5
    0 2
    0 3

## Simple File Format

With the option --simple-file-format, the format used is the following:
It is a list of pairs (i,j) or triplets (i,j,w).
Labels i and j are considered as strings so you can use "2" or "b" as labels.


    1 10
    1 34
    # Lines starting with a # are ignored
    2 12
    a 12 # Labels can be string
    0 12 # Labels can start at 0
    a b 0.5 # third value will be weight
    # by default weight is 1

Option --ignore-weights ignore the weights of the files and use 1 by default.


## Launch

Launch

    ./clustiflor.cli.exe cluster data.edges

Or with the simple file format : 

    ./clustiflor.cli.exe cluster data.simple.edges --simple-file-format

Clusters are written in data.edges.clusters with one line per cluster.

## Biclustering or One-sided Biclustering

### Matrix format

For example the following bipartite graph has 4 lines and 3 columns:

    0 1 0
    1 1 1
    1 1 1
    0 0 1

You can replace the weights by any number.

### Simple Bipartite Graph Format

List of edges (row, column) with string labels and possible weight.
If weight is missing, the default weight is 1.
Comment lines start with a #.
Remark that it is possible that a column and a row can have the same label.

    a b 0.5
    # comment
    1 2
    0 0

### Bicluster 

The bicluster algorithm uses the Matrix Format

Example:

    ./clustiflor-cli.exe bicluster data.matrix

Biclusters are printed in file data.matrix.BiMarkov3.results

Result format is:

    Biclusters
    2 3    # the number of columns, the number of rows in the bicluster
    1 2    # the list of the columns (columns indices are in 1 to n)
    4 5 6  # the list of the rows (rows indices are in n+1 to n+m)


### One-sided Bicluster

One-sided biclsuter will use Simple Bipartite Graph Format by default.
If you want to use the Matrix Format, use option --matrix-format.

Example:

    ./clustiflor-cli.exe onesided-bicluster data.matrix --matrix-format

Results are printed in .biclusters

Columns indices are from 0 to n-1 and rows indices are from n to n+m-1


## Options for every algorithm

    -s, --split-threshold <SPLIT_THRESHOLD>
            Split threshold: [default: 1] in [0,1]
        --size-sensitivity <SIZE_SENSITIVITY>
            Size sensitivity: [default: 1] should be >= 1
        --samples-size <SAMPLES_SIZE>
            Samples size [default: 10] should be >= 1
    -i, --ignore-weights
            Ignore weights (all weights are set to 1)
        --simple-file-format
            Simple file format
    -v, --verbose <VERBOSE>
            Verbose: 0 means print nothing, 1 prints details for each iteration, 2 prints a lot of details [default: 0]       
    -h, --help
            Print help
    -V, --version
            Print version




## Help

Use

    ./clustiflor.cli.exe --help

