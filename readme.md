# Clustiflor

Clustering of graphs and biclustering of bipartite graphs algorithms with overlaps.



## Usage


## Clustering

Edge data files should be in the following format:

    1 10
    1 34
    # Lines starting with a # are ignored
    2 12
    a 12 # Labels can be string
    0 12 # Labels can start at 0
    a b 0.5 # third value will be weight
    # by default weight is 1
    # use --ignore-weights to ignore weights

Hyperparameters:

    --split-threshold=1 (float >= 1., default value: 1)
    --samples-size=10 (integer, default value: 10)
    --ignore-weights to ignore weights when loading a file


Launch

    cargo run --release data.txt
    cargo run --release data.txt --split-threshold=2



## Asymetrical Biclustering

Data files should be in the following format:

    1 10
    1 34
    # Lines starting with a # are ignored
    2 11
    2 12
    a 12 # labels can be string
    0 12 # labels can start at 0
    3 22 1 # default weight is 1
    4 33 2.5 # this edge has weight 2.5

Only labels on the left can be split.


Hyperparameters:

    --cost=1 should be in [0,1] (default value: 1)
    --split=1 should be >= 1 (default value: 1)
    --power=3 should be an integer (default value: 3)


Launch

    cargo run --release data.adj " "