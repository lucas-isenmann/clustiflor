# Clustiflor

Clustering and biclustering algorithm with overlaps.



## Usage

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

    --cost=1 should be in [0,1]
    --split=1 should be positive
    --power=3 should be an integer