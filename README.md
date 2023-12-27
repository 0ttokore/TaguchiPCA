
# TaguchiPCA

TaguchiPCA performs clustering of big non-vectorized data (in
particular, nucleotide sequences) using data dimensionality reduction
methods (PCA, PCoA). The package also provides the ability to create an
imitation model of nucleotide sequences with specified parameters, such
as: number of sequences, their length, number of clusters formed from
these sequences.

## Installation

You can install the development version of TaguchiPCA from
[GitHub](https://github.com/0ttokore/TaguchiPCA) with:

``` r
install.packages("devtools")
devtools::install_github("0ttokore/TaguchiPCA")
```

## Imitation model description

    The function creates a list of nucleotide sequences based on the first generated
    sequence with specified parameters. Imitation model consists of two 
    parts.generate_positions() function returns nucleotide positions that will be 
    regenerated, while the generate_sequences() function performs generation of 
    sequences. As a result, the obtained sequences will have some similarity to each
    other, leading them to be clustered together.

## TaguchiClustering description

    The function takes a list of sequences and a distance metric, based on which the
    given sequences will be vectorized. Then the function computes and returns the 
    PC_num number of principal coordinates and also performs a transition to a new 
    feature space.
