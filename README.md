# Phydelity 

_Inferring putative transmission clusters from phylogenetic trees_

## Latest updates 

* 25-Jul-2019: V2.1 - Fixed int64 type declaration for cross-platform compatibility; tested on Windows and Mac
* 10-May-2019: V2.0 - Improved algorithm yielding clusters with higher purity and lower probability of misclassification (see Manuscript for more details).  
* 6-Dec-2018: Fixed bug which was overly strict when cleaning up clusters that violated within-cluster limit.  

## Overview

Phydelity, a redesign of [PhyCLIP](http://github.com/alvinxhan/PhyCLIP), is a statistically-principled and phylogeny-informed tool capable of identifying putative transmission clusters in pathogen phylogenies without the introduction of arbitrary distance thresholds.  

Minimally, Phydelity only requires a phylogeny (in NEWICK format) as input.  

Phydelity infers the within-cluster divergence of putative transmission clusters by first determining the pairwise patristic distance distribution of closely-related tips. This distance distribution comprises of the pairwise distances of sequence _j_ and its closest _k_-neighbouring tips, where the closest _k_-neighbours includes sequence _j_. The user can **optionally** input the desired _k_ parameter **_OR_** allow Phydelity to **automatically scale _k_ to yield the supremum distribution with the lowest overall divergence**.  

To cite Phydelity:  
Alvin X Han, Edyth Parker, Sebastian Maurer-Stroh, Colin A Russell, Inferring putative transmission clusters with Phydelity, Virus Evolution, Volume 5, Issue 2, July 2019, vez039, https://doi.org/10.1093/ve/vez039

## Quickstart (for users of an academic institution only)
1. **Installation**  

* Install dependencies using [Anaconda](http://www.anaconda.com/download/).
  * Phydelity is written in **Python 2** (Users using Python 3 for base Conda environment can build a separate Python 2 environment. 
```
$ conda install -c etetoolkit ete3
$ conda install -c anaconda cython
$ conda install numpy scipy 
```
* For Windows users, you may also need to compile C extensions using Microsoft Visual C++ Compiler for Python. You can install that from [http://www.microsoft.com/en-us/download/details.aspx?id=44266](http://www.microsoft.com/en-us/download/details.aspx?id=44266).

* Install Gurobi (linear programming solver) using Anaconda as well. 
```
$ conda config --add channels http://conda.anaconda.org/gurobi
$ conda install gurobi
```

* Obtain Gurobi license.  
  Method I: If your machine is connected to the internet via a recognized academic domain (e.g. '.edu' addresss) 
    * Register a FREE account via http://www.gurobi.com/registration/academic-license-reg 
    * Login and access https://www.gurobi.com/downloads/end-user-license-agreement-academic/ to download license.  
    * Follow instructions on page to retrieve your license or check [here](https://github.com/alvinxhan/PhyCLIP/wiki/II.-Installation) for more details.   
  
  Method II: If your machine is **NOT** connected via an academic domain but you can verify that you are an academic user.  
    * To be added. 

* Download and install Phydelity. 
```
$ cd /path/to/Phydelity-master/
$ python setup.py install 
$ python setup.py clean --all 
```

2. **Run Phydelity**  
Minimal input command: 

```
$ phydelity.py --tree </path/to/treefile.newick>
```

See **Full options** below for other analysis options. 

3. **Outputs**
* cluster\_phydelity\_k<\d+>\_sol<\d+>\_<treefname>.txt - Tab-delimited text file of tip names and corresponding cluster ID. 
* tree\_phydelity\_k<\d+>\_sol<\d+>\_<treefname>.txt - Figtree-formatted NEXUS tree file with cluster annotations. 
* pdftree\_phydelity\_k<\d+>\_sol<\d+>\_<treefname>.txt - **Optional** PDF tree file if --pdf_tree is called. 

## Full options 
```
usage: phydelity.py [-h] [-t TREE] [--k K] [--outgroup OUTGROUP]
                    [--collapse_zero_branch_length]
                    [--equivalent_zero_length EQUIVALENT_ZERO_LENGTH]
                    [--solver_verbose {0,1}] [--solver_check] [--pdf_tree]

Phydelity v1.0

optional arguments:
  -h, --help            show this help message and exit

Required:
  -t TREE, --tree TREE  Input phylogenetic tree in NEWICK format.

Analysis options:
  --k K                 Custom k neighbours (optional).
  --outgroup OUTGROUP   Taxon (name as appeared in tree) to be set as outgroup
                        OR type 'midpoint' for midpoint-rooting.
  --collapse_zero_branch_length
                        Collapse internal nodes with zero branch length of
                        tree before running Phydelity.
  --equivalent_zero_length EQUIVALENT_ZERO_LENGTH
                        Maximum branch length to be rounded to zero if the
                        --collapse_zero_branch_length flag is passed (default
                        = 1e-06).

Solver options:
  --solver_verbose {0,1}
                        Gurobi solver verbose (default: 0)
  --solver_check        Check if Gurobi is installed.

Output options:
  --pdf_tree            PDF tree output annotated with cluster results (X
                        server required).

```
