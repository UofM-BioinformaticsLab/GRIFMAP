# GRIFMAP
GRIFMAP (GRoup I Fungal Mitogenome Accurate Prediction)

## Installation instructions

GRIFMAP runs HMMER and Infernal, so both must be installed. 

### hmmer

Run `apt install hmmer` (or the appropriate equivalent for your system) to install *hmmer*.

See: http://hmmer.org/documentation.html for more details or to compile from source. 

### Infernal 

Run `apt-get install infernal infernal-doc` (or the appropriate equivalent for your system) to install *Infernal*.

See: http://eddylab.org/infernal/ for more details or to compile from source.

## Running GRIFMAP

Begin by cloning this repository and making it your working directory. 

1) Copy the AugmentedModel files from 'Group I Intron CM' to the main directory: `cp Group I Intron CM/* .`
 
    Alternatively, copy the single file *AugmentedModels.cm* to the main directory and run `cmpress AugmentedModels.cm`

2) Copy the *geneDatabasewoCountModel3* files from 'Group I Intron HMMER/' to the main directory: `cp Group I Intron HMMER/geneDatabasewoCountModel3* .`

3) To test using our example sequence, run `python GRIFMAP.py Example/WINM117-smallRegion.fa`

### Main output table of the test example
target name    query name      from    to    E-value    score  Possible Twin-tron found
---------------  ------------  ------  ----  ---------  -------  --------------------------
cox1             Graphilbum         3   734   4.1e-233    770.2  -
CMIB4S23         Graphilbum       735  1763   3.5e-22     117.6  No
cox1             Graphilbum      1764  1899   6.3e-31      99.9  -
CMIB4S23         Graphilbum      1900  3030   5.8e-18      96.9  No
cox1             Graphilbum      3031  3292   1.3e-73     241.4  -
CMIB4S23         Graphilbum      3293  5429   3.5e-13      73.5  No
cox1             Graphilbum      5430  5891   1.6e-110    363.8  -
