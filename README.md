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

3) To test using our example sequence, run `python GRIFMAP.py Example/WINM57-cox1Region.txt`
