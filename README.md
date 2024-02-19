# HPC_Drug version 2
## Overview
A middleware python tool for computational drug discovery on HPC architectures.
It automates the setup of hamiltonian replica exchange HREM and non equilibrium alchemical transformations for protein-ligand systems. Its primary use is to calculate the binding free energy between a protein and a ligand with the virtual double system single box method [vDSSB](https://pubs.acs.org/doi/10.1021/acs.jctc.0c00634).

At the moment it only works for [Gromacs](https://www.gromacs.org/) (with [Plumed](https://www.plumed.org/)) but it can easily be adapted for any MD program that is supported by [Parmed](https://github.com/ParmEd/ParmEd). 

The latest version of the software can be found on GitHub [here](https://github.com/MauriceKarrenbrock/HPC_Drug), and the documentation can be found [here](https://mauricekarrenbrock.github.io/HPC_Drug/)

## What it can do and how to use it
HPC_Drug offers some out of the box scripts to set up a vDSSB by starting from a protein-ligand PDB file and a ligand SDF file (for example downloaded from the PDB databank).
Of course you can also create some custom ones.

To use them activate the conda environment you created
```
conda activate HPC_Drug
```
copy the script(s) that you want in the root directory
```
cp scripts/<SCRIPT NAME>.py .
```

At this point, for every script, you can read the documentation with the `-h` or `--help` option.
```
python <SCRIPT NAME>.py -h
```

Briefly what you can do:

* Once you have a PDB of the protein-ligand system and an SDF file of the ligand (for example downloaded from the PDB databank) you can add missing atoms and residues, protonate, rename the residues, and solvate the system with `prepare_pdb.py`

* Parametrize the system with any forcefield supported by [Openmm](https://github.com/openmm/openmm)/[Openmmforcefields](https://github.com/openmm/openmmforcefields) with `parametrize.py`

* Set up an hamiltonian replica exchange (solute tempering) for the protein-ligand system and for the ligand alone with `setup_hrem.py`

* Starting from the HREM create the input for the non equilibrium alchemical transformations, both for the bound (protein-ligand) and the unbound (only ligand) systems, with `fsdam_input.py`

* Postprocess the results of the non equilibrium alchemical transformations in order to get a binding free energy value and the various volume and charge corrections with `fsdam_postprocessing.py`

Of course this scripts will never be able to cover every use case, therefore advanced users might need to develop custom scripts and use HPC_Drug, but more so the other python packages that HPC_Drug depends on and that you can find on [my github](https://github.com/MauriceKarrenbrock?tab=repositories) as libraries (for experienced python developers it should be straightforward in most of the cases).

## How to set up the program
The easy way is to create a conda environment by using the environment.yml file. Because of the fact that this program has a lot of dependecies it might happen that the environment created in this way will not always be up to date.

The command to give on the terminal is
```
conda env create -f environment.yml
```
This command will create an environment called `HPC_Drug`. While creating the environment because of the large number of dependencies conda might complain about conflicts or take ages to create it, so you might have to play aroud a bit.

At this point you only have to copy one of the scripts from the `scripts` directory in the root directory of `HPC_Drug` (or you can create a custom one) and run it!

## Differences with version 1
The main difference is that version 2 dropped the support for the [Orac](http://www1.chim.unifi.it/orac/) MD program and doesn't use [Primadorac](http://www1.chim.unifi.it/orac/primadorac/) to parametrize the ligands anymore; instead, through the [Openmmforcefields](https://github.com/openmm/openmmforcefields) interface, it allows to parametrize the ligands with the Gaff forcefields through antechamber or with the [Openff](https://github.com/openforcefield/openff-forcefields) or [Smirnoff99Frost](https://github.com/openforcefield/smirnoff99Frosst/) ones through the [openff-toolkit](https://github.com/openforcefield/openff-toolkit).
Another difference is that the version one was more automathic but also less flexible than this version.

If for any reason you want to use the version 1 (that is the version referenced in the papers DOI: [10.1021/acs.jctc.0c00634](https://pubs.acs.org/doi/10.1021/acs.jctc.0c00634), and DOI: [10.1021/acs.jcim.1c00909](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00909)) switch to the branch version v1.0

## Issues
If you find a bug please open an issue on [GitHub](https://github.com/MauriceKarrenbrock/HPC_Drug/issues)

## Contributions
If you would like to contribute please first open an issue in order to discuss about it and then feel free to make a pull request, but please follow good coding practices and write some tests for the code you wrote.

## License
[GNU Affero General Public License (agpl v3)](https://www.gnu.org/licenses/agpl-3.0.en.html)

See the `LICENSE` file

## Citation
Please cite these papers:

1. Addressing Suboptimal Poses in Nonequilibrium Alchemical Calculations
Maurice Karrenbrock, Valerio Rizzi, Piero Procacci, Francesco Luigi Gervasio
The Journal of Physical Chemistry B 2024
DOI: [10.1021/acs.jpcb.3c06516](https://pubs.acs.org/doi/10.1021/acs.jpcb.3c06516)

1. Virtual Double-System Single-Box: A Nonequilibrium Alchemical Technique for Absolute Binding Free Energy Calculations: Application to Ligands of the SARS-CoV-2 Main Protease
Marina Macchiagodena, Marco Pagliai, Maurice Karrenbrock, Guido Guarnieri, Francesco Iannone, and Piero Procacci
Journal of Chemical Theory and Computation 2020 16 (11), 7160-7172
DOI: [10.1021/acs.jctc.0c00634](https://pubs.acs.org/doi/10.1021/acs.jctc.0c00634)

1. Virtual Double-System Single-Box for Absolute Dissociation Free Energy Calculations in GROMACS
Marina Macchiagodena, Maurice Karrenbrock, Marco Pagliai, and Piero Procacci
Journal of Chemical Information and Modeling 2021 61 (11), 5320-5326
DOI: [10.1021/acs.jcim.1c00909](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00909) 
