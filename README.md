# TreeTop

TreeTop is an algorithm for single-cell data analysis to identify and give confidences scores to branch points in biological processes with possibly multi-level branching hierarchies. We demonstrate branch point identification for processes with varying topologies, including T cell maturation, B cell differentiation and hematopoiesis. Our analyses are consistent with recent experimental studies suggesting a shallow hierarchy of differentiation events in hematopoiesis, rather than the classical multi-level hierarchy.

TreeTop is described in the paper _Tree-ensemble analysis tests for presence of multifurcations in single cell data_ (Will Macnair, Laura De Vargas Roditi, Stefan Ganscha, Manfred Claassen).

## Installation

TreeTop is implemented in MATLAB. To install TreeTop, copy the files to somewhere suitable for storing MATLAB packages (either via `git clone https://github.com/wmacnair/TreeTop.git`, or via downloading and extracting a zip file).

Once you have some unpacked files, change into the directory and run the function `install_treetop.m`. This will add the relevant folders to the path, and check that various necessary functions work.

It is possible that you will need to compile some files in which case your platform will need to have a C compiler installed.

## Usage

Once you have successfully set up TreeTop, we recommend following the `treetop_tutorial` file for a guided tour of TreeTop, and a description of how to interpret the outputs. This is found in the folder `./TreeTop/`.

To replicate the results shown in the paper, please also download the data in the repository [TreeTop_data](https://github.com/wmacnair/TreeTop_data), and follow the instructions at the end of `treetop_tutorial`.
