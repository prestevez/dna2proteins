[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.34023.svg)](http://dx.doi.org/10.5281/zenodo.34023)
# dna2proteins
This Python script was produced as part of the course *Introduction to Scientific Programming in Python* of the UCL Graduate School.

More information on the course can be found on its [home page](http://www.cs.ucl.ac.uk/scipython/).

## The script

The script puts together a collection of functions that essentially import a fastafile containing sequences of DNA and produce a fastafile with the most likely protein sequence for each DNA sequence.

The steps in the script are roughly the following:  
1. Reads in the fasta file  
2. Stores the sequences in a dictionary  
3. Generates the six possible frames for each sequence (+1, +2, +3 and -1, -2, -3)  
4. Swaps the DNA sequences for protein sequences  
5. Finds the longest protein sequence between an open and close marker  
6. Stores the longest protein sequence for each DNA sequence in a dictionary  
7. Can save the protein sequences on a fasta file or print the sequences on the terminal  

### Usage

The script is quite simple. It contains three options that can be passed from the command line:
- -h prints a very simple help
- -i (--ifile) must be followed by the fasta file
- -o (--ofile) must be followed by the name where the protein sequences will be stored
- -p is an option that allows printing the protein sequences on the terminal

To use the script enter the following in the terminal:

`$ python dna2proteins.py -i sequences.fa -o proteins.fa -p`

And substitute `sequences.fa` and `proteins.fa` for the appropriate filenames and paths.

### Credits

The code for this script was developed jointly by:
- Erin Vehstedt
- Johanna Fischer
- Maragatham Kumar
- Andrés Calderón
- Marya Koleva
- Patricio R. Estévez Soto
- With the guidance and help of Fabian Zimmer
