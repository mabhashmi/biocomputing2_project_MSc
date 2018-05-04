# MSc Bioinformatics Biocomputing II - Project work

========================================================================================

## Aim and Scope

The aim of this project was to create a genome browser to explore the content of a given chromosome using Genbank data.

This program allows the identification of the coding regions for a given gene, the alignment of its coding DNA sequence to its protein sequence, as well as the identification of restriction enzyme sites including those with a site within a coding region. In addition, the codon usage frequencies of a given gene can be calculated and compared to the codon usage across all coding regions.

========================================================================================

## Execution of the script

The python script attached to this documentation is not meant to be used in isolation, but in conjunction with 2 other scripts. The first of these scripts is coding a data access tier allowing the extraction of appropriate pieces of data from an SQL database filled with the details of a GenBank file. The other script is coding a web-based graphical user interface. When it is fully deployed on a live system, the script described in this documentation imports pieces of Genbank data from the data access tier, process the data, and provides the web interface with a set of variables containing different features of a given gene.

The user of this script is not supposed to edit the main python file. A configuration file is provided for the configuration of the database access and different options such as the list of restriction sites used by the genome browser.

However, this script can be tested on a local machine. If the python interpreter is able to import the dependencies correctly, the program will behave like if it is installed on a live system and it will only act as a set of functions. However, if the database layer is not accessible, the python interpreter will produce an error message and ask if the user wants to run some dummy data locally to verify the integrity of the script. The output of this test will be stored as an HTML file in the same directory as the main script.

========================================================================================

## Functions

### get_ gene_features()
	
This function is importing a list containing the gene name, its location, and its accession number. This function can be used by the webpage developer to display the whole list of genes as a table. If the argument is an empty string, the list will be populated by dictionaries containing the data for all the genes.

```
print(get_ gene_features((dummy))
[['TMEM1'], {'TMEM1': 'AB001523'}, {'TMEM1': '21q22.3'}, {'TMEM1': 'BAA21099.1'}]
```

### get_genbank_data()
	
This function is importing the arguments used by other functions from the database for a given query. The access to the database information is provided through three functions. 

The function get_dictionary() imports a dictionary containing the coordinates of DNA coding sequences. Since python uses a 0-based indexing, we subtracted one to each of the dictionary values to place the start and end of each codon in the correct reading frame. The function get_fasta() imports a string containing the fasta sequence of the query gene. The content of this string is then checked for non-nucleotide characters. We assume that the presence of errors in the sequence can have dramatic downstream effect for the user. In this case, the function exit is called to exit the interpreter. The function get_aa_sequence() will import the complete aminoacid sequence for the protein coded by the query gene.

A cds variable containing a list of the coding DNA sequences is then created to be used in the subsequent functions. This variable contains uppercase characters by opposition to the fasta variable containing lowercase characters, which is a useful feature to distinguish the two types of sequences.

All these 4 variables (the fasta, the dictionary of cds coordinates, the aminoacid sequence, and the list of cds sequences) are returned as a dictionary.
	
```	
print(get_genbank_data('TMEM1'))
{'fasta':'aagcttgtcggaattcgggctgctaacttacacttcagaggcctgtgtcccaaaggcctggctgcgtttgccgtgctgtgcgaggacctgtgtacacaggcaggtgtgcgcctgcccgaattcgagtagctcttgtgtagtggtgaaatgctgcaggcatctgtttattaaaat', 'cds_coordinates_dict': [{'start': 30, 'end': 61, 'cds_id': 1}, {'start': 90, 'end': 125, 'cds_id': 2}], 'aa_sequence': 'HFRGLCPKGLCTQAGVRLPERE', 'cds': ['CACTTCAGAGGCCTGTGTCCCAAAGGCCTGG', 'TGTACACAGGCAGGTGTGCGCCTGCCCGAATTCGA']}
```	

### get_coding_sequence()

The function get_coding_sequence() is very simple. It only aggregates the elements in the list of DNA coding sequence (cds variable) in order to return the full coding sequence that can be aligned with the protein sequence.

```
print(get_coding_sequence('TMEM1'))
'CACTTCAGAGGCCTGTGTCCCAAAGGCCTGGTGTACACAGGCAGGTGTGCGCCTGCCCGAATTCGA'
```

### highlight_cds ()

The function highlight_cds() return a the fasta sequence corresponding to the query gene, with the addition of marks which can be used by the web page developer in order to highlight the DNA coding regions. The nature of these marks can be changed as they are provided as a set of variables, open_cds_tag and close_cds_tag, respectively. In the example data given to run the script locally, the marks are HTML <span> marks containing a 'background-color:#FFFF00' argument changing the sequence background colour between the opening and closing marks.

```
print(highlight_cds(dummy))
'aagcttgtcggaattcgggctgctaactta<span style='background-color:#FFFF00'>CACTTCAGAGGCCTGTGTCCCAAAGGCCTGG</span>ctgcgtttgccgtgctgtgcgaggacctg<span style='background-color:#FFFF00'>TGTACACAGGCAGGTGTGCGCCTGCCCGAATTCGA</span>gtagctcttgtgtagtggtgaaatgctgcaggcatctgtttattaaaat'
```

### align_sequences()

The function align_sequences() provides an aminoacid sequence that can be aligned with the coding DNA sequence if the font used has a fixed-width (e.g. monospace). The reason is that a codon can be represented as three characters. Thus, a protein sequence with two-spaces inserted between each aminoacid residue can be aligned to the coding DNA sequence.

```
print(align_sequences(dummy))
‘H  F  R  G  L  C  P  K  G  L  C  T  Q  A  G  V  R  L  P  E  R  E’
```

### get_codon_usage()

The function get_codon_usage() provides the front end with a dictionary containing codon usage for a given gene expressed as a frequency. This function operates using a sliding window of 3 residues on the DNA coding sequence. It is used to populate a predefined dictionary containing the 64 combinations of nucleotides of the codon code. 

The codon usage of a given gene can be compared to the codon usage of the whole chromosome. The codon usage in the whole chromosome has already been calculated and stored in the database. 

If the user wants to calculate again the codon usage for the whole chromosome, the user can edit the script by removing the two hashtags at the start of the following statements  (#gene_name = "ALL_THE_GENES", #print(get_codon_usage(gene_name)) and execute the script from the command line. 

```
print(get_codon_usage(dummy))
‘{'UUU': 0.0, 'UUC': 4.55, 'UUA': 0.0, 'UUG': 0.0, 'CUU': 0.0, 'CUC': 0.0, 'CUA': 0.0, 'CUG': 9.09, 'AUU': 4.55, 'AUC': 0.0, 'AUA': 0.0, 'AUG': 0.0, 'GUU': 0.0, 'GUC': 0.0, 'GUA': 0.0, 'GUG': 4.55, 'UAU': 0.0, 'UAC': 4.55, 'UAA': 0.0, 'UAG': 0.0, 'CAU': 0.0, 'CAC': 4.55, 'CAA': 0.0, 'CAG': 0.0, 'AAU': 0.0, 'AAC': 0.0, 'AAA': 4.55, 'AAG': 0.0, 'GAU': 0.0, 'GAC': 0.0, 'GAA': 0.0, 'GAG': 0.0, 'UCU': 0.0, 'UCC': 0.0, 'UCA': 0.0, 'UCG': 0.0, 'CCU': 4.55, 'CCC': 4.55, 'CCA': 0.0, 'CCG': 0.0, 'ACU': 0.0, 'ACC': 0.0, 'ACA': 4.55, 'ACG': 0.0, 'GCU': 0.0, 'GCC': 4.55, 'GCA': 0.0, 'GCG': 4.55, 'UGU': 9.09, 'UGC': 0.0, 'UGA': 0.0, 'UGG': 0.0, 'CGU': 0.0, 'CGC': 0.0, 'CGA': 9.09, 'CGG': 0.0, 'AGU': 0.0, 'AGC': 0.0, 'AGA': 4.55, 'AGG': 4.55, 'GGU': 0.0, 'GGC': 13.64, 'GGA': 0.0, 'GGG': 0.0}
```

### highlight_re()

The function highlight_re() is used to display the position of a defined set of restriction sites. This function works in a similar way than highlight_cds(). The front end manager can change the marks added to the fasta string by editing a set of variables, open_re_tag and close_re_tag, respectively. Note that the output of the function highlight(cds) is used as an input, meaning that the DNA coding sequences will be highlighted as defined by highlight(cds). This will allow the user to identify restriction enzymes cutting in the middle of a coding region.

In the example given below, the EcoRI restriction sites (GAATTC) are surrounded by am HTML mark providing an argument for a green background colour. 

```
print(highlight_re(dummy))
aagcttgtcg<span style='background-color:#23ff01'>GAATTC</span>gggctgctaactta<span style='background-color:#FFFF00'>CACTTCAGAGGCCTGTGTCCCAAAGGCCTGG</span>ctgcgtttgccgtgctgtgcgaggacctg<span style='background-color:#FFFF00'>TGTACACAGGCAGGTGTGCGCCTGCCC<span style='background-color:#23ff01'>GAATTC</span>GA</span>gtagctcttgtgtagtggtgaaatgctgcaggcatctgtttattaaaat
```

### get_exon_restriction_sites() and get_restriction_sites ()

These two function are used to obtain the coordinates of the restriction sites. The function get_restriction_site_list() provides a dictionary containing the coordinates of the restriction sites in the fasta sequence. The function get_restriction_site_exon_list() provides the coordinates of the restriction sites in the DNA coding sequence when a restriction enzyme is cutting in the coding region.

```
print(get_restriction_sites(dummy))
{'EcoRI': ['10', '117']}
print(get_exon_restriction_sites(dummy))
{'EcoRI': ['58']}
```

========================================================================================

## Author

Robin Mesnage          robin.mesnage [at] kcl.ac.uk

========================================================================================

## License

This project is licensed under the GNU Public Licence (GPL V3)

========================================================================================

## Version 1.0

========================================================================================

## Acknowledgments 

I thank Rebecca Brooker, Chelsea Wescott and Laurentino Quiroga Moreno for their invaluable help in creating the genome browser presented in this project


