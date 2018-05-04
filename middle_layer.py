#!/usr/bin/env python3

"""
Program:    Mesnage_MSc_Coursework_middle_layer
File:       middle_layer.py

Version:    V1.0
Date:       03.05.18
Function:   Take requests from the genome browser interface.
            Extract data from a database containing pieces of data from Genbank file.
            Perform processing of the data, including any calculations that need to be done.

Copyright:  (c) Robin Mesnage, Birkbeck, 2018
Author:     Robin Mesnage
Contact:    robin.mesnage@kcl.ac.uk
            
--------------------------------------------------------------------------

This program is released under the GNU Public Licence (GPL V3)

--------------------------------------------------------------------------
Description:
============
Provide the following functions to process pieces of data from a GenBank file:
           - get_ gene_features: Get the features of a gene (name, accession number, location, protein product)
           - get_genbank_data: Get the genbank data from the database. Generate the coding DNA sequence from the DNA.
           - highlight_cds: Identify where the coding regions are in the genomic DNA.
           - align_sequences: Provide a DNA sequence and an aminoacid sequence that can be aligned together.
           - get_codon_usage: Create a table with codon usage frequencies expressed as a percentage
           - highlight_re: Identify where restriction sites are in the genomic DNA
           - get_restriction_site_exon_list: Identify if restriction enzymes cut in the coding sequence and provide the genomic coordinates

--------------------------------------------------------------------------
Usage:
======
The set of functions described above can be called to process Genbank data using the query gene name as an argument.

This script uses the output of functions from the database layer, used to extract data from the SQL database containing the genbank information

The following functions can be called by importing the database layer module for a given gene (e.g. TMEM1)

            get_dictionary('TMEM1') -> [{'cds_id': 1, 'start': 41, 'end': 64}, {'cds_id': 2, 'start': 70, 'end': 95}]
            get_fasta('TMEM1') -> 'atgatagtacagatcgatccgatatac'
            get_aa_sequence('TMEM1') -> 'EEPLPPVIYTMENKPIVTCAGDQNLFTSVYPTLSQQ'

--------------------------------------------------------------------------

Revision History:
=================
V1.0   
"""
    

import re


# Access the database layer
# If the database import is failing, some dummy data will be loaded

try:
    mode = 'live_version'                     
    from data_access_layer_v0_2 import *
    from config_db import * 

except ImportError:
    mode = 'local_test'
    gene_name = 'dummy'
    print ("Impossible to import the database layer or the configuration file, please check the depedencies. Importing dummy data")

    # Create dummy functions containing dummy data
    def get_dictionary(gene_name):
        cds_coordinates_dict = [{'start': 31, 'end': 61, 'cds_id': 1}, {'start': 91, 'end': 125, 'cds_id': 2}]
        return (cds_coordinates_dict)

    def get_fasta(gene_name):
        fasta = "aagcttgtcggaattcgggctgctaacttacacttcagaggcctgtgtcccaaaggcctggctgcgtttgccgtgctgtgcgaggacctgtgtacacaggcaggtgtgcgcctgcccgaattcgagtagctcttgtgtagtggtgaaatgctgcaggcatctgtttattaaaat"
        return (fasta)

    def get_aa_sequence(gene_name):
        aa_sequence = "HFRGLCPKGLCTQAGVRLPERE"
        return (aa_sequence)

    def data_access(gene_name):
        genbank_data = "[['TMEM1'], {'TMEM1': 'AB001523'}, {'TMEM1': '21q22.3'}, {'TMEM1': 'BAA21099.1'}]"

    d_RE = {}
    d_RE['EcoRI'] = 'GAATTC'
    open_cds_tag = "<span style='background-color:#FFFF00'>"
    close_cds_tag = "</span>"
    open_re_tag = ["<span style='background-color:#23ff01'>"]
    close_re_tag = "</span>"


def get_gene_features(gene_name):

    """ Provide genbank information for a given gene.
        If the argument is an empty string, the list will be populated by dictionaries containing all the data available

                Args:
                    gene_name: query gene name

                Returns:
                    Return the information for a given query as a list

                Output example:
                    [['TMEM1'], {'TMEM1': 'AB001523'}, {'TMEM1': '21q22.3'}, {'TMEM1': 'BAA21099.1'}]


    """

    query_information         = data_access(gene_name)
    return(query_information)


def get_genbank_data(gene_name):

    """ This function is importing the arguments used by other functions, from the database, for a given query gene.

            Args:
                gene_name: A gene name from the Genbank data

            Returns:
                cds_coordinates_dict: a dictionary containing the DNA coding sequence coordinates
                fasta: the gene fasta sequence
                aa_sequence: the aminoacid sequence
                cds: A list of the DNA coding sequences

            Output examples:
                cds_coordinates_dict = [{'cds_id': 1, 'start': 41, 'end': 64}, {'cds_id': 2, 'start': 70, 'end': 95}]
                fasta = 'atgatacccaaaggcctggctgcgtttgccgtacagatcgatatgtatagcaccgatatac'
                aa_sequence = 'EEPLPPVIYTMENKPIVTCAGDQNLFTSVYPTLSQQ'
                cds = ['CCCAAAGGCCTGGCTGCGTTTGCC', 'ATGTATAGCA']

    """


    # Import the genbank variables corresponding to the query gene
    cds_coordinates_dict = get_dictionary(gene_name)
    fasta                = get_fasta(gene_name)
    aa_sequence          = get_aa_sequence(gene_name)
    
    # Check if the fasta sequence only contain DNA nucleotides   
    for letter in fasta:
        if letter not in "atcg":
            print("invalid input, sequence could only contain the following letters (a,t,c,g)")
            exit()

    # Place the start and the end of each codon in the correct reading frame
    for i in range(len(cds_coordinates_dict)):
        cds_coordinates_dict[i]['start'] = cds_coordinates_dict[i]['start'] - 1

    # Create a list containing the different exons
    cds = []
    for i in range(len(cds_coordinates_dict)):
        exon = fasta[cds_coordinates_dict[i]['start']:cds_coordinates_dict[i]['end']]
        cds.append(exon.upper())

    return {'fasta': fasta, 'cds_coordinates_dict': cds_coordinates_dict, 'aa_sequence': aa_sequence, 'cds': cds}


def get_coding_sequence(gene_name):

    """ Generate the coding DNA sequence from the DNA

                Args:
                    gene_name: A gene name from the Genbank data

                Returns:
                    Returns a string containing the full coding DNA sequence of the query gene

                Output example:
                    coding_DNA = 'CCCAAAGGCCTGGCTGCGTTTGCCATGTATAGCA'

    """

    genbank_data        = get_genbank_data(gene_name)
    cds                  = genbank_data['cds']

    # Join the different cds regions to get the full coding DNA sequence
    coding_DNA           = "".join(cds)

    return(coding_DNA)


def highlight_cds(gene_name):

    """ Identify where the coding regions are in the genomic DNA by highlighting them with HTML marks.
        The HTML marks are <span> marks containing a background colour argument used to change the background colour of the regions spanned by the marks.

                Args:
                    gene_name: query gene name

                Returns:
                    Returns a fasta sequence with opening and closing HTML marks surrounding each DNA coding region

                Output example:
                    marked_fasta = cttcagaggcctgtgt<span style='background-color:#FFFF00'>CCCAAAGGCCTGGCTGCGTTTGCC</span>gtgctgtgcgag

    """
    
    genbank_data        = get_genbank_data(gene_name)
    fasta                = genbank_data['fasta']
    cds                  = genbank_data['cds']

    # Create a formatted fasta with a colour code for the coding sequences
    marked_fasta = fasta
    for i in range(len(cds)):
        regex = re.compile(r'((' + cds[i] + ')(\s*(' + cds[i] + '))*)', re.I)
        marked_fasta = regex.sub(r'' + open_cds_tag + '' + cds[i] + '' + close_cds_tag + '', marked_fasta)
    
    return(marked_fasta)


def align_sequences(gene_name):

    """ Provide a DNA sequence and an aminoacid sequence that can be aligned if the font used has a fixed-width (e.g. monospace)

                Args:
                    gene_name: query gene name

                Returns:
                    Return the protein sequence in a format that can be aligned to the DNA coding sequence

                Output example:
                    protein_seq = 'E  E  P  L  P  P  V  I  Y  T  M  E  N  K  P  I'
    """

    genbank_data        = get_genbank_data(gene_name)
    aa_sequence          = genbank_data['aa_sequence']

    # Introduce two white spaces between each aminoacid letter to align the protein sequence to the DNA coding sequence
    protein_seq          = "  ".join(aa_sequence)

    return(protein_seq)


def get_codon_usage(gene_name):

    """ Create a table with codon usage frequencies expressed as a percentage

                    Args:
                        gene_name: query gene name. Alternatively, the query "ALL_THE_GENES" will calculate codon usage across the whole chromosome

                    Returns:
                        a dictionary containing codon usage frequencies

                    Output example:
                        CodonsDict = {'UUU': 2.4, 'UUC': 5.7, 'UUA': 0, 'UUG': 1.1, 'CUU': 1.3}
    """
    if gene_name == "ALL_THE_GENES":
        list_genes = get_gene_names()
        full_cds = []
        for gene_name in list_genes:
            genbank_data = get_genbank_data(gene_name)
            cds = genbank_data['cds']
            coding_dna = "".join(cds)
            full_cds.append(coding_dna)
        cds = "".join(full_cds)

    else:
        genbank_data        = get_genbank_data(gene_name)
        cds                  = genbank_data['cds']
        cds                  = "".join(cds)

    mrna                 = ''

    for i in range(len(cds)):
        mrna += cds[i].replace('T', 'U')

    CodonsDict = {
        'UUU': 0, 'UUC': 0, 'UUA': 0, 'UUG': 0, 'CUU': 0,
        'CUC': 0, 'CUA': 0, 'CUG': 0, 'AUU': 0, 'AUC': 0,
        'AUA': 0, 'AUG': 0, 'GUU': 0, 'GUC': 0, 'GUA': 0,
        'GUG': 0, 'UAU': 0, 'UAC': 0, 'UAA': 0, 'UAG': 0,
        'CAU': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAU': 0,
        'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAU': 0, 'GAC': 0,
        'GAA': 0, 'GAG': 0, 'UCU': 0, 'UCC': 0, 'UCA': 0,
        'UCG': 0, 'CCU': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
        'ACU': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCU': 0,
        'GCC': 0, 'GCA': 0, 'GCG': 0, 'UGU': 0, 'UGC': 0,
        'UGA': 0, 'UGG': 0, 'CGU': 0, 'CGC': 0, 'CGA': 0,
        'CGG': 0, 'AGU': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
        'GGU': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}

    for n in range(0, len(mrna), 3):
        CodonsDict[mrna[n:n + 3]] += 1

    # Calculate the codon frequency based on the number of occurrences
    def f(x):
        return x / (((len(cds) / 3)) / 100)

    CodonsDict = dict((k, round(f(v),2)) for k, v in CodonsDict.items())

    return(CodonsDict)


def highlight_re(gene_name):

    """ Identify where restriction sites  are in the genomic DNA

                Args:
                    gene_name: query gene name

                Returns:
                    A fasta sequence containing HTML marks that can be used to tisplay restriction sites on the fasta sequence having highlighted DNA coding sequences

                Output example for restriction site EcoRI (GAATTC):
                    marked_fasta = cttcagaggcctgtgt<span style='background-color:#ff0202'>GAATTC</span>gtgctgtgcgag
    """

    # Get the coding regions of the genes (cds)
    marked_fasta = highlight_cds(gene_name)

    RElist          = list(d_RE.values())

    for i in range(len(RElist)):
        regex = re.compile(r'((' + RElist[i] + ')(\s*(' + RElist[i] + '))*)', re.I)
        marked_fasta = regex.sub(r'' + open_re_tag[i] + '' + RElist[i] + '' + close_re_tag + '', marked_fasta)

    return (marked_fasta)


def get_restriction_sites(gene_name):

    """ Identify the coordinates of restriction sites in the genomic DNA

                Args:
                    gene_name: query gene name

                Returns:
                    Return the coordinates of the restriction sites for different enzymes in a dictionary

                Output example:
                    dict_RE_sites = {'EcoRI': [10, 50], 'BamHI': [79, 96]}
    """

    genbank_data = get_genbank_data(gene_name)
    fasta = genbank_data['fasta']
    fasta = fasta.upper()

    # Search the fasta file for restriction sites, list the coordinates of these sites in a dictionary
    dict_RE_sites = {}
    for key in d_RE:
        matches   = re.finditer(r'(' + d_RE[key] + ')', fasta)
        positions = []
        for m in matches:
            pos   = m.start()
            positions.append(str(pos))
            name  = key
            dict_RE_sites[name] = positions

    return (dict_RE_sites)


def get_exon_restriction_sites(gene_name):

    """ Identify if restriction enzymes cut in the coding sequence and provide the genomic coordinates

                Args:
                    gene_name: query gene name

                Returns:
                    Return the coordinates of the restriction sites for different enzymes cutting in the CDS

                Output example:
                    dict_RE_sites_exons = {'EcoRI': 10}


    """

    genbank_data        = get_genbank_data(gene_name)

    # Get the coding regions of the genes (cds)
    cds = genbank_data['cds']
    cds = "".join(cds)

    # Search the coding sequence for restriction sites, list the coordinates of these sites in a dictionary
    dict_RE_sites_exons = {}
    for key in d_RE:
        matches = re.finditer(r'(' + d_RE[key] + ')', cds)
        positions_exon = []
        for m in matches:
            pos = m.start()
            positions_exon.append(str(pos))
            name = key
            dict_RE_sites_exons[name] = positions_exon

    return (dict_RE_sites_exons)


# Activate this option if you want to calculate codon usage across the whole chomosome
#gene_name = "ALL_THE_GENES"
#print(get_codon_usage(gene_name))


'''
==========================================================================

  TESTING THE CODE
  
    Description. 
        This section allows the user to test the code by generating an HTML page based on the dummy data
        The dummy data is loaded only if the database layer or the config file import has failed.

==========================================================================
'''


if mode == 'local_test':
    choice = input('Do you want to execute the functions on dummy data to produce a test html page? Y or N :')
    if choice == "Y":
        # Use the business layer to process the genbank data
        highlighted_cds = highlight_cds(gene_name)
        coding_sequence = get_coding_sequence(gene_name)
        alignment = align_sequences(gene_name)
        codon_usage = get_codon_usage(gene_name)
        highlighted_re = highlight_re(gene_name)
        dict_RE_sites = get_restriction_sites(gene_name)
        dict_RE_sites_exons = get_exon_restriction_sites(gene_name)

        # Write an HTML variable for a CSS-formatted webpage
        html = "<html>\n"
        html += "<head>\n"
        html += "<style type='text/css'>"
        html += "<!--"
        html += "body { background: #FFFFFF; }"
        html += ".highlight_cds_head {padding: 5px; font: bold 10pt Helvetica, Arial, sans-serif; background: #336699; color: #FFFFFF; margin: 10px 0px 0px 0px;}"
        html += ".highlight_cds {background: #ffffff; color: #000000; font: 8pt Helvetica; margin: 0px 0px 0px 0px; word-break: break-word;}"
        html += ".alignment_head {padding: 5px; font: bold 10pt Helvetica, Arial, sans-serif; background: #336699; color: #FFFFFF; margin: 10px 0px 0px 0px;}"
        html += ".alignment {background: #ffffff; color: #000000; font-family: monospace; margin: 0px 0px 0px 0px; white-space: pre;vertical-align: bottom;}"
        html += ".text {background: #ffffff; color: #000000; font: 8pt Helvetica; margin: 0px 0px 0px 0px;padding: 5px;}"
        html += "table {border-collapse: collapse; border: 1px solid black;}"
        html += "th, td {text-align: left;padding: 2px;border-bottom: 1px solid #ddd;font: 8pt Helvetica;}"
        html += "th {background-color: #B5CDE6; color: black;}"
        html += "-->"
        html += "-->"
        html += "</style>"
        html += "<title>This is the gene list</title>\n"
        html += "<h1>Sequence analysis results</h1>\n"
        html += "<body>\n"

        html_codon = "<p class='alignment_head'> Create a table with codon usage frequencies expressed as a percentage</p>\n"
        html_codon += """<table><tr><th>Codon</th><th>Usage Frequency (%)</th></tr>"""
        for codon in codon_usage:
            html_codon += "<tr><td>{}</td>".format(codon)
            html_codon += "<td>{}</td>".format(codon_usage[codon])
            html_codon += "</tr>"
        html_codon += "</table>"

        html += "<p class='highlight_cds_head'> Identify where the coding regions are in the genomic DNA of " + gene_name  + "</p>\n"
        html += "<p class='highlight_cds'>"+ highlighted_cds +"</p>\n"
        html += "<p class='highlight_cds_head'> Generate the coding DNA sequence from the DNA of " + gene_name + "</p>\n"
        html += "<p class='highlight_cds'>"+ coding_sequence +"</p>"
        html += "<p class='alignment_head'> Align the protein sequence translation with the DNA coding sequence of " + gene_name  + "</p>\n"
        html += "<p class='alignment'>"+ coding_sequence +"</p>"
        html += "<p class='alignment'>"+ alignment +"</p>\n"
        html += "<p class='highlight_cds_head'> Identify where the restriction sites are in the genomic DNA of " + gene_name  + "</p>\n"
        html += "<p class='highlight_cds'>"+ highlighted_re +"</p>\n"
        html += "<p class='alignment_head'> Provide a list of available REs to the front end for the gene " + gene_name  + "</p>\n"
        for k, v in dict_RE_sites.items():
                html += "<p class='text'> - The enzyme "+ k +" has restriction sites at positions: "+ str(v) +"</p>"
        html += "<p class='alignment_head'> Identify whether an RE has restriction sites within the coding region of " + gene_name  + "</p>\n"
        for k, v in dict_RE_sites_exons.items():
                html += "<p class='text'> - The enzyme "+ k +" has restriction sites in the coding regions at positions: "+ str(v) +"</p>"

        html += html_codon

        html += "</body>\n"

        web_page = open("filename.html", "w")
        web_page.write(html)
        web_page.close()

    else:
        exit()


