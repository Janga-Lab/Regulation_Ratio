# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 19:00:07 2020

@author: Xander
"""

import argparse, numpy as np, os, pandas as pd
os.path.join(os.path.dirname(__file__))

def Get_Chromosme_Length(row):
    """
    Returns the length of a chromosome based on the given row.

    Parameters:
    - row: A pandas DataFrame row containing the 'Chromosome' column.

    Returns:
    - The length of the chromosome.

    Note:
    - If the 'Chromosome' value is not found in the function, it returns pd.NA.
    """

    # Check the value of the 'Chromosome' column and return the corresponding length
    if (row['Chromosome'] == '1'):
        return 248956422
    elif (row['Chromosome'] == '2'):
        return 242193529
    elif (row['Chromosome'] == '3'):
        return 198295559
    elif (row['Chromosome'] == '4'):
        return 190214555
    elif (row['Chromosome'] == '5'):
        return 181538259
    elif (row['Chromosome'] == '6'):
        return 170805979
    elif (row['Chromosome'] == '7'):
        return 159345973
    elif (row['Chromosome'] == '8'):
        return 145138636
    elif (row['Chromosome'] == '9'):
        return 138394717
    elif (row['Chromosome'] == '10'):
        return 133797422
    elif (row['Chromosome'] == '11'):
        return 135086622
    elif (row['Chromosome'] == '12'):
        return 133275309
    elif (row['Chromosome'] == '13'):
        return 114364328
    elif (row['Chromosome'] == '14'):
        return 107043718
    elif (row['Chromosome'] == '15'):
        return 101991189
    elif (row['Chromosome'] == '16'):
        return 90338345
    elif (row['Chromosome'] == '17'):
        return 83257441
    elif (row['Chromosome'] == '18'):
        return 80373285
    elif (row['Chromosome'] == '19'):
        return 58617616
    elif (row['Chromosome'] == '20'):
        return 64444167
    elif (row['Chromosome'] == '21'):
        return 46709983
    elif (row['Chromosome'] == '22'):
        return 50818468
    elif (row['Chromosome'] == 'X'):
        return 156040895
    elif (row['Chromosome'] == 'Y'):
        return 57227415
    elif (row['Chromosome'] == 'MT'):
        return 16569
    else:
        return pd.NA

def Make_Gene_Index(file='Gene_Information.txt'):
    """
    This function reads a gene information file and creates a pandas DataFrame with additional columns for gene length and chromosome length.

    Parameters:
    - file (str): The file path of the gene information file. Defaults to 'Gene_Information.txt'.

    Returns:
    - pandas DataFrame: A DataFrame containing the gene information with additional columns for gene length and chromosome length.

    Note:
    - The input file should be a tab-separated file with columns: 'Gene Name', 'Gene Symbol', 'Chromosome', 'Gene Start', 'Gene End', 'Strand'.
    - The function drops rows with missing values in the 'Chromosome' column.
    - The function calculates the gene length by subtracting 'Gene Start' from 'Gene End'.
    - The function calculates the chromosome length using the 'Get_Chromosme_Length' function.
    - The function calculates the gene location by dividing the midpoint of the gene by the chromosome length.
    """
    gene_Index = pd.read_csv(file, sep='\t', header=0)
    gene_Index.columns = ['Gene Name', 'Gene Symbol', 'Chromosome', 'Gene Start', 'Gene End', 'Strand']
    gene_Index['Gene Start'] = gene_Index['Gene Start'].astype(int)
    gene_Index['Gene End'] = gene_Index['Gene End'].astype(int)
    gene_Index['Chromosome'] = gene_Index['Chromosome'].apply(lambda x: str(x))
    gene_Index['Gene Length'] = gene_Index['Gene End'] - gene_Index['Gene Start']
    gene_Index['Chromosome Length'] = gene_Index.apply(lambda x: Get_Chromosme_Length(x), axis=1)
    gene_Index['Gene Location'] = ((gene_Index['Gene End'] + gene_Index['Gene Start']) / 2) / gene_Index['Chromosome Length']
    gene_Index = gene_Index.dropna(axis='index')
    gene_Index.index = range(len(gene_Index.index))
    return gene_Index

def Get_Peaks(file_Name):
    """
    This function reads a file containing genomic peak data and returns a dictionary 
    where keys are chromosome names and values are lists of peak start and end positions.

    Parameters:
    - file_Name (str): The name of the file containing genomic peak data. The file should be tab-separated and contain at least three columns: chromosome name, peak start position, and peak end position.

    Returns:
    - dict: A dictionary where keys are chromosome names and values are lists of peak start and end positions.

    Note:
    - The function skips lines that start with 'chr' and lines that have less than two columns.
    - The function converts peak start and end positions to integers.
    """
    gene_Peaks = {}
    with open (file_Name) as this_File:
        this_Text = this_File.readlines()
        this_File.close()
    for line in this_Text:
        line = line.split('\t')
        if (len(line) > 1) and (line[0]!= 'chr'):
            if (line[0] not in gene_Peaks.keys()):
                gene_Peaks[line[0]] = []
            gene_Peaks[line[0]].append([int(line[1]),int(line[2])])
    return gene_Peaks

def Get_Categories(row):
    """
    This function categorizes genes based on their regulation ratio.

    Parameters:
    - row (pandas Series): A row from a pandas DataFrame containing gene information.
        It should have a column named 'Regulation Ratio'.

    Returns:
    - str: A string representing the category of the gene.
        It can be either 'Predominantly Post-Transcriptional',
        'Predominantly Transcriptional', or 'Neutral'.

    Note:
    - The function assumes that the input DataFrame is a result of the 'Get_Correlation' function.
    - The function uses the 'Regulation Ratio' column to categorize genes.
    """
    if (row['Regulation Ratio'] > 1):
        return 'Predominantly Post-Transcriptional'
    elif (row['Regulation Ratio'] < 1):
        return 'Predominantly Transcriptional'
    else:
        return 'Neutral'

def Get_Correlation(atac_File, pop_File, upstream = 5000, downstream = 500):
    """
    This function calculates the correlation between ATAC-seq and POP-seq data.

    Parameters:
    - atac_File (str): The file path of the ATAC-seq data.
    - pop_File (str): The file path of the POP-seq data.
    - upstream (int, optional): The upstream distance for peak detection. Defaults to 5000.
    - downstream (int, optional): The downstream distance for peak detection. Defaults to 500.

    Returns:
    - pandas DataFrame: A DataFrame containing the correlation data.

    The function calculates the number of ATAC-seq and POP-seq peaks for each gene,
    and then calculates the correlation between the two.
    """
    gene_Index = Make_Gene_Index()
    gene_Index.index = range(len(gene_Index.index))
    atac_Peaks = Get_Peaks(atac_File)
    pop_Peaks = Get_Peaks(pop_File)
    atac_Hits = []
    atac_Indices = {}
    pop_Hits = []
    pop_Indices = {}
    for i in gene_Index.index:
        line = gene_Index.iloc[i]
        this_Atac = []
        this_Pop = []
        atac_Count = 0
        if (line['Chromosome'] in atac_Peaks.keys()):
            for j in (atac_Peaks[line['Chromosome']]):
                if (line['Strand'] == 1):
                    if (j[0] > (line['Gene Start'] - upstream)) and (j[0] < line['Gene End'] + downstream):
                        atac_Count += 1
                        this_Atac.append((j[0],j[1]))
                elif (line['Strand'] == -1):
                    if (j[0] > (line['Gene Start'] - downstream) and (j[0] < line['Gene End'] + upstream)):
                        atac_Count += 1
                        this_Atac.append((j[0], j[1]))
                if (j[0] > line['Gene End']):
                    break
        atac_Indices[line['Gene Name']] = this_Atac
        atac_Hits.append(atac_Count)
        if (line['Chromosome'] in pop_Peaks.keys()):
            pop_Count = 0
            for j in (pop_Peaks[line['Chromosome']]):
                if (j[0] > line['Gene Start']) and (j[0] < line['Gene End']):
                    pop_Count += 1
                    this_Pop.append((j[0],j[1]))
                if (j[0] > line['Gene End']):
                    break
        pop_Indices[line['Gene Name']] = this_Pop
        pop_Hits.append(pop_Count)
    atac_Hits = []
    pop_Hits = []
    for i in gene_Index['Gene Name']:
        atac_Hits.append(len(atac_Indices[i]))
        pop_Hits.append(len(pop_Indices[i]))
    gene_Index['Gene Sites'] = atac_Hits
    gene_Index['Transcript Sites'] = pop_Hits
    gene_Index['Gene Regulation'] = gene_Index['Gene Sites'] / gene_Index['Gene Length']
    gene_Index['Transcript Regulation'] = gene_Index['Transcript Sites'] / gene_Index['Gene Length']
    gene_Index['Regulation Ratio'] = gene_Index['Transcript Sites'] / gene_Index['Gene Sites']
    gene_Index['Gene Regulation Group'] = gene_Index.apply(lambda x: Get_Categories(x), axis=1)
    return gene_Index

def Get_RNA_Index(rna_File):
    """
    This function reads a RNA-seq data file and creates a pandas DataFrame with additional column for transcript abundance.

    Parameters:
    - rna_File (str): The file path of the RNA-seq data file. It should be a CSV file with at least three columns: transcript identifiers, and two columns representing the abundance of the transcripts in two conditions.

    Returns:
    - pandas DataFrame: A DataFrame containing the RNA-seq data with additional column for transcript abundance. The transcript abundance is calculated as the average of the abundance values in the two conditions.

    Note:
    - The input file should be a CSV file with a header.
    - The function drops the columns representing the abundance of the transcripts in the two conditions.
    """
    rna_Index = pd.read_csv(rna_File)
    rna_Cols = rna_Index.columns
    rna_Index['Transcript Abundance'] = (rna_Index[rna_Cols[1]] + rna_Index[rna_Cols[2]]) / 2
    rna_Index = rna_Index.drop([rna_Cols[1], rna_Cols[2]], axis = 1)
    return rna_Index

def Get_Protein_Index(protein_File, cell_Type):
    """
    This function reads a protein expression data file, filters the data based on the given cell type,
    and calculates the average protein expression for each gene.

    Parameters:
    - protein_File (str): The file path of the protein expression data file. It should be a CSV file.
    - cell_Type (str): The name of the cell type for which the protein expression data needs to be filtered.

    Returns:
    - pandas DataFrame: A DataFrame containing the gene names, gene symbols, protein names, and average protein expressions.

    Note:
    - The input file should have a column named 'cell_line_display_name' to filter the data based on the cell type.
    - The protein expression values are assumed to be in columns starting from the 6th column.
    - The function calculates the average protein expression for each gene by taking the mean of the protein expression values for that gene.
    """
    prot_Data = pd.read_csv(protein_File)
    prot_Row = prot_Data[prot_Data['cell_line_display_name'] == cell_Type.upper()]
    prot_Index = pd.read_csv('Protein_Information.txt', sep = '\t')
    prot_Index['Protein Expression'] = np.nan
    if (prot_Row.size > 0):
        for i in prot_Row.columns[6:]:
            gene_Name = i.split(' ')[0]
            prot_Name = i.split(' ')[1][1:-1]
            if gene_Name in list(prot_Index['Gene name']):
                if prot_Name in list(prot_Index[prot_Index['Gene name'] == gene_Name]['UniProtKB Gene Name ID']):
                    prot_Index.loc[(prot_Index['Gene name'] == gene_Name) & (prot_Index['UniProtKB Gene Name ID'] == prot_Name),'Protein Expression'] = 2 ** float(prot_Row[i])
    prot_Index.columns = ['Gene Name', 'Gene Symbol', 'Protein Name', 'Protein Expression']
    prot_Index = prot_Index.groupby('Gene Name').agg({'Protein Expression':'mean'})
    prot_Index['Gene Name'] = prot_Index.index
    prot_Index.index = range(len(prot_Index))
    prot_Index = prot_Index.drop_duplicates()
    return prot_Index

def Count_Isoform(iso_File):
    """
    Count the number of isoforms for each gene in the given GTF file.

    Parameters:
    iso_File (str): The file path of the GTF file containing the transcript information.

    Returns:
    iso_List (list): A list of integers representing the number of isoforms for each gene.

    The function reads the GTF file line by line, extracts the transcript information,
    and counts the number of isoforms for each gene. The function returns a list of integers
    representing the number of isoforms for each gene.
    """
    gene_Index = Make_Gene_Index()
    gene_Index.index = range(len(gene_Index.index))
    with open (iso_File) as this_File:
        this_Text = this_File.readlines()
        this_File.close()
    iso_Counts = {}
    for i in gene_Index['Gene Name']:
        iso_Counts[i] = 0
    for i in range(len(this_Text)):
        line = this_Text[i].strip().split('\t')
        if (len(line) ==  9):
            if (line[2] == 'transcript'):
                trx_Info = line[8].split('"')
                if (len(trx_Info) == 9):
                    if (trx_Info[5] in iso_Counts.keys()):
                        iso_Counts[trx_Info[5]] += 1
    iso_List = []
    for i in iso_Counts:
        iso_List.append(iso_Counts[i])
    return iso_List

def Make_Profile(cell_type, atac_file, pop_file, iso_file = None, rna_file = None, prot_file = "DepMap_Proteomics.csv"):
    """
    This function generates a profile for a given cell line based on ATAC-seq, RNA-seq, and proteomics data.

    Parameters:
    - atac_file (str): The file path of the protien-DNA interaction data (ATAC-seq data).
    - pop_file (str): The file path of the protein-RNA interaction data (POP-seq data).
    - iso_file (str, optional): The file path of the isoform data. Defaults to None.
    - rna_file (str, optional): The file path of the RNA-seq data. Defaults to None.
    - prot_file (str, optional): The file path of the proteomics data. Defaults to "DepMap_Proteomics.csv".

    Returns:
    - None. However, it writes a CSV file named "<cell_type>_Piranha_Profile.csv" in a directory named "<cell_type>".

    Note:
    - The function uses the functions Get_Correlation(), Count_Isoform(), Get_RNA_Index(), and Get_Protein_Index() to generate the profile.
    """
    iso_list = []
    if iso_file:
        iso_list = Count_Isoform(iso_file)

    gene_index = Get_Correlation(atac_file, pop_file, upstream = arg_Vals.upstream, downstream = arg_Vals.downstream)
    #gene_index.to_csv(cell_type + "_Total_Gene_Profile.csv")
    gene_index['Isoform Count'] = iso_list

    comm_index = gene_index[(gene_index['Gene Sites'] > 0) & (gene_index['Transcript Sites'] > 0)]
    comm_index.index = range(len(comm_index.index))
    comm_index['Chromosome'] = comm_index['Chromosome'].replace({'X':'23','Y':'24','MT':'0'})
    comm_index['Gene Location'] = ((comm_index['Gene End'] + comm_index['Gene Start']) / 2) / comm_index['Chromosome Length']

    if rna_file:
        rna_index = Get_RNA_Index(rna_file)
        comm_index = comm_index.merge(rna_index, on='Gene Name')

    if prot_file:
        prot_index = Get_Protein_Index(prot_file, cell_type)
        comm_index = comm_index.merge(prot_index, on='Gene Name')

    comm_index = comm_index.drop_duplicates()
    if (arg_Vals.upstream != None):
        comm_index.to_csv(cell_type + "_Piranha_Profile_" + str(arg_Vals.upstream) + ".csv")
    else:
        comm_index.to_csv(cell_type + "_Piranha_Profile.csv")

# Parse command line arguments
arg_Parser = argparse.ArgumentParser()
arg_Parser.add_argument('cell_type', type = str, help = 'a cell line in the database')
arg_Parser.add_argument('atac_file', type = str, help = 'ATAC-seq data file path')
arg_Parser.add_argument('pop_file', type = str, help = 'POP-seq data file path')
arg_Parser.add_argument('--iso_file', type = str, help = 'Isoform data file path (optional)')
arg_Parser.add_argument('--rna_file', type = str, help = 'RNA-seq data file path (optional)')
arg_Parser.add_argument('--prot_file', type = str, help = 'Proteomics data file path (optional)')
arg_Parser.add_argument('--upstream', type = int, help = 'number of bases upstream of gene boundary to consider when evaluating ATAC peaks (optional)')
arg_Parser.add_argument('--downstream', type = int, help = 'number of bases downstream of gene boundary to consider when evaluating ATAC peaks (optional)')

arg_Vals = arg_Parser.parse_args()

# Call the Make_Profile function with command line arguments
Make_Profile(arg_Vals.cell_type, arg_Vals.atac_file, arg_Vals.pop_file, iso_file = arg_Vals.iso_file, rna_file = arg_Vals.rna_file, prot_file = arg_Vals.prot_file)

'''
prot_index = Get_Protein_Index('DepMap_Proteomics.csv', 'HepG2')

test_index = pd.read_csv('HepG2_Tot_Tot_Piranha_Profile.csv', index_col = 0)
test_index = test_index.drop('Protein Expression', axis = 1)

test_index = test_index.merge(prot_index, on='Gene Name')
test_index.to_csv('HepG2_Tot_Tot_Piranha_Profile.csv')

hmm_index = pd.read_csv('HepG2_Tot_Tot_HMM_Piranha_Profile.csv', index_col = 0)
hmm_index = hmm_index.drop('Protein Expression', axis = 1)

hmm_index = hmm_index.merge(prot_index, on='Gene Name')
hmm_index.to_csv('HepG2_Tot_Tot_HMM_Piranha_Profile.csv')

prot_index = Get_Protein_Index('DepMap_Proteomics.csv', 'K562')

test_index = pd.read_csv('K562_Tot_Tot_Piranha_Profile.csv', index_col = 0)
test_index = test_index.drop('Protein Expression', axis = 1)

test_index = test_index.merge(prot_index, on='Gene Name')
test_index.to_csv('K562_Tot_Tot_Piranha_Profile.csv')

hmm_index = pd.read_csv('K562_Tot_Tot_HMM_Piranha_Profile.csv', index_col = 0)
hmm_index = hmm_index.drop('Protein Expression', axis = 1)

hmm_index = hmm_index.merge(prot_index, on='Gene Name')
hmm_index.to_csv('K562_Tot_Tot_HMM_Piranha_Profile.csv')
'''