#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 11:46:34 2022

@author: Xander
"""

import argparse, matplotlib.pyplot as plt, numpy as np, os, pandas as pd, scipy, seaborn as sns, statannot, cliffs_delta
os.path.join(os.path.dirname(__file__))

def Print_Gene_Lists(cell_Type, data_Frame, report, criteria, directory=os.getcwd()):
    if (report in data_Frame.columns) and (criteria in data_Frame.columns):
        if (directory != os.getcwd()):
            os.mkdir(directory)
            os.chdir('./' + directory)
        gene_Categories = list(set(data_Frame[criteria]))
        for i in gene_Categories:
            gene_List = data_Frame[data_Frame[criteria] == i][report]
            gene_List = pd.Series(list(set(gene_List)))
            gene_List.to_csv(cell_Type + '_' + 'Gene_List_' + str(i) + '.txt', header=False, index=False)
    else:
        print ('Error Message')

def Plot_Covariate(cell_Type, data_Frame, report, criteria, log_Axes=False):
    report_Title = {'Gene Length':'Gene Length', 'Isoform Count':'Number of Isoforms', 'Chromosome':'Distribution of Chromosome Numbers', 'Gene Location':'Gene Location', 'Transcript Abundance':'Abundance of Gene Product', 'Protein Expression':'Abundance of Protein Product'}
    data_Frame = data_Frame[data_Frame[report] != 0]
    cat_List = Get_Subsets(data_Frame, report, criteria)
    if (len(cat_List) == 2):
        try:
            p_Value = scipy.stats.mannwhitneyu(*cat_List).pvalue
        except ValueError:
            print ('No numeric values observed for selected feature.')
            return None
        criteria_List = sorted(set(data_Frame[criteria]), reverse = True)
    else:
        p_Value = scipy.stats.kruskal(*cat_List).pvalue
        if (p_Value is np.nan):
            print ('No numeric values observed for selected feature.')
            return None
        criteria_List = ['Predominantly Transcriptional', 'Balanced', 'Predominantly Post-Transcriptional']
    if (len(cat_List) == 2):
        if (np.median(cat_List[0]) >= np.median(cat_List[1])):
            print('Ratio of Medians:' , str(np.median(cat_List[0]) / np.median(cat_List[1])))
        else:
            print('Ratio of Medians:' , str(np.median(cat_List[1]) / np.median(cat_List[0])))
    else:
        print(data_Frame.groupby(criteria)[report].median())
    sns.set_style('white')
    fig = plt.figure(dpi = 900)
    ax = plt.subplot(111)
    if (log_Axes):
        data_Frame[report] = np.log10(data_Frame[report])
    if (type(data_Frame[criteria].unique()[0]) == str):
        if (len(cat_List) == 2):
            if (p_Value > 0):
                sns.violinplot(x = criteria, y = report, data = data_Frame, order = criteria_List, ax = ax).set(title = output_file + ' ' + report + ' by Regulation Group\nMann-Whitney p-Value: ' + str('{:.2e}'.format(float(p_Value))))
            else:
                sns.violinplot(x = criteria, y = report, data = data_Frame, order = criteria_List, ax = ax).set(title = output_file + ' ' + report + ' by Regulation Group\nMann-Whitney p-Value: < 1e-300')
        else:
            if (p_Value > 0):
                sns.violinplot(x = criteria, y = report, data = data_Frame, order = criteria_List, ax = ax).set(title = output_file + ' ' + report + ' by Regulation Group\nKruskal p-Value: ' + str('{:.2e}'.format(float(p_Value))))
            else:
                sns.violinplot(x = criteria, y = report, data = data_Frame, order = criteria_List, ax = ax).set(title = output_file + ' ' + report + ' by Regulation Group\nKruskal p-Value: < 1e-300')
    else:
        sns.violinplot(x = criteria, y = report, data = data_Frame, palette = sns.color_palette("Set2"), ax = ax)
    if (log_Axes):
        from matplotlib import ticker as mticker
        ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
        ymin, ymax = ax.get_ylim()
        tick_range = np.arange(np.floor(ymin), ymax)
        ax.yaxis.set_ticks(tick_range)
        ax.yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)], minor=True)
    if (report == 'Transcript Abundance'):
        plt.ylabel('Transcript Abundance (FPKM)')
    box_Pairs = Get_Pairs(criteria_List)
    statannot.add_stat_annotation(ax, data = data_Frame, x = criteria, y = report, order = criteria_List, box_pairs = box_Pairs, test = 'Mann-Whitney', text_format = 'star', pvalue_thresholds = [[1, ''], [0.05, '*'], [10 ** -5, '**'], [10 ** -10, '***'], [10 ** -50, '****'], [10 ** -100, '*****']], fontsize = 6)
    text_String = 'Significance level\n'
    text_String += '*' + ' : p < 0.05\n'
    text_String += '**' + ' : p < 1e-5\n'
    text_String += '***' + ' : p < 1e-10\n'
    text_String += '****' + ' : p < 1e-50\n'
    text_String += '*****' + ' : p < 1e-100\n'
    props = dict(boxstyle = 'round', facecolor = 'white', alpha = 0.5)
    ax.text(1.025, 0, text_String, transform = ax.transAxes, bbox = props, fontsize = 12)
    plt.savefig(cell_Type + '_' + report + '.png', bbox_inches = 'tight')
    plt.close(fig)
      
def Get_Subsets(data_Frame, report, criteria):
    if (report in data_Frame.columns) and (criteria in data_Frame.columns):
        set_Title = set(data_Frame[criteria])
        subset_List = []
        for i in set_Title:
            subset_List.append(data_Frame[(data_Frame[criteria] == i) & (data_Frame[report] != 0)][report].dropna())
        return subset_List
    else:
        print('Error')

def Find_Group_Indices(data_Set, n_Groups):
    n_Tiles = []
    for i in range(1, n_Groups):
        n_Tiles.append(1 / n_Groups * i)
    bin_Vals = scipy.stats.mstats.mquantiles(data_Set, prob = n_Tiles)
    bin_Vals = np.insert(bin_Vals, 0, 0)
    bin_Vals = np.append(bin_Vals, max(data_Set) + 0.01)
    return bin_Vals

def Find_Group_Lists(data_Set, group_Index):
    group_Lists = []
    for i in range(len(group_Index) - 1):
        this_Group = data_Set.loc[(data_Set['Regulation Ratio'] > group_Index[i]) & (data_Set['Regulation Ratio'] <= group_Index[i + 1])]
        group_Lists.append(this_Group['Gene Name'])
    return group_Lists

def Print_Lists(group_List, cell_Line, group_Category):
    for i in range(len(group_List)):
        group_List[i].to_csv(path_or_buf = './' + cell_Line + '/' + cell_Line + '_' + max(group_Category['Regulation Group']) + '_' + str(i + 1) + '.csv', header = False, index = False)

def Print_Genes(gene_List, cell_Type, file_Title):
    with open ('./' + cell_Type + '/' + cell_Type + '_' + file_Title + '.txt', 'w') as this_File:
        for i in gene_List:
            this_File.write(i + '\n')
        this_File.close()

def Make_Gene_Index(file='Gene_Information.txt'):
    gene_Index = pd.read_csv(file, sep='\t', header=0)
    gene_Index.columns = ['Gene Name', 'Gene Symbol', 'Chromosome', 'Gene Start', 'Gene End', 'Strand']
    gene_Index['Gene Start'] = gene_Index['Gene Start'].astype(int)
    gene_Index['Gene End'] = gene_Index['Gene End'].astype(int)
    gene_Index['Chromosome'] = gene_Index['Chromosome'].apply(lambda x: str(x))
    gene_Index['Gene Length'] = gene_Index['Gene End'] - gene_Index['Gene Start']
    gene_Index = gene_Index.dropna(axis='index')
    gene_Index.index = range(len(gene_Index.index))
    return gene_Index

def Merge_Profiles(file_List):
    if (type(file_List) == list):
        file_Merge = pd.DataFrame()
        for i in file_List:
            cell_Type = i.split('/')[1]
            this_Frame = pd.read_csv(i, index_col= 0)
            this_Frame['Cell Line'] = cell_Type
            file_Merge = pd.concat([file_Merge, this_Frame], axis = 0)
        return file_Merge
    else:
        print ('Error')

def Categorical_Graph(data_Frame, category, group = None, threshold = 1, title = ''):
    cat_List = list(set(data_Frame[category]))
    if (group != None):
        sig_Results = Transcript_Type_Enrichment_Test(data_Frame, group)
        group_List = ['Predominantly Transcriptional','Balanced', 'Predominantly Post-Transcriptional']
        group_Frame = pd.DataFrame()
        group_Frame[category] = cat_List
        for i in range(len(group_List)):
            value_Counts = data_Frame[data_Frame[group] == group_List[i]][category].value_counts()
            value_Counts = value_Counts / sum(value_Counts) * 100
            group_Frame = group_Frame.join(value_Counts, on = category, rsuffix = '_' + str(i))
        group_Frame.columns = [category, *group_List]
        group_Frame = group_Frame.melt(id_vars = category, value_vars = group_List)
        group_Frame.columns = [category, group, 'Percentage']
        group_Frame = pd.merge(group_Frame, sig_Results, on = [category, group], how = 'left')
        group_Frame = group_Frame.dropna()
    else:
        group_Frame = pd.DataFrame()
        group_Frame[category] = cat_List
        value_Counts = data_Frame[category].value_counts()
        value_Counts = value_Counts / sum(value_Counts) * 100
        group_Frame = group_Frame.join(value_Counts, on = category)
        group_Frame.columns = [category, 'Percentage']
    group_Frame[category] = group_Frame[category].str.replace('_', ' ')
    
    group_Frame = group_Frame[group_Frame['Percentage'] > threshold]

    for i in pd.unique(group_Frame[category]):
        for j in pd.unique(group_Frame[group]):
            if (i not in group_Frame.loc[group_Frame[group] == j][category].tolist()):
                group_Frame.loc[max(group_Frame.index) + 1] = [i, j, 0, 1, 0, '']
    group_Frame = group_Frame.sort_values(by = [category])
    print(group_Frame[['Transcript type', 'Gene Regulation Group']])
    print(group_Frame[['Transcript type', 'P-Value']])
    sns.set_style('white')
    fig = plt.figure(dpi = 900)
    ax = plt.subplot(111)
    if (group != None):
        sns.barplot(data = group_Frame, x = 'Percentage', y = category, hue = group, ax = ax, hue_order = ['Predominantly Transcriptional','Balanced', 'Predominantly Post-Transcriptional']).set(title = title + ' Gene Type Enrichment')
    else:
        sns.barplot(data = group_Frame, x = 'Percentage', y = category, ax = ax).set(title = title + ' Gene Type Enrichment')
    if (group != None):
        counter = 0
        for container in ax.containers:
            ax.bar_label(container, labels = group_Frame.loc[group_Frame[group] == group_List[counter], 'Significance'], fontsize = 6)
            print(group_Frame.loc[group_Frame[group] == group_List[counter]])
            counter += 1
    plt.legend(bbox_to_anchor=(1,0.6))
    if (group != None):
        text_String = 'Significance level\n'
        text_String += '*' + ' : p < 0.05\n'
        text_String += '**' + ' : p < 1e-5\n'
        text_String += '***' + ' : p < 1e-10\n'
        text_String += '****' + ' : p < 1e-50\n'
        text_String += '*****' + ' : p < 1e-100\n'
        props = dict(boxstyle = 'round', facecolor = 'white', alpha = 0.5)
        ax.text(1.025, 0, text_String, transform = ax.transAxes, bbox = props, fontsize = 12)
    plt.savefig(title + '_Categorical.png', bbox_inches = 'tight')
    plt.close()

def Transcript_Type_Enrichment_Test(data_Frame, group):
    test_Results = pd.DataFrame(columns = [group, 'Transcript type', 'P-Value', 'Odds Ratio'])
    cat_List = list(set(data_Frame[group]))
    background_Data = pd.read_csv('Gene_Type_Data.txt', sep = '\t')['Transcript type'].value_counts()
    for i in cat_List:
        this_Group = data_Frame[data_Frame[group] == i]
        this_Data = this_Group['Transcript type'].value_counts()
        for j in this_Data.index:
            test_Result = scipy.stats.fisher_exact([[this_Data[j], background_Data[j]], [this_Data.sum() - this_Data[j], background_Data.sum() - background_Data[j]]], alternative = 'two-sided')
            test_Results.loc[len(test_Results)] = [i, j, test_Result.pvalue, test_Result.statistic]
    test_Results['Significance'] = ''
    test_Results.loc[test_Results['P-Value'] < 0.05, 'Significance'] = '*'
    test_Results.loc[test_Results['P-Value'] < 10 ** -5, 'Significance'] = '**'
    test_Results.loc[test_Results['P-Value'] < 10 ** -10, 'Significance'] = '***'
    test_Results.loc[test_Results['P-Value'] < 10 ** -50, 'Significance'] = '****'
    test_Results.loc[test_Results['P-Value'] < 10 ** -100, 'Significance'] = '*****'
    test_Results.to_csv(output_file + " Transcript Type Enrichment.csv")
    return test_Results

def Get_Pairs(this_List):
    pairs = []
    for i in range(len(this_List)):
        for j in range(i + 1, len(this_List)):
            pairs.append((this_List[i], this_List[j]))
    return pairs

# Parse command line arguments
arg_Parser = argparse.ArgumentParser()
arg_Parser.add_argument('profile', type = str, help = 'a profile file to analyze')
arg_Parser.add_argument('output', type = str, help = 'a name for the output file')

arg_Vals = arg_Parser.parse_args()

gene_Index = pd.read_csv(arg_Vals.profile, index_col= 0)

gene_Index.loc[gene_Index['Gene Regulation Group'] == "Neutral", "Gene Regulation Group"] = "Balanced"

output_file = arg_Vals.output

print(output_file + '\t\t\t\t\t ' + str(len(gene_Index)))
print(gene_Index['Gene Regulation Group'].value_counts())
    
post_Index = gene_Index[gene_Index['Gene Regulation Group'] == 'Predominantly Post-Transcriptional']
neutral_Index = gene_Index[gene_Index['Gene Regulation Group'] == 'Balanced']
trans_Index = gene_Index[gene_Index['Gene Regulation Group'] == 'Predominantly Transcriptional']
    
color_palette = ['#66c2a5', '#fc8d62']

sns.set_palette(color_palette)

p_Value = scipy.stats.mannwhitneyu(gene_Index['Gene Regulation'], gene_Index['Transcript Regulation']).pvalue
diff_Data = pd.melt(gene_Index[['Gene Regulation','Transcript Regulation']])
diff_Data.replace('Gene Regulation', 'Transcriptional Regulation', inplace = True)
diff_Data.replace('Transcript Regulation', 'Post-Transcriptional Regulation', inplace = True)
diff_Data.columns = ['Regulatory Density', 'Peaks per Gene Length in kb']
fig = plt.figure(dpi = 900)
ax = plt.subplot(111)
if (p_Value > 0):
    sns.boxplot(data = diff_Data, x = 'Regulatory Density', y = 'Peaks per Gene Length in kb', ax = ax).set(title =  output_file + ' Gene Regulation Densities' +'\nMann-Whitney p-Value: ' + str('{:.2e}'.format(float(p_Value))))
else:
    sns.boxplot(data = diff_Data, x = 'Regulatory Density', y = 'Peaks per Gene Length in kb', ax = ax).set(title = output_file + ' Gene Regulation Densities' + '\nMann-Whitney p-Value: < 1e-300')
ax.set_yscale('log')
plt.legend(bbox_to_anchor=(1,0.5), loc='center left')
plt.savefig(output_file + '_Regulation_Density.png', bbox_inches = 'tight')
plt.close(fig)

color_palette = ['#66c2a5', '#e78ac3', '#fc8d62']

sns.set_palette(color_palette)

r_Value = scipy.stats.linregress(gene_Index['Gene Regulation'], gene_Index['Transcript Regulation']).rvalue
r_Sqd = round(r_Value ** 2, 3)
fig = plt.figure(dpi = 900)
ax = plt.subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')
sns.scatterplot(data=gene_Index, x='Gene Regulation', y='Transcript Regulation', marker = '.', hue='Gene Regulation Group', hue_order=['Predominantly Transcriptional','Balanced', 'Predominantly Post-Transcriptional']).set(title = output_file + ' Regulation Ratio' + '\nPearson R-squared: ' + str(r_Sqd))
plt.xlabel('-Log10 DNA accessible regions (ATAC-seq sites) / kb')
plt.ylabel('-Log10 Protein RNA binding sites (POP-seq sites) / kb')
plt.yticks(ticks = [10 ** -6, 10 ** -5, 10 ** -4, 10 ** -3, 10 ** -2, 10 ** -1], labels = [6, 5, 4, 3, 2, 1])
plt.xticks(ticks = [10 ** -6, 10 ** -5, 10 ** -4, 10 ** -3, 10 ** -2, 10 ** -1], labels = [6, 5, 4, 3, 2, 1])
plt.legend(bbox_to_anchor=(1,0.5), loc='center left')
plt.savefig(output_file + '_Regulation_Ratio.png', bbox_inches = 'tight')
#plt.show()
plt.close(fig)

print("Trx vs Post", gene_Index[gene_Index['Gene Regulation Group'] == 'Predominantly Transcriptional']['Gene Regulation'].median() / gene_Index[gene_Index['Gene Regulation Group'] == 'Predominantly Post-Transcriptional']['Transcript Regulation'].median())
print("Post vs Trx", gene_Index[gene_Index['Gene Regulation Group'] == 'Predominantly Post-Transcriptional']['Transcript Regulation'].median() / gene_Index[gene_Index['Gene Regulation Group'] == 'Predominantly Transcriptional']['Gene Regulation'].median())

type_Frame = pd.read_csv('Gene_Type_Data.txt', sep = '\t')
type_Frame = type_Frame.drop_duplicates()
merge_Frame = pd.merge(gene_Index, type_Frame, on = 'Gene Name', how = 'left')

#Gene Length vs Gene Sites
r_Value, p_Value = scipy.stats.spearmanr(merge_Frame['Gene Length'], merge_Frame['Gene Sites'])
r_Sqd = round(r_Value ** 2, 3)

fig = plt.figure(dpi = 900)
ax = plt.subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')
if (p_Value > 0):
    sns.scatterplot(data=merge_Frame, x='Gene Length', y='Gene Sites', marker = '.', hue='Gene Regulation Group', hue_order=['Predominantly Transcriptional', 'Balanced', 'Predominantly Post-Transcriptional']).set(title = output_file +'\nSpearman R-squared: ' + str(r_Sqd) + '\np-value: ' + str('{:.2e}'.format(float(p_Value))))
else:
    sns.scatterplot(data=merge_Frame, x='Gene Length', y='Gene Sites', marker = '.', hue='Gene Regulation Group', hue_order=['Predominantly Transcriptional', 'Balanced', 'Predominantly Post-Transcriptional']).set(title = output_file + '\nSpearman R-squared: ' + str(r_Sqd) + '\np-value: < 1e-300')
plt.xlabel('Gene Length (bp)')
plt.ylabel('ATAC-Seq Peaks')
plt.legend(bbox_to_anchor=(1,0.5), loc='center left')
plt.savefig(output_file + '_Gene_Length_vs_Gene_Peaks.png', bbox_inches = 'tight')
plt.close(fig)

#Gene Length vs Transcript Sites
r_Value, p_Value = scipy.stats.spearmanr(merge_Frame['Gene Length'], merge_Frame['Transcript Sites'])
r_Sqd = round(r_Value ** 2, 3)

fig = plt.figure(dpi = 900)
ax = plt.subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')
if (p_Value > 0):
    sns.scatterplot(data=merge_Frame, x='Gene Length', y='Transcript Sites', marker = '.', hue='Gene Regulation Group', hue_order=['Predominantly Transcriptional', 'Balanced', 'Predominantly Post-Transcriptional']).set(title = output_file + '\nSpearman R-squared: ' + str(r_Sqd) + '\np-value: ' + str('{:.2e}'.format(float(p_Value))))
else:
    sns.scatterplot(data=merge_Frame, x='Gene Length', y='Transcript Sites', marker = '.', hue='Gene Regulation Group', hue_order=['Predominantly Transcriptional', 'Balanced', 'Predominantly Post-Transcriptional']).set(title = output_file + '\nSpearman R-squared: ' + str(r_Sqd) + '\np-value: < 1e-300')
plt.xlabel('Gene Length (bp)')
plt.ylabel('POP-Seq Peaks')
plt.legend(bbox_to_anchor=(1,0.5), loc='center left')
plt.savefig(output_file + '_Gene_Length_vs_Transcript_Sites.png', bbox_inches = 'tight')
plt.close(fig)

#Transcript Abundance vs Transcript Sites
r_Value, p_Value = scipy.stats.spearmanr(merge_Frame['Transcript Abundance'], merge_Frame['Transcript Sites'])
r_Sqd = round(r_Value ** 2, 3)

fig = plt.figure(dpi = 900)
ax = plt.subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')
if (p_Value > 0):
    sns.scatterplot(data=merge_Frame, x='Transcript Abundance', y='Transcript Sites', marker = '.', hue='Gene Regulation Group', hue_order=['Predominantly Transcriptional', 'Balanced', 'Predominantly Post-Transcriptional']).set(title = output_file + '\nSpearman R-squared: ' + str(r_Sqd) + '\np-value: ' + str('{:.2e}'.format(float(p_Value))))
else:
    sns.scatterplot(data=merge_Frame, x='Transcript Abundance', y='Transcript Sites', marker = '.', hue='Gene Regulation Group', hue_order=['Predominantly Transcriptional', 'Balanced', 'Predominantly Post-Transcriptional']).set(title = output_file + '\nSpearman R-squared: ' + str(r_Sqd) + '\np-value: < 1e-300')
plt.xlabel('Transcipt Abundance (FPKM)')
plt.ylabel('POP-Seq Peaks')
plt.legend(bbox_to_anchor=(1,0.5), loc='center left')
plt.savefig(output_file + '_Transcript_Abundance_vs_Transcript_Sites.png', bbox_inches = 'tight')
#plt.show()
plt.close(fig)

Categorical_Graph(merge_Frame, 'Transcript type', group = 'Gene Regulation Group', title = str(output_file))
#diff_Index = gene_Index[gene_Index['Gene Regulation Group'] != 'Balanced']
#Print_Gene_Lists(cell_line, gene_Index, 'Gene Name', 'Gene Regulation Group')

#Gene Length Comparison
Plot_Covariate(output_file, gene_Index, 'Gene Length', 'Gene Regulation Group', log_Axes=True)

print(scipy.stats.mannwhitneyu(trans_Index['Gene Length'], neutral_Index['Gene Length']))
print("Effect Size:", cliffs_delta.cliffs_delta(trans_Index['Gene Length'], neutral_Index['Gene Length']))
print(scipy.stats.mannwhitneyu(neutral_Index['Gene Length'], post_Index['Gene Length']))
print("Effect Size:", cliffs_delta.cliffs_delta(neutral_Index['Gene Length'], post_Index['Gene Length']))
print(scipy.stats.mannwhitneyu(trans_Index['Gene Length'], post_Index['Gene Length']))
print("Effect Size:", cliffs_delta.cliffs_delta(trans_Index['Gene Length'], post_Index['Gene Length']))
#Plot_Covariate(output_file, diff_Index, 'Gene Length', 'Gene Regulation Group', log_Axes=True)
   
#Isoform Count Comparison
Plot_Covariate(output_file, gene_Index, 'Isoform Count', 'Gene Regulation Group')
print(scipy.stats.mannwhitneyu(trans_Index['Isoform Count'], neutral_Index['Isoform Count']))
print("Effect Size:", cliffs_delta.cliffs_delta(trans_Index['Isoform Count'], neutral_Index['Isoform Count']))
print(scipy.stats.mannwhitneyu(neutral_Index['Isoform Count'], post_Index['Isoform Count']))
print("Effect Size:", cliffs_delta.cliffs_delta(neutral_Index['Isoform Count'], post_Index['Isoform Count']))
print(scipy.stats.mannwhitneyu(trans_Index['Isoform Count'], post_Index['Isoform Count']))
print("Effect Size:", cliffs_delta.cliffs_delta(trans_Index['Isoform Count'], post_Index['Isoform Count']))
#Plot_Covariate(output_file, diff_Index, 'Isoform Count', 'Gene Regulation Group')

#Gene Position
#Plot_Covariate(output_file, gene_Index, 'Chromosome', 'Gene Regulation Group')
#print(scipy.stats.mannwhitneyu(trans_Index['Chromosome'], neutral_Index['Chromosome']))
#print(scipy.stats.mannwhitneyu(neutral_Index['Chromosome'], post_Index['Chromosome']))
#print(scipy.stats.mannwhitneyu(trans_Index['Chromosome'], post_Index['Chromosome']))
#Plot_Covariate(output_file, diff_Index, 'Chromosome', 'Gene Regulation Group')
    
#Gene Location
#Plot_Covariate(output_file, gene_Index, 'Gene Location', 'Gene Regulation Group')
#print(scipy.stats.mannwhitneyu(trans_Index['Gene Location'], neutral_Index['Gene Location']))
#print(scipy.stats.mannwhitneyu(neutral_Index['Gene Location'], post_Index['Gene Location']))
#print(scipy.stats.mannwhitneyu(trans_Index['Gene Location'], post_Index['Gene Location']))
#Plot_Covariate(output_file, diff_Index, 'Gene Location', 'Gene Regulation Group')
    
#Transcript Abundance
Plot_Covariate(output_file, gene_Index, 'Transcript Abundance', 'Gene Regulation Group', log_Axes=True)
print(scipy.stats.mannwhitneyu(trans_Index['Transcript Abundance'], neutral_Index['Transcript Abundance']))
print("Effect Size:", cliffs_delta.cliffs_delta(trans_Index['Transcript Abundance'], neutral_Index['Transcript Abundance']))
print(scipy.stats.mannwhitneyu(neutral_Index['Transcript Abundance'], post_Index['Transcript Abundance']))
print("Effect Size:", cliffs_delta.cliffs_delta(neutral_Index['Transcript Abundance'], post_Index['Transcript Abundance']))
print(scipy.stats.mannwhitneyu(trans_Index['Transcript Abundance'], post_Index['Transcript Abundance']))
print("Effect Size:", cliffs_delta.cliffs_delta(trans_Index['Transcript Abundance'], post_Index['Transcript Abundance']))
#Plot_Covariate(output_file, diff_Index, 'Transcript Abundance', 'Gene Regulation Group', log_Axes=True)
    
#Protein Abundance
Plot_Covariate(output_file, gene_Index, 'Protein Expression', 'Gene Regulation Group', log_Axes=True)
print(scipy.stats.mannwhitneyu(trans_Index['Protein Expression'], neutral_Index['Protein Expression']))
print("Effect Size:", cliffs_delta.cliffs_delta(trans_Index['Protein Expression'], neutral_Index['Protein Expression']))
print(scipy.stats.mannwhitneyu(neutral_Index['Protein Expression'], post_Index['Protein Expression']))
print("Effect Size:", cliffs_delta.cliffs_delta(neutral_Index['Protein Expression'], post_Index['Protein Expression']))
print(scipy.stats.mannwhitneyu(trans_Index['Protein Expression'], post_Index['Protein Expression']))
print("Effect Size:", cliffs_delta.cliffs_delta(trans_Index['Protein Expression'], post_Index['Protein Expression']))
#Plot_Covariate(output_file, diff_Index, 'Protein Expression', 'Gene Regulation Group', log_Axes=True)
