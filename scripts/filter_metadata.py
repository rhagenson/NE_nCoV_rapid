# -*- coding: utf-8 -*-

# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Release date: 2020-03-24
# Last update: 2021-07-07


import pycountry_convert as pyCountry
import pycountry
import argparse
from Bio import SeqIO
import pandas as pd
from epiweeks import Week
import time

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter nextstrain metadata files re-formmating and exporting only selected lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genomes", required=True, help="FASTA file genomes to be used")
    parser.add_argument("--metadata1", required=True, help="Metadata file from NextStrain")
    parser.add_argument("--metadata2", required=False, help="Custom lab metadata file")
    parser.add_argument("--filter", required=False, nargs='+', type=str,  help="List of filters for tagged rows in lab metadata")
    parser.add_argument("--output1", required=True, help="Filtered metadata file")
    parser.add_argument("--output2", required=True, help="Reformatted, final FASTA file")
    args = parser.parse_args()

    genomes = args.genomes
    metadata1 = args.metadata1
    metadata2 = args.metadata2
    filterby = args.filter
    output1 = args.output1
    output2 = args.output2

    # path =
    # genomes = path + 'pre-analyses/temp_sequences.fasta'
    # metadata1 = path + 'pre-analyses/metadata_nextstrain.tsv'
    # metadata2 = path + 'pre-analyses/SC2_Project_NE_DHHS.xlsx'
    # filterby =
    # output1 = path + 'pre-analyses/metadata_filtered.tsv'
    # output2 = path + 'pre-analyses/sequences.fasta'
    
    # temporal boundaries
    today = time.strftime('%Y-%m-%d', time.gmtime())
    min_date = '2019-12-15'
    
    variants = {'B.1.1.7': 'Alpha (B.1.1.7)',
                'B.1.351': 'Beta (B.1.351)',
                'B.1.351.2': 'Beta (B.1.351.2)',
                'B.1.351.3': 'Beta (B.1.351.3)',
                'P.1': 'Gamma (P.1)',
                'P.1.1': 'Gamma (P.1.1)',
                'P.1.2': 'Gamma (P.1.2)',
                'B.1.617.2': 'Delta (B.1.617.2)',
                'AY.1': 'Delta (AY.1)',
                'AY.2': 'Delta (AY.2)',
                'AY.3': 'Delta (AY.3)',
                'B.1.525': 'Eta (B.1.525)',
                'B.1.526': 'Iota (B.1.526)',
                'B.1.617.1': 'Kappa (B.1.617.1)',
                'C.37': 'Lambda (C.37)',
                'B.1.427': 'Epsilon (B.1.427/B.1.429)',
                'B.1.429': 'Epsilon (B.1.427/B.1.429)',
                'P.2': 'Zeta (P.2)'
                }
    
    
    # get ISO alpha3 country codes
    isos = {}
    def get_iso(country):
        global isos
        if country not in isos.keys():
            try:
                isoCode = pyCountry.country_name_to_country_alpha3(country, cn_name_format="default")
                isos[country] = isoCode
            except:
                try:
                    isoCode = pycountry.countries.search_fuzzy(country)[0].alpha_3
                    isos[country] = isoCode
                except:
                    isos[country] = ''
        return isos[country]

    # create epiweek column
    def get_epiweeks(date):
        date = pd.to_datetime(date)
        epiweek = str(Week.fromdate(date, system="cdc"))  # get epiweeks
        epiweek = epiweek[:4] + '_' + 'EW' + epiweek[-2:]
        return epiweek
    
    # add state code
    us_state_abbrev = {
        'Alabama': 'AL',
        'Alaska': 'AK',
        'American Samoa': 'AS',
        'Arizona': 'AZ',
        'Arkansas': 'AR',
        'California': 'CA',
        'Colorado': 'CO',
        'Connecticut': 'CT',
        'Delaware': 'DE',
        'District of Columbia': 'DC',
        'Washington DC': 'DC',
        'Florida': 'FL',
        'Georgia': 'GA',
        'Guam': 'GU',
        'Hawaii': 'HI',
        'Idaho': 'ID',
        'Illinois': 'IL',
        'Indiana': 'IN',
        'Iowa': 'IA',
        'Kansas': 'KS',
        'Kentucky': 'KY',
        'Louisiana': 'LA',
        'Maine': 'ME',
        'Maryland': 'MD',
        'Massachusetts': 'MA',
        'Michigan': 'MI',
        'Minnesota': 'MN',
        'Mississippi': 'MS',
        'Missouri': 'MO',
        'Montana': 'MT',
        'Nebraska': 'NE',
        'Nevada': 'NV',
        'New Hampshire': 'NH',
        'New Jersey': 'NJ',
        'New Mexico': 'NM',
        'New York': 'NY',
        'North Carolina': 'NC',
        'North Dakota': 'ND',
        'Northern Mariana Islands': 'MP',
        'Ohio': 'OH',
        'Oklahoma': 'OK',
        'Oregon': 'OR',
        'Pennsylvania': 'PA',
        'Puerto Rico': 'PR',
        'Rhode Island': 'RI',
        'South Carolina': 'SC',
        'South Dakota': 'SD',
        'Tennessee': 'TN',
        'Texas': 'TX',
        'Utah': 'UT',
        'Vermont': 'VT',
        'Virgin Islands': 'VI',
        'Virginia': 'VA',
        'Washington': 'WA',
        'West Virginia': 'WV',
        'Wisconsin': 'WI',
        'Wyoming': 'WY'
    }


    # nextstrain metadata
    dfN = pd.read_csv(metadata1, encoding='utf-8', sep='\t', dtype='str')
    dfN.insert(5, 'iso', '')
    dfN['category'] = ''
    dfN['batch'] = ''
    dfN['sequencing_date'] = ''
    dfN['group'] = ''
    dfN['Cluster_ID'] = ''
    dfN.fillna('', inplace=True)

     # add tag of variant category
    def variant_category(lineage):
        var_category = 'Other variants'
        for name in variants.keys():
            if lineage == name:
                var_category = variants[lineage]
        return var_category

    dfN['category'] = dfN['pango_lineage'].apply(lambda x: variant_category(x))
    
    list_columns = dfN.columns.values  # list of column in the original metadata file

    # Lab genomes metadata
    dfE = pd.read_excel(metadata2, index_col=None, header=0, sheet_name=0,
                        # 'sheet_name' must be changed to match the Excel sheet name
                        converters={'sample': str, 'sample_id': str, 'collection-date': str, 'category': str, 'batch': str, 'group': str, 'Cluster_ID': str, 'Filter': str})  # this need to be tailored to your lab's naming system
    dfE.fillna('', inplace=True)
    
    dfE = dfE.rename(columns={'sample_id': 'id', 'collection-date': 'date', 'lab': 'originating_lab', 'Filter': 'filter' })
    dfE['epiweek'] = ''
    
    # exclude rows with no ID
    if 'id' in dfE.columns.to_list():
        dfE = dfE[~dfE['id'].isin([''])]

    lab_sequences = dfE['id'].tolist()
    # exclude unwanted lab metadata row
    if len(filterby) > 0:
        print('\nFiltering metadata by category: ' + ', '.join(filterby) + '\n')
    dfL = pd.DataFrame(columns=dfE.columns.to_list())
    for value in filterby:
        dfF = dfE[dfE['filter'].isin([value])]  # batch inclusion of specific rows
        dfL = pd.concat([dfL, dfF]) # add filtered rows to dataframe with lab metadata

    # list of relevant genomes sequenced
    keep_only = dfL['id'].tolist()
    excluded = [id for id in lab_sequences if id not in keep_only]

    # create a dict of existing sequences
    sequences = {}
    for fasta in SeqIO.parse(open(genomes), 'fasta'):  # as fasta:
        id, seq = fasta.description, fasta.seq
        if id not in sequences.keys() and id not in excluded:
            sequences[id] = str(seq)

    # add inexistent columns
    for col in list_columns:
        if col not in dfL.columns:
            dfL[col] = ''


    # output dataframe
    outputDF = pd.DataFrame(columns=list_columns)
    found = []
    lab_label = {}
    metadata_issues = {}
    # process metadata from excel sheet
    for idx, row in dfL.iterrows():
        id = dfL.loc[idx, 'id']
        if id in sequences:
            dict_row = {}
            for col in list_columns:
                dict_row[col] = ''
                if col in row:
                    dict_row[col] = dfL.loc[idx, col]  # add values to dictionary
                    
             # check for missing geodata
            geodata = ['country'] # column
            for level in geodata:
                if len(dict_row[level]) < 1:
                    if id not in metadata_issues:
                        metadata_issues[id] = [level]
                    else:
                        metadata_issues[id].append(level)

            if dict_row['location'] in ['', None]:
                dict_row['location'] = dfL.loc[idx, 'location']

            collection_date = ''
            if len(str(dict_row['date'])) > 1:
                collection_date = dict_row['date'].split(' ')[0].replace('.', '-').replace('/', '-')
                dict_row['date'] = collection_date
                # check is date is appropriate: not from the 'future', not older than 'min_date'
                if pd.to_datetime(today) < pd.to_datetime(collection_date) or pd.to_datetime(min_date) > pd.to_datetime(collection_date):
                    if id not in metadata_issues:
                        metadata_issues[id] = ['date']
                    else:
                        metadata_issues[id].append('date')
            else: # missing date
                if id not in metadata_issues:
                    metadata_issues[id] = ['date']
                else:
                    metadata_issues[id].append('date')
            
            # fix exposure
            columns_exposure = ['country_exposure', 'division_exposure']
            for level_exposure in columns_exposure:
                level = level_exposure.split('_')[0]
                dict_row[level_exposure] = dfL.loc[idx, level_exposure]
                if dict_row[level_exposure] in ['', None]:
                    if level_exposure == 'country_exposure':
                        dict_row[level_exposure] = dict_row[level]
                    else:
                        if dict_row['country_exposure'] != dfL.loc[idx, 'country']:
                            dict_row[level_exposure] = dict_row['country_exposure']
                        else:
                            dict_row[level_exposure] = dict_row[level]        
                    
            if row['state'] == '':
                code = 'un'  # change this line to match the acronym of the most likely state of origin if the 'State' field is unknown
            else:
                code = row['state']
                
            strain = code + '/' + row['sample'] + '/' + id # new strain name

            dict_row['strain'] = strain
            dict_row['iso'] = get_iso(dict_row['country'])
            dict_row['originating_lab'] = dfL.loc[idx, 'originating_lab']
            dict_row['submitting_lab'] = 'Nebraska DHHS Sequencing Initiative'
            dict_row['authors'] = row['group']
            dict_row['batch'] = 'Batch' + str('0' * (3 - len(row['batch']))) + row['batch']
            dict_row['sequencing_date'] = row['sequencing-collection-date'].strftime('%Y-%m-%d')

            # add lineage
            lineage = ''
            if dfL.loc[idx, 'pango_lineage'] != '':
                lineage = dfL.loc[idx, 'pango_lineage']
            dict_row['pango_lineage'] = lineage

            # variant classication (VOI, VOC, VHC)
            dict_row['category'] = variant_category(lineage)
            
            # assign epiweek
            if len(dict_row['date']) > 0:
                dict_row['epiweek'] = get_epiweeks(collection_date)
            else:
                dict_row['epiweek'] = ''

            # record sequence and metadata as found
            found.append(strain)
            if id not in metadata_issues.keys():
                lab_label[id] = strain
                outputDF = outputDF.append(dict_row, ignore_index=True)
                
    # process metadata from TSV
    dfN = dfN[dfN['strain'].isin(sequences.keys())]
    for idx, row in dfN.iterrows():
        strain = dfN.loc[idx, 'strain']
        if strain in sequences:
            if strain in outputDF['strain'].to_list():
                continue
            dict_row = {}
            date = ''
            for col in list_columns:
                if col == 'date':
                    date = dfN.loc[idx, col]
                dict_row[col] = ''
                if col in row:
                    dict_row[col] = dfN.loc[idx, col]

            # fix exposure
            columns_exposure = ['country_exposure', 'division_exposure']
            for level_exposure in columns_exposure:
                level = level_exposure.split('_')[0]
                if dict_row[level_exposure] in ['', None]:
                    dict_row[level_exposure] = dict_row[level]

            dict_row['iso'] = get_iso(dict_row['country'])
            dict_row['epiweek'] = get_epiweeks(date)
            found.append(strain)

            outputDF = outputDF.append(dict_row, ignore_index=True)


    # write new metadata files
    outputDF = outputDF.drop(columns=['region'])
    outputDF.to_csv(output1, sep='\t', index=False)


    # write sequence file
    exported = []
    with open(output2, 'w') as outfile2:
        # export new metadata lines
        for id, sequence in sequences.items():
            if id in lab_label and id not in metadata_issues.keys(): # export lab generated sequences
                if lab_label[id] not in exported:
                    entry = '>' + lab_label[id] + '\n' + sequence + '\n'
                    outfile2.write(entry)
                    print('* Exporting newly sequenced genome and metadata for ' + id)
                    exported.append(lab_label[id])
            else:  # export publicly available sequences
                if id not in exported and id in outputDF['strain'].tolist():
                    entry = '>' + id + '\n' + sequence + '\n'
                    outfile2.write(entry)
                    exported.append(id)

    if len(metadata_issues) > 0:
        print('\n\n### WARNINGS!\n')
        print('\nPlease check for metadata issues related to these samples and column (which will be otherwise ignored)\n')
        for id, columns in metadata_issues.items():
            print('\t- ' + id + ' (issues found at: ' + ', '.join(columns) + ')')



    # # write fasta file
    # exported = []
    # print('\n### Exporting genomes and metadata')
    # print('\t Exporting all selected lab sequences, publicly available genomes and metadata')
    # with open(output2, 'w') as outfile2:
    #     # export new fasta entries
    #     for id, sequence in sequences.items():
    #         if len(id) < 5:
    #             if lab_label[id] not in exported:
    #                 entry = '>' + lab_label[id] + '\n' + sequence + '\n'
    #                 outfile2.write(entry)
    #                 print('\t\t* Newly sequenced genome and metadata: ' + id)
    #                 exported.append(lab_label[id])
    #         elif len(id) > 4:
    #             if id not in exported:
    #                 entry = '>' + id + '\n' + sequence + '\n'
    #                 outfile2.write(entry)
    #                 exported.append(id)

print('\nMetadata file successfully reformatted and exported!\n')
