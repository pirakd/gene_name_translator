import pandas as pd
import pickle as pl
import tqdm
from os import path


class GeneTranslator:
    def __init__(self, verbosity=True):
        self.root_folder = path.dirname(path.realpath(__file__))
        self.raw_data_file = 'hgnc_complete_set.txt'
        self.dictionary_file_name = 'gene_dictionary.pl'
        self.raw_data_file_path = path.join(self.root_folder, self.raw_data_file)
        self.dictionary_file_path = path.join(self.root_folder, self.dictionary_file_name)
        self.old_entrez_id_file_path = path.join(self.root_folder, 'old_entrez.txt')
        self.dictionary = None
        self.verbosity = verbosity
        self.load_old_enterz_ids = True

    def translate(self, query, query_type, return_type):
        keys_not_found = list()
        targets_not_found = list()
        result_dict = dict()
        dictionary = self.dictionary[query_type]
        for q in query:
            if q in dictionary:
                result = dictionary[q][return_type]
                if result is None:
                    targets_not_found.append(q)
                result_dict[q] = dictionary[q][return_type]
            else:
                keys_not_found.append(q)

        if len(keys_not_found):
            print('{} genes were not found ({}))'.format(len(keys_not_found), keys_not_found))
        if len(targets_not_found):
            print('{} translations are missing ({})'.format(len(targets_not_found), targets_not_found))

        return result_dict

    def load_dictionary(self):
        assert path.isfile(self.dictionary_file_path), \
            'Gene dictionary is missing! call GeneTranslator.generate_dictionaries first'

        with open(self.dictionary_file_path, 'rb') as f:
            self.dictionary = pl.load(f)
        if self.verbosity:
            print('available query types : {}'.format([x for x in self.dictionary.keys()]))

    def generate_dictionaries(self, keys=None):

        if keys is None:
            keys = ['symbol', 'entrez_id', 'uniprot', 'ensemble_gene_id']
        dictionaries = dict()
        for key in keys:
            dictionaries[key] = self.generate_dictionary(key)

        if self.load_old_enterz_ids:
            dictionaries = self.load_old_entrez_ids_from_file(dictionaries)

        with open(self.dictionary_file_path, 'wb') as f:
            pl.dump(dictionaries, f)


        return dictionaries

    def generate_dictionary(self, key):
        name_df = pd.read_csv(self.raw_data_file_path, delimiter='\t', index_col=False, low_memory=False)

        dictionary = dict()
        row_values = dict()

        for idx, row in tqdm.tqdm(name_df.iterrows(), desc='Generating dictionary with {} keys'.format(key),
                                  total=len(name_df)):
            row_values['symbol'] = row.symbol
            row_values['uniprot'] = row.uniprot_ids if not pd.isnull(row.uniprot_ids) else None
            row_values['entrez_id'] = int(row.entrez_id) if not pd.isnull(row.entrez_id) else None
            row_values['ensemble_gene_id'] = row.ensembl_gene_id if not pd.isnull(row.ensembl_gene_id) else None
            if not pd.isnull(row_values[key]):
                symbol = row.symbol
                aliases = []
                prev_symbols = []
                if not pd.isnull(row.alias_symbol):
                    aliases = row.alias_symbol.split('|')
                if not pd.isnull(row.prev_symbol):
                    prev_symbols = row.prev_symbol.split('|')

                aliases = list(set([symbol] + aliases + prev_symbols))

                query_dict = dict(entrez_id=row_values['entrez_id'], uniprot=row_values['uniprot'],
                                  ensembl_gene_id=row_values['ensemble_gene_id'],
                                  symbol=row_values['symbol'], symbol_aliases=aliases)

                if key == 'symbol':
                    for alias in aliases:
                        if alias == symbol or alias not in dictionary:
                            dictionary[alias] = query_dict
                else:
                    dictionary[row_values[key]] = query_dict
        return dictionary

    def load_old_entrez_ids_from_file(self, dictionaries):
        with open(self.old_entrez_id_file_path, 'r') as f:
            lines = f.readlines()
        for line in lines[1:]:
            line = line.strip()
            official_name, old_name = line.split('\t')
            dictionaries['entrez_id'][int(old_name)] = dictionaries['entrez_id'][int(official_name)]

        return dictionaries