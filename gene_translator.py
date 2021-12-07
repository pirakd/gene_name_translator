import pandas as pd
import pickle as pl
import tqdm
from os import path
from collections import defaultdict


class GeneTranslator:
    def __init__(self, verbosity=True):
        self.root_folder = path.dirname(path.realpath(__file__))
        self.raw_data_file = 'Homo_sapiens.gene_info'
        self.dictionary_file_name = 'gene_dictionary.pl'
        self.old_names_file = 'gene_history.txt'
        self.raw_data_file_path = path.join(self.root_folder, self.raw_data_file)
        self.dictionary_file_path = path.join(self.root_folder, self.dictionary_file_name)
        self.old_names_file_path = path.join(self.root_folder, self.old_names_file)
        self.dictionary = None
        self.old_names_mapping = None
        self.verbosity = verbosity
        self.use_old_names = True
        self.new_to_old_keys_mapping = {'symbol': 'obsolete_symbol',
                               'entrez_id': 'obsolete_entrez_id'}

    def translate(self, query, query_type, return_type):
        keys_not_found = list()
        targets_not_found = list()
        result_dict = dict()
        dictionary = self.dictionary[query_type]
        i = 0
        for q in query:
            if q in dictionary:
                result = dictionary[q][return_type]
                if result is None:
                    targets_not_found.append(q)
            else:
                result = self.query_old_name(q, query_type, return_type)
                i+=1
                if result is None:
                    keys_not_found.append(q)

            result_dict[q] = result

        if len(keys_not_found):
            print('{} genes were not found ({}))'.format(len(keys_not_found), keys_not_found))
        if len(targets_not_found):
            print('{} translations are missing ({})'.format(len(targets_not_found), targets_not_found))

        return result_dict

    def load_dictionary(self):
        assert path.isfile(self.dictionary_file_path), \
            'Gene dictionary is missing! call GeneTranslator.generate_dictionaries first'

        with open(self.dictionary_file_path, 'rb') as f:
            save_dict = pl.load(f)
        self.dictionary = save_dict['dictionary']
        self.old_names_mapping = save_dict['old_names_mapping']
        if self.verbosity:
            print('available query types : {}'.format([x for x in self.dictionary.keys()]))

    def generate_dictionaries(self, keys=None):

        keys = ['symbol', 'entrez_id']
        dictionaries = dict()
        for key in keys:
            dictionaries[key] = self.generate_dictionary(key)

        old_names_map = self.load_old_names(dictionaries)

        save_dict = {'dictionary': dictionaries,
                     'old_names_mapping': dict(old_names_map)}
        with open(self.dictionary_file_path, 'wb') as f:
            pl.dump(save_dict, f)

        return dictionaries

    def generate_dictionary(self, key):
        name_df = pd.read_csv(self.raw_data_file_path, delimiter='\t', index_col=False, low_memory=False)

        dictionary = dict()
        row_values = dict()

        for idx, row in tqdm.tqdm(name_df.iterrows(), desc='Generating dictionary with {} keys'.format(key),
                                  total=len(name_df)):
            row_values['symbol'] = row.symbol
            row_values['entrez_id'] = int(row.entrez_id) if not pd.isnull(row.entrez_id) else None
            if not pd.isnull(row_values[key]):
                symbol = row.symbol
                aliases = []
                if not pd.isnull(row.alias_symbol):
                    aliases = row.alias_symbol.split('|')

                aliases = list(set([symbol] + aliases))

                query_dict = dict(entrez_id=row_values['entrez_id'],
                                  symbol=row_values['symbol'], symbol_aliases=aliases)

                if key == 'symbol':
                    for alias in aliases:
                        if alias == symbol or alias not in dictionary:
                            dictionary[alias] = query_dict
                else:
                    dictionary[row_values[key]] = query_dict
        return dictionary

    def load_old_names(self, dictionaries):
        keys = ['entrez_id', 'obsolete_entrez_id',	'obsolete_symbol']

        old_names_dictionary = {key:dict() for key in keys}

        with open(self.old_names_file_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                entrez_id, old_entrez_id , old_symbol = line.strip().split('\t')
                entrez_id, old_entrez_id = int(entrez_id) if entrez_id != '-' else None, int(old_entrez_id),
                if entrez_id == '-':
                    entrez_id = None
                else:
                    if entrez_id not in old_names_dictionary[keys[0]]:
                        old_names_dictionary[keys[0]][entrez_id] = defaultdict(list)
                    old_names_dictionary[keys[0]][entrez_id][keys[1]].append(old_entrez_id)
                    old_names_dictionary[keys[0]][entrez_id][keys[2]].append(old_symbol)

                if old_entrez_id not in old_names_dictionary[keys[1]]:
                    old_names_dictionary[keys[1]][old_entrez_id] = defaultdict(list)
                if old_symbol not in old_names_dictionary[keys[2]]:
                    old_names_dictionary[keys[2]][old_symbol] = defaultdict(list)

                # if there is a mapping to an official gene
                if entrez_id:
                    old_names_dictionary[keys[1]][old_entrez_id][keys[0]].append(entrez_id)
                    old_names_dictionary[keys[2]][old_symbol][keys[0]].append(entrez_id)
                else:
                    old_names_dictionary[keys[1]][old_entrez_id][keys[2]].append(old_symbol)
                    old_names_dictionary[keys[2]][old_symbol][keys[1]].append(old_entrez_id)

        return old_names_dictionary

    def query_old_name(self, q, query_type, return_type):
        query = self.old_names_mapping[self.new_to_old_keys_mapping[query_type]].get(q, None)
        if query is not None:
            if 'entrez_id' in query:
                result = self.dictionary['entrez_id'].get(query['entrez_id'][0], None)
                if result:
                    return result[return_type]
            return query[self.new_to_old_keys_mapping[return_type]][0]

        if query_type == 'entrez_id':
            query = self.old_names_mapping['entrez_id'].get(q, None)
            if query:
                return query[self.new_to_old_keys_mapping[return_type]][0]
        return None
