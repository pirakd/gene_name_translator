import pandas as pd
import pickle as pl
from tqdm import tqdm
from os import path, remove
import numpy as np
from collections import defaultdict
import gzip
import urllib.request
import zlib


class GeneTranslator:
    def __init__(self, verbosity=True):
        self.root_folder = path.dirname(path.realpath(__file__))
        self.raw_data_file = 'Homo_sapiens.gene_info'
        self.dictionary_filename = 'gene_dictionary.pl'
        self.old_names_file = 'gene_history'
        self.raw_data_file_path = path.join(self.root_folder, self.raw_data_file)
        self.dictionary_file_path = path.join(self.root_folder, self.dictionary_filename)
        self.old_names_file_path = path.join(self.root_folder, self.old_names_file)
        self.dictionary = None
        self.old_names_mapping = None
        self.verbosity = verbosity
        self.use_old_names = True
        self.new_to_old_keys_mapping = {'symbol': 'Discontinued_Symbol',
                                        'entrez_id': 'Discontinued_GeneID'}
        self.keys_to_data_files_headers = {'symbol': 'Symbol',
                                           'entrez_id': 'GeneID',
                                           'alias_symbol': 'Synonyms'}

    def translate(self, query, query_type, return_type):
        single_query_flag = False
        if not isinstance(query, (tuple, list, np.ndarray, set)):
            query = [query]
            single_query_flag = True
        keys_not_found = list()
        targets_not_found = list()
        result_dict = dict()
        dictionary = self.dictionary[query_type]
        for q in query:
            if q in dictionary:
                result = dictionary[q][return_type]
                if result is None:
                    targets_not_found.append(q)
            elif isinstance(q, str) and q.upper() in dictionary:
                result = dictionary[q.upper()][return_type]
                if result is None:
                    targets_not_found.append(q)
            elif isinstance(q, str) and q.lower() in dictionary:
                result = dictionary[q.lower()][return_type]
                if result is None:
                        targets_not_found.append(q)
            else:
                result = self.query_old_name(q, query_type, return_type)
                if result is None:
                    keys_not_found.append(q)
                    continue
            result_dict[q] = result

        if self.verbosity:
            if len(keys_not_found):
                print('{} genes were not found ({}))'.format(len(keys_not_found), keys_not_found))
            if len(targets_not_found):
                print('{} translations are missing ({})'.format(len(targets_not_found), targets_not_found))

        if single_query_flag:
            if len(result_dict):
                return result_dict[query[0]]
            else:
                return None
        return result_dict

    def load_dictionary(self):
        if not path.isfile(self.dictionary_file_path):
            print('Gene dictionary is missing! Generating dictionary, this might take a few minutes...')
            self.init_translator()

        with open(self.dictionary_file_path, 'rb') as f:
            save_dict = pl.load(f)
        self.dictionary = save_dict['dictionary']
        self.old_names_mapping = save_dict['old_names_mapping']
        if self.verbosity:
            print('available query types : {}'.format([x for x in self.dictionary.keys()]))

    def init_translator(self):
        self._download_files_()
        self._generate_dictionaries_()
        remove(self.raw_data_file_path)
        remove(self.old_names_file_path)

    def _generate_dictionaries_(self):
        keys = ['symbol', 'entrez_id']
        dictionaries = dict()
        for key in keys:
            dictionaries[key] = self._generate_dictionary_(key)

        old_names_map = self.load_old_names()

        save_dict = {'dictionary': dictionaries,
                     'old_names_mapping': dict(old_names_map)}
        with open(self.dictionary_file_path, 'wb') as f:
            pl.dump(save_dict, f)

        return dictionaries

    def _generate_dictionary_(self, key):
        name_df = pd.read_csv(self.raw_data_file_path, delimiter='\t', index_col=False, low_memory=False)

        dictionary = dict()
        row_values = dict()

        for idx, row in tqdm(name_df.iterrows(), desc='Generating dictionary with {} keys'.format(key),
                             total=len(name_df)):
            row_values['symbol'] = row[self.keys_to_data_files_headers['symbol']]
            row_values['entrez_id'] = int(row[self.keys_to_data_files_headers['entrez_id']]) \
                if not pd.isnull(self.keys_to_data_files_headers['entrez_id']) else None
            row_values['alias_symbol'] = row[self.keys_to_data_files_headers['alias_symbol']]
            if not pd.isnull(row_values[key]):
                aliases = []
                if not pd.isnull(row_values['alias_symbol']):
                    aliases = row_values['alias_symbol'].split('|')

                aliases = list(set([row_values['symbol']] + aliases))

                query_dict = dict(entrez_id=row_values['entrez_id'],
                                  symbol=row_values['symbol'], symbol_aliases=aliases)

                if key == 'symbol':
                    for alias in aliases:
                        if alias == row_values['symbol'] or alias not in dictionary:
                            dictionary[alias] = query_dict
                else:
                    dictionary[row_values[key]] = query_dict
        return dictionary

    def load_old_names(self):
        keys = ['entrez_id', self.new_to_old_keys_mapping['entrez_id'], self.new_to_old_keys_mapping['symbol']]
        old_names_dictionary = {key: dict() for key in keys}

        with open(self.old_names_file_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                entrez_id, old_entrez_id, old_symbol = line.strip().split('\t')[1:-1]
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
            result = query[self.new_to_old_keys_mapping[return_type]]
            if result:
                return result[0]
            else:
                return None

        if query_type == 'entrez_id':
            query = self.old_names_mapping['entrez_id'].get(q, None)
            if query:
                return query[self.new_to_old_keys_mapping[return_type]][0]
        return None

    def _download_files_(self):

        obsolete_gene_gz_file_path = self.old_names_file_path + '.gz'
        url = 'https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'

        # download entrez-symbol dict
        print('Downloading files from the web:')
        print('Downloading entrez id <-> symbol dict')
        request = urllib.request.urlopen(url).read()
        lines = str(zlib.decompress(request, 16 + zlib.MAX_WBITS), 'utf-8').split('\n')[:-1]
        lines = [x + '\n' for x in lines]
        if path.isfile(self.raw_data_file_path):
            remove(self.raw_data_file_path)
        with open(path.join(self.raw_data_file_path), 'w') as f:
            f.writelines(lines)

        # download old names dict
        url = 'https://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz'
        request_var = urllib.request.urlopen(url)
        size_in_bytes = request_var.length
        block_size_in_bytes = 1024 * (10 ** 3)
        progress_bar = tqdm(total=size_in_bytes, desc='Downloading obsolete gene names file', unit='iB',
                            unit_scale=True)
        with open(obsolete_gene_gz_file_path, 'wb') as f:
            while True:
                buf1 = request_var.read(block_size_in_bytes)
                if not buf1:
                    break
                f.write(buf1)
                progress_bar.update(block_size_in_bytes)
        progress_bar.close()

        print('Filtering unnecessary organisms from obsolete genes file')
        history_lines = []
        with gzip.open(obsolete_gene_gz_file_path, 'rb') as f:
            history_lines.append(str(f.readline(), 'utf-8'))
        with gzip.open(obsolete_gene_gz_file_path, 'rb') as f:
            for line in f:
                if line.startswith(b'9606'):
                    history_lines.append(str(line, 'utf-8'))
        remove(obsolete_gene_gz_file_path)
        with open(self.old_names_file_path, 'w') as f:
            f.writelines(history_lines, )

        print('Finished downloading up-to-date files')
