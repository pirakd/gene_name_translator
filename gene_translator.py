import pandas as pd
import pickle as pl
import tqdm


class DictionaryEntry:
    def __init__(self, entrez_id, uniprot, ensembl_gene_id, symbol, symbol_aliases):
        self.entrez_id = entrez_id
        self.uniprot = uniprot
        self.ensembl_gene_id = ensembl_gene_id
        self.symbol = symbol
        self.symbol_aliases = symbol_aliases


class GeneTranslator:
    def __init__(self):
        self.raw_data_folder = 'hgnc_complete_set.txt'
        self.dictionary_file_name = 'gene_dictionary.pl'
        self.dictionary = None

    def translate(self, query, query_type, return_type):
        result_dict = dict()
        dictionary = self.dictionary[query_type]
        for q in query:
            result_dict[q] = getattr(dictionary[q], return_type)
        return result_dict

    def load_dictionary(self):
        with open(self.dictionary_file_name, 'rb') as f:
            self.dictionary = pl.load(f)

    def generate_dictionaries(self, keys=None):

        if keys is None:
            keys = ['symbol', 'entrez_id', 'uniprot']
        dictionaries = dict()
        for key in keys:
            dictionaries[key] = self.generate_dictionary(key)
        with open(self.dictionary_file_name, 'wb') as f:
            pl.dump(dictionaries, f)
        return dictionaries

    def generate_dictionary(self, key):
        name_df = pd.read_csv(self.raw_data_folder, delimiter='\t', index_col=False)

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
                if not pd.isnull(row.alias_symbol):
                    aliases = row.alias_symbol.split('|')
                    aliases = [symbol] + aliases
                else:
                    aliases = [symbol]
                entrez_id = int(row.entrez_id) if not pd.isnull(row.entrez_id) else None
                uniprot = row.uniprot_ids if not pd.isnull(row.uniprot_ids) else None
                ensembl_gene_id = row.ensembl_gene_id if not pd.isnull(row.ensembl_gene_id) else None

                if key == 'symbol':
                    for alias in aliases:
                        dictionary[alias] = DictionaryEntry(entrez_id=entrez_id, uniprot=uniprot,
                                                            ensembl_gene_id=ensembl_gene_id, symbol=symbol,
                                                            symbol_aliases=aliases)
                else:
                    dictionary[row_values[key]] = DictionaryEntry(entrez_id=entrez_id, uniprot=uniprot,
                                                                  ensembl_gene_id=ensembl_gene_id, symbol=symbol,
                                                                  symbol_aliases=aliases)
        return dictionary
