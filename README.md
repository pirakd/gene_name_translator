
# (Human) Gene Name Translator

-- See working example in [example.py](https://github.com/pirakd/gene_name_translator/blob/main/example.py)

-- Default supported formats: 'symbol', 'entrez_id', 'uniprot', 'ensemble_gene_id'


### Further Use 
-- Additional query types can be added by calling: 
```python 
GeneTranslator.generate_dictionaries(keys=[list_of_new_keys])
```
-- Available query types can be found in the [raw data file](https://github.com/pirakd/gene_name_translator/blob/main/hgnc_complete_set.txt) headers

### Data Source 

All data is taken from https://www.genenames.org
