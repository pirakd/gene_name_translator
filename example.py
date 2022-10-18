from gene_translator import GeneTranslator

if __name__ == '__main__':

    translator = GeneTranslator()

    # On first run translator will initiate itself
    translator.load_dictionary()

    entrez_to_symbol = translator.translate([339457, 46], 'entrez_id', 'symbol')
    symbol_to_entrez = translator.translate([entrez_to_symbol[339457], entrez_to_symbol[46]], 'symbol', 'entrez_id')
    print(entrez_to_symbol)
    print(symbol_to_entrez)

    # 339457 is an obsolete entrez_id that was replaced by another id.
    # translating it would return the up to date symbol.
    # translating the updated symbol of 339457 will result in an updated entrez_id.
    # this means that 339457 != symbol_to_entrez[entrez_to_symbol[339457]]
