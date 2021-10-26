from gene_translator import GeneTranslator


if __name__ == '__main__':
    translator = GeneTranslator()
    # this should be done only once after clone
    translator.generate_dictionaries()
    translator.load_dictionary()
    a = translator.translate([23, 1], 'entrez_id', 'symbol')
    b = translator.translate([a[23], a[1], 46], 'symbol', 'entrez_id')
    print(a)
    print(b)
