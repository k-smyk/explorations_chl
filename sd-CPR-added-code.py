def read_custom_csv(fname):
    sounds = [x.strip() for x in open('pmi_model/sounds41.txt').readlines()]
    d = defaultdict(lambda: defaultdict())
    cogd = defaultdict(lambda: defaultdict())
    with open(fname, 'r') as f:
        header = ['index', 'concept', 'language', 'word', 'ID', 'cc']
        reader = csv.DictReader(f, fieldnames=header)
        next(reader)  # skip the header row
        asjp_idx = header.index("word")
        cog_idx = header.index("ID")
        print("Reading indexes ", asjp_idx, cog_idx)
        ID = 1
        for row in reader:
            print(row)
            concept = row['concept']
            word = row['word']
            lang = row['language']
            id = int(row['ID'])
            cc = row['cc']  # =cogidx
            word = word.replace("~", "").replace(" ", "").replace("%", "").replace("*", "").replace("$", "").replace(
                "\"", "").replace("K", "").replace("D", "d")
            for i in word:
                if i not in sounds:
                    print("Sound %s not in ASJP alphabet" % (i))
            if len(word) < 1: continue
            d[concept][ID] = word
            cogd[concept][ID] = cc
            ID += 1
    return d, cogd

# python3 pmi_CRP.py -i albanoRomanceCC.csv -o albanoRomanceCC_CPR.csv --sample
