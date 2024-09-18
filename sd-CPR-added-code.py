# def read_custom_csv(fname):
#     sounds = [x.strip() for x in open('pmi_model/sounds41.txt').readlines()]
#     d = defaultdict(lambda: defaultdict())
#     cogd = defaultdict(lambda: defaultdict())
#     with open(fname, 'r') as f:
#         header = ['index', 'concept', 'language', 'word', 'ID', 'cc']
#         reader = csv.DictReader(f, fieldnames=header)
#         next(reader)  # skip the header row
#         asjp_idx = header.index("word")
#         cog_idx = header.index("ID")
#         print("Reading indexes ", asjp_idx, cog_idx)
#         ID = 1
#         for row in reader:
#             print(row)
#             concept = row['concept']
#             word = row['word']
#             lang = row['language']
#             id = int(row['ID'])
#             cc = row['cc']  # =cogidx
#             word = word.replace("~", "").replace(" ", "").replace("%", "").replace("*", "").replace("$", "").replace(
#                 "\"", "").replace("K", "").replace("D", "d")
#             for i in word:
#                 if i not in sounds:
#                     print("Sound %s not in ASJP alphabet" % (i))
#             if len(word) < 1: continue
#             d[concept][ID] = word
#             cogd[concept][ID] = cc
#             ID += 1
#     return d, cogd

# python3 pmi_CRP.py -i albanoRomanceCC.csv -o albanoRomanceCC_CPR.csv --sample

ccData = pd.DataFrame()
th = .5
for c in concepts:
    cData = data[data.concept==c].copy()
    cPairs = synpairs[synpairs.concept==c]
    cIDs = cData.ID.values
    simMtx = zeros((len(cIDs),len(cIDs)))
    simMtx[pd.match(cPairs.ID1.values,cIDs),
           pd.match(cPairs.ID2.values,cIDs)] = cPairs.prediction.values
    simMtx[pd.match(cPairs.ID2.values,cIDs),
           pd.match(cPairs.ID1.values,cIDs)] = cPairs.prediction.values
    simMtx[simMtx<th]=0
    G = igraph.Graph.Weighted_Adjacency(list(simMtx))
    clusters = G.community_label_propagation(weights='weight')
    ccDict = {cIDs[x]:i for i,cl in enumerate(clusters)
              for x in cl}
    cData['cc'] = [c+':'+str(ccDict[i]) for i in cData.ID.values]
    ccData = pd.concat([ccData,cData])


taxa = ccData.language.unique()

ccMtx = pd.DataFrame(index=taxa)
for c in concepts:
    cData = ccData[ccData.concept==c]
    cMtx = pd.crosstab(cData.language,cData.cc)
    cMtx[cMtx>1]=1
    cMtx = cMtx.reindex(taxa,fill_value='-')
    ccMtx = pd.concat([ccMtx,cMtx],axis=1)



ccMtx.to_csv('albanoRomanceCCbin.csv')

nexCharOutput(ccMtx.values,ccMtx.index,'albanoRomanceCC.nex')


ccData.to_csv('albanoRomanceCC.csv',index='False')
