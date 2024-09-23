def read_custom_csv(fname):
    sounds = [x.strip() for x in open('pmi_model/sounds41.txt').readlines()]
    d = defaultdict(lambda: defaultdict())
    cogd = defaultdict(lambda: defaultdict())
    langd = defaultdict(lambda: defaultdict())
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
            langd[concept][ID] = lang
            ID += 1
    return d, cogd, langd
#
# python3 pmi_CRP.py -i albanoRomanceCC.csv -o albanoRomanceCC_CPR.csv --sample

all_socher_result = []
seen_entries = set()
def socherCRP(pair_dist, gold_dict, gloss, lang_dict):
    """Socher et al. (2011) sampling based CRP. Modified to save word-to-cluster mapping to a file."""
    alpha = float(args.calpha)
    items_list = list(pair_dist.keys())
    assert len(items_list) >= 1
    n_items = len(items_list)

    cluster_idx = [[x] for x in items_list]
    n_iter = 1
    bcubed_prec, bcubed_recall, bcubed_fscore, ari_vec = [], [], [], []
    g = igraph.Graph(n_items, directed=True)

    for v in range(n_items):
        g.vs[v]["name"] = items_list[v].split("::")[0]

    for n_iter in range(max_iter):
        n_single_clusters = 0
        for i, item in enumerate(items_list):
            sim_vec = [0.0] * n_items

            for j, v in enumerate(items_list):
                if i == j:
                    sim_vec[j] = alpha * pair_dist[v][v]
                    continue

                if g.are_connected(i, j):
                    g.delete_edges([(i, j)])

                sitbehind = g.subcomponent(j, mode="in")
                sSim = 0.0
                for s in sitbehind:
                    sSim += pair_dist[item][items_list[s]]
                sim_vec[j] = sSim
            insert_index = np.argmax(sim_vec)
            g.add_edges([(i, insert_index)])

        predicted_labels, gold_labels = [], []
        cluster_count = {}

        for k_idx, k in enumerate(g.clusters(mode="weak")):
            for x in k:
                word = items_list[x].split("::")[0]
                gloss_index = int(items_list[x].split("::")[1])
                gloss_name = gold_dict[gloss_index]
                lang_name = lang_dict[gloss_name.split(":")[0]][gloss_index]

                if word not in cluster_count:
                    cluster_count[word] = 1
                else:
                    cluster_count[word] += 1

                cluster_name = f"{gloss_name[:-2]}:{cluster_count[word]}"
                predicted_labels.append(cluster_name)
                gold_labels.append(gloss_name)
                entry_identifier = (word, gloss_name, lang_name)

                if entry_identifier not in seen_entries:
                    seen_entries.add(entry_identifier)  # Add to the set to track uniqueness
                    all_socher_result.append((len(all_socher_result), gloss_name, lang_name, word, cluster_name))

            # if (word, gloss_name, lang_name) not in all_socher_result:
            #     all_socher_result.append((len(all_socher_result), gloss_name, lang_name, word, cluster_name))

        p, r, f_score = utils.b_cubed(gold_labels, predicted_labels)
        bcubed_prec.append(float(p))
        bcubed_recall.append(float(r))
        bcubed_fscore.append(float(f_score))
        ari = metrics.adjusted_rand_score(gold_labels, predicted_labels)
        ari_vec.append(ari)

        print(gloss, "Socher", str(bcubed_prec[-1]), str(bcubed_recall[-1]), str(bcubed_fscore[-1]), str(ari_vec[-1]),
              str(len(set(predicted_labels))), str(len(set(gold_labels))), alpha, sep="\t", file=fh)

        if args.sample:
            alpha = sample_alpha(n_single_clusters, len(set(predicted_labels)), alpha)
    return



# ccData = pd.DataFrame()
# th = .5
# for c in concepts:
#     cData = data[data.concept==c].copy()
#     cPairs = synpairs[synpairs.concept==c]
#     cIDs = cData.ID.values
#     simMtx = zeros((len(cIDs),len(cIDs)))
#     simMtx[pd.match(cPairs.ID1.values,cIDs),
#            pd.match(cPairs.ID2.values,cIDs)] = cPairs.prediction.values
#     simMtx[pd.match(cPairs.ID2.values,cIDs),
#            pd.match(cPairs.ID1.values,cIDs)] = cPairs.prediction.values
#     simMtx[simMtx<th]=0
#     G = igraph.Graph.Weighted_Adjacency(list(simMtx))
#     clusters = G.community_label_propagation(weights='weight')
#     ccDict = {cIDs[x]:i for i,cl in enumerate(clusters)
#               for x in cl}
#     cData['cc'] = [c+':'+str(ccDict[i]) for i in cData.ID.values]
#     ccData = pd.concat([ccData,cData])
#
#
# taxa = ccData.language.unique()
#
# ccMtx = pd.DataFrame(index=taxa)
# for c in concepts:
#     cData = ccData[ccData.concept==c]
#     cMtx = pd.crosstab(cData.language,cData.cc)
#     cMtx[cMtx>1]=1
#     cMtx = cMtx.reindex(taxa,fill_value='-')
#     ccMtx = pd.concat([ccMtx,cMtx],axis=1)
#
#
#
# ccMtx.to_csv('albanoRomanceCCbin.csv')
#
# nexCharOutput(ccMtx.values,ccMtx.index,'albanoRomanceCC.nex')
#
#
# ccData.to_csv('albanoRomanceCC.csv',index='False')


###

fname = args.infile
#d, cogd = utils.readCSV(fname)
d, cogd, langd = utils.read_custom_csv(fname)
cluster_CRP(d, cogd, fname, langd)

fh.close()
