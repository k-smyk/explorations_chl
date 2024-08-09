# class CorrespondenceModel:
# 	def __init__(self, symbol_table):
# 		self.symbol_table = symbol_table  # PhoneticSymbolTable object
# 		self.scores = {}
#
# 	def set_score(self, symbol1_id, symbol2_id, score):
# 		symbol_pair_id = self.symbol_table.get_size() * symbol1_id + symbol2_id
# 		self.scores[symbol_pair_id] = score
#
# 	def get_score(self, symbol1_id, symbol2_id):
# 		symbol_pair_id = self.symbol_table.get_size() * symbol1_id + symbol2_id
# 		return self.scores.get(symbol_pair_id, 0.0)
#
# 	def get_score_by_symbols(self, symbol1, symbol2):
# 		symbol1_id = self.symbol_table.to_int(symbol1)
# 		symbol2_id = self.symbol_table.to_int(symbol2)
# 		return self.get_score(symbol1_id, symbol2_id)

from phonetic_symbol_table import PhoneticSymbolTable


class CorrespondenceModel:
    def __init__(self, symbol_table):
        self.symbol_table = symbol_table #PhoneticSymbolTable
        self.scores = {}
        #self.db_path = None

    def set_score(self, symbol1_id, symbol2_id, score):
        symbol_pair_id = self.get_symbol_pair_id(symbol1_id, symbol2_id)
        self.scores[symbol_pair_id] = score

    def get_score(self, symbol1_id, symbol2_id):
        symbol_pair_id = self.get_symbol_pair_id(symbol1_id, symbol2_id)
        return self.scores.get(symbol_pair_id, 0.0)

    def get_score_by_symbols(self, symbol1, symbol2):
        return self.get_score(self.symbol_table.to_int(symbol1), self.symbol_table.to_int(symbol2))

    def get_symbol_pair_id(self, symbol1_id, symbol2_id):
        return self.symbol_table.get_size() * symbol1_id + symbol2_id

    def compute_similarity_components(self, word_a, word_b):
        """Compute the similarity score components for distance normalization."""
        len_a = len(word_a)
        len_b = len(word_b)

        sc_ab = self.calculate_similarity_score(word_a, word_b)
        sc_aa = self.calculate_similarity_score(word_a, word_a)
        sc_bb = self.calculate_similarity_score(word_b, word_b)

        return sc_ab, sc_aa, sc_bb, len_a, len_b

    def normalize_score(self, sc_ab, sc_aa, sc_bb, len_a, len_b):
        """Normalize the score based on the lengths of the sequences."""
        normalized_distance = 0.5 * (sc_ab / ((sc_aa / len_a) + (sc_bb / len_b)))
        return normalized_distance

    def compute_distance(self, word_a, word_b):
        """Compute the distance between two words."""
        sc_ab, sc_aa, sc_bb, len_a, len_b = self.compute_similarity_components(word_a, word_b)
        return self.normalize_score(sc_ab, sc_aa, sc_bb, len_a, len_b)

    def calculate_similarity_score(self, word_a, word_b):
        """Calculate similarity score between two words based on stored scores."""
        score = 0.0
        for symbol1, symbol2 in zip(word_a, word_b):
            score += self.get_score_by_symbols(symbol1, symbol2)
        return score

    def __str__(self):
        output = []
        for i in range(self.symbol_table.get_size()):
            row = [self.symbol_table.to_symbol(i)]
            for j in range(self.symbol_table.get_size()):
                row.append(f"{self.get_score(i, j):.4f}")
            output.append("\t".join(row))
        return "\n".join(output)

    # def set_db_path(self, db_path):
    #     self.db_path = db_path
    #
    # def get_db_path(self):
    #     return self.db_path




import math
import random
from collections import defaultdict, namedtuple

WordForm = namedtuple('WordForm', ['concept', 'language', 'word'])
data = [
    WordForm("night", "ALBANIAN_GHEG", "nat3"),
    WordForm("night", "ALBANIAN", "nat3"),
    WordForm("night", "ALBANIAN_TOSK", "nat3"),
    WordForm("full", "ALBANIAN_GHEG", "cit3"),
    WordForm("full", "ALBANIAN", "embuSurplotme"),
    WordForm("full", "ALBANIAN_TOSK", "plot3"),
    WordForm("new", "ALBANIAN_GHEG", "ri"),
    WordForm("new", "ALBANIAN", "iri"),
    WordForm("new", "ALBANIAN_TOSK", "ri"),
    WordForm("name", "ALBANIAN_GHEG", "em3n"),
    WordForm("name", "ALBANIAN", "em3r")
]
def organize_data_by_concept(data):
    concept_dict = defaultdict(lambda: defaultdict(list))
    for entry in data:
        concept_dict[entry.concept][entry.language].append(entry.word)
    return concept_dict
concept_dict = organize_data_by_concept(data)


class CorrespondenceModelInference:
    def __init__(self, concept_dict, symbol_table, num_global_reestimations=3, verbose=False):
        self.concept_dict = concept_dict
        self.symbol_table = symbol_table
        self.num_global_reestimations = num_global_reestimations
        self.verbose = verbose

    def infer_global_correspondence_model(self):
        random_pair_dist = self._simulate_non_cognates()
        cognate_pair_dist = self._find_cognate_candidates()
        global_corr = self._calculate_pmi(random_pair_dist, cognate_pair_dist)
        global_corr = self._reestimate_model(global_corr)
        return global_corr

    def _simulate_non_cognates(self):
        random_pair_dist = defaultdict(float)
        languages = list({lang for concepts in self.concept_dict.values() for lang in concepts})

        for _ in range(1000):  # to change?
            lang1, lang2 = random.sample(languages, 2)
            word1 = random.choice(self.concept_dict[random.choice(list(self.concept_dict))][lang1])
            word2 = random.choice(self.concept_dict[random.choice(list(self.concept_dict))][lang2])
            alignment = self._align_forms(word1, word2)

            for pos in range(len(alignment)):
                symbol_pair_id = self._get_symbol_pair_id(alignment, pos)
                random_pair_dist[symbol_pair_id] += 1.0

        return random_pair_dist

    def _find_cognate_candidates(self):
        cognate_pair_dist = defaultdict(float)
        num_pairs = 0
        num_cognate_pairs = 0

        for concept, languages in self.concept_dict.items():
            lang_pairs = [(lang1, lang2) for lang1 in languages for lang2 in languages if lang1 != lang2]
            for lang1, lang2 in lang_pairs:
                for word1 in languages[lang1]:
                    for word2 in languages[lang2]:
                        alignment = self._align_forms(word1, word2)
                        num_pairs += 1

                        if alignment.normalized_distance_score <= 0.35:
                            for pos in range(len(alignment)):
                                symbol_pair_id = self._get_symbol_pair_id(alignment, pos)
                                cognate_pair_dist[symbol_pair_id] += 1.0
                            num_cognate_pairs += 1

        if self.verbose:
            print(f"Aligned {num_pairs} form pairs, {num_cognate_pairs} look like cognates.")

        return cognate_pair_dist

    def _calculate_pmi(self, random_pair_dist, cognate_pair_dist):
        global_corr = CorrespondenceModel(self.symbol_table)

        for symbol_pair_id in range(self.symbol_table.get_size() ** 2):
            cognate_prob = cognate_pair_dist[symbol_pair_id]
            random_prob = random_pair_dist[symbol_pair_id]
            if random_prob > 0:  # avoid division by zero
                pmi_score = math.log(cognate_prob / random_prob)
                global_corr.set_score(symbol_pair_id // self.symbol_table.get_size(),
                                      symbol_pair_id % self.symbol_table.get_size(), pmi_score)

        return global_corr

    def _reestimate_model(self, global_corr):
        for _ in range(self.num_global_reestimations):
            cognate_pair_dist = self._find_cognate_candidates()
            global_corr = self._calculate_pmi(global_corr, cognate_pair_dist)

        return global_corr

    def _align_forms(self, form1, form2):
        pass

    def _compute_distance(self, form1, form2):
        pass

    def _get_symbol_pair_id(self, alignment, pos):
        symbol1_id = self.symbol_table.to_int(alignment[0][pos])
        symbol2_id = self.symbol_table.to_int(alignment[1][pos])
        return self.symbol_table.get_size() * symbol1_id + symbol2_id


