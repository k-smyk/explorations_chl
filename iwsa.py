import numpy as np
class PhoneticSymbolTable:
    """
    From:
    https://github.com/jdellert/iwsa/blob/master/src/main/java/de/jdellert/iwsa/sequence/PhoneticSymbolTable.java
    """
    BOUNDARY_ID = 0
    BOUNDARY_SYMBOL = "#"
    EMPTY_ID = 1
    EMPTY_SYMBOL = "-"
    UNKNOWN_SYMBOL = "?"

    def __init__(self, symbols=None):
        self.id_to_symbol = []
        self.symbol_to_id = {}
        self.next_id = 0

        # initialize with boundary and empty symbols
        self.id_to_symbol.append(self.BOUNDARY_SYMBOL)
        self.symbol_to_id[self.BOUNDARY_SYMBOL] = self.BOUNDARY_ID
        self.id_to_symbol.append(self.EMPTY_SYMBOL)
        self.symbol_to_id[self.EMPTY_SYMBOL] = self.EMPTY_ID
        self.next_id = 2

        # define additional symbols if provided
        if symbols is not None:
            self.define_symbols(symbols)

    def define_symbols(self, symbols):
        for symbol in symbols:
            self.define_symbol(symbol)

    def define_symbol(self, symbol):
        if symbol in self.symbol_to_id:
            return self.symbol_to_id[symbol]
        else:
            self.id_to_symbol.append(symbol)
            self.symbol_to_id[symbol] = self.next_id
            self.next_id += 1
            return self.next_id - 1

    def contains(self, symbol):
        return symbol in self.symbol_to_id

    def to_int(self, symbol):
        return self.symbol_to_id.get(symbol, -1)

    def to_symbol(self, id):
        return self.UNKNOWN_SYMBOL if id < 0 else self.id_to_symbol[id]

    def encode(self, segments):
        return [self.to_int(segment) for segment in segments]

    def decode(self, segment_ids):
        return [self.to_symbol(segment_id) for segment_id in segment_ids]

    def get_defined_symbols(self):
        return set(self.symbol_to_id.keys())

    def get_size(self):
        return len(self.id_to_symbol)

    def to_symbol_pair(self, symbol_pair_id):
        size = self.get_size()
        return f"({self.to_symbol(symbol_pair_id // size)},{self.to_symbol(symbol_pair_id % size)})"

    def __str__(self):
        line1 = "\t".join(map(str, range(len(self.id_to_symbol))))
        line2 = "\t".join(self.id_to_symbol)
        return f"{line1}\n{line2}"

    def __iter__(self):
        return self.PhoneticSymbolTableIterator(self)

    class PhoneticSymbolTableIterator:
        def __init__(self, table):
            self.table = table
            self.index = 1

        def __iter__(self):
            return self

        def __next__(self):
            if self.index + 1 < len(self.table.id_to_symbol):
                self.index += 1
                return self.table.id_to_symbol[self.index]
            else:
                raise StopIteration


class CorrespondenceModel:
    """
    From:
    https://github.com/jdellert/iwsa/blob/master/src/main/java/de/jdellert/iwsa/corrmodel/CorrespondenceModel.java
    """
    def __init__(self, symbol_table):
        self.symbol_table = symbol_table
        self.scores = {}

    def set_score(self, symbol1_id, symbol2_id, score):
        symbol_pair_id = self.symbol_table.get_size() * symbol1_id + symbol2_id
        self.scores[symbol_pair_id] = score

    def get_score(self, symbol1_id, symbol2_id):
        symbol_pair_id = self.symbol_table.get_size() * symbol1_id + symbol2_id
        return self.scores.get(symbol_pair_id, 0.0)


# def get_correspondence_score(glo_corr_model, loc_corr_model, s1, s2, i, j):
#     if i == -1 or j == -1:
#         return -1  # gap penalty
#     symbol1_id = s1[i]
#     symbol2_id = s2[j]
#     return glo_corr_model.get_score(symbol1_id, symbol2_id)


class PhoneticString:
    def __init__(self, segments):
        self.segments = segments

    def get_length(self):
        return len(self.segments)

    def __str__(self):
        return "[" + " ".join(map(str, self.segments)) + "]"

    def to_string(self, symbol_table):
        return " ".join(symbol_table.decode(self.segments))

    def to_untokenized_string(self, symbol_table):
        return "".join(symbol_table.decode(self.segments))

    def copy_without_gaps(self):
        reduced_segments = [segment for segment in self.segments if segment > 1]
        return PhoneticString(reduced_segments)

    def __eq__(self, other):
        if isinstance(other, PhoneticString):
            return self.segments == other.segments
        return False

    def __hash__(self):
        return hash(tuple(self.segments))


class PhoneticStringAlignment:
    def __init__(self, str1=None, str2=None, alignment_score=0, normalized_distance_score=0):
        self.str1 = str1
        self.str2 = str2
        self.alignment_score = alignment_score
        self.normalized_distance_score = normalized_distance_score

    def __str__(self):
        return f"Alignment A:\n{self.str1}\nAlignment B:\n{self.str2}\nAlignment Score: {self.alignment_score}\nNormalized Distance Score: {self.normalized_distance_score}"

    def to_string(self, symbol_table):
        return self.str1.to_string(symbol_table) + "\n" + self.str2.to_string(symbol_table)

    def get_length(self):
        return len(self.str1.segments)

    def get_symbol_pair_id_at_pos(self, pos, symbol_table):
        return self.str1.segments[pos] * symbol_table.get_size() + self.str2.segments[pos]

    def get_symbol1_id_at_pos(self, pos):
        return self.str1.segments[pos]

    def get_symbol2_id_at_pos(self, pos):
        return self.str2.segments[pos]

    def get_symbol_pairs(self, symbol_table):
        symbol_pairs = []
        for i in range(self.get_length()):
            symbol_pairs.append([symbol_table.to_symbol(self.str1.segments[i]), symbol_table.to_symbol(self.str2.segments[i])])
        return symbol_pairs

    def __eq__(self, other):
        if isinstance(other, PhoneticStringAlignment):
            return (self.alignment_score == other.alignment_score and
                    self.normalized_distance_score == other.normalized_distance_score and
                    self.str1 == other.str1 and
                    self.str2 == other.str2)
        return False

    def __hash__(self):
        return (5 * (hash(self.alignment_score) +
                     19 * (hash(self.normalized_distance_score) +
                           37 * (hash(self.str1) + 13 * hash(self.str2)))))


class InformationWeightedSequenceAlignment:
    """
    Information-Weighted Sequence Alignment, from:
    https://github.com/jdellert/iwsa/blob/master/src/main/java/de/jdellert/iwsa/align/iwsa.java
    """
    GAP_SYMBOL = "-"
    NEW_DISTANCE_TRANSFORMATION = True

    @staticmethod
    def iwsa(str1, str2, glo_corr_model, loc_corr_model, self_sim_model1, self_sim_model2, info_model1,
             info_model2):
        GAP_SYMBOL = 1
        m = str1.get_length() + 1
        n = str2.get_length() + 1

        mtx = np.zeros((m, n))
        a_subst = np.zeros((m, n), dtype=int)
        b_subst = np.zeros((m, n), dtype=int)

        for i in range(1, m):
            mtx[i][0] = mtx[i - 1][0] + glo_corr_model.get_score(str1.segments[i - 1],
                                                                 GAP_SYMBOL) * info_model1.get_mean_info_score(
                str1.segments[i - 1], str1.segments[i - 1])
            a_subst[i][0] = str1.segments[i - 1]
            b_subst[i][0] = 1  # corresponds to gap symbol

        for j in range(1, n):
            mtx[0][j] = mtx[0][j - 1] + glo_corr_model.get_score(GAP_SYMBOL,
                                                                 str2.segments[
                                                                     j - 1]) * info_model2.get_mean_info_score(
                str2.segments[j - 1], str2.segments[j - 1])
            a_subst[0][j] = 1  # corr to gap symbol
            b_subst[0][j] = str2.segments[j - 1]

        for i in range(1, m):
            for j in range(1, n):
                match_value = mtx[i - 1][j - 1] + glo_corr_model.get_score(str1.segments[i - 1], str2.segments[
                    j - 1]) * info_model1.get_mean_info_score(str1.segments[i - 1], str2.segments[j - 1])
                insertion_value = mtx[i][j - 1] + glo_corr_model.get_score(GAP_SYMBOL,
                    str2.segments[j - 1]) * info_model2.get_mean_info_score(str2.segments[j - 1],
                                                                            str2.segments[j - 1])
                deletion_value = mtx[i - 1][j] + glo_corr_model.get_score(str1.segments[i - 1],
                                                                          GAP_SYMBOL) * info_model1.get_mean_info_score(
                    str1.segments[i - 1], str1.segments[i - 1])

                mtx[i][j] = max(match_value, max(insertion_value, deletion_value))

                if mtx[i][j] == match_value:
                    a_subst[i][j] = str1.segments[i - 1]
                    b_subst[i][j] = str2.segments[j - 1]
                elif mtx[i][j] == insertion_value:
                    a_subst[i][j] = 1
                    b_subst[i][j] = str2.segments[j - 1]
                else:
                    a_subst[i][j] = str1.segments[i - 1]
                    b_subst[i][j] = 1

        i = m - 1
        j = n - 1
        result1 = []
        result2 = []
        while i > 0 or j > 0:
            a_part = a_subst[i][j]
            b_part = b_subst[i][j]
            result1.insert(0, a_part)
            result2.insert(0, b_part)
            if a_part != 1:
                i -= 1
            if b_part != 1:
                j -= 1
            if a_part == 1 and b_part == 1:
                i -= 1
                j -= 1
            if i < 0 or j < 0:
                break

        similarity_score = mtx[m - 1][n - 1]
        str1_self_similarity = 0.0
        for i in range(str1.get_length()):
            str1_self_similarity += InformationWeightedSequenceAlignment.get_correspondence_score(glo_corr_model, self_sim_model1, str1.segments[i], str1.segments[i]) * InformationWeightedSequenceAlignment.get_mean_info_score(str1, str1, i, i, info_model1, info_model1)

        str2_self_similarity = 0.0
        for j in range(str2.get_length()):
            str2_self_similarity += InformationWeightedSequenceAlignment.get_correspondence_score(glo_corr_model, self_sim_model2, str2.segments[j], str2.segments[j]) * InformationWeightedSequenceAlignment.get_mean_info_score(str2, str2, j, j, info_model2, info_model2)

        if InformationWeightedSequenceAlignment.NEW_DISTANCE_TRANSFORMATION:
            similarity_score /= len(result1)
            str1_self_similarity /= m - 1
            str2_self_similarity /= n - 1

        normalized_distance_score = 1 - (2 * similarity_score) / (str1_self_similarity + str2_self_similarity)

        alignment = PhoneticStringAlignment(
            str1=PhoneticString(result1),
            str2=PhoneticString(result2),
            alignment_score=int(similarity_score),
            normalized_distance_score=int(normalized_distance_score)
        )

        return alignment

    @staticmethod
    def get_info_score(str, pos, info_model):
        return info_model.information_content(str.segments, pos)

    @staticmethod
    def get_mean_info_score(str1, str2, pos1, pos2, info_model1, info_model2):
        info_content1 = info_model1.information_content(str1.segments, pos1)
        info_content2 = info_model2.information_content(str2.segments, pos2)
        return np.sqrt((info_content1 ** 2 + info_content2 ** 2) / 2)

    @staticmethod
    def get_correspondence_score(glo_corr_model, loc_corr_model, ci, cj):
        score = loc_corr_model.get_score_or_null(ci, cj)
        if score is None:
            score = glo_corr_model.get_score_or_null(ci, cj)
        if score is None:
            score = 0.0
        return score

    @staticmethod
    def combined_info_scores_for_alignment(alignment, info_model1, info_model2):
        pos1 = -1
        pos2 = -1

        str1_reduced = alignment.str1.copy_without_gaps()
        str2_reduced = alignment.str2.copy_without_gaps()

        info_scores = np.zeros(alignment.get_length())

        for pos in range(alignment.get_length()):
            symb1 = alignment.str1.segments[pos]
            symb2 = alignment.str2.segments[pos]

            if symb1 > 1:
                pos1 += 1
            if symb2 > 1:
                pos2 += 1

            if symb1 == 1:
                info_scores[pos] = InformationWeightedSequenceAlignment.get_mean_info_score(str2_reduced, str2_reduced, pos2, pos2, info_model2, info_model2)
            elif symb2 == 1:
                info_scores[pos] = InformationWeightedSequenceAlignment.get_mean_info_score(str1_reduced, str1_reduced, pos1, pos1, info_model1, info_model1)
            else:
                info_scores[pos] = InformationWeightedSequenceAlignment.get_mean_info_score(str1_reduced, str2_reduced, pos1, pos2, info_model1, info_model2)

        return info_scores

