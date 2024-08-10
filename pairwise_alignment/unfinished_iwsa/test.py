import numpy as np
from info_score import InformationModel
from phonetic_symbol_table import PhoneticSymbolTable
from corr_model import CorrespondenceModel

class InformationWeightedSequenceAlignment:
	def __init__(self, str1, str2, corr_model, info_model1, info_model2):
		self.str1 = str1  # PhoneticString object
		self.str2 = str2  # PhoneticString object
		self.corr_model = corr_model  # Unified CorrespondenceModel object
		self.info_model1 = info_model1  # InformationModel object for the first language
		self.info_model2 = info_model2  # InformationModel object for the second language

	@staticmethod
	def get_correspondence_score(corr_model, seg1, seg2):
		return corr_model.get_score(seg1, seg2)

	@staticmethod
	def get_mean_info_score(info_model1, info_model2, i1, i2):
		return (info_model1.get_score(i1) + info_model2.get_score(i2)) / 2.0

	def construct_alignment(self):
		m = len(self.str1) + 1
		n = len(self.str2) + 1

		mtx = np.zeros((m, n))
		a_subst = np.zeros((m, n), dtype=int)
		b_subst = np.zeros((m, n), dtype=int)

		mtx[0][0] = 0

		for i in range(1, m):
			mtx[i][0] = (mtx[i - 1][0] +
						 self.get_correspondence_score(self.corr_model, self.str1[i - 1], 1) *
						 self.get_mean_info_score(self.info_model1, self.info_model2, i - 1, i - 1))
			a_subst[i][0] = self.str1[i - 1]
			b_subst[i][0] = 1  # gap symbol

		for j in range(1, n):
			mtx[0][j] = (mtx[0][j - 1] +
						 self.get_correspondence_score(self.corr_model, 1, self.str2[j - 1]) *
						 self.get_mean_info_score(self.info_model1, self.info_model2, j - 1, j - 1))
			a_subst[0][j] = 1  # gap symbol
			b_subst[0][j] = self.str2[j - 1]

		for i in range(1, m):
			for j in range(1, n):
				match_score = (self.get_correspondence_score(self.corr_model, self.str1[i - 1], self.str2[j - 1]) *
							   self.get_mean_info_score(self.info_model1, self.info_model2, i - 1, j - 1))
				insert_score = (self.get_correspondence_score(self.corr_model, self.str1[i - 1], 1) *
								self.get_mean_info_score(self.info_model1, self.info_model2, i - 1, j))
				delete_score = (self.get_correspondence_score(self.corr_model, 1, self.str2[j - 1]) *
								self.get_mean_info_score(self.info_model1, self.info_model2, i, j - 1))

				mtx[i][j] = max(
					mtx[i - 1][j - 1] + match_score,  # Match/mismatch
					mtx[i - 1][j] + delete_score,  # Deletion
					mtx[i][j - 1] + insert_score  # Insertion
				)

				if mtx[i][j] == mtx[i - 1][j - 1] + match_score:
					a_subst[i][j] = self.str1[i - 1]
					b_subst[i][j] = self.str2[j - 1]
				elif mtx[i][j] == mtx[i - 1][j] + delete_score:
					a_subst[i][j] = self.str1[i - 1]
					b_subst[i][j] = 1  # gap
				else:
					a_subst[i][j] = 1  # gap
					b_subst[i][j] = self.str2[j - 1]
		return mtx, a_subst, b_subst


class PhoneticString:
	def __init__(self, segments):
		self.segments = segments

	def __len__(self):
		return len(self.segments)

	def __getitem__(self, index):
		return self.segments[index]


str1 = PhoneticString([2, 3, 4])
str2 = PhoneticString([3, 4, 5])
corr_model = CorrespondenceModel()
info_model1 = InformationModel()
str2 = PhoneticSymbolTable() #?
info_model2 = InformationModel()

iwsa = InformationWeightedSequenceAlignment(str1, str2, corr_model, info_model1, info_model2)
alignment_mtx, a_subst, b_subst = iwsa.construct_alignment()

print(alignment_mtx, a_subst, b_subst)
