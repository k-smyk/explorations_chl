# Adapted from https://github.com/jdellert/iwsa/blob/master/src/main/java/de/jdellert/iwsa/infomodel/InformationModel.java
import math
from collections import defaultdict
from phonetic_symbol_table import PhoneticSymbolTable


class InformationModel:
	VERBOSE = False

	def __init__(self, symbol_table: PhoneticSymbolTable):
		self.symbol_table = symbol_table
		self.num_symbols = symbol_table.get_size()
		self.num_trigrams = (self.num_symbols - 1) ** 3
		self.num_gappy_bigrams = (self.num_symbols - 1) ** 2
		self.observation_counts = defaultdict(int)
		self.observation_counts_sum = 0.0
		self.smoothing_mass_ratio = 0.2

	def set_smoothing_mass_ratio(self, ratio: float):
		self.smoothing_mass_ratio = ratio

	def add_trigram_observation(self, a: int, b: int, c: int):
		self.store_observation(self.trigram_id(a, b, c))  # ABC
		self.store_gappy_observation(self.trigram_id(a, b, 1))  # AB_
		self.store_gappy_observation(self.trigram_id(a, 1, c))  # A_C
		self.store_gappy_observation(self.trigram_id(1, b, c))  # _BC

	def trigram_count(self, a: int, b: int, c: int) -> int:
		return self.observation_counts[self.trigram_id(a, b, c)]

	def smoothed_trigram_count(self, a: int, b: int, c: int) -> float:
		return self.trigram_count(a, b, c) + (
					(self.smoothing_mass_ratio * self.observation_counts_sum) / self.num_trigrams)

	def smoothed_trigram_count_with_bigram(self, a: int, b: int, c: int, bigram_count: float) -> float:
		return self.trigram_count(a, b, c) + ((self.smoothing_mass_ratio * bigram_count) / (self.num_symbols - 1))

	def smoothed_gappy_bigram_count(self, a: int, b: int, c: int) -> float:
		return (1 + self.smoothing_mass_ratio) * self.trigram_count(a, b, c) + (
					(self.smoothing_mass_ratio * self.observation_counts_sum) / self.num_gappy_bigrams)

	def information_content(self, s: list, i: int) -> float:
		a = s[i - 2] if i > 1 else 0
		b = s[i - 1] if i > 0 else 0
		c = s[i]
		d = s[i + 1] if i < len(s) - 1 else 0
		e = s[i + 2] if i < len(s) - 2 else 0
		return self.information_content_abcde(a, b, c, d, e)

	def information_content_abcde(self, a: int, b: int, c: int, d: int, e: int) -> float:
		if self.VERBOSE:
			print(f"c({a:3d},{b:3d},{c:3d},{d:3d},{e:3d}):  ", end="")

		ab_count = self.smoothed_gappy_bigram_count(a, b, 1)
		bd_count = self.smoothed_gappy_bigram_count(b, 1, d)
		de_count = self.smoothed_gappy_bigram_count(1, d, e)
		abc_count = self.smoothed_trigram_count(a, b, c)
		bcd_count = self.smoothed_trigram_count(b, c, d)
		cde_count = self.smoothed_trigram_count(c, d, e)

		obs_prob = (abc_count + bcd_count + cde_count) / (ab_count + bd_count + de_count)

		if self.VERBOSE:
			print(f"{obs_prob:.3f} = ", end="")
			print(f"({round(abc_count)}+{round(bcd_count)}+{round(cde_count)})/", end="")
			print(f"({round(ab_count)}+{round(bd_count)}+{round(de_count)})")

		return -math.log(obs_prob)

	def store_observation(self, trigram_id: int):
		self.observation_counts[trigram_id] += 1
		self.observation_counts_sum += 1.0

	def store_gappy_observation(self, trigram_id: int):
		self.observation_counts[trigram_id] += 1

	def trigram_id(self, symbol_a: int, symbol_b: int, symbol_c: int) -> int:
		return symbol_a * self.num_symbols * self.num_symbols + symbol_b * self.num_symbols + symbol_c

	def print_counts(self, output_file, smoothing: bool):
		with open(output_file, "w") as out:
			for a in range(self.symbol_table.get_size()):
				if a == 1:
					continue
				for b in range(self.symbol_table.get_size()):
					if b == 1:
						continue
					ab_count = self.trigram_count(a, b, 1)
					if smoothing:
						ab_count = self.smoothed_gappy_bigram_count(a, b, 1)
					if ab_count >= 1.0:
						ab_string = self.symbol_table.to_symbol(a) + self.symbol_table.to_symbol(b)
						out.write(f"{ab_string}-:{ab_count}")
						for c in range(self.symbol_table.get_size()):
							if c == 1:
								continue
							abc_count = self.trigram_count(a, b, c)
							if smoothing:
								abc_count = self.smoothed_trigram_count(a, b, c)
							if abc_count > 0.002:
								out.write(f" {ab_string}{self.symbol_table.to_symbol(c)}:{abc_count}")
						out.write("\n")

			for a in range(self.symbol_table.get_size()):
				if a == 1:
					continue
				for c in range(self.symbol_table.get_size()):
					if c == 1:
						continue
					ac_count = self.trigram_count(a, 1, c)
					if smoothing:
						ac_count = self.smoothed_gappy_bigram_count(a, 1, c)
					if ac_count >= 1.0:
						a_string = self.symbol_table.to_symbol(a)
						c_string = self.symbol_table.to_symbol(c)
						out.write(f"{a_string}-{c_string}:{ac_count}")
						for b in range(self.symbol_table.get_size()):
							if b == 1:
								continue
							abc_count = self.trigram_count(a, b, c)
							if smoothing:
								abc_count = self.smoothed_trigram_count(a, b, c)
							if abc_count > 0.002:
								out.write(f" {a_string}{self.symbol_table.to_symbol(b)}{c_string}:{abc_count}")
						out.write("\n")

			for b in range(self.symbol_table.get_size()):
				if b == 1:
					continue
				for c in range(self.symbol_table.get_size()):
					if c == 1:
						continue
					bc_count = self.trigram_count(1, b, c)
					if smoothing:
						bc_count = self.smoothed_gappy_bigram_count(1, b, c)
					if bc_count >= 1.0:
						bc_string = self.symbol_table.to_symbol(b) + self.symbol_table.to_symbol(c)
						out.write(f"-{bc_string}:{bc_count}")
						for a in range(self.symbol_table.get_size()):
							if a == 1:
								continue
							abc_count = self.trigram_count(a, b, c)
							if smoothing:
								abc_count = self.smoothed_trigram_count(a, b, c)
							if abc_count > 0.002:
								out.write(f" {self.symbol_table.to_symbol(a)}{bc_string}:{abc_count}")
						out.write("\n")
			out.close()
