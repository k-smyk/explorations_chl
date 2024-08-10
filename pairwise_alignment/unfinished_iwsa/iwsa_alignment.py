from asjp import ipa2asjp, asjp2ipa

from phonetic_symbol_table import PhoneticSymbolTable
from corr_model import CorrespondenceModel, CorrespondenceModelInference

phonetic_symbols = asjp2ipa({"p", "b", "f", "v", "m", "w", "8", "t", "d",
                    "s", "z", "c", "n", "r", "l", "S", "Z", "C", "j",
                    "T", "5", "y", "k", "g", "x", "N", "q", "X", "h", "7", "L", "4", "G", "!",
                    "i", "e", "E", "3", "a", "u", "o"})

symbol_table = PhoneticSymbolTable(symbols=phonetic_symbols)
print(symbol_table)

for symbol in phonetic_symbols:
    print(f"Symbol: {symbol}, ID: {symbol_table.to_int(symbol)}")

segments = ["a", "b", "t"]
encoded = symbol_table.encode(segments)
print(f"Encoded: {encoded}")

decoded = symbol_table.decode(encoded)
print(f"Decoded: {decoded}")

correspondence_model = CorrespondenceModel(symbol_table=symbol_table)


# inference = CorrespondenceModelInference(concept_dict, symbol_table, verbose=True)
# global_corr_model = inference.infer_global_correspondence_model()

print(global_corr_model)


# distance = correspondence_model.compute_distance("bat", "bat")
# print(f"Distance between 'bat' and 'bat': {distance:.4f}")
#
# sc_ab, sc_aa, sc_bb, len_a, len_b = correspondence_model.compute_similarity_components("bat", "bad")
# normalized_distance = correspondence_model.normalize_score(sc_ab, sc_aa, sc_bb, len_a, len_b)
# print(f"Normalized distance: {normalized_distance:.4f}")
#
# distance = correspondence_model.compute_distance("bat", "bad")
# print(f"Distance between 'bat' and 'bad': {distance:.4f}")

# print("Correspondence Matrix:")
# print(correspondence_model)
