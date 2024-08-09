# Adapted from https://github.com/jdellert/iwsa/blob/master/src/main/java/de/jdellert/iwsa/sequence/PhoneticSymbolTable.java
import itertools
from collections import defaultdict
from typing import List, Dict, Set, Iterator

class PhoneticSymbolTable:
    VERBOSE = False

    BOUNDARY_ID = 0
    BOUNDARY_SYMBOL = "#"
    EMPTY_ID = 1
    EMPTY_SYMBOL = "-"
    UNKNOWN_SYMBOL = "?"

    def __init__(self, symbols: Set[str] = None):
        self.id_to_symbol: List[str] = [""] * 2
        self.symbol_to_id: Dict[str, int] = {}
        self.next_id: int = 0

        self.id_to_symbol[self.BOUNDARY_ID] = self.BOUNDARY_SYMBOL
        self.id_to_symbol[self.EMPTY_ID] = self.EMPTY_SYMBOL
        self.symbol_to_id[self.BOUNDARY_SYMBOL] = self.BOUNDARY_ID
        self.symbol_to_id[self.EMPTY_SYMBOL] = self.EMPTY_ID
        self.next_id = 2

        if symbols:
            self.define_symbols(symbols)

    def define_symbols(self, symbols: Set[str]):
        for symbol in symbols:
            self.define_symbol(symbol)

    def define_symbol(self, symbol: str) -> int:
        if symbol in self.symbol_to_id:
            if self.VERBOSE:
                print(f"WARNING from PhoneticSymbolTable: symbol {symbol} is already defined, returning the existing ID.")
            return self.symbol_to_id[symbol]
        else:
            self.id_to_symbol.append(symbol)
            self.symbol_to_id[symbol] = self.next_id
            self.next_id += 1
            return self.next_id - 1

    def contains(self, symbol: str) -> bool:
        return symbol in self.symbol_to_id

    def to_int(self, symbol: str) -> int:
        return self.symbol_to_id.get(symbol, -1)

    def to_symbol(self, id: int) -> str:
        return self.UNKNOWN_SYMBOL if id < 0 else self.id_to_symbol[id]

    def encode(self, segments: List[str]) -> List[int]:
        return [self.to_int(segment) for segment in segments]

    def decode(self, segment_ids: List[int]) -> List[str]:
        return [self.to_symbol(segment_id) for segment_id in segment_ids]

    def get_defined_symbols(self) -> Set[str]:
        return set(self.symbol_to_id.keys())

    def get_size(self) -> int:
        return len(self.id_to_symbol)

    def to_symbol_pair(self, symbol_pair_id: int) -> str:
        size = self.get_size()
        return f"({self.to_symbol(symbol_pair_id // size)},{self.to_symbol(symbol_pair_id % size)})"

    def __str__(self) -> str:
        line1 = "\t".join(map(str, range(len(self.id_to_symbol))))
        line2 = "\t".join(self.id_to_symbol)
        return f"{line1}\n{line2}"

    def __iter__(self) -> Iterator[str]:
        return self.PhoneticSymbolTableIterator(self)

    class PhoneticSymbolTableIterator:
        def __init__(self, pst):
            self.pst = pst
            self.i = 1

        def __iter__(self):
            return self

        def __next__(self) -> str:
            if self.i + 1 >= len(self.pst.id_to_symbol):
                raise StopIteration
            self.i += 1
            return self.pst.id_to_symbol[self.i]
