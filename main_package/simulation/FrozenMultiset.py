from collections import Counter


class FrozenMultiset:
    """
    Immutable type to be used as a dict key.
    Serves as a unique identifier for an interaction term (e.g. lambda_A_A_C_C).
    """
    def __init__(self, *items: str, on: bool = False):
        # Initial On/off state of the interaction term.
        # TODO: This is kludge. Possibly move self._on somewhere else.
        self._on = on
        # tuple is immutable
        # Needs to be sorted to compare tuple with other tuples
        self.items: tuple[str, ...] = tuple(sorted(items))

    def __contains__(self, val, n: int):
        # Overrides keyword "in"
        # returns True when val is in the multiset n times, else returns False
        return val in self.items and self.items.count(val) == n

    def __len__(self):
        # returns the number of elements in the multiset
        return len(self.items)

    def __eq__(self, other):
        """Necessary to implement equality check to use this class as a dict key."""

        # Check if both tuples contain the same elements
        return hasattr(other, 'items') and self.items == other.items

    def __hash__(self):
        """Necessary to make this class hashable to use it as a dict key."""
        return hash(self.items)

    def __repr__(self):
        """Called for printing this object to the console. Useful for debugging."""
        return self.get_name()

    def __str__(self):
        """Called for printing this object to the console. Useful for debugging."""
        return self.get_name()

    def to_dict(self) -> dict[str, int]:
        """Returns as a dict of name:frequency pairs, e.g. {A: 3, B: 1, C: 2}"""
        return dict(Counter(self.items))

    def get_items(self) -> list[str]:
        """Returns all items, e.g. [A, A, A, B, C, C]"""
        return list(self.items)

    def distinct_items(self) -> list[str]:
        """Returns list of distinct items, e.g. [A, B, C] instead of [A, A, A, B, C, C]"""
        return list(Counter(self.items).keys())

    def frequencies(self) -> list[int]:
        """Returns number of times a string appears as a list, e.g. [3, 1, 2] for [A, A, A, B, C, C]"""
        return list(Counter(self.items).values())

    def get_name(self) -> str:
        """Make a string representing this object."""
        return '_'.join(self.items)

    def get_on(self) -> bool:
        return self._on

    def set_on(self, on: bool) -> None:
        self._on = on