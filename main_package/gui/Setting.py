from enum import Enum, auto, IntEnum
from typing import Union


class TwoWayDict(dict):
    def __init__(self, dictionary: dict):
        super().__init__()
        for key, value in dictionary.items():
            self.__setitem__(key, value)

    def __setitem__(self, key, value):
        # Remove any previous connections with these values
        if key in self:
            del self[key]
        if value in self:
            del self[value]
        dict.__setitem__(self, key, value)
        dict.__setitem__(self, value, key)

    def __delitem__(self, key):
        dict.__delitem__(self, self[key])
        dict.__delitem__(self, key)

    def __len__(self):
        """Returns the number of connections"""
        return dict.__len__(self) // 2


# Have to use IntEnum here. Using just Enum results in strange behaviour: Two enums with the same value are not equal.
# See https://stackoverflow.com/questions/28125055/enum-in-python-doesnt-work-as-expected
class Setting(IntEnum):
    DEFAULT = auto()
    INTERMEDIATE = auto()
    ODD = auto()
    SILLY1 = auto()
    SILLY2 = auto()
    SILLY3 = auto()

    # NOTE: Type hint 'Setting' must be in quotes
    @classmethod
    def lookup(cls, setting_or_string: Union['Setting', str]) -> Union['Setting', str]:
        strings = {cls.DEFAULT: 'Default',
                   cls.INTERMEDIATE: 'Intermediate',
                   cls.ODD: 'Odd',
                   cls.SILLY1: 'Silly 1',
                   cls.SILLY2: 'Silly 2',
                   cls.SILLY3: 'Silly 3'}
        two_way_dict = TwoWayDict(strings)
        return two_way_dict[setting_or_string]
