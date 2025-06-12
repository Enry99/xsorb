"""
Base module for adsorption data classes.
"""

from abc import ABC, abstractmethod

class JsonableBase(ABC):
    """
    Base class for objects that can be serialized to and from JSON,
    to work with the ASE database
    """

    @abstractmethod
    def todict(self) -> dict:
        """Convert class instance to a dictionary, used by ase json and db.
        None values are not allowed in the dictionary and must be removed.
        """

    @classmethod
    @abstractmethod
    def fromdict(cls, dct: dict) -> 'JsonableBase':
        """
        Create an instance of the class from a dictionary.
        Used by xsorb to reconstruct objects after reading from JSON or database.
        """
