
from abc import ABC, abstractmethod
import argparse


class CLICommandBase(ABC):
    '''
    Abstract Base class for all CLI commands.
    '''

    @staticmethod
    @abstractmethod
    def help():
        ...

    @staticmethod
    @abstractmethod
    def add_arguments(parser : argparse.ArgumentParser):
        ...

    @staticmethod
    @abstractmethod
    def run(args : argparse.Namespace):
        ...

    @staticmethod
    @abstractmethod
    def bind_function(parser: argparse.ArgumentParser):
        ...
