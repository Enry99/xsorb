'''
Small module to define an abstract base class for CLI commands.
This is useful to ensure that all CLI commands have the same interface.

It also defines two helper functions to check if a value is a positive integer or float.
'''

from abc import ABC, abstractmethod
import argparse


class CLICommandBase(ABC):
    '''
    Abstract Base class for all CLI commands.
    '''

    @staticmethod
    @abstractmethod
    def add_arguments(parser : argparse.ArgumentParser):
        '''
        Add arguments to the parser for this command.
        This method should be implemented by all subclasses.
        '''

    @staticmethod
    @abstractmethod
    def run(args : argparse.Namespace):
        '''
        Run the command with the given arguments.
        This method should be implemented by all subclasses.
        '''

    @staticmethod
    @abstractmethod
    def bind_function(parser: argparse.ArgumentParser):
        '''
        Bind the run method to the parser.
        This method should be implemented by all subclasses.
        '''


def nonnegative_int(value):
    '''
    Check if a value is a non-negative integer.
    '''
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError(f"{value} is negative")
    return ivalue

def positive_float(value):
    '''
    Check if a value is a positive float.
    '''
    fvalue = float(value)
    if fvalue < 1e-10:
        raise argparse.ArgumentTypeError(f"{value} is not positive")
    return fvalue

class CustomFormatter(argparse.RawDescriptionHelpFormatter,
                      argparse.ArgumentDefaultsHelpFormatter):
    '''
    Combine the formatters to have a more informative help message.
    '''
