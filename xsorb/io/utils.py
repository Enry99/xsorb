#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Enrico Pedretti

'''
General utility functions for the I/O module
'''

from __future__ import annotations
import time
import math

from ase.io import read, write

from xsorb.ase_custom import write_xyz_custom, AtomsCustom


def ase_custom_read(filename, **kwargs):
    """
    Modified version of ase.io.read that returns an ase_custom.AtomsCustom object
    """
    atoms_or_atoms_list = read(filename=filename, **kwargs)
    if isinstance(atoms_or_atoms_list, list):
        return [AtomsCustom(at) for at in atoms_or_atoms_list]
    else:
        return AtomsCustom(atoms_or_atoms_list)


def overwrite_question(file_path : str) -> str:
    '''
    Asks the user if they want to overwrite the file at the specified path.

    Args:
    - file_path: str, path to the file

    Returns:
    - str, 'y' if the user wants to overwrite, 'n' otherwise,
        'yall' if the user wants to overwrite all files, 'nall' otherwise
    '''

    while True:
        answer = input(f'{file_path} already exists. Overwrite? '\
                       '("y" = yes to this one, "yall" = yes to all,'\
                         ' "n" = no to this one, "nall" = no to all): ')
        if answer.lower() in ['y', 'n', 'yall', 'nall']:
            return answer
        else:
            print('Value not recognized. Try again.')


def continue_even_if_not_all_completed_question() -> bool:
    '''
    Asks the user if they want to continue even if not all calculations are completed.

    Returns:
    - bool, True if the user wants to continue, False otherwise
    '''

    while True:
        answer = input('Not all calculations are completed. '\
                       'Continue anyway with those that are present? (y/n): ')
        if answer.lower() in ['y', 'n']:
            return answer.lower() == 'y'
        else:
            print('Value not recognized. Try again.')



def progressbar(it, prefix="", size=50):
    '''
    Simple progress bar for any iterable, similar to tqdm, but without
    the need of external dependencies.
    It is a generator function, so it will yield the items from the iterable.
    Adapted from https://stackoverflow.com/a/3160819

    Parameters
    ----------
    it : iterable
        The iterable to be iterated over.
    prefix : str
        Prefix to be shown in the progress bar.
    size : int
        Size of the progress bar.
    '''
    count = len(it)
    spacing = int(math.log10(count+0.1)) + 1
    start = time.time() # time estimate start

    def show(j):
        x = int(size*j/count)

        if j != 0:
            # time estimate calculation and string
            time_per_iter = (time.time() - start) / j
            remaining = time_per_iter * (count - j) #time/iter * remaining iterations
            mins, sec = divmod(remaining, 60) # limited to minutes
            time_str = f"{int(mins):02}:{sec:03.1f} ({time_per_iter:3.2f}s/it)"
        else:
            time_str = "?"

        print(f"{prefix}[{'|'*x}{('.'*(size-x))}] {j:>{spacing}}/{count} ETA: {time_str}",
              end='\r',flush=True)

    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    print("\n", flush=True)
