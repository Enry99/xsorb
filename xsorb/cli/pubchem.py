'''
CLI parser for command: pubchem
'''

import argparse

from xsorb.cli.command import CLICommandBase


class CLICommand(CLICommandBase):
    """download molecular conformers from PubChem, saving the to molecule_conformers.xyz
Example:
 $ xsorb pubchem name 1-Hexene
 $ xsorb pubchem cid 1234
    """

    @staticmethod
    def add_arguments(parser : argparse.ArgumentParser):
        parser.add_argument('identifier_type',
                        type=str,
                        choices=['name', 'cid'],
                        help='''specify whether the identifier is a name or a PubChem CID''')
        parser.add_argument('identifier',
                        type=str,
                        help='''the molecule name or PubChem CID to search for''')


    @staticmethod
    def run(args : argparse.Namespace):
        from xsorb.io.pubchem import download_pubchem_conformers # pylint: disable=import-outside-toplevel
        download_pubchem_conformers(args.identifier_type, args.identifier)


    @staticmethod
    def bind_function(parser: argparse.ArgumentParser):
        parser.set_defaults(func=CLICommand.run)
