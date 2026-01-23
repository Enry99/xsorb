'''
Module to download molecular conformers from PubChem.
'''

import logging

from ase.data.pubchem import pubchem_atoms_conformer_search
from ase.io import write

from xsorb.io.filenames import CONFORMERS_FILENAME


def download_pubchem_conformers(identifier_type: str, identifier: str) -> None:
    """Download molecular conformers from PubChem and save them to molecule_conformers.xyz

    Args:
        identifier_type (str): 'name' or 'cid' to specify the type of identifier
        identifier (str): the molecule name or PubChem CID to search for
    """
    if identifier_type == 'name':
        conformers = pubchem_atoms_conformer_search(name=identifier)
    elif identifier_type == 'cid':
        conformers = pubchem_atoms_conformer_search(cid=int(identifier))
    else:
        raise ValueError("identifier_type must be 'name' or 'cid'")

    if not conformers:
        raise ValueError(f'No conformers found for {identifier_type}: {identifier}')

    # Write conformers to an XYZ file
    write(CONFORMERS_FILENAME, conformers)
    logging.info(f'Saved {len(conformers)} conformers to {CONFORMERS_FILENAME}') # pylint: disable=logging-fstring-interpolation
