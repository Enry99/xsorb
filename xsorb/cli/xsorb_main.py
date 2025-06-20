#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import logging

from xsorb import __version__
from xsorb.cli.xsorb_parser import build_xsorb_parser
from xsorb.cli.logger import setup_logging

def main():
    '''
    Main entry point for the xsorb command line interface.
    '''

    print(f'xsorb version {__version__}')

    if len(sys.argv) == 1:
        print("No command provided. The program will now terminate.")
        return

    #parse the command line arguments
    parser = build_xsorb_parser()
    args = parser.parse_args()

    setup_logging()

    #run the command
    try:
        args.func(args)
    except Exception as e: #pylint: disable=broad-except
        if args.traceback:
            raise
        logging.critical('Error: %s', e)
        return

if __name__ == '__main__':
    sys.exit(main())
