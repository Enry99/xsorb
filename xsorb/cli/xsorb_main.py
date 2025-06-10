#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

from xsorb import __version__
from xsorb.cli.xsorb_parser import build_xsorb_parser

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

    #run the command
    try:
        args.func(args)
    except Exception as e: #pylint: disable=broad-except
        if args.traceback:
            raise
        print(f'Error: {e}')
        return

if __name__ == '__main__':
    sys.exit(main())
