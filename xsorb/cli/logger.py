'''
Module for configuring logging in the xsorb CLI application.
'''

import logging

def setup_logging():
    '''
    Configures the logging for the xsorb CLI application.
    Sets up logging to both file (xsorb.log) and the console.
    The console output is formatted to show only the messages.
    '''

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter("%(message)s"))

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler("xsorb.log"),
            console_handler
        ]
    )
