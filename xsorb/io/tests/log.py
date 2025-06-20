import logging

consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logging.Formatter("%(message)s"))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("debug.log"),
        consoleHandler
    ]
)


logging.info("Logging initialized")

logging.warning("This is a warning message")