import logging
import sys


class ColoredFormatter(logging.Formatter):
    blue = "\n\033[94m"
    yellow = "\033[93m"
    red = "\033[91m"
    reset = "\033[0m"
    format = "%(levelname)s: %(message)s"

    FORMATS = {
        logging.INFO: blue + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


def setup_logging(logfile=".python.log") -> None:
    """
    Setup logging configuration of the script
    """
    # a basic config to save logs to metadata.log
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
        filename=logfile,
        filemode="w",
    )

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.WARNING)
    # tell the handler to use colored format
    console.setFormatter(ColoredFormatter())
    # add the handler to the root logger
    logging.getLogger("").addHandler(console)
