import logging
import os

import coloredlogs


class Logger:
    """
    A singleton class used to contain the `Logger` object used for all modules in the package.
    """
    _logger = None

    def __new__(cls, verbose: bool = False):
        """
        Creates a new `Logger` object for the class if one does not already exist,
        or if the current process is different from the process that originally created the `Logger`.
        Otherwise, the existing `Logger` is used.

        By default, the logging level is `logging.INFO`.

        Args:
            verbose: The level of logging: (False = `INFO`, True = `DEBUG`).
        """
        pid = os.getpid()
        if cls._logger is None or cls._logger[0] != pid:
            cls._logger = (pid, super().__new__(cls))
            cls._logger[1]._setup_logger(verbose)
        return cls._logger[1]

    @classmethod
    def _setup_logger(cls, verbose: bool = False):
        cls.logger = logging.getLogger(__name__)
        process_str = f" (process {cls._logger[0]}) " if verbose else " "
        log_format = f"%(asctime)s{process_str}[%(levelname)s] %(message)s"
        log_level = logging.DEBUG if verbose else logging.INFO
        coloredlogs.install(fmt=log_format, level=log_level, logger=cls.logger)

    @classmethod
    def reinitialize_logger(cls, verbose: bool = False):
        """
        Reinitialize the logger for child processes in a multiprocessing context.

        Args:
            verbose: The level of logging: (False = `INFO`, True = `DEBUG`).
        """
        cls._logger = None
        cls.__new__(cls, verbose)

    @classmethod
    def set_logger(cls, verbose: bool):
        """
        Mutator method for updating the `Logger` object.

        Args:
            verbose: The level of logging: (False = `INFO`, True = `DEBUG`).
        """
        cls._setup_logger(verbose)

    @classmethod
    def get_logger(cls) -> logging.Logger:
        """
        Accessor method for retrieving the Logger object.

        Returns:
            The `Logger` object contained by the class.

        """

        return cls.logger


logger = Logger().get_logger()
