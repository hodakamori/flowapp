import logging
import os
from datetime import datetime
from logging.handlers import RotatingFileHandler
from typing import Optional, Type


class BlockExecutionLogger:
    _instance: Optional["BlockExecutionLogger"] = None
    _logger: Optional[logging.Logger] = None

    def __new__(cls: Type["BlockExecutionLogger"]) -> "BlockExecutionLogger":
        if cls._instance is None:
            cls._instance = super(BlockExecutionLogger, cls).__new__(cls)
            cls._setup_logger()
        return cls._instance

    @classmethod
    def _setup_logger(cls) -> None:
        if cls._logger is not None:
            return

        logger = logging.getLogger("flowapp")
        logger.setLevel(logging.INFO)

        log_dir = "logs"
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)

        current_date = datetime.now().strftime("%Y%m%d")
        log_file = os.path.join(log_dir, f"flowapp_{current_date}.log")

        file_handler = RotatingFileHandler(
            log_file,
            maxBytes=10 * 1024 * 1024,  # 10MB
            backupCount=5,
            encoding="utf-8",
        )
        file_handler.setLevel(logging.INFO)

        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)

        formatter = logging.Formatter(
            "[%(asctime)s] %(levelname)s [%(filename)s:%(lineno)d] %(message)s"
        )
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)

        logger.addHandler(file_handler)
        logger.addHandler(console_handler)

        cls._logger = logger

    @classmethod
    def get_logger(cls) -> logging.Logger:
        if cls._logger is None:
            cls._setup_logger()
        assert cls._logger is not None
        return cls._logger
