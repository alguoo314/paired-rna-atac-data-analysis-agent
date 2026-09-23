import logging

from multiome_agent import config
from multiome_agent.logging_utils import get_logger


def test_returns_logger():
    logger = get_logger("test.logger.a")
    assert isinstance(logger, logging.Logger)


def test_respects_log_level():
    logger = get_logger("test.logger.b")
    assert logger.level == logging.getLevelName(config.LOG_LEVEL)


def test_no_duplicate_handlers_on_repeat_calls():
    name = "test.logger.c"
    logger1 = get_logger(name)
    n_handlers = len(logger1.handlers)
    logger2 = get_logger(name)
    assert logger1 is logger2
    assert len(logger2.handlers) == n_handlers
