"""General-purpose console logging.

Not to be confused with the agent's reasoning trail (structured log of its
plan, tool calls, and justifications, per CLAUDE.md's "everything is logged"
principle) — that's a separate structured writer built alongside the
step-5 agent loop, since it needs to be machine-readable, not just console output.
"""

import logging

from multiome_agent.config import LOG_LEVEL


def get_logger(name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(LOG_LEVEL)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s %(levelname)-8s %(name)s: %(message)s")
        )
        logger.addHandler(handler)
        logger.propagate = False
    return logger
