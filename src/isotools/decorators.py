import functools
import logging

logger = logging.getLogger('isotools')

# warn unsave functions


def deprecated(func):
    """Warns about use of deprecated function"""
    @functools.wraps(func)
    def wrapper_depreciated(*args, **kwargs):
        logger.warning(f"Calling deprecated function {func.__name__}")
        value = func(*args, **kwargs)
        return value
    return wrapper_depreciated


def experimental(func):
    """Informs about use of untested functionality"""
    @functools.wraps(func)
    def wrapper_experimental(*args, **kwargs):
        logger.warning(f"Calling {func.__name__}, which is untested/experimental")
        value = func(*args, **kwargs)
        return value
    return wrapper_experimental

# helpers for debugging


def debug(func):
    """Print the function signature and return value"""
    @functools.wraps(func)
    def wrapper_debug(*args, **kwargs):
        args_repr = [repr(a) for a in args]
        kwargs_repr = [f"{k}={v!r}" for k, v in kwargs.items()]
        signature = ", ".join(args_repr + kwargs_repr)
        logger.info(f"Calling {func.__name__}({signature})")
        value = func(*args, **kwargs)
        logger.info(f"{func.__name__!r} returned {value!r}")
        return value
    return wrapper_debug


def traceback(func):
    """In case of exception, print the arguments"""
    @functools.wraps(func)
    def wrapper_try(*args, **kwargs):
        try:
            value = func(*args, **kwargs)
            return value
        except Exception as e:
            args_repr = [repr(a) for a in args]
            kwargs_repr = [f"{k}={v!r}" for k, v in kwargs.items()]
            signature = ", ".join(args_repr + kwargs_repr)
            logger.info(f"Exception during call {func.__name__}({signature})")
            raise e
    return wrapper_try
