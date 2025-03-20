from importlib import import_module
from os import environ


def load_code(community_code, symbol):
    """Load the given symbol from the given package.

    This tries to import the given symbol (the name of a class usually) from the
    amuse-<community_code> package, and return it. If the package is not installed, it
    inspects the environment and produces a helpful error explaining to the user how to
    install it.

    Args:
        community_code: Name of the community code to load
        symbol: Name of the symbol to return
    """
    try:
        interface = import_module(f"amuse_{community_code}.interface")
    except ImportError:
        msg = f"Error: {community_code} is not installed, so we cannot use it."
        if "CONDA_PREFIX" in environ:
            env = environ["CONDA_DEFAULT_ENV"]
            msg += (
                    f" To install {community_code}, go to the terminal, make sure that"
                    f" the '{env}' conda environment is active and that you are in the"
                    " amuse source directory, then install the code using './setup"
                    f" install amuse-{community_code}' and restart your script or"
                    " reload the kernel.")

        if "VIRTUAL_ENV" in environ:
            env = environ["VIRTUAL_ENV"]
            msg += (
                    f" To install {community_code}, go to the terminal, make sure that"
                    f" the '{env}' virtual environment is active and that you are in"
                    " the amuse source directory, then install the code using './setup"
                    f" install amuse-{community_code}' and restart your script or"
                    " reload the kernel.")
        raise ImportError(msg) from None

    try:
        return vars(interface)[symbol]
    except KeyError:
        raise ImportError(
                f"cannot import '{symbol}' from 'amuse.community.{community_code}'"
                ) from None
