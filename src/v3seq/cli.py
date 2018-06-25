"""
Module that contains the command line app.

Why does this file exist, and why not put this in __main__?

  You might be tempted to import things from __main__ later, but that will cause
  problems: the code will get executed twice:

  - When you run `python -mv3seq` python will execute
    ``__main__.py`` as a script. That means there won't be any
    ``v3seq.__main__`` in ``sys.modules``.
  - When you import __main__ it will get executed again (as a module) because
    there's no ``v3seq.__main__`` in ``sys.modules``.

  Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""
import argparse
import sys

from pkg_resources import (get_distribution, DistributionNotFound)

try:
    __version__ = get_distribution('v3seq').version
except DistributionNotFound:
    __version__ = 'unknown'

files_to_remove = []

# parse command line
parser = argparse.ArgumentParser()
# First define all option groups
group1 = parser.add_argument_group('Input files', 'Required input')
group1.add_argument("-f", "--fastq", default="", type=str, dest="f",
                    help="input reads in fastq format")
group1.add_argument("-k", "--keep", action="store_true",
                    help="keep intermediate files <default: %(default)s>",
                    default=False)
group1.add_argument('-v', '--version', action='version',
                    version=__version__)


# exit so that log file is not written
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()


def main(args=None):
    """What the main does."""
    import logging
    import logging.handlers

    args = parser.parse_args(args=args)

    log_format = '%(levelname)s %(asctime)s %(filename)s: %(funcName)s() %(lineno)d: \t%(message)s'
    logging.basicConfig(filename='v3seq.log', level=logging.INFO, format=log_format, datefmt='%Y/%m/%d %H:%M:%S')
    logging.info(' '.join(sys.argv))

    from v3seq import run
    run.main(args.f)
