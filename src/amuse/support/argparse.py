import argparse


def new_argument_parser(
    **kwargs
):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_amuse_arguments(parser, **kwargs)
    return parser


def add_amuse_arguments(
    parser,
    literature=True,
    units=False,
):
    if literature:
        # Note: at the moment these are dummy arguments, with
        # amuse.support.literature parsing sys.argv on its own.
        # They are therefore usable as command line arguments, but the values
        # here are not used further.
        parser.add_argument(
            '--no-report-references',
            dest='no_report_references',
            action='store_true',
            default=False,
            help="don't report literature references to stdout",
        )
        parser.add_argument(
            '--bibtex',
            dest='create_bibtex',
            action='store_true',
            default=False,
            help="create bibtex literature file",
        )
    if units:
        # These need to be parsed to be usable.
        parser.add_argument(
            '--length',
            dest='unit_length',
            type=str,
            default='parsec',
            help="length unit",
        )
        parser.add_argument(
            '--mass',
            dest='unit_mass',
            type=str,
            default='MSun',
            help="mass unit",
        )
        parser.add_argument(
            '--time',
            dest='unit_time',
            type=str,
            default='Myr',
            help="time unit",
        )
        parser.add_argument(
            '--speed',
            dest='unit_speed',
            type=str,
            default='kms',
            help="speed unit",
        )
    return parser
