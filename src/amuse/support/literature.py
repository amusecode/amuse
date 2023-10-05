import sys
try:
    from docutils import core
except ValueError:
    import os
    import locale
    os.environ['LC_CTYPE'] = 'C'
    os.environ['LANG'] = 'C'
    from docutils import core

import atexit
import shutil
import traceback
import importlib
from os.path import exists

from collections import namedtuple
from docutils import nodes
from amuse.support import exceptions
try:
    from amuse._version import version as amuse_version
except ImportError:
    amuse_version = "unknown version"

import amuse


ClassWithLiteratureReferences = namedtuple(
    "ClassWithLiteratureReferences", 
    "name_of_class_with_refs literature_references_of_class"
)
LiteratureReference = namedtuple(
    "LiteratureReference",
    "id footnote"
)


class TrackLiteratureReferences:
    """
    .. [#] DOI:10.5281/zenodo.1435860
    .. [#] ADS:2018araa.book.....P (Portegies Zwart, S. & McMillan, S.L.W., 2018)
    .. [#] ADS:2013CoPhC.183..456P ** (Portegies Zwart, S. et al., 2013)
    .. [#] ADS:2013A&A...557A..84P ** (Pelupessy, F. I. et al., 2013)
    .. [#] ADS:2009NewA...14..369P (Portegies Zwart, S. et al., 2009)
    """
    INSTANCE = None

    def __init__(self):
        self.registered_classes = set([])
        self.must_show_literature_references_atexit = True
        self.original_excepthook = None

    @classmethod
    def default(cls):
        if cls.INSTANCE is None:
            cls.INSTANCE = cls()
            cls.INSTANCE.register()
        return cls.INSTANCE

    def register(self):
        self.original_excepthook = sys.excepthook
        sys.excepthook = self.exception_hook
        if "--bibtex" in sys.argv:
            atexit.register(self.atexit_bibtex_hook)
        else:
            atexit.register(self.atexit_hook)

    @classmethod
    def suppress_output(cls):
        cls.default().must_show_literature_references_atexit = False

    def register_class(self, cls):
        self.registered_classes.add(cls)

    def exception_hook(self, *arguments):
        # print "exception", arguments, self.original_excepthook
        self.must_show_literature_references_atexit = False
        lines = traceback.format_exception(*arguments)
        # print ''.join(lines)

        self.original_excepthook(*arguments)

    def atexit_bibtex_hook(self):
        if self.original_excepthook is not None:
            sys.excepthook = self.original_excepthook
        self.original_excepthook = None

        if (
            self.must_show_literature_references_atexit
            and "--no-report-references" not in sys.argv
        ):
            texstring = self.all_literature_references_texstring()
            if texstring:
                tex_filename = f"bib-{sys.argv[0]}.tex"
                bib_filename = f"bib-{sys.argv[0]}.bib"
                filenumber = 0
                while exists(tex_filename) or exists(bib_filename):
                    filenumber += 1
                    tex_filename = f"bib-{sys.argv[0]}-{filenumber}.tex"
                    bib_filename = f"bib-{sys.argv[0]}-{filenumber}.bib"

                terminal_message = f"""

In this session you have used the modules below.
Please use the {tex_filename} and {bib_filename} files to include the relevant
citations.

"""
                with open(tex_filename, 'w') as tex_out:
                    tex_out.write(
                        f"{texstring}"
                    )
                shutil.copyfile(amuse.get_data('AMUSE.bib'), bib_filename)
                print(terminal_message)
                print(self.all_literature_references_string())

    def atexit_hook(self):
        if self.original_excepthook is not None:
            sys.excepthook = self.original_excepthook
        self.original_excepthook = None

        if (
            self.must_show_literature_references_atexit
            and "--no-report-references" not in sys.argv
        ):
            string = self.all_literature_references_string()
            if string:
                prefix = """

In this session you have used the AMUSE modules below. 
Please cite any relevant articles:

"""
                print(prefix + string)

    def get_literature_dict_of_class(self, cls):
        """
        get the name and bibkeys from the class, as a dict.
        """
        result = {}
        for current_class in cls.__mro__:
            docstring_in = current_class.__doc__
            if docstring_in:
                if hasattr(current_class, "version"):
                    version = current_class.version()
                else:
                    version = amuse_version
                name = current_class.__name__
                if name.endswith("Interface"):
                    name = "AMUSE-" + name[:-9]
                objectname = f"{name} ({version})"
                doctree = core.publish_doctree(source=docstring_in)
                ref_keys = list(doctree.ids.keys())
                natsort(ref_keys)
                ref_values = [doctree.ids[key] for key in ref_keys]
                for ival in ref_values:
                    if isinstance(ival, nodes.footnote):
                        line = ival.rawsource.split()[0]
                        if (
                            line.startswith('ADS:')
                            or line.startswith('DOI:')
                        ):
                            if objectname in result.keys():
                                result[objectname] += [line[4:]]
                            else:
                                result[objectname] = [line[4:]]
        return result

    def get_literature_list_of_class(self, cls):
        """
        filter the refs from the docstring, if there are no refs nothing is
        appended
        """

        result = []
        for current_class in cls.__mro__:
            docstring_in = current_class.__doc__
            if docstring_in:
                if hasattr(current_class, "version"):
                    version = current_class.version()
                else:
                    version = amuse_version
                name = current_class.__name__
                if name.endswith("Interface"):
                    name = "AMUSE-" + name[:-9]
                objectname = f"{name} ({version})"
                doctree = core.publish_doctree(source=docstring_in)
                ref_keys = list(doctree.ids.keys())
                natsort(ref_keys)
                ref_values = [doctree.ids[key] for key in ref_keys]
                literature_references_of_class = []
                for ikey, ival in zip(ref_keys, ref_values):
                    if isinstance(ival, nodes.footnote):
                        literature_references_of_class.append(
                            LiteratureReference(ikey, ival.rawsource)
                        )
                filled = len(literature_references_of_class) > 0
                if filled:
                    result.append(
                        ClassWithLiteratureReferences(
                            objectname, 
                            literature_references_of_class
                        )
                    )
        return result

    def get_literature_dict(self):
        result = {}
        for x in self.registered_classes:
            if sys.version_info.major == 3 and sys.version_info.minor >= 9:
                result |= self.get_literature_dict_of_class(x)
            else:
                result = {**result, **self.get_literature_dict_of_class(x)}
        return result

    def get_literature_list(self):
        result = []
        for x in self.registered_classes:
            result.extend(self.get_literature_list_of_class(x))
        return result

    def all_literature_references_string(self):
        lines = []
        for s in self.get_literature_list():
            lines.append(
                f'\n\t"{s.name_of_class_with_refs}"'
            )
            for literature_reference_of_class_item in s.literature_references_of_class:
                lines.append(
                    f'\t\t{literature_reference_of_class_item.footnote}'
                )

        lines.append(f'\n\t"AMUSE ({amuse_version})"')
        amuse_list = self.get_literature_list_of_class(type(self))
        for x in amuse_list:
            for literature_reference_of_class_item in x.literature_references_of_class:
                lines.append(
                    f'\t\t{literature_reference_of_class_item.footnote}'
                )

        return "\n".join(lines)

    def all_literature_references_texstring(self):
        result = 'In this article, we used the following AMUSE modules: '
        result += f'AMUSE-framework {amuse_version} \\citep{{'
        amuse_lib = self.get_literature_dict_of_class(type(self))
        for name in amuse_lib.keys():
            for j, key in enumerate(amuse_lib[name]):
                result += f'{key}'
                if j != len(amuse_lib[name]) - 1:
                    result += ', '
        result += '}'

        lib = self.get_literature_dict()
        for name in lib.keys():
            result += (
                f', {name} \\citep{{'
            )
            for j, key in enumerate(lib[name]):
                result += f'{key}'
                if j != len(lib[name]) - 1:
                    result += ', '
            result += '}'
        result += '.\n'

        return result

    def names_of_classes_with_references(self):
        return [x.name_of_class_with_refs for x in self.get_literature_list()]


def literature_references():
    return TrackLiteratureReferences.default().all_literature_references_string()


class LiteratureReferencesMixIn(object):
    def __init__(self):
        self.register_use()

    @classmethod
    def version(cls):
        try:
            version = importlib.import_module(
                '.._version',
                cls.__module__
            ).version
        except (ImportError, ValueError):
            try:
                from amuse.version import version
                version = f"framework {version}"
            except ImportError:
                version = "unknown version"
        return version

    @classmethod
    def print_literature_references(cls):
        print("You are currently using the following codes, which contain literature references")
        print(TrackLiteratureReferences.default().all_literature_references_string())

    @classmethod
    def export2html(cls):
        pass

    @classmethod
    def export2bibtex(cls):
        pass

    @classmethod
    def names_of_classes_with_references(cls):
        return TrackLiteratureReferences.default().names_of_classes_with_references()

    @classmethod
    def all_literature_references_string(cls):
        return TrackLiteratureReferences.default().all_literature_references_string()

    @classmethod
    def register_use(cls):
        TrackLiteratureReferences.default().register_class(cls)


# ------------------------------------------------------------------------------
# from natsort.py: Natural string sorting by Seo Sanghyeon and Connelly Barnes.
# ------------------------------------------------------------------------------

def try_int(s):
    "Convert to integer if possible."
    try:
        return int(s)
    except ValueError:
        return s


def natsort_key(s):
    "Used internally to get a tuple by which s is sorted."
    import re
    return list(map(try_int, re.findall(r'(\d+|\D+)', s)))


def natsort(seq):
    "In-place natural string sort."
    seq.sort(key=natsort_key)
