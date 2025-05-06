.. _documenting:

======================
Writing documentation
======================

Getting started
===============

The documentation for AMUSE is generated from ReStructured Text using the Sphinx_
documentation generation tool.

It's easiest to install the required packages using Conda:

.. code-block:: bash

    conda install jupyter sphinx make numpy sphinx_rtd_theme matplotlib rst2pdf

.. _Sphinx: http://www.sphinx-doc.org/en/master/

The documentation sources are found in the :file:`doc/` directory in the trunk.
To build the AMUSE documentation in html format, cd into :file:`doc/` and
do::

  make html

You can also pass a ``pdflatex`` flag to ``make`` to build a pdf, or pass no
arguments to show help information.

The output produced by Sphinx can be configured by editing the :file:`conf.py`
file located in the :file:`doc/`.


Organization of the AMUSE documentation
==========================================

The actual ReStructured Text files are kept in :file:`doc/install`,
:file:`doc/developing`, :file:`doc/tutorial` and :file:`doc/interactive_tutorial`. The
main entry point is :file:`doc/index.rst`. The documentation suite is built as a single
document in order to make the most effective use of cross referencing, we want to make
navigating the AMUSE documentation as easy as possible.

Additional files can be added to the various sections by including their base
file name (the .rst extension is not necessary) in the table of contents.
