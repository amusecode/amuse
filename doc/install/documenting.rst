.. _documenting:

======================
Writing documentation
======================

Getting started
===============

The documentation for AMUSE is generated from ReStructured Text 
using the Sphinx_ documentation generation tool. Sphinx version 1.0 
or later is required. You might still run into problems, so most 
developers work from the sphinx source repository (Mercurial based) 
because it is a rapidly evolving project:

.. code-block:: bash

  > hg clone http://bitbucket.org/birkenfeld/sphinx/
  > cd sphinx
  > python setup.py install

.. _Sphinx: http://www.sphinx-doc.org/en/master/

The documentation sources are found in the :file:`doc/` directory in the trunk.
To build the AMUSE documentation in html format, cd into :file:`doc/` and
do::

  make html

You can also pass a ``pdflatex`` flag to make to build a pdf, or pass no
arguments to show help information.

The output produced by Sphinx can be configured by editing the :file:`conf.py`
file located in the :file:`doc/`.


Organization of the AMUSE documentation
==========================================

The actual ReStructured Text files are kept in :file:`doc/install`,
:file:`doc/design`, :file:`doc/tutorial`. The main entry point is
:file:`doc/index.rst`. The documentation suite is
built as a single document in order to make the most effective use of cross
referencing, we want to make navigating the AMUSE documentation as easy as
possible.

Additional files can be added to the various sections by including their base
file name (the .rst extension is not necessary) in the table of contents.
