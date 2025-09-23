.. _releasing:

================
Making a release
================

Every once in a while, we make a new release of AMUSE, when we feel that enough new
features have accumulated or when important fixes have been made. AMUSE is quite mature
by now, and uses a time-based versioning scheme with the `CalVer <https://calver.org/>`_
pattern vYYYY.MM.MICRO.

Making a release takes three steps: collecting changes, tagging the release, and making
a release on GitHub.


Collecting changes
==================

Every release should come with a description of changes, so that the users know what to
expect. GitHub can make these automatically, but it's nice to check by hand to make sure
nothing and in particular no one (contributor) goes missing.

All the changes are in the git log, so this is a matter of running ``git log`` with your
favourite options and making a list of changes and contributors.


Tagging the release
===================

Each release needs to be tagged in git to make sure we know exactly which version of the
software is being released, and so that the build system knows which version it's
building and installing and can set the metadata correctly.

Git knows two different kinds of tags, annotated tags and lightweight tags. Annotated
tags are intended for releases, and because of the way our build system is set up, they
*must* be used. Tagging the current main version of AMUSE is done using:

.. code-block:: bash

    git switch main
    git pull
    git tag -a -m 'Release v2025.9.0' v2025.9.0


Or whatever the current year and month is. The third number should be used to
distinguish multiple releases within a single month.

The ``-a`` in the final command makes an annotated tag, ``-m`` adds a (required) message
to it, and ``v2025.9.0`` is the name of the tag.

Make sure that you follow exactly this pattern (with a ``v`` at the start and periods,
not underscores or dashes), as the build system relies on it.

If you run

.. code-block:: bash

    git cat-file -t v2025.9.0


it should show ``tag``, not ``commit``.

If the tag has been made correctly, then you can push it to GitHub using

.. code-block:: bash

    git push tag v2025.9.0


Making a GitHub release
=======================

Finally, we will make a GitHub release, which adds the release to GitHub's Releases page
with a nice description, causes GitHub to make a tar.gz archive available for the
release, and archives the release on Zenodo.

Go to `the AMUSE releases page <https://github.com/amusecode/amuse/releases>`_, and
click "Draft a new release" at the top. Select the tag you just made, add a title for
the release (generally just "v2025.9.0" or whatever the version is), and then click
"Generate release notes".

Check your list of changes and contributors from above against what GitHub generated,
and edit as needed. Consider whether the release notes make sense from the perspective
of a user of AMUSE; the developers already know what's going on. In particular, you may
want to organise the changes a bit, starting with new user-visible features,
user-visible fixes, then more technical improvements that don't affect the users' every
day usage.

Once you have something that will help the users quickly see what's new and whether they
should upgrade, you can set the new release as the latest one, don't make it a
pre-release (see below), and don't create a discussion, and when all is ready click
"Publish release" to publish it.


Pre-releases
============

Pre-releases can be a good way of inviting people to test some big changes, giving you a
chance to fix any issues before making the final release. Pre-releases should be tagged
as above, with a version of the form v2025.5.0.rc1 (rc stands for release candidate).
Then a GitHub release can be made and marked as a pre-release.


