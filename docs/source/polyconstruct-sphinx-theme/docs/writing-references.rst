==================
Writing references
==================

If you want to add references to your documentation,
you can use the `sphinxcontrib-bibtex`_ extension.

Please refer to the `sphinxcontrib-bibtex`_ documentation
for full details. However, some notes and tips are provided
below for compatibility with the ``mdanalysis-sphinx-theme``.

Basic use
=========

Setting up the reference file
-----------------------------

Generally references should be provided as a BibTeX ``.bib`` file.
In MDAnalysis, this is traditionally called ``references.bib`` and
is placed in the same directory as the ``conf.py``.
For an example, see the `MDAnalysis reference file`_.

The extension and the ``references.bib`` file should then added to
the ``conf.py`` file as below:

.. code-block:: python

   extensions = [
       'sphinxcontrib.bibtex',
   ]

   bibtex_bibfiles = ['references.bib']


Citing references in text
-------------------------

References can be cited in text using ``:cite:`` or ``:footcite:``.
Which one to use depends on if you want to display a **"local"**
or a **"global"** bibliography. A **global** bibliography
is a single reference list that is rendered on its own page, with the
following directive:

.. code-block:: rst

   .. bibliography::


A ``bibliography`` directive picks up all references cited using
the ``:cite:`` directive across *all* the documentation.


A **local** bibliography is one that is rendered on a particular
page with other documentation. This is done using the ``:footcite:``
and with the following directive:

.. code-block:: rst

   .. footbibliography::


A ``footbibliography`` directive only displays the references
cited on that particular page for display.

.. important::

    If you are trying to render a ``.. bibliography::`` on a
    page with other documentation, it **must** be the *last*
    item on the page. Otherwise, it will displace any other
    content that follows. In general, it is safest to use a
    ``.. footbibliography::`` if there is any other content on the page.
    Please see the `mdanalysis-sphinx-theme Issue #25`_ for more details.


.. _`sphinxcontrib-bibtex`: https://sphinxcontrib-bibtex.readthedocs.io/en/latest/
.. _`MDAnalysis reference file`: https://github.com/MDAnalysis/mdanalysis/blob/develop/package/doc/sphinx/source/references.bib
.. _`mdanalysis-sphinx-theme Issue #25`: https://github.com/MDAnalysis/mdanalysis-sphinx-theme/issues/25