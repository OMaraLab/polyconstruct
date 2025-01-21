.. _customization:

=============
Customization
=============

There are two methods to alter the theme.
The first, and simplest, uses the options exposed through ``html_theme_options`` in ``conf.py``.
This site's options are:

.. literalinclude:: conf.py
   :language: python
   :lines: 104-145
   :lineno-start: 104
   :linenos:

Many of these settings are provided in this site as examples; for a minimalistic example, see :ref:`quickstart`.

Configuration Options
=====================

``mda_official``
   Whether to use official MDAnalysis styling.
   If True, this adds a footer with the `MDAnalysis privacy policy`_
   and other common features. For example, if a logo and a favicon
   are not specified, the MDAnalysis logo and favicon are used.

``sidebar_logo_background``
   The colour of the background of the logo in the sidebar.
   This can be a hex string (e.g. ``'#ffffff'``) or any of
   the named colours (``mdanalysis-orange``, ``mdanalysis-code-orange``,
   ``white``, ``dark-gray``)

``mobile_navbar_background``
   The colour of the top navigation header in mobile format.
   This can be a hex string (e.g. ``'#ffffff'``) or any of
   the named colours (``mdanalysis-orange``, ``mdanalysis-code-orange``,
   ``white``, ``dark-gray``)

``extra_nav_options``
   This is a dictionary of additional navigation links.
   They appear *before* the rest of the tables of contents
   in the sidebar. The keys are the labels of the links,
   and the values are the URLs.



.. _`MDAnalysis privacy policy`: https://www.mdanalysis.org/pages/privacy/

Sidebars
========
You must set ``html_sidebars`` in order for the side bar to appear.
There are four in the complete set.

.. code-block:: python

   html_sidebars = {
       "**": ["globaltoc.html", "localtoc.html", "searchbox.html"]
   }


You can exclude any to hide a specific sidebar.
For example, if this is changed to

.. code-block:: python

   html_sidebars = {
       "**": ["globaltoc.html"]
   }

then only the global ToC would appear on all pages (``**`` is a glob pattern).

Customizing the layout
======================

You can customize the theme by overriding Jinja template blocks.
For example, 'layout.html' contains several blocks that can be overridden or extended.

Place a 'layout.html' file in your project's '/_templates' directory.

.. code-block:: bash

    mkdir source/_templates
    touch source/_templates/layout.html

Then, configure your 'conf.py':

.. code-block:: python

    templates_path = ['_templates']

Finally, edit your override file ``source/_templates/layout.html``:

.. code-block:: jinja

    {# Import the theme's layout. #}
    {% extends '!layout.html' %}

    {%- block extrahead %}
    {# Add custom things to the head HTML tag #}
    {# Call the parent block #}
    {{ super() }}
    {%- endblock %}

New Blocks
==========
The theme has a small number of new blocks to simplify some types of
customization:

``footerrel``
   Previous and next in the footer.
``fonticon``
   Block that contains the icon font. You should probably call ``{{ super() }}`` at the end of the block to include the default icon fonts as well. (Font Awesome and Academicons)

