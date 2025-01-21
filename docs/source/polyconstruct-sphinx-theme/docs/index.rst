.. _quickstart:

==================================
MDAnalysis Sphinx theme Quickstart
==================================

This theme provides a responsive theme for Sphinx documentation by the Open Force Field Initiative.
It is inspired by `Material for Sphinx <https://bashtage.github.io/sphinx-material/>`_ and `Material for MkDocs <https://squidfunk.github.io/mkdocs-material/>`_, but has been rewritten with the `Bulma <https://bulma.io>`_ CSS framework to remove any JavaScript dependencies.

Getting Started
---------------
Install from git

.. code-block:: shell-session

   $ pip install git+https://github.com/mdanalysis/mdanalysis-sphinx-theme.git@main

Or, add to your ReadTheDocs environment.yml

.. code-block:: yaml

    dependencies:
        - pip
        - <conda dependency>
        - <conda dependency>
        - <conda dependency>
        # --- snip --- #
        - pip:
            - git+https://github.com/mdanalysis/mdanalysis-sphinx-theme.git@main

Update your ``conf.py`` with the required changes:

.. code-block:: python

    extensions.append("mdanalysis_sphinx_theme")
    html_theme = "mdanalysis_sphinx_theme"
    html_sidebars = {"**": ["globaltoc.html", "localtoc.html", "searchbox.html"]}


There are a lot more ways to customize this theme. See :ref:`Customization`
or ``theme.conf`` for more details.

.. code-block:: python

    # Enable the theme itself
    extensions.append("mdanalysis_sphinx_theme")
    html_theme = "mdanalysis_sphinx_theme"

    # (Optional) Logo in PNG format.
    # If not provided and mda_official is False, no logo will be added.
    # If not provided and mda_official is True, the MDAnalysis logo will be used.
    html_logo = "_static/logo/placeholder_logo.png"

    # (Optional) favicon.
    # If not provided and mda_official is False, will default to a placeholder favicon.
    # If not provided and mda_official is True, will default to the MDAnalysis favicon.
    html_favicon = "_static/logo/placeholder_favicon.svg"

    # Theme options are theme-specific and customize the look and feel of a
    # theme further.
    html_theme_options = {
        # ===== mdanalysis-sphinx-theme options =====
        # whether to apply official MDAnalysis styling
        # e.g. using the official MDAnalysis logo and favicon
        # and using the MDAnalysis privacy policy
        "mda_official": False,
        # The background colour of the logo area in the navigation bar
        "sidebar_logo_background": "white",
        # The background colour of the top navigation bar on mobile
        "mobile_navbar_background": "dark-gray",
        # Extra navigation links to show on the sidebar, before the table of contents
        "extra_nav_links": extra_nav_links,

        # ===== inherited options =====
        # For more details, please see
        # https://sphinx-rtd-theme.readthedocs.io/en/stable/
        # Default values are shown below
        
        # Only display logo and not the project name on sidebar
        "logo_only": True,
        # Display the version number on the sidebar
        "display_version": True,
        # Where to display "next" and "previous" buttons
        "prev_next_buttons_location": "bottom",
        # Add an icon next to external links
        "style_external_links": False,
        # If enabled, navigation entries are not expandable
        "collapse_navigation": True,
        # If enabled, the navigation bar scrolls with the main page
        "sticky_navigation": True,
        # Maximum depth of the contents tree -- set to -1 for unlimited
        "navigation_depth": 4,
        # Whether to include hidden toctrees in the navigation bar
        "includehidden": True,
        # Whether to hide page subheadings from navigation
        "titles_only": False,

    }

    # Custom sidebar templates, must be a dictionary that maps document names
    # to template names.
    html_sidebars = {
        # By default, show everything
        "**": ["globaltoc.html", "localtoc.html", "searchbox.html"]
    }


.. toctree::
    :caption: Basic Use
    :maxdepth: 1

    customization
    writing-references
    specimen
    additional_samples
    subpage/index


.. toctree::
    :caption: Other Examples and Uses
    :maxdepth: 1

    pymethod
    numpydoc_example
    autodoc_pydantic_example
    mdanalysis
    notebook.ipynb
    rst-cheatsheet/rst-cheatsheet
    primer

.. toctree::
    :caption: Changes and License
    :maxdepth: 1

    change-log
    license



Index
~~~~~
:ref:`genindex`
