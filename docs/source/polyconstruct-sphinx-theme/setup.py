# -*- coding: utf-8 -*-

import re
from setuptools import setup
from pathlib import Path
import shutil
import warnings


def dynamic_author_list():
    """Generate __authors__ from AUTHORS
    This function generates authors.py that contains the list of the
    authors from the AUTHORS file. This avoids having that list maintained in
    several places. Note that AUTHORS is sorted chronologically while we want
    __authors__ in authors.py to be sorted alphabetically.
    The authors are written in AUTHORS as bullet points under the
    "Chronological list of authors" title.
    """
    authors = []
    current_dir = Path(__file__).parent.resolve()
    author_file = current_dir / "AUTHORS.md"

    with open(author_file, "r") as f:
        lines = [x.strip() for x in f.readlines()]

    pattern = r"^- ([\w '.]+) @[\w_-]+$"
    for line in lines:
        match = re.match(pattern, line)
        if match:
            authors.append(match.group(1))

    authors.sort(key=lambda name: name.split()[-1])
    authors.remove("Lily Wang")
    prior_authors = [
        "Matthew Harrigan",
        "Robert T. McGibbon",
        "Christian Schwantes",
        "Martin K. Scherer",
        "Joshua A. Mitchell",
        "Simon Boothroyd",
    ]
    for name in prior_authors:
        authors.remove(name)

    authors = ["Lily Wang"] + authors

    # Write the authors.py file.
    out_path = current_dir / "polyconstruct_sphinx_theme" / "authors.py"
    author_lines = "\n".join([f'    u"{name}",' for name in authors])
    template = f"""\
#-*- coding:utf-8 -*-
# This file is generated from the AUTHORS file during the installation process.
# Do not edit it as your changes will be overwritten.
__authors__ = [
{author_lines}
]
"""

    with out_path.open("w") as f:
        f.write(template)


def copy_mda_assets():
    static_dir = Path('.') / 'polyconstruct_sphinx_theme' / 'static'

    # outdir for logos and icons
    static_logos = static_dir / 'logo'

    # main logo to appear in nav bar
    mda_logo_base_file = 'polyconstruct_logo.png'
    mda_logo = static_dir / 'logo' / mda_logo_base_file

    mda_favicon = static_dir / 'logo' / 'placeholder_favicon.svg'

    if not mda_logo.exists():
        raise FileNotFoundError("Could not find the MDAnalysis logo. "
                                "Check that the branding submodule is "
                                "initialized.")

    if not mda_favicon.exists():
        raise FileNotFoundError("Could not find the MDAnalysis icon. "
                                "Check that the branding submodule is "
                                "initialized.")

    if not static_dir.exists():
        raise FileNotFoundError("Could not find the theme static directory.")

    shutil.copy(mda_logo, static_logos / 'mda_logo.png')
    shutil.copy(mda_favicon, static_logos / 'mda_favicon.ico')


if __name__ == '__main__':
    try:
        dynamic_author_list()
    except (OSError, IOError):
        warnings.warn('Cannot write the list of authors.')

    copy_mda_assets()
    setup()

