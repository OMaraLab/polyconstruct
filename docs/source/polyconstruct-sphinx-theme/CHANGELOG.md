# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<!--
The rules for this file:
  * entries are sorted newest-first.
  * summarize sets of changes - don't reproduce every git log comment here.
  * don't ever delete anything.
  * keep the format consistent:
    * do not use tabs but use spaces for formatting
    * 79 char width
    * YYYY-MM-DD date format (following ISO 8601)
  * accompany each entry with github issue/PR number (Issue #xyz)
-->

## [Unreleased]

### Authors
IAlibay

### Reviewers
<!-- GitHub usernames of reviewers of this release -->

### Added
<!-- New added features -->

### Fixed
<!-- Bug fixes -->

### Changed
- Modernized packaging to be PEP518 compliant.

### Removed
<!-- Removed features -->




## v1.3.1

### Authors
<!-- GitHub usernames of contributors to this release -->
- lilyminium
- IAlibay

### Reviewers
<!-- GitHub usernames of reviewers of this release -->
- lilyminium

### Fixed
- Changed PyPI deployment to use OpenID Connect (PR #72, #77)

### Removed
<!-- Removed features -->
- Upper Sphinx<7 pin (PR #74)



## v1.3.0

### Authors
<!-- GitHub usernames of contributors to this release -->
- lilyminium

### Reviewers
<!-- GitHub usernames of reviewers of this release -->
- IAlibay

### Added
<!-- New added features -->
- Added authors to documentation footer (Issue #43, PR #64)
- Added navigation options (Issue #67, PR #68)

### Fixed
<!-- Bug fixes -->
- Fixed background and navbar issue so they can now be
  specified separately (Issue #66, PR #68)


## v1.2.0

### Authors
<!-- GitHub usernames of contributors to this release -->
- lilyminium

### Reviewers
- IAlibay
- RMeli

### Fixed
<!-- Bug fixes -->
- Fixed mobile rendering by changing the bar default background to gray (Issue #59, PR #60, PR #61)

## v1.1.0

### Authors
<!-- GitHub usernames of contributors to this release -->
- ianmkenney
- ialibay
- lilyminium

### Fixed
<!-- Bug fixes -->
- Fixed documentation and added bibtex docs (PR #53, Issues #10, #26)
- Fixed non-inclusion of versions javascript (PR #55)
- Fixed warning for package directories (Issue #52, PR #56)
- Fix long descriptions issue (Issue #49, PR #51)

### Changed
<!-- Changes in existing functionality -->
- Allow favicon override when using mda_official (PR #37, Issue #35)

## v1.0.1
### Authors
<!-- GitHub usernames of contributors to this release -->
- lilyminium
- IAlibay

### Added
- Atom license details (PR #47)
- Pin to sphinx_rtd_theme (PR #38)

## v1.0.0

### Authors
<!-- GitHub usernames of contributors to this release -->
- lilyminium
- IAlibay
### Fixed
<!-- Bug fixes -->
- Updated path to deploy action (PR #31)
- Remove pytest in deploy action (PR #39)

### Changed
<!-- Changes in existing functionality -->
- Customised theme further (PR #32)

### Removed
<!-- Removed features -->
- Removed deployment from CI check (PR #38)


## v0.1.2

### Authors
<!-- GitHub usernames of contributors to this release -->
- ianmkenney
- IAlibay

### Added
<!-- New added features -->
- Deployment to PyPI on release (PR #27)

### Fixed
<!-- Bug fixes -->
- Fixed path resolution to logo and favicon (PR #22)

### Changed
<!-- Changes in existing functionality -->
- Unpinned Sphinx in dependencies (PR #29)


## v0.1.1

### Authors
<!-- GitHub usernames of contributors to this release -->
- lilyminium

### Fixed
- Fixed metadata to allow for PyPI release

## v0.1.0

### Authors
<!-- GitHub usernames of contributors to this release -->
- lilyminium
- ianmkenney

### Added
<!-- New added features -->
- Meta files (PR #15)
- Conditional content (PR #6)
- ReadTheDocs config (PR #1, PR #2)


### Changed
<!-- Changes in existing functionality -->
- Customised to MDAnalysis colours (PR #4, #5, #9, #13)
