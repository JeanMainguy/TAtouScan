# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.1.1] – 2025-04-29

### Fixed

- Use of a list of attributes to map protein id and gene in the GFF.

## [v0.1.0] – 2025-04-29

### Added

- Support for parsing genome annotations in **GFF** format and protein sequences in **FAA** format.
- TA annotation of CDS using **pyhmmer** and curated **HMM profiles**.
- Grouping of CDS hits into putative **toxin-antitoxin systems** based on genomic proximity.
- Export of annotated CDS hits and metadata to a **TSV file**.
- Command-line interface (**CLI**) with arguments for input files, thresholds, and output configuration.
