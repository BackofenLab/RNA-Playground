# Freiburg RNA algorithms for teaching

This repository contains the Javascript and Html sources for various
algorithms related to RNA structure and RNA-RNA interaction prediction.

###We are currently supporting:

- nested RNA structure prediction
  - nested structure counting
  - maximal base pair structure prediction via the Nussinov algorithm with different recursions
  - base pair probability and unpaired probability computation via a variant of the McCaskill algorithm
  - maximum expected accuracy (MEA) structure prediction
- RNA-RNA interaction prediction
  - co-folding approach via Nussinov-variant
  - maximal intermolecular base pair interaction prediction
  - accessibility-incorporating interaction prediction

The algorithms were developed in order to enable an example-driven learning and teaching of
RNA structure related algorithms. To reduce the level of complexity,
all algorithms use a simple Nussinov-like energy scoring scheme, i.e.
the energy of an RNA structure is directly related to its number 
of base pairs without further distinction.


For local usage, just download or clone the repository content and open the
`index.html` in a recent browser.
