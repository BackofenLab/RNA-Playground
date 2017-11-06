# Freiburg RNA algorithms for teaching

This repository contains the Javascript and Html sources for various
algorithms related to RNA structure and RNA-RNA interaction prediction.
Furthermore, we provide interactive implementations for general
sequence alignment approaches.

All approaches are [available in the Freiburg RNA Teaching webserver](http://rna.informatik.uni-freiburg.de/Teaching/).

### We are currently supporting:

- **Nested RNA structure prediction**
  - nested structure [counting](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=counting)
  - maximal base pair structure prediction via the [Nussinov](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Nussinov) algorithm with different recursions
  - base pair probability and unpaired probability computation via a variant of the [McCaskill](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=McCaskill) algorithm
  - maximum expected accuracy ([MEA](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=MEA)) structure prediction
- **RNA-RNA interaction prediction**
  - [co-folding](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=co-folding) approach via Nussinov-variant
  - maximal [hybrid-only](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=hybrid-only) intermolecular base pair interaction prediction
  - [accessibility](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=accessibility)-incorporating interaction prediction
- **General sequence alignment approaches**
  - [Needleman-Wunsch](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Needleman-Wunsch) (global, linear gap costs)
  - [Gotoh](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Gotoh) (global, affine gap costs)
  - [Waterman-Smith-Beyer](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Waterman-Smith-Beyer) (global, arbitrary gap costs)
  - [Smith-Waterman](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Smith-Waterman) (local, linear gap costs)
  
Each algorithm comes with a short introduction with links to related
literature, the according dynamic programming recursions 
and a short explanation of the algorithms.

The algorithms were developed in order to enable an example-driven learning and teaching of
RNA structure related algorithms. To reduce the level of complexity,
all algorithms use a simple Nussinov-like energy scoring scheme, i.e.
the energy of an RNA structure is directly related to its number 
of base pairs without further distinction.

For all approaches an exhaustive enumeration of optimal solutions is provided. 
Some interfaces also allow for suboptimal solution enumeration.
Predicted structures/interactions are graphically visualized and where possible
an interactive traceback visualization is enabled. 

For local usage, just download or clone the repository content and open the
`index.html` (structure/interaction prediction) or 
`alignment.html` (sequence alignment approaches) 
in a recent browser.



