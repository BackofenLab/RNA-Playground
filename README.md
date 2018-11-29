# Freiburg RNA algorithms for teaching

This repository contains the Javascript and Html sources for various
algorithms related to RNA structure and RNA-RNA interaction prediction.
Furthermore, we provide interactive implementations for general
sequence alignment approaches.

All approaches are [available in the Freiburg RNA Teaching webserver](http://rna.informatik.uni-freiburg.de/Teaching/).

### We are currently supporting:

- **Nested RNA structure prediction**
  - [Counting](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=counting) nested structures (see [`NussinovDPAlgorithm_structuresCount`](js/nussinovmatrix.js))
  - [Nussinov](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Nussinov)'s maximal base pair structure prediction algorithm using different recursions (see [`NussinovDPAlgorithm_Unique`, `..Ambiguous`, `..Ambiguous2`, `..MostAmbiguous`](js/nussinovmatrix.js))
  - [McCaskill](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=McCaskill)'s base pair probability and unpaired probability computation (see [`NussinovDPAlgorithm_McCaskill`](js/nussinovmatrix.js))
  - [MEA](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=MEA): maximum expected accuracy structure prediction (see [`DPAlgorithm_MEA`](js/nussinovmatrix.js))
- **RNA-RNA interaction prediction**
  - [co-folding](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=co-folding) approach via Nussinov-variant (see [`DPAlgorithm_coFold`](js/nussinovmatrix.js))
  - [hybrid-only](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=hybrid-only) intermolecular base pair maximization (see [`DPAlgorithm_hybrid`](js/nussinovmatrix4d.js))
  - [accessibility](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=accessibility)-incorporating interaction prediction (see [`DPAlgorithm_rnaup`](js/nussinovmatrix4d.js))
- **General sequence alignment approaches**
  - **pairwise alignment**
    - [Needleman-Wunsch](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Needleman-Wunsch) (global, linear gap costs, quadratic space)  (see [js](js/needleman_wunsch.js))
    - [Hirschberg](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Hirschberg) (global, linear gap costs, linear space)  (see [js](js/hirschberg.js))
    - [Gotoh](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Gotoh) (global, affine gap costs)  (see [js](js/ 	gotoh.js))
    - [Waterman-Smith-Beyer](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Waterman-Smith-Beyer) (global, arbitrary gap costs)  (see [js](js/waterman_smith_beyer.js))
    - [Smith-Waterman](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Smith-Waterman) (local, linear gap costs)  (see [js](smith_waterman.js))
    - [Arslan-Egecioglu-Pevzner](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Arslan-Egecioglu-Pevzner) (local, linear gap costs, normalized score)  (see [js](js/arslan_egecioglu_pevzner.js))
    - [Gotoh (Local)](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Gotoh%20(Local)) (local, affine gap costs)  (see [js](js/gotoh_local.js))
  - **multiple alignment**
    - [Feng-Doolittle (Gotoh-based)](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Feng-Doolittle) (progressive, affine gap costs)  (see [js](js/feng_doolittle.js))
    - [Iterative Refinement (Feng-Doolittle-based)](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Iterative%20Refinement) (progressive, affine gap costs, iterative refinement)  (see [js](js/iterative_refinement.js))
    - [Notredame-Higgins-Heringa (t-coffee)](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Notredame-Higgins-Heringa) (progressive, data-derived scoring)  (see [js](js/notredame_higgins_heringa.js))
- **Agglomerative Clustering**  (see [js](js/agglomerative_clustering.js))
  - [Complete Linkage](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Agglomerative%20Clustering) (Furthest Neighbour)
  - [Neighbour-Joining](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Agglomerative%20Clustering) (Nearest Neighbour on clusters)
  - [Single Linkage](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Agglomerative%20Clustering)  (Nearest Neighbour on elements between two different clusters)
  - [Unweighted Pair Group Method with Arithmetic Mean (UPGMA)](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Agglomerative%20Clustering) (Group Average or Average Linkage)
  - [Weighted Pair Group Method with Arithmetic Mean (WPGMA)](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Agglomerative%20Clustering) (Simple Average)


  
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



