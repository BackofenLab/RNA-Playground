/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

/**
 * Used to store all constants, structs and lists used in the program.
 * Namespaces, HTML class-names and events like "mouse-over" etc. aren't stored in this file
 * because it would slow down the development process.
 * Multiple constants of same "logic" (for example with similar naming) stored in structures.
 */

/* Hint: alphabetically ordered (also within definitions) and structures (i.e. {}) after constants. */

/*---- CONSTANTS ----*/
var END_SO_ON = "..";
var EPSILON = 0.000000001;  // some very low number to test against

var FENG_DOOLITTLE_CONSTANT = 0.001;  // Hint: it has not to be used 0.001, but it is the paper constant!

var HIRSCHBERG_LOWER_NODE = 2;  // name for a node of the computation tree
var HIRSCHBERG_UPPER_NODE = 1;  // name for a node of the computation tree

var MAX_NUMBER_ITERATIONS = 5;  // number of iterations in algorithm with convergence
var MAX_NUMBER_TRACEBACKS = 10;  // stores the number of tracebacks after which an alignment algorithm stops to compute
var MAX_TRACE_FUNCTION_ARG_LEN = 10;  // tells the allowed length of an argument to avoid a string which goes over the page border

var REUPDATE_TIMEOUT_MS = 100;  // time in ms after which new LaTeX-Code is reinterpreted or outputs updated
var REACTION_TIME_HIGHLIGHT = REUPDATE_TIMEOUT_MS + 50;  // to highlight tracebacks only after outputs have definitely been updated

var SMITH_WATERMAN_STOP = "0";

/*---- STRUCTS ----*/
/**
 * Stores the implemented algorithm names.
 */
var ALGORITHMS = {  // contains a list of algorithms with an own javascript-file (javascript names without extension)
    AGGLOMERATIVE_CLUSTERING: "agglomerative_clustering",
    ARSLAN_EGECIOGLU_PEVZNER: "arslan_egecioglu_pevzner",
    FENG_DOOLITTLE: "feng_doolittle",
    GOTOH: "gotoh",
    GOTOH_LOCAL: "gotoh_local",
    HIRSCHBERG: "hirschberg",
    ITERATIVE_REFINMENT: "iterative_refinement",
    NEEDLEMAN_WUNSCH: "needleman_wunsch",
    NEIGHBOUR_JOINING: "neighbour_joining",
    NONE: "none",  // this is not an algorithm :)
    NOTREDAME_HIGGINS_HERINGA: "notredame_higgins_heringa",
    SMITH_WATERMAN: "smith_waterman",
    WATERMAN_SMITH_BEYER: "waterman_smith_beyer"
};

/**
 * Stores the types of calculation.
 */
var ALIGNMENT_TYPES = {  // contains a list of alignment types
    DISTANCE: "distance",
    SIMILARITY: "similarity"
};

/**
 * Stores HTML-code for arrows.
 */
var ARROWS = {  // Hint: inner quotes have to be this here: " " or it won't work!
    LEFT: '<div class="arrows_l"></div>',
    LEFT_NAME: ".arrows_l",
    TOP: '<div class="arrows_t"></div>',
    TOP_NAME: ".arrows_t",
    DIAGONAL: '<div class="arrows_d"></div>',
    DIAGONAL_NAME: ".arrows_d"
};

var CELL_PERCENT = {
    LINE: 0.2,  // position of a "between-table"-arrow in the Gotoh algorithm
    LINE_2: 0.1,  // position of another "between-table"-arrow in the Gotoh algorithm (leads to more clarity)
    LINE_HEAD_PENETRATION: 0.1  // tells how much a line-head of a long "in-table"-arrow penetrates into a cell
};

/**
 * Defines allowed and disallowed input characters with regular expressions (regex).
 */
var CHARACTER = {
    BASE: /[a-zA-Z]/i,
    BASES: /^[a-zA-Z-]+$/,
    CSV_SYMBOLS: /^[0-9;\n\s]+$/,
    NON_BASES: /[^a-zA-Z-]+/g,  // g to replace globally
    NON_CSV_SYMBOLS: /[^0-9;\n\s]+/g,
    NUMBER: /[0-9]/,
    NUMBERS: /[-+]?[0-9]+\.[0-9]*/
};

/**
 * Stores the displayed clustering algorithm names.
 */
var CLUSTERING_ALGORITHMS = {
    COMPLETE_LINKAGE: "Complete-Linkage",
    NEIGHBOUR_JOINING: "Neighbour-Joining", /* Hint: Change name also in the Agglomerative-Clustering.html or it won't work anymore! */
    SINGLE_LINKAGE: "Single-Linkage",
    UPGMA: "Unweighted PGMA",
    WPGMA: "Weighted PGMA"
};

/**
 * Allowed file extensions.
 */
var FILE_EXTENSIONS = {
    HYPERTEXT_MARKUP_LANGUAGE: ".html",
    JAVASCRIPT: ".js"
};

/**
 * Stores all errors which can be displayed to the user.
 */
var ERRORS = {
    DIFFERENT_NUMBER_OF_COLUMNS_AND_ROWS: "The number of columns and rows should be equal!",
    WRONG_NUMBER_OF_COLUMNS_IN_ROW: "wrong number of columns in row: "
};

/**
 * Allowed max values for inputs.
 */
var INPUT = {
    ALIGNMENTS_MAX: 10,
    ALIGNMENTS_MIN: 1,
    LENGTH_MAX: 1000,
    LENGTH_MIN: 0,
    MAX: 10,  // abs: absolute value
    MIN: -10
};

var ITERATIVE_REFINEMENT_ORDERS = {
    INPUT_ORDER: "input-sequences order",
    LEFT_FIRST: "guide-tree left-to-right traversal",
    RIGHT_FIRST: "guide-tree right-to-left traversal"
};

var ITERATIVE_REFINEMENT_STRATEGIES = {
    MIN_DISTANCE_PAIR: "Minimum Distance Pair", /* Hint: Change name also in the Iterative-Refinement.html or it won't work anymore! */
    ONE_VS_ALL: "One-vs-All", /* Hint: Change name also in the Iterative-Refinement.html or it won't work anymore! */
    PAIRWISE_BEST_PAIR: "Pairwise Best Pair"
};

/**
 * Stores key codes from keyboard.
 */
var KEY_CODES = {
    ENTER: 13
};

/**
 * Stores LaTeX-code for displayed formulas.
 */
var LATEX = {
    ALIGNED_PLUS: "& + &",
    ALPHA: "\\alpha",
    BEGIN_CASES: "\\begin{cases}",
    BETA: "\\beta",
    CLOSE: "}",
    CURLY_BRACKET_LEFT: "{",
    CURLY_BRACKET_RIGHT: "}",
    DOT: "\\cdot",
    END_CASES: "\\end{cases}",
    FACTOR: " k",
    LN: "\\ln",
    MATH_REGION: "$",
    MAX: "\\max",
    MIN: "\\min",
    NEGATIVE_INFINITY: "-$\\infty$",
    NEW_LINE: "\\\\",
    POSITIVE_INFINITY: "$\\infty$",
    POW2: "^2",
    SPACE: "\\phantom{-}",
    SPACE_SMALL: "\\;",
    SUBORDINATE: "_",
    SUM: "\\sum",
    SUPERSCRIPT: "^",
    TEXT_START: "\\texttt{",

    FORMULA: {
        /* contains subparts of a formula */
        CURRENT: "D_{i,j}",
        CURRENT_BACKWARD: "D'_{i,j}",
        CURRENT_P: "P_{i,j}",
        CURRENT_Q: "Q_{i,j}",
        D: "D \\,",
        D_BIG: /D/g,
        D_BIG_UNDERSCORE: /D_/g,
        D_PRIME: "D'",
        D_PURE: "D",
        D_STAR: "M",
        D_PRIME_UNDERSCORE: "D'_",
        DELETION: "b_j = -",
        DIAGONAL: "D_{i-1,j-1}",
        DPM: "\\text{Trace}",
        I_IS: "i = ",
        I_MINUS_ONE: /i-1/g,
        I_PLUS_ONE: "i+1",
        INSERTION: "a_i = -",
        J_MINUS_ONE: /j-1/g,
        J_PLUS_ONE: "i+1",
        GAP: "g(k)",
        LEFT: "D_{i,j-1}",
        LEFT_Q: "Q_{i,j-1}",
        MATCH: "a_i = b_j",
        MAXIMIZE_HORIZONTAL: "\\displaystyle \\max_{1 \\leq k \\leq j} \\{D_{i,j-k} & + & g(k) \\}",
        MAXIMIZE_VERTICAL: "\\displaystyle \\max_{1 \\leq k \\leq j} \\{D_{i-k,j} & + & g(k) \\}",
        MINIMIZE_HORIZONTAL: "\\displaystyle \\min_{1 \\leq k \\leq j} \\{D_{i,j-k} & + & g(k) \\}",
        MINIMIZE_VERTICAL: "\\displaystyle \\min_{1 \\leq k \\leq j} \\{D_{i-k,j} & + & g(k) \\}",
        MISMATCH: "a_i \\neq b_j",
        SEQ_A_I: /a_i/g,
        SEQ_A_I_PLUS_1: "a_{i+1}",
        SEQ_B_J: /b_j/g,
        SEQ_B_J_PLUS_1: "b_{j+1}",
        S_BIG: "S",
        TOP: "D_{i-1,j}",
        TOP_P: "P_{i-1,j}",
        ZERO: "0"
    },

    FORMULAS: {
        /* contains complete formulas */
        COMPLETE_LINKAGE_DELTA_1: "$$\\delta_m(k, ij) = \\max_{s \\in k,\\, t \\in ij} D_{s,t}$$",
        COMPLETE_LINKAGE_DELTA_2: "$$\\delta_t(i, ij) = \\delta_t(j, ij) = \\frac{D_{min}}{2}$$",
        NEIGHBOUR_JOINING_DELTA_1:
            "$$\\delta_m(k, ij) = \\frac{D_{k,i} + D_{k,j} - D_{i,j}}{2}$$",
        NEIGHBOUR_JOINING_DELTA_2:
            "$$\\delta_t(i, ij) = \\frac{D_{i,j} - \\Delta_{i,j}}{2}, \\quad \\delta_t(j, ij) = \\frac{D_{i,j} + \\Delta_{i,j}}{2} = D_{i,j} - \\delta_t(i, ij)$$",
        SINGLE_LINKAGE_DELTA_1: "$$\\delta_m(k, ij) = \\min_{s \\in k,\\, t \\in ij} D_{s,t}$$",
        SINGLE_LINKAGE_DELTA_2: "$$\\delta_t(i, ij) = \\delta_t(j, ij) = \\frac{D_{min}}{2}$$",
        UPGMA_DELTA_1:
        "$$\\delta_m(k, ij) = \\frac{D_{k,i} \\cdot |i| + D_{k,j} \\cdot |j|}{|i| + |j|} " +
        "= \\frac{1}{|k| |ij|} \\sum_{s \\in k} \\sum_{t \\in ij} D_{s,t}$$",
        UPGMA_DELTA_2:
            "$$\\delta_t(i, ij) = \\delta_t(j, ij) = \\frac{D_{min}}{2}$$",
        WPGMA_DELTA_1:
            "$$\\delta_m(k, ij) = \\frac{D_{k,i} + D_{k,j}}{2}$$",
        WPGMA_DELTA_2:
            "$$\\delta_t(i, ij) = \\delta_t(j, ij) = \\frac{D_{min}}{2}$$"
    },

    RECURSION: {
        GOTOH:
        "\\begin{cases}" +
        "D_{i-1,j-1}	& + & s(a_i,b_j)    \\\\" +
        "P_{i,j}                            \\\\" +
        "Q_{i,j}                                " +
        "\\end{cases}",

        GOTOH_LOCAL:
        "\\begin{cases}" +
        "S_{i-1,j-1}	& + & s(a_i,b_j)    \\\\" +
        "P_{i,j}                            \\\\" +
        "Q_{i,j}                            \\\\" +
        "0                                      " +
        "\\end{cases}",

        GOTOH_P:
        "\\begin{cases}" +
        "D_{i-1,j}      & + & g(1)		    \\\\" +
        "P_{i-1,j}      & + & \\beta" +
        "\\end{cases}",

        GOTOH_Q:
        "\\begin{cases}" +
        "D_{i,j-1}      & + & g(1)		    \\\\" +
        "Q_{i,j-1}      & + & \\beta" +
        "\\end{cases}",

        HIRSCHBERG_BACKWARD:
        "\\begin{cases}" +
        "D'_{i+1,j+1}   & + &  s(a_i,b_j)   \\\\" +
        "D'_{i+1,j}     & + &  \\gamma      \\\\" +
        "D'_{i,j+1}     & + &  \\gamma" +
        "\\end{cases}",

        HIRSCHBERG_FORWARD:
        "\\begin{cases}" +
        "D_{i-1,j-1}    & + & s(a_i,b_j)    \\\\" +
        "D_{i-1,j}      & + & \\gamma       \\\\" +
        "D_{i,j-1} 		& + & \\gamma" +
        "\\end{cases}",

        HIRSCHBERG_INITIALIZATION:
            "D'_{n,m} = 0, \\quad D'_{n,j}= (m-j)\\gamma, \\quad D'_{i,m}= (n-i)\\gamma",

        NEEDLEMAN_WUNSCH:
        "\\begin{cases}" +
        "D_{i-1,j-1}    & + & s(a_i,b_j)    \\\\" +
        "D_{i-1,j}      & + & s(a_i,-)      \\\\" +
        "D_{i,j-1} 		& + & s(-,b_j)" +
        "\\end{cases}",

        SMITH_WATERMAN:
        "\\begin{cases}" +
        "S_{i-1,j-1}    & + & s(a_i,b_j)    \\\\" +
        "S_{i-1,j}      & + & s(a_i,-)      \\\\" +
        "S_{i,j-1} 		& + & s(-,b_j)      \\\\" +
        "0" +
        "\\end{cases}",

        SMITH_WATERMAN_MODIFIED:
        "\\begin{cases}" +
        "S_{i-1,j-1}    & + & s^r(a_i,b_j)    \\\\" +
        "S_{i-1,j}      & + & s^r(a_i,-)      \\\\" +
        "S_{i,j-1} 		& + & s^r(-,b_j)      \\\\" +
        "0" +
        "\\end{cases}",

        WATERMAN_SMITH_BEYER_MIN:
        "\\begin{cases}" +
        "\\displaystyle     \\min_{1 \\leq k \\leq j} \\{   D_{i,j-k}   &   +   &   g(k)          \\}   \\\\" +
        "                                                   D_{i-1,j-1} &   +   &   s(a_i,b_j)          \\\\[5pt]" +
        "\\displaystyle     \\min_{1 \\leq k \\leq i} \\{   D_{i-k,j}   &   +   &   g(k)          \\}" +
        "\\end{cases}",

        WATERMAN_SMITH_BEYER_MAX:
        "\\begin{cases}" +
        "\\displaystyle     \\max_{1 \\leq k \\leq j} \\{   D_{i,j-k}   &   +   &   g(k)          \\}   \\\\" +
        "                                                   D_{i-1,j-1} &   +   &   s(a_i,b_j)          \\\\[5pt]" +
        "\\displaystyle     \\max_{1 \\leq k \\leq i} \\{   D_{i-k,j}   &   +   &   g(k)          \\}" +
        "\\end{cases}"
    },

    SUB_FORMULAS: {
        ARSLAN_EGECIOGLE_PEVZNER_SCORING:
        "s^r(a_i,b_j) = " +
        "\\begin{cases}" +
        "s(a_i,b_j) - 2 \\lambda^{r-1}    &  a_i, b_j \\in \\Sigma     \\\\" +
        "s(a_i,b_j) - \\lambda^{r-1}      &  a_i =\\_ \\vee b_j = \\_ " +
        "\\end{cases}",

        GOTOH_GAP_FUNCTION: "g(k) = \\alpha + \\beta \\cdot k"
    }
};

/**
 * Stores strings, numbers to access matrices used in the program.
 */
var MATRICES = {
    DEFAULT: "D",
    DEFAULT_NUMBER: 1,
    HORIZONTAL: "Q",
    HORIZONTAL_NUMBER: 2,
    ITERATION_NUMBER_1: -1,
    ITERATION_NUMBER_2: -2,
    ITERATION_NUMBER_3: -3,
    ITERATION_NUMBER_4: -4,
    ITERATION_NUMBER_5: -5,
    VERTICAL: "P",
    VERTICAL_NUMBER: 0
};

/**
 * Stores the type of moves you can do from a matrix to another matrix in affine algorithms like Gotoh.
 */
var MOVE = {
    HORIZONTAL: "horizontal",
    P_TO_X: "pToX",
    Q_TO_X: "qToX",
    VERTICAL: "vertical",
    X_TO_P: "xToP",
    X_TO_Q: "xToQ"
};

/**
 * Symbols which are used to be for example globally replaced.
 */
var MULTI_SYMBOLS = {
    BRACKET_LEFT: /\(/g,
    BRACKET_RIGHT: /\)/g,
    COMMA: /,/g,
    DELIMITER: /-/g,
    GAP: /_/g,
    G_LITTLE_SPECIAL: /ğ/g,
    MATCH: /\*/g,
    MISMATCH: /\|/g,
    NONE: /#/g,
    NUMBERS: /[0-9]/g,
    SEPARATORS: /;/g,
    SPACE: / /g,
    STRINGS: /[^0-9]/g
};

/**
 * Stores used class paths and some library paths.
 * Hint: Reimports are needed when different classes using objects of same class.
 * Else the objects of this classes are not correctly initialized.
 */
var PATHS = {
    ALIGNMENT_INTERFACE: "js/interfaces/alignment_interface.js",
    INPUT_PROCESSOR: "js/post_processing/input_processor.js",
    INTERFACE: "js/interfaces/interface.js",
    LINEAR_ALIGNMENT_INTERFACE: "js/interfaces/linear_alignment_interface.js", /* only if needed, else can be commented out */
    MULTI_SEQUENCE_INTERFACE: "js/interfaces/multi_sequence_interface.js", /* only if needed, else can be commented out */
    SUBADDITIVE_ALIGNMENT_INTERFACE: "js/interfaces/subadditive_alignment_interface.js", /* only if needed, else can be commented out */
    VISUALIZER: "js/post_processing/visualizer.js",

    LIBS: {
        KNOCKOUT: "js/libs/knockout-3.4.2.js"
    },

    MAIN: {
        PAGES: "",
        SCRIPTS: "js/"
    }
};

var PHYLOGENETIC_TREE = {
    SVG_CANVAS: '<div id=\"phylogenetic_tree\"></div>',
    SVG_CANVAS_NAME: "phylogenetic_tree",
    SVG_DIMENSION_FACTOR: 40,
    SVG_WIDTH: 450
};

/**
 * Stores HTML sub-tags.
 */
var SUB = {
    END_TAG: "</sub>",
    END_TAGS: /<\/sub>/g,
    START_TAG: "<sub>",
    START_TAGS: /<sub>/g
};

/**
 * Stores subadditive function names.
 */
var SUBADDITIVE_FUNCTIONS = {
    AFFINE: "affine",
    LOGARITHMIC: "logarithmic",
    QUADRATIC: "quadratic"
};

/**
 * Stores SVG-parameters of the SVG-overlay above the page.
 */
var SVG = {
    FLOW_LONG_ARROW_COLOR: "black",
    NAME_SPACE: "http://www.w3.org/2000/svg",
    PX: "px",
    STROKE_DASHARRAY: "3,6",
    TRACEBACK_LONG_ARROW_COLOR: "black",

    MARKER: {
        ID_TRACEBACK: "triangle_traceback",
        ID_FLOW: "triangle_flow",
        ORIENT: "auto",
        URL_TRACEBACK: "url(#triangle_traceback)",  // dependant on SVG.MARKER.ID_TRACEBACK
        URL_FLOW: "url(#triangle_flow)",  // dependant on SVG.MARKER.ID_FLOW
        VIEW_BOX: "0 0 8 8",

        BOUNDS: {
            HEIGHT: "4",
            REF_X: "0",  // relative marker coordinate
            REF_Y: "4",  // relative marker coordinate
            WIDTH: "4"
        }
    },

    TRIANGLE: {
        D: "M 0 0 L 0 8 L 8 4 z"  // M: move to, L: line to
    }
};

/**
 * Stores non-LaTeX symbols used in the program.
 */
var SYMBOLS = {  // contains all non-LaTeX symbols used in the project
    ALIGN: "~",
    AND: "&",
    BRACKET_LEFT: "(",
    BRACKET_RIGHT: ")",
    COLON: ":",
    COMMA: ",",
    DUMMY: "/",  // have to be a non-letter
    EMPTY: "",
    EQUAL: "=",
    GAP: "_",  // Hint: Do not change without a good reason! Can break the interface code!
    G_LITTLE: "g",
    HYPHEN: "-",  // Hint: Do not change without a good reason! Can break the loader!
    INFINITY: "∞",
    MINUS: "-",  // Hint: Do not change without a good reason! Can break the interface code! Regular expressions have to be changed.
    NEGATIVE_INFINITY: "-∞",
    NEW_LINE: "\n",
    NONE: "#",  // Hint: Do not change without a good reason! Can break the interface code!
    PLUS: "+",
    SEMICOLON: ";",
    SEPARATOR: "_",
    SPACE: " ",
    STAR: "*",
    VERTICAL_BAR: "|"
};

/**
 * Stores the download-name and the encoding of the table you want to download.
 */
var TABLE = {
    DOWNLOAD_NAME: "table.csv",
    TEXT_FILE_ENCODING: "text/csv;charset=utf-8"
};

var TOGGLE_LINK_TEXT = {
    HIDE: "hide intermediate steps",
    SHOW: "show intermediate steps"
};

/*---- DEFAULTS ----*/
/**
 * Stores the default parameters for alignment algorithms.
 */
var ALIGNMENT_DEFAULTS = {
    CALCULATION: "similarity",
    CALCULATION_HIRSCHBERG: "distance",
    SEQUENCE_1: "AATCG",  // Hint: UPPERCASE letters!
    SEQUENCE_2: "AACG",  // Hint: UPPERCASE letters!

    FUNCTION: {
        GAP: -2,
        MATCH: 1,
        MISMATCH: -1
    }
};

/**
 * Stores the default parameters for clustering algorithms.
 */
var HIERARCHICAL_CLUSTERING_DEFAULTS = {
    APPROACHES: [
        CLUSTERING_ALGORITHMS.COMPLETE_LINKAGE,
        CLUSTERING_ALGORITHMS.NEIGHBOUR_JOINING,
        CLUSTERING_ALGORITHMS.SINGLE_LINKAGE,
        CLUSTERING_ALGORITHMS.UPGMA,
        CLUSTERING_ALGORITHMS.WPGMA
    ],
    CSV_TABLE: /* input from https://en.wikipedia.org/wiki/Neighbor_joining */
    " 0 ;  5 ;  9 ;  9 ;  8" + "\n" +
    "   ;  0 ; 10 ; 10 ;  9" + "\n" +
    "   ;    ;  0 ;  8 ;  7" + "\n" +
    "   ;    ;    ;  0 ;  3" + "\n" +
    "   ;    ;    ;    ;  0",
    FORMULAS: [
        LATEX.FORMULAS.COMPLETE_LINKAGE_DELTA_1,
        LATEX.FORMULAS.NEIGHBOUR_JOINING_DELTA_1,
        LATEX.FORMULAS.SINGLE_LINKAGE_DELTA_1,
        LATEX.FORMULAS.UPGMA_DELTA_1,
        LATEX.FORMULAS.WPGMA_DELTA_1
    ],

    SUB_FORMULAS: [
        LATEX.FORMULAS.COMPLETE_LINKAGE_DELTA_2,
        LATEX.FORMULAS.NEIGHBOUR_JOINING_DELTA_2,
        LATEX.FORMULAS.SINGLE_LINKAGE_DELTA_2,
        LATEX.FORMULAS.UPGMA_DELTA_2,
        LATEX.FORMULAS.WPGMA_DELTA_2],

    STANDARD_APPROACH: [CLUSTERING_ALGORITHMS.UPGMA]
};

var ITERATIVE_SEQUENCE_DEFAULTS = {
    /* sequences from lecture */
    APPROACHES: [
        ITERATIVE_REFINEMENT_STRATEGIES.MIN_DISTANCE_PAIR,
        ITERATIVE_REFINEMENT_STRATEGIES.ONE_VS_ALL,
        ITERATIVE_REFINEMENT_STRATEGIES.PAIRWISE_BEST_PAIR
    ],

    ORDERS: [
        ITERATIVE_REFINEMENT_ORDERS.INPUT_ORDER,
        ITERATIVE_REFINEMENT_ORDERS.LEFT_FIRST,
        ITERATIVE_REFINEMENT_ORDERS.RIGHT_FIRST
    ],

    SEQUENCES: [/* input from T-Coffee paper does not work */
        "ACGT",
        "AT",
        "GCT",
        "GC"],
    SEQUENCES_COPY: [
        "ACGT",
        "AT",
        "GCT",
        "GC"], /* some bugfix for Knockout problem */

    STANDARD_APPROACH: [ITERATIVE_REFINEMENT_STRATEGIES.MIN_DISTANCE_PAIR],
    STANDARD_ORDER: [ITERATIVE_REFINEMENT_ORDERS.LEFT_FIRST],

    FUNCTION: {
        BASE_COSTS: -1,
        ENLARGEMENT: -3,
        MATCH: 1,
        MISMATCH: 0
    }
};

var MULTI_SEQUENCE_DEFAULTS = {
    /* example from paper */
    CALCULATION: "similarity",
    GLOBAL_ALIGNMENTS_PER_SEQUENCE: 1,
    LOCAL_ALIGNMENTS_PER_SEQUENCE: 2,

    SEQUENCES: [/* input from T-Coffee paper */
        "GARFIELD-THE-LAST-FAT-CAT",
        "GARFIELD-THE-FAST-CAT",
        "GARFIELD-THE-VERY-FAST-CAT",
        "THE-FAT-CAT"],
    SEQUENCES_COPY: [
        "GARFIELD-THE-LAST-FAT-CAT",
        "GARFIELD-THE-FAST-CAT",
        "GARFIELD-THE-VERY-FAST-CAT",
        "THE-FAT-CAT"], /* some bugfix for Knockout problem */

    USE_LOCAL_LIBRARY: false,

    FUNCTION: {
        BASE_COSTS: -3,
        BASE_COSTS_LOCAL: -3,

        ENLARGEMENT: -2,
        ENLARGEMENT_LOCAL: -3,

        MATCH: 1,
        MATCH_LOCAL: 2,

        MISMATCH: -1,
        MISMATCH_LOCAL: -2
    }
};

/**
 * Stores the default parameters for normalized alignment algorithms.
 */
var NORMALIZED_ALIGNMENT_DEFAULTS = {
    CALCULATION: "similarity",
    LENGTH: 10,
    SEQUENCE_1: "CTTGACCATG",  // Hint: UPPERCASE letters!
    SEQUENCE_2: "GCATTAGCCGG",  // Hint: UPPERCASE letters!

    FUNCTION: {
        GAP: -2,
        MATCH: 3,
        MISMATCH: -1
    }
};

/**
 * MathJax and other plugins can disrupt the rendering of the overlay and in these cases the overlay has to be redrawn.
 */
var REACTION_TIME_REDRAWING = {
    FIRST: 1000,
    SECOND: 4000,  // same as above, but for slow computers and tablets
    LONG_ARROWS: {
        FIRST: 500,
        SECOND: 4000  // same as above, but for slow computers and tablets
    }
};

/**
 * Stores the default parameters for subadditive alignment algorithms.
 */
var SUBADDITIVE_ALIGNMENT_DEFAULTS = {
    CALCULATION: "similarity",
    GAP_FUNCTION: "affine",
    SEQUENCE_1: "CG",  // Hint: UPPERCASE letters!
    SEQUENCE_2: "CCGA",  // Hint: UPPERCASE letters!

    FUNCTION: {
        BASE_COSTS: -3,
        ENLARGEMENT: -1,
        MATCH: 1,
        MISMATCH: -1
    }
};

/*---- LISTS ----*/
var CLUSTER_NAMES =
    ["a", "b", "c", "d", "e", "f", "g", "h",
        "i", "j", "k", "l", "m", "n", "o", "p",
        "q", "r", "s", "t", "u", "v", "w", "x",
        "y", "z"];   // Hint: There have to be at least one element in the list!

var EMPTY_ALIGNMENT = [SYMBOLS.EMPTY, SYMBOLS.EMPTY, SYMBOLS.EMPTY];

var GLOBAL_ALGORITHMS = [ALGORITHMS.GOTOH, ALGORITHMS.HIRSCHBERG, ALGORITHMS.NEEDLEMAN_WUNSCH, ALGORITHMS.WATERMAN_SMITH_BEYER];
var LOCAL_ALGORITHMS = [ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER, ALGORITHMS.GOTOH_LOCAL, ALGORITHMS.SMITH_WATERMAN];
var MULTI_SEQUENCE_ALGORITHMS = [ALGORITHMS.FENG_DOOLITTLE, ALGORITHMS.ITERATIVE_REFINMENT, ALGORITHMS.NOTREDAME_HIGGINS_HERINGA];

/**
 * Algorithms which use exactly three tables for their computations.
 * One table for the horizontal gaps,
 * one for the vertical gaps and one for matches and mismatches.
 */
var MULTI_TABLE_ALGORITHMS = [ALGORITHMS.GOTOH, ALGORITHMS.GOTOH_LOCAL];

/**
 * Algorithms with more than one-step arrows in the traceback!
 * It is used to automatically activate SVG-redrawing support
 * when the browser window is resized
 * or if a table slider is used.
 */
var SVG_ARROW_ALGORITHMS = [ALGORITHMS.GOTOH, ALGORITHMS.GOTOH_LOCAL, ALGORITHMS.HIRSCHBERG, ALGORITHMS.WATERMAN_SMITH_BEYER];

/**
 * Algorithms for which an initial table highlighting should be activated.
 * It is used to automatically activate initial highlighting
 * of for example traceback paths.
 * Hint: HIRSCHBERG and AEP are not in the list,
 * because they work with multiple tables without any "between-table" dependencies.
 */
var TABLE_INITIAL_HIGHLIGHT_ALGORITHMS =
    [ALGORITHMS.GOTOH, ALGORITHMS.GOTOH_LOCAL,
        ALGORITHMS.NEEDLEMAN_WUNSCH, ALGORITHMS.SMITH_WATERMAN, ALGORITHMS.WATERMAN_SMITH_BEYER];

/**
 * Algorithms which displaying a phylogenetic tree.
 * It is used to activate initial tree drawing (redrawing have to be done manually, to save resources).
 */
var TREE_ALGORITHMS = [ALGORITHMS.AGGLOMERATIVE_CLUSTERING, ALGORITHMS.FENG_DOOLITTLE,
    ALGORITHMS.ITERATIVE_REFINMENT, ALGORITHMS.NOTREDAME_HIGGINS_HERINGA];
