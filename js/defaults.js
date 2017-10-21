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

/* HINT: alphabetically ordered (also within definitions) and structures after constants */

// constants
var EPSILON = 0.000000001;  // some very low number to test against
var MAX_NUMBER_TRACEBACKS = 10;  // stores the number of tracebacks after which an alignment algorithm stops to compute

var MAX_NUMBER_ITERATIONS = 5;  // number of iterations in algorithm with convergence
var REUPDATE_TIMEOUT_MS = 100;  // time in ms after which new LaTeX-Code is reinterpreted or outputs updated
var REACTION_TIME_HIGHLIGHT = REUPDATE_TIMEOUT_MS + 50;  // to highlight tracebacks only after outputs have been updated

var SMITH_WATERMAN_STOP = "0";

// structs
/**
 * Stores the implemented algorithm names.
 */
var ALGORITHMS = {  // contains a list of all implemented algorithms (javascript names without extension)
	ARSLAN_EGECIOGLU_PEVZNER: "arslan_egecioglu_pevzner",
    GOTOH: "gotoh",
    GOTOH_LOCAL: "gotoh_local",
    NEEDLEMAN_WUNSCH: "needleman_wunsch",
    NONE: "none",
    SMITH_WATERMAN: "smith_waterman",
    WATERMAN_SMITH_BEYER: "waterman_smith_beyer"
};

/**
 * Stores the default parameters for alignment algorithms.
 */
var ALIGNMENT_DEFAULTS = {
    CALCULATION: "similarity",
    SEQUENCE_1: "AACG",  // hint: UPPERCASE letters!
    SEQUENCE_2: "AATCG",  // hint: UPPERCASE letters!

    FUNCTION: {
        GAP: -2,
        MATCH: 1,
        MISMATCH: -1
    }
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
var ARROWS = {  // HINT: inner quotes have to be this here: " " or it won't work!
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
 * Defines allowed and disallowed input characters.
 */
var CHARACTER = {
    BASE: /[a-zA-Z]/i,
    BASES: /^[a-zA-Z-]+$/,
    NON_BASES: /[^a-zA-Z-]+/g,  // g to replace globally
    NUMBERS: /[-+]?[0-9]+\.[0-9]*/
};

/**
 * Allowed file extensions.
 */
var FILE_EXTENSIONS = {
    HYPERTEXT_MARKUP_LANGUAGE: ".html",
    JAVASCRIPT: ".js"
};

/**
 * Allowed max values for inputs.
 */
var INPUT = {
    LENGTH_MIN: 0,
    LENGTH_MAX: 1000,
    MAX: 10,  // abs: absolute value
    MIN: -10
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

    FORMULA: {
        CURRENT: "D_{i,j}",
        CURRENT_P: "P_{i,j}",
        CURRENT_Q: "Q_{i,j}",
        DELETION: "b_j = -",
        DIAGONAL: "D_{i-1,j-1}",
        INSERTION: "a_i = -",
        GAP: "g(k)",
        LEFT: "D_{i,j-1}",
        LEFT_Q: "Q_{i,j-1}",
        MATCH: "a_i = b_j",
        MAXIMIZE_HORIZONTAL: "\\displaystyle \\max_{1 \\leq k \\leq j} \\{D_{i,j-k} & + & g(k) \\}",
        MAXIMIZE_VERTICAL: "\\displaystyle \\max_{1 \\leq k \\leq j} \\{D_{i-k,j} & + & g(k) \\}",
        MINIMIZE_HORIZONTAL: "\\displaystyle \\min_{1 \\leq k \\leq j} \\{D_{i,j-k} & + & g(k) \\}",
        MINIMIZE_VERTICAL: "\\displaystyle \\min_{1 \\leq k \\leq j} \\{D_{i-k,j} & + & g(k) \\}",
        MISMATCH: "a_i \\neq b_j",
        TOP: "D_{i-1,j}",
        TOP_P: "P_{i-1,j}",
        ZERO: "0"
    },

    RECURSION: {
        GOTOH:
        "\\begin{cases}"                            +
        "D_{i-1,j-1}	& + & s(a_i,b_j)    \\\\"   +
        "P_{i,j}                            \\\\"   +
        "Q_{i,j}                                "   +
        "\\end{cases}",

        GOTOH_LOCAL:
        "\\begin{cases}"                            +
        "S_{i-1,j-1}	& + & s(a_i,b_j)    \\\\"   +
        "P_{i,j}                            \\\\"   +
        "Q_{i,j}                            \\\\"   +
        "0                                      "   +
        "\\end{cases}",

        GOTOH_P:
        "\\begin{cases}"                            +
        "D_{i-1,j}      & + & g(1)		    \\\\"   +
        "P_{i-1,j}      & + & \\beta"              +
        "\\end{cases}",

        GOTOH_Q:
        "\\begin{cases}"                            +
        "D_{i,j-1}      & + & g(1)		    \\\\"   +
        "Q_{i,j-1}      & + & \\beta"              +
        "\\end{cases}",

        NEEDLEMAN_WUNSCH:
        "\\begin{cases}"                            +
        "D_{i-1,j-1}    & + & s(a_i,b_j)    \\\\"   +
        "D_{i-1,j}      & + & s(a_i,-)      \\\\"   +
        "D_{i,j-1} 		& + & s(-,b_j)"             +
        "\\end{cases}",

        SMITH_WATERMAN:
        "\\begin{cases}"                            +
        "S_{i-1,j-1}    & + & s(a_i,b_j)    \\\\"   +
        "S_{i-1,j}      & + & s(a_i,-)      \\\\"   +
        "S_{i,j-1} 		& + & s(-,b_j)      \\\\"   +
        "0"                                         +
        "\\end{cases}",

      SMITH_WATERMAN_MODIFIED:
        "\\begin{cases}"                            +
        "S_{i-1,j-1}    & + & s^r(a_i,b_j)    \\\\"	+
        "S_{i-1,j}      & + & s^r(a_i,-)      \\\\"	+
        "S_{i,j-1} 		& + & s^r(-,b_j)      \\\\"	+
        "0"                                         +
        "\\end{cases}",
        WATERMAN_SMITH_BEYER_MIN:
        "\\begin{cases}"                                                                                                +
        "\\displaystyle     \\min_{1 \\leq k \\leq j} \\{   D_{i,j-k}   &   +   &   g(k)          \\}   \\\\"           +
        "                                                   D_{i-1,j-1} &   +   &   s(a_i,b_j)          \\\\[5pt]"      +
        "\\displaystyle     \\min_{1 \\leq k \\leq i} \\{   D_{i-k,j}   &   +   &   g(k)          \\}"                  +
        "\\end{cases}",

        WATERMAN_SMITH_BEYER_MAX:
        "\\begin{cases}"                                                                                                +
        "\\displaystyle     \\max_{1 \\leq k \\leq j} \\{   D_{i,j-k}   &   +   &   g(k)          \\}   \\\\"           +
        "                                                   D_{i-1,j-1} &   +   &   s(a_i,b_j)          \\\\[5pt]"      +
        "\\displaystyle     \\max_{1 \\leq k \\leq i} \\{   D_{i-k,j}   &   +   &   g(k)          \\}"                  +
        "\\end{cases}"
    },

    SUB_FORMULAS: {
        ARSLAN_EGECIOGLE_PEVZNER_SCORING:
        "s^r(a_i,b_j) = "                                                    +
        "\\begin{cases}"                                                    +
        "s(a_i,b_j) - 2 \\lambda^{r-1}    &  a_i, b_j \\in \\Sigma     \\\\"      +
        "s(a_i,b_j) - \\lambda^{r-1}      &  a_i =\\_ \\vee b_j = \\_ "           +
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
    D_BIG: /D/g,
    DELIMITER: /-/g,
    G_LITTLE_SPECIAL: /ğ/g,
    SPACE: / /g
};

var NORMALIZED_ALIGNMENT_DEFAULTS = {
    CALCULATION: "similarity",
    LENGTH: 10,
    SEQUENCE_1: "GCATTUGCCUU",  // hint: UPPERCASE letters!
    SEQUENCE_2: "CTTGACCATU",  // hint: UPPERCASE letters!

    FUNCTION: {
        GAP: -2,
        MATCH: 3,
        MISMATCH: -1
    }
};

/**
 * Stores used class paths and some library paths.
 * Hint: Reimports are needed when different classes using objects of same class.
 * Else the objects of this classes are not correctly initialized.
 */
var PATHS = {
    ALIGNMENT: "js/bases/alignment.js",
    ALIGNMENT_INTERFACE: "js/interfaces/alignment_interface.js",
    GOTOH: "js/gotoh.js",
    INPUT_PROCESSOR: "js/post_processing/input_processor.js",
    LINEAR_ALIGNMENT_INTERFACE: "js/interfaces/linear_alignment_interface.js",
    SMITH_WATERMAN: "js/smith_waterman.js",
    SUBADDITIVE_ALIGNMENT_INTERFACE: "js/interfaces/subadditive_alignment_interface.js",
    VISUALIZER: "js/post_processing/visualizer.js",

    LIBS: {
        KNOCKOUT: "js/libs/knockout-3.4.2.js",
        MATH_JAX: "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-MML-AM_CHTML"
    },

    MAIN: {
        PAGES: "",
        SCRIPTS: "js/"
    }
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
 * Stores the default parameters for subadditive alignment algorithms.
 */
var SUBADDITIVE_ALIGNMENT_DEFAULTS = {
    CALCULATION: "similarity",
    GAP_FUNCTION: "affine",
    SEQUENCE_1: "CCGA",  // hint: UPPERCASE letters!
    SEQUENCE_2: "CG",  // hint: UPPERCASE letters!

    FUNCTION: {
        BASE_COSTS: -3,
        ENLARGEMENT: -1,
        MATCH: 1,
        MISMATCH: -1
    }
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
 * Stores non-LaTeX symbols used in the program.
 */
var SYMBOLS = {  // contains all non-LaTeX symbols used in the project
    AND: "&",
    BRACKET_LEFT: "(",
    BRACKET_RIGHT: ")",
    COMMA: ",",
    DUMMY: "/",  // have to be a non-letter
    EMPTY: "",
    EQUAL: "=",
    GAP: "_",
    G_LITTLE: "g",
    INFINITY: "∞",
    NEGATIVE_INFINITY: "-∞",
    NEW_LINE: "\n",
    PLUS: "+",
    S_BIG: "S",
    SEPARATOR: "_",
    SPACE: " ",
    STAR: "*",
    VERTICAL_BAR: "|"
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

    MARKER : {
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

    TRIANGLE : {
        D: "M 0 0 L 0 8 L 8 4 z"  // M: move to, L: line to
    }
};

/**
 * Stores the download-name and the encoding of the table you want to download.
 */
var TABLE = {
    DOWNLOAD_NAME: "table.csv",
    TEXT_FILE_ENCODING: "text/csv;charset=utf-8"
};

// lists (used for check ups)
var GLOBAL_ALGORITHMS = [ALGORITHMS.GOTOH, ALGORITHMS.NEEDLEMAN_WUNSCH, ALGORITHMS.WATERMAN_SMITH_BEYER];
var LOCAL_ALGORITHMS = [ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER, ALGORITHMS.GOTOH_LOCAL, ALGORITHMS.SMITH_WATERMAN];

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
var SVG_ARROW_ALGORITHMS = [ALGORITHMS.GOTOH, ALGORITHMS.GOTOH_LOCAL, ALGORITHMS.WATERMAN_SMITH_BEYER];
