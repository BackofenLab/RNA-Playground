/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

/**
 * Used to store all constants used in the program.
 * Namespaces, HTML class-names and events like "mouse-over" etc. aren't stored in this file
 * because it would slow down the development process.
 * Multiple constants of same "logic" (for example with similar naming) stored in structures.
 */

// constants
var CELL_PERCENT = 0.2;  // can change position of a "between-table"-arrow in the Gotoh algorithm

var MAX_NUMBER_TRACEBACKS = 10;  // stores the number of tracebacks after which an alignment algorithm stops to compute

var REUPDATE_TIMEOUT_MS = 100;  // time in ms after which new LaTeX-Code is reinterpreted or outputs updated
var REACTION_TIME_HIGHLIGHT = REUPDATE_TIMEOUT_MS + 50;  // to highlight tracebacks only after outputs have been updated

var UNIT_TEST_WEBTITLE = "Console Runner";  // title of the Unit-test site

// structs
/**
 * Stores the default parameters for affine alignment algorithms.
 */
var AFFINE_ALIGNMENT_DEFAULTS = {
    CALCULATION: "similarity",
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
 * Stores the implemented algorithms.
 */
var ALGORITHMS = {  // contains a list of all implemented algorithms
    GOTOH: "gotoh",
    NEEDLEMAN_WUNSCH: "needleman_wunsch",
    SMITH_WATERMAN: "smith_waterman"
};

/**
 * Stores the default parameters for alignment algorithms.
 */
var ALIGNMENT_DEFAULTS = {
    CALCULATION: "similarity",
    SEQUENCE_1: "AACG",  // hint: UPPERCASE letters!
    SEQUENCE_2: "AATCG",  // hint: UPPERCASE letters!

    FUNCTION: {
        DELETION: -2,
        INSERTION: -2,
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

/**
 * Allowed input characters.
 */
var CHARACTER = {
    BASE: /[a-zA-Z]/i,
    BASES: /^[a-zA-Z]+$/,
    NON_BASES: /[^a-zA-Z]+/g  // g to replace globally
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
    ABS_MAX: 9,  // abs: absolute value
    ABS_MIN: -9
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
    BEGIN_CASES: "\\begin{cases}",
    END_CASES: "\\end{cases}",
    NEGATIVE_INFINITY: "-$\\infty$",    // LaTeX
    POSITIVE_INFINITY: "$\\infty$",     // LaTeX
    MATH_REGION: "$",
    MAX: "\\max",
    MIN: "\\min",
    MULTIPLIER: "k \\cdot",
    NEW_LINE: "\\\\",
    SPACE: "\\phantom{-}",

    FORMULA: {
        CURRENT: "D_{i,j}",
        CURRENT_P: "P_{i,j}",
        CURRENT_Q: "Q_{i,j}",
        DELETION: "b_j = -",
        DIAGONAL: "D_{i-1,j-1}",
        INSERTION: "a_i = -",
        LEFT: "D_{i,j-1}",
        LEFT_Q: "Q_{i,j-1}",
        MATCH: "a_i = b_j",
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

        GOTOH_P:
        "\\begin{cases}"                            +
        "D_{i-1,j}      & + & g(1)		    \\\\"   +
        "P_{i-1,j}      & + & \\delta"              +
  	    "\\end{cases}",

        GOTOH_Q:
        "\\begin{cases}"                            +
        "D_{i,j-1}      & + & g(1)		    \\\\"   +
        "Q_{i,j-1}      & + & \\delta"              +
        "\\end{cases}",

        GOTOH_GAP_FUNCTION:
        "g(k) = k \\cdot \\delta + \\alpha",

        NEEDLEMAN_WUNSCH:
        "\\begin{cases}"                            +
        "D_{i-1,j-1}    & + & s(a_i,b_j)    \\\\"   +
        "D_{i-1,j}      & + & s(a_i,-)      \\\\"   +
        "D_{i,j-1} 		& + & s(-,b_j)"             +
        "\\end{cases}",

        SMITH_WATERMAN:
        "\\begin{cases}"                            +
        "D_{i-1,j-1}    & + & s(a_i,b_j)    \\\\"   +
        "D_{i-1,j}      & + & s(a_i,-)      \\\\"   +
        "D_{i,j-1} 		& + & s(-,b_j)      \\\\"   +
        "0"                                         +
        "\\end{cases}"
    }
};

/**
 * Stores strings for matrix types used in the program.
 */
var MATRICES = {
    DEFAULT: "D",
    HORIZONTAL: "Q",
    VERTICAL: "P"
};

/**
 * Stores the type of moves you can do from a matrix to another matrix in affine algorithms like Gotoh.
 */
var MOVE = {
    P_TO_X: "pToX",
    Q_TO_X: "qToX",
    X_TO_P: "xToP",
    X_TO_Q: "xToQ"
};

/**
 * Symbols which are used to be for example globally replaced.
 */
var MULTI_SYMBOLS = {
    DELIMITER: /-/g,
    SPACE: / /g
};

/**
 * Stores all class paths and some library paths.
 */
var PATHS = {
    AFFINE_ALIGNMENT_INTERFACE: "js/interfaces/affine_alignment_interface.js",
    ALIGNMENT: "js/procedures/bases/alignment.js",
    ALIGNMENT_INTERFACE: "js/interfaces/alignment_interface.js",
    BACKTRACKING: "js/procedures/backtracking.js",
    INPUT_PROCESSOR: "js/post_processing/input_processor.js",
    VISUALIZER: "js/post_processing/visualizer.js",

    LIBS: {
        KNOCKOUT: "js/libs/knockout-3.4.2.js"
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
 * Stores non-LaTeX symbols used in the program.
 */
var SYMBOLS = {  // contains all strings used in the project
    AND: "&",
    BRACKET_LEFT: "(",
    BRACKET_RIGHT: ")",
    COMMA: ",",
    DUMMY: "/",  // have to be a non-letter
    EMPTY: "",
    EQUAL: "=",
    GAP: "_",
    INFINITY: "∞",
    NEGATIVE_INFINITY: "-∞",
    NEW_LINE: "\r\n",
    PLUS: "+",
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
            REF_X: "0",
            REF_Y: "4",
            WIDTH: "4"
        },
    },

    TRIANGLE : {
        D: "M 0 0 L 0 8 L 8 4 z"
    }
};

/**
 * Stores the download-name and the encoding of the table you want to download.
 */
var TABLE = {
    DOWNLOAD_NAME: "table.csv",
    TEXT_FILE_ENCODING: "text/csv;charset=utf-8"
};