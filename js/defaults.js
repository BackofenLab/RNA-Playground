/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

// constants
var CELL_PERCENT = 0.2;  // can change position of a "between-table"-arrow in the Gotoh algorithm

var MAX_NUMBER_TRACEBACKS = 10;

var REUPDATE_TIMEOUT_MS = 100;  // time in ms after which new LaTeX-Code is reinterpreted or outputs updated
var REACTION_TIME_HIGHLIGHT = REUPDATE_TIMEOUT_MS + 50;  // to highlight tracebacks only after outputs have been updated

var UNIT_TEST_WEBTITLE = "Console Runner";  // title of the Unit-test site

// structs
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

var ALGORITHMS = {  // contains a list of all implemented algorithms
    GOTOH: "gotoh",
    NEEDLEMAN_WUNSCH: "needleman_wunsch",
    SMITH_WATERMAN: "smith_waterman"
};

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

var ALIGNMENT_TYPES = {  // contains a list of alignment types
    DISTANCE: "distance",
    SIMILARITY: "similarity"
};

var ARROWS = {  // HINT: inner quotes have to be this here: " " !
    LEFT: '<div class="arrows_l"></div>',
    LEFT_NAME: ".arrows_l",
    TOP: '<div class="arrows_t"></div>',
    TOP_NAME: ".arrows_t",
    DIAGONAL: '<div class="arrows_d"></div>',
    DIAGONAL_NAME: ".arrows_d"
};

var CHARACTER = {
    BASE: /[a-zA-Z]/i,
    BASES: /^[a-zA-Z]+$/,
    NON_BASES: /[^a-zA-Z]+/g  // g to replace globally
};

var FILE_EXTENSIONS = {
    HYPERTEXT_MARKUP_LANGUAGE: ".html",
    JAVASCRIPT: ".js"
};

var INPUT = {
    ABS_MAX: 9,  // abs: absolute value
    ABS_MIN: -9
};

var KEY_CODES = {
    ENTER: 13
};

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

var MATRICES = {
    DEFAULT: "D",
    HORIZONTAL: "Q",
    VERTICAL: "P"
};

var MOVE = {
    P_TO_X: "pToX",
    Q_TO_X: "qToX",
    X_TO_P: "xToP",
    X_TO_Q: "xToQ"
};

var MULTI_SYMBOLS = {
    DELIMITER: /-/g,
    SPACE: / /g
};

var PATHS = {
    AFFINE_ALIGNMENT_INTERFACE: "js/interfaces/affine_alignment_interface.js",
    ALIGNMENT: "js/procedures/bases/alignment.js",
    ALIGNMENT_INTERFACE: "js/interfaces/alignment_interface.js",
    BACKTRACKING: "js/procedures/backtracking.js",
    INPUT_PROCESSOR: "js/post_processing/input_processor.js",
    VISUALIZER: "js/post_processing/visualizer.js",

    HTML: {
        AFFINE_ALIGNMENT_INTERFACE: "interfaces/affine_alignment_interface.html",
        ALIGNMENT_INTERFACE: "interfaces/alignment_interface.html",
    },

    LIBS: {
        KNOCKOUT: "js/libs/knockout-3.4.2.js"
    },

    MAIN: {
        PAGES: "",
        SCRIPTS: "js/"
    }
};

var SUB = {
    END_TAG: "</sub>",
    END_TAGS: /<\/sub>/g,
    START_TAG: "<sub>",
    START_TAGS: /<sub>/g
};

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

var TABLE = {
    DOWNLOAD_NAME: "table.csv",
    TEXT_FILE_ENCODING: "text/csv;charset=utf-8"
};