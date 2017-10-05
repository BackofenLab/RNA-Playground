/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

// constants
var CELL_PERCENT = 0.2;  // can change position of a "between-table"-arrow in the Gotoh algorithm

var DOUBLE_INFINITIES = /∞∞/g;

var MATH_JAX_TAGS = /(<[^>]*>)|(&amp;#x221E;)|(" role="presentation" style="font-size: 101%; position: relative;")|(\\infty)|(>)/g;

var NEW_LINE = "\r\n";

var REUPDATE_TIMEOUT_MS = 100;  // time in ms after which new LaTeX-Code is reinterpreted or outputs updated
var REACTION_TIME_HIGHLIGHT = REUPDATE_TIMEOUT_MS + 50;  // to highlight tracebacks only after outputs have been updated

var UNIT_TEST_WEBTITLE = "Console Runner";  // title of the Unit-test site

// structs
var AFFINE_ALIGNMENT_DEFAULTS = {
    CALCULATION: "similarity",
    SEQUENCE_1: "CCGA",
    SEQUENCE_2: "CG",

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
    SEQUENCE_1: "AACG",
    SEQUENCE_2: "AATCG",

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

var MATRICES = {
    DEFAULT: "X",
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

var SUB = {
    END_TAG: "</sub>",
    END_TAGS: /<\/sub>/g,
    START_TAG: "<sub>",
    START_TAGS: /<sub>/g
};

var SYMBOLS = {  // contains all strings used in the project
    DUMMY: "X",
    EMPTY: "",
    GAP: "_",
    INFINITY: "∞",
    NEGATIVE_INFINITY: "-$\\infty$",  // LaTeX
    POSITIVE_INFINITY: "$\\infty$",  // LaTeX
    SEPARATOR: "_",
    SPACE: " ",
    STAR: "*",
    VERTICAL_BAR: "|"
};

var SVG = {
    NAME_SPACE: "http://www.w3.org/2000/svg",
    OBJECT_COLOR: "black",

    MARKER : {
        ID: "triangle",
        ORIENT: "auto",
        URL: "url(#triangle)",  // dependant on SVG.MARKER.ID
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