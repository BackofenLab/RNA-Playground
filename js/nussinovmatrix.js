/**
 * @file Main file containing backend algorithms for RNA-algorithms-JS project.
 Main Items contains in theis file are:
 -Matrix class
 -Nussinov algorithms
 -Traceback algorithms
 * @authors "Martin Mann", "Syed Mohsin Ali".
 */
"use strict";

/**
 * Utility class that covers RNA specific functions.
 *
 */
var RnaUtil = {

    /**
     * checks base pair complementary of two nucleotides
     * @param {string} nt1 first nucleotide
     * @param {string} nt2 second nucleotide
     * @returns {boolean} true if complementary; false otherwise
     */
    areComplementary: function (nt1, nt2) {
        //console.log("areComp:", nt1, nt2);
        var complementary =
            (nt1 === "A" && nt2 === "U") || (nt1 === "U" && nt2 === "A") ||
            (nt1 === "G" && nt2 === "C") || (nt1 === "C" && nt2 === "G") ||
            (nt1 === "G" && nt2 === "U") || (nt1 === "U" && nt2 === "G");
        return complementary;
    },

    /**
     * checks whether or not the sequence is a valid RNA sequence
     * @param sequence
     */
    isRnaSequence: function (sequence) {
        var isValid =
            // check if sequence given
            (sequence !== null)
            // check RNA alphabet
            && sequence.match("^[ACGU]+$");
        return isValid;
    }


};


/**
 * Ancestor information for a certain traceback
 */
var NussinovCellTrace = {

// list of parent cells
    parents: null,

// list of base pairs added
    bps: null,

    /**
     * initializes the object
     * @param {object} parents the parents to set
     * @param {object} bps the base pairs to set
     * @returns {NussinovCellTrace} this for call chaining
     */
    init: function (parents, bps) {
        this.parents = null;
        this.bps = null;
        // input check
        if ((parents !== null && bps === null) || (parents === null && bps !== null)) {
            console.log("ERROR : NussinovCellTrace.init : only one value null");
            return this;
        }
        // store sane data
        this.parents = parents;
        this.bps = bps;
        return this;
    }


};

/**
 * data stored within a cell of a nussinov matrix
 */
var NussinovCell = {

// row
    i: -1,

// column
    j: -1,

// value
    value: null,

// traces for the current value
    traces: null,

    /**
     * inits a cell with the given data and sets traces to an empty list
     * @param i the row of the cell
     * @param j the column of the cell
     * @param value the value of the cell
     * @return this : cell access for chaining
     */
    init: function (i, j, value) {
        // init data
        this.i = i;
        this.j = j;
        this.value = value;
        this.traces = [];
        // this access for chaining
        return this;
    }

};

/**
 * Encodes a full/partial traceback for a given cell
 */
var Traceback = {

    /** structure in dot-bracket-notation */
    structure: "",
    /** list of cell traces of the form [source, parent1, parent2, ...] */
    traces: [],
    //partialStructure: [],

};

/**
 * Nussinov matrix object
 */
var NussinovMatrix = {

        /**
         * Access to the sequence for this matrix
         */
        sequence: null,

        /**
         * Access name of recursion used
         */
        name: null,

        /**
         * Minimal loop length within computation
         */
        minLoopLength: 0,

        /**
         * cells of the matrix
         */
        cells: [],

        /**
         * The latex representation of the formula computing the matrix.
         */
        latex_representation: "$$",

        /**
         * The table description
         */
        descripion: "default description",

        /**
         * Tracebacks allowance
         */
        allowTraceBack: true,


        /**
         * initialize a matrix of dim = n+1 with indices 0..n, where n is the
         * length of the provided sequence
         * @param {string} sequence the RNA sequence (not null or empty)
         * @returns {NussinovMatrix} this
         */
        init: function
            (sequence, name) { //initialize matrix

            // reset data
            this.sequence = null;
            this.name = null;
            this.cells = [];

            // check input
            if (sequence == null || sequence === "" || name == null) {
                console.log("Matrix init failed for sequence (", sequence, ")");
                return this;
            }

            // store sequence
            this.sequence = sequence;
            this.name = name;

            // create matrix cells
            var n = this.sequence.length;
            for (var i = 0; i <= n; i++) {
                this.cells[i] = [];
                for (var j = 0; j <= n; j++) {
                    // create new cell and initialize
                    console.log("initializing matrix", this.name, this.sequence);
                    if (this.name === "unique" || this.name === "ambiguous" || this.name === "Ambiguous2" || this.name === "NussinovFold") {

                        this.cells[i][j] = Object.create(NussinovCell).init(i, j, 0);
                    }
                    else {
                        this.cells[i][j] = Object.create(NussinovCell).init(i, j, null);
                    }
                }
                ;
            }
            ;

            return this;
        }
        ,


        computeCell: function(i, j) {
           // updateCell(i, j);
        },
        /**
         * access whole object at cell location (i,j) in matrix
         * @param {int} i row #.
         * @param {int} j column #.
         * @returns {NussinovCell} cell object or null if not available
         */
        getCell: function (i, j) {
            // check border cases
            if (i < 0 || j < 0 || i >= this.getDim() || j >= this.getDim() || i > j + 1) {
                return null;
            }
            if (this.cells[i][j] === null) {
                //this.cells[i][j] = computeCell(i, j);
            }
            return this.cells[i][j];
        }
        ,

        /**
         * Gives the dimension of the quadratic matrix
         * @returns {int} the dimension (#rows and columns)
         */
        getDim: function () {
            if (this.cells === null) {
                return 0;
            }
            return this.cells.length;
        }
        ,

        /**
         * returns traces of cell (i,j), i.e. ancestor cells and the basepairs they form
         * @param {int} i row #.
         * @param {int} j column #.
         * @returns {object} ancestor's object Eg. {parents:[],bPs:[]} or null if not available
         */
        getTraces: function (i, j) { //get traceback info for each cell
            // access cell at location (i,j) in the matrix
            var cell = this.getCell(i, j);
            // check if valid cell
            if (cell === null) {
                return null;
            }
            return cell.traces;
        }
        ,

        computeValue: function (i, j) {
            // Computes the M[i, j] by accessing cells using getValue function for memoization.
            
            // base case
            return 0;
            
            // computation
        },
        /**
         * access value of cell at location (i,j) in matrix, and compute it if it's not computed.
         * @param {int} i row #.
         * @param {int} j column #.
         * @returns {int} cell value or 0 if invalid cell
         */
        getValue: function (i, j) {
            // access cell at location (i,j) in the matrix
            var cell = this.getCell(i, j);
            if (cell === null) {
                return null;
            }
            // check if invalid cell
            if (cell.value === null) {
                cell.value = this.computeValue(i, j);
            }
            // get cell value
            return cell.value;
        }
        ,

        /**
         * Updates the ancestor list of a given cell if the curVal is higher or
         * equal to the current value within the cell.
         * If the value is equal, curAncestor is added to the list.
         * If the value is smaller than curVal, curAncestor will be set to be the
         * only list entry.
         * @param i the row of the cell to update
         * @param j the column of the cell to update
         * @param curAncestor
         */
        updateCell: function (i, j, curAncestor) {
            // get cell to update
            var curCell = this.getCell(i, j);
            // check if something to update
            if (curCell === null) {
                return;
            }

            // init value with number of additional base pairs
            var curVal = curAncestor.bps.length;
            // add scores of ancestor cells
            for (var x = 0; x < curAncestor.parents.length; x++) {
                curVal += this.getValue(curAncestor.parents[x][0], curAncestor.parents[x][1]);
            }
            // check if we have to update
            if (curCell.value <= curVal) {
                // check for new maximal value
                if (curCell.value < curVal) {
                    // reset ancestor list
                    curCell.traces = [];
                    // store new maximum
                    curCell.value = curVal;
                }
                ;
                // store this ancestor
                curCell.traces.push(curAncestor);
            }
            ;

        }
        ,

        /**
         * Fills the matrix according to the recursion.
         *
         * NOTE: this function has to be overwritten by subclasses
         *
         * @param {string} sequence the RNA sequence to compute the matrix for
         * @param {int} minLoopLength the minimal loop length to be used for computation
         *
         * @returns {NussinovMatrix} this for call chaining
         */
        computeMatrix: function (input) {
            console.log("WARNING: computeMatrix() not implemented in NussinovMatrix superclass; overwrite in subclass!");

            // resize and init matrix
            this.init(sequence);

            // set minimal loop length
            this.minLoopLength = input.loopLength;

            return this;
        }
        ,

        /**
         * Access to the recursion in LaTeX encoding that is used in the computeMatrix implementation
         *
         * NOTE: this function has to be overwritten by subclasses
         *
         * @returns {string} latex encoding of the recursion
         */
        getRecursionInLatex: function () {
            console.log("WARNING: getRecursionInLatex() not implemented in NussinovMatrix superclass; overwrite in subclass!");
            return "";
        }
        ,

        /**
         * Returns a description for the implemented recursion
         *
         * NOTE: this function has to be overwritten by subclasses
         *
         * @returns {string} description of the recursion
         */
        getRecursionDescription: function () {
            console.log("WARNING: getRecursionDescription() not implemented in NussinovMatrix superclass; overwrite in subclass!");
            return "";
        }
        ,

        /**
         * creates a string representation of the matrix
         * @returns {string} matrix as string
         */
        toString: function () {
            var str = this.minLoopLength + "    ";
            for (var i = 0; i + 1 < this.getDim(); i++) {
                str += this.sequence[i] + "  ";
            }
            str += "\n";
            for (var i = 1; i < this.getDim(); i++) {
                // print sequence
                if (i === 0) {
                    str += "  ";
                } else {
                    str += this.sequence[i - 1] + " ";
                }
                // print values
                for (var j = 0; j < this.getDim(); j++) {
                    if (j !== 0) {
                        str += ", ";
                    }
                    // get value
                    str += this.getValue(i, j);
                }
                ;
                str += "\n";
            }
            ;
            return str;
        }
        ,


        /**
         * countBasepairs(for wuchty)
         */
        countBasepairs: function (bps, sigma) {
            var NSprime = bps.length;
            for (var s in sigma) {
                var i = sigma[s][0];
                var j = sigma[s][1];

                NSprime += this.getValue(i, j);
            }
            return NSprime;
        },




        conv_str: function (x, length) {
            var str = "";
            for (var l = 0; l < length; l++) {
                str += ".";
            }
            var linked = this.sequence.indexOf("XXX") > -1;

            //console.log(str, length);
            for (var i in x) {
                str = str.substr(0, x[i][0] - 1) + "(" + str.substr(x[i][0], str.length - 1);
                str = str.substr(0, x[i][1] - 1) + ")" + str.substr(x[i][1], str.length - 1);
            }
            return str;
        },

    }
    ;

/**
 * Dynamic programming algorithm.
 *
 * DP Algorithms will work by memoization usually. If so, there's a function getValue for each of the tables, that
 * that computes a value and memoize it if it's not computed(default null cells), and returns the memoized value.
 *
 * The computation is done through the computeValue in each of the tables, this function should be overriden for each
 * table depending on how an entry in the table is computed. The function should include the base case, and should use
 * the other tables by the getValue function and not by accessing the tables directly, in order to ensure the
 * correctness of the memoization.
 *
 * The tables should be usable after invoking the computeMatrix method, this method should be overriden to set the tables
 * parameters  and compute all the dynamic programming. Unless there's a special way need to do this, it can usually 
 * (with memoization) be done by invoking getValue for all the entries of all the tables.
 *
 *
 * @type {{Description: Algorithm description,
 *        Tables: Array of tables for all the recursive formulas,
 *        defaultPars: {}, Parameters that can be fed to the algorithm, when needed. // To Be Removed
 *        computeMatrix: Compute All the Tables
 *        getRecursionInLatex: Get Latex String describing the recursive equations for all the tables.}}
 */
var DPAlgorithm = {
    Description: "Algorithm",

    Tables: [], // create new array

    defaultPars: {},

    computeMatrix: function (args_dict) {
    },

    // Return an aligned latex array that contains the latex formula of each table, (seperated with empty lines).
    getRecursionInLatex: function () {
        var formula = " \\begin{array} ";
        for (var i = 0; i < this.Tables.length; ++i) {
            formula += " \\\\ \\\\ " + this.Tables[i].latex_representation;
        }
        formula += " \\end{array} ";
        //console.log(formula);
        return formula;
    },

};


var NussinovDPAlgorithm_Ambiguous = Object.create(DPAlgorithm);

NussinovDPAlgorithm_Ambiguous.Description = "Ambiguous recursion";
NussinovDPAlgorithm_Ambiguous.Tables = new Array();
NussinovDPAlgorithm_Ambiguous.Tables.push(Object.create(NussinovMatrix));
NussinovDPAlgorithm_Ambiguous.Tables[0].latex_representation = "D(i,j) = \\max \\begin{cases} D(i+1,j) & S_i \\text{ unpaired} \\\\ D(i,j-1) & S_j \\text{ unpaired} \\\\ D(i+1,j-1)+1 &  S_i,S_j \\text{ compl. base pair and } i+ l< j \\\\ \\max_{i< k< (j-1)} D(i,k)+D(k+1,j) & \\text{decomposition} \\end{cases}";


NussinovDPAlgorithm_Ambiguous.computeMatrix = function (input) {
// resize and initialize matrix
    this.Tables[0].init(input.sequence(), "ambiguous");
// store minimal loop length
    var minLL = parseInt(input.loopLength());
    this.Tables[0].minLoopLength = minLL;
    //console.log("computing ambiguos matrix");
// fill matrix by diagonals
// iterate over all substructure spans that can have a base pair
    for (var span = this.Tables[0].minLoopLength; span < this.Tables[0].getDim(); span++) {
        // iterate over all rows
        for (var i = 1; i < this.Tables[0].getDim() - this.Tables[0].minLoopLength; i++) {
            // get column for current span
            var j = i + span;

            // i unpaired
            this.Tables[0].updateCell(i, j, Object.create(NussinovCellTrace).init([[i + 1, j]], []));

            // j unpaired
            this.Tables[0].updateCell(i, j, Object.create(NussinovCellTrace).init([[i, j - 1]], []));

            // check (i,j) base pair
            if ((j - i > input.loopLength()) && RnaUtil.areComplementary(input.sequence()[i - 1], input.sequence()[j - 1])) {
                // get value for base pair
                this.Tables[0].updateCell(i, j, Object.create(NussinovCellTrace).init([[i + 1, j - 1]], [[i, j]]));
            }
            ;

            // check decomposition into substructures (minLength==2)
            for (var k = i + 1; k < (j - 1); k++) {
                // get decomposition value
                this.Tables[0].updateCell(i, j, Object.create(NussinovCellTrace).init([[i, k], [k + 1, j]], []));
            }
            ;
         }
        ;
    }
    ;

    return this.Tables;
};

NussinovDPAlgorithm_Ambiguous.Tables[0].getSubstructures = function (sigma, P, traces, delta, maxLengthR) {
    var Nmax = this.getCell(1, this.sequence.length).value;
    var R = [];
    var ij = sigma.pop();
    //console.log(ij);

    // check for sane interval
    // if i>j dont continue
    if (ij[0] >= ij[1]) {
        //console.log("ij[0] > ij[1]", ij[0], ij[1]);
        var S_prime = {};

        var sigma_prime = JSON.stringify(sigma);
        sigma_prime = JSON.parse(sigma_prime);

        var tmp_P = JSON.stringify(P);
        tmp_P = JSON.parse(tmp_P);

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);

        S_prime.sigma = sigma_prime;
        S_prime.P = tmp_P;
        S_prime.traces = tmp_traces;

        R.push(S_prime);
        //console.log("returning R:", JSON.stringify(R));
        return R;
    }

    // if (i,j) == (i+1,j-1) + bp(ij)
    {
        if (ij[1] - ij[0] > this.minLoopLength) {
            //console.log(this.sequence);
            //console.log(this.sequence[ij[0] - 1], this.sequence[ij[1] - 1]);
            if (RnaUtil.areComplementary(this.sequence[ij[0] - 1], this.sequence[ij[1] - 1])) {
                var sigma_prime = JSON.stringify(sigma);
                sigma_prime = JSON.parse(sigma_prime);
                sigma_prime.push([ij[0] + 1, ij[1] - 1]);

                var tmp_P = JSON.stringify(P);
                tmp_P = JSON.parse(tmp_P);
                tmp_P.push([ij[0], ij[1]]);

                var tmp_traces = JSON.stringify(traces);
                tmp_traces = JSON.parse(tmp_traces);

                var NSprime = this.countBasepairs(tmp_P, sigma_prime);

                if (NSprime >= Nmax - delta) {
                    var S_prime = {};
                    S_prime.sigma = sigma_prime;
                    S_prime.P = tmp_P;
                    tmp_traces.unshift([ij, [[ij[0] + 1, ij[1] - 1]]]);
                    S_prime.traces = tmp_traces;
                    //console.log("i+1,j-1:", JSON.stringify(S_prime));
                    // push to the front to keep base pair most prominent to refine
                    R.unshift(S_prime);
                }
            }
        }

        // check if enough structures found so far
        if (R.length >= maxLengthR) {
            //console.log("returning R:", JSON.stringify(R));
            return R;
        }
    }

    // if (i,j) == (i+1,j)
    {
        var sigma_prime = JSON.stringify(sigma);
        sigma_prime = JSON.parse(sigma_prime);
        sigma_prime.unshift([ij[0] + 1, ij[1]]);

        var tmp_P = JSON.stringify(P);
        tmp_P = JSON.parse(tmp_P);

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);

        var NSprime = this.countBasepairs(P, sigma_prime);

        if (NSprime >= Nmax - delta) {
            var S_prime = {};
            S_prime.sigma = sigma_prime;
            S_prime.P = tmp_P;
            tmp_traces.unshift([ij, [[ij[0] + 1, ij[1]]]]);
            S_prime.traces = tmp_traces;
            //console.log("i+1,j:", JSON.stringify(S_prime));
            // push to the front to keep base pair most prominent to refine
            R.unshift(S_prime);
        }

        // check if enough structures found so far
        if (R.length >= maxLengthR) {
            //console.log("returning R:", JSON.stringify(R));
            return R;
        }
    }

    // if (i,j) == (i,j-1)
    {
        var sigma_prime = JSON.stringify(sigma);
        sigma_prime = JSON.parse(sigma_prime);
        sigma_prime.unshift([ij[0], ij[1] - 1]);

        var tmp_P = JSON.stringify(P);
        tmp_P = JSON.parse(tmp_P);

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);

        var NSprime = this.countBasepairs(P, sigma_prime);

        if (NSprime >= Nmax - delta) {
            var S_prime = {};
            S_prime.sigma = sigma_prime;
            S_prime.P = tmp_P;
            tmp_traces.unshift([ij, [[ij[0], ij[1] - 1]]]);
            S_prime.traces = tmp_traces;
            //console.log("i,j-1:", JSON.stringify(S_prime));
            // push to the front to keep base pair most prominent to refine
            R.unshift(S_prime);
        }

        // check if enough structures found so far
        if (R.length >= maxLengthR) {
            //console.log("returning R:", JSON.stringify(R));
            return R;
        }
    }

    // if (i,j) == (i,l) + (l+1, j)
    for (var l = ij[0] + 1; l < ij[1] - 1; l++) {

        var sigma_prime = JSON.stringify(sigma);
        sigma_prime = JSON.parse(sigma_prime);
        sigma_prime.push([ij[0], l]);
        sigma_prime.push([l + 1, ij[1]]);

        var tmp_P = JSON.stringify(P);
        tmp_P = JSON.parse(tmp_P);

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);

        var NSprime = this.countBasepairs(tmp_P, sigma_prime);

        if (NSprime >= Nmax - delta) {

            var S_prime = {};
            S_prime.sigma = sigma_prime;
            S_prime.P = tmp_P;
            tmp_traces.unshift([ij, [[ij[0], l], [l + 1, ij[1]]]]);
            S_prime.traces = tmp_traces;
            //console.log("ilj:", JSON.stringify(S_prime));
            // push to the front to keep base pair most prominent to refine
            R.unshift(S_prime);
        }

        // check if enough structures found so far
        if (R.length >= maxLengthR) {
            //console.log("returning R:", JSON.stringify(R));
            return R;
        }

    }

    console.log("returning R:", JSON.stringify(R));
    return R;
}
;


var NussinovDPAlgorithm_Unique = Object.create(DPAlgorithm);

NussinovDPAlgorithm_Unique.Description = "Recursion by Nussinov et al. (1978) with unique decomposition";
NussinovDPAlgorithm_Unique.Tables = new Array();
NussinovDPAlgorithm_Unique.Tables.push(Object.create(NussinovMatrix));
NussinovDPAlgorithm_Unique.Tables[0].latex_representation = "D(i,j) = \\max \\begin{cases} D(i,j-1) & S_j \\text{ unpaired} \\\\ \\max_{i\\leq k< (j-l)} D(i,k-1)+D(k+1,j-1)+1 & S_k,S_j \\text{ compl. base pair} \\end{cases}";

NussinovDPAlgorithm_Unique.computeMatrix = function (input) {

    console.log('nussiUnique', input.sequence());
// resize and initialize matrix
    this.Tables[0].init(input.sequence(), "unique");
// store minimal loop length
    var minLL = parseInt(input.loopLength());
    this.Tables[0].minLoopLength = minLL;

// fill matrix by diagonals
// iterate over all substructure spans that can have a base pair
    for (var span = this.Tables[0].minLoopLength; span < this.Tables[0].getDim(); span++) {
        // iterate over all rows
        for (var i = 0; i < this.Tables[0].getDim() - this.Tables[0].minLoopLength ; i++) {
            // get column for current span
            var j = i + span;

            // j unpaired
            this.Tables[0].updateCell(i, j, Object.create(NussinovCellTrace).init([[i, j - 1]], []));

            // check base pair based decomposition : (k,j) base pair
            for (var k = i; k + this.Tables[0].minLoopLength < j; k++) {
                // check if sequence positions are compatible
                if (RnaUtil.areComplementary(input.sequence()[k - 1], input.sequence()[j - 1])) {
                    this.Tables[0].updateCell(i, j, Object.create(NussinovCellTrace).init([[i, k - 1], [k + 1, j - 1]], [[k, j]]));
                }
                ;
            }
            ;
        }
        ;
    }
    ;

    return this.Tables;
};


NussinovDPAlgorithm_Unique.Tables[0].getSubstructures = function (sigma, P, traces, delta, maxLengthR) {
    var Nmax = this.getCell(1, this.sequence.length).value;
    var R = [];
    var ij = sigma.pop();

    // if i>j dont countinue
    if (ij[0] > ij[1]) {
        //console.log("ij[0] > ij[1]", ij[0], ij[1]);
        var S_prime = {};
        S_prime.sigma = sigma;
        S_prime.P = P;
        S_prime.traces = traces;
        R.push(S_prime);
        //console.log("returning R:", JSON.stringify(R));
        return R;
    }

    // if (i,j) == (i,l-1) + (l+1, j-1) + 1
    for (var l = ij[0]; l <= ij[1] - 1; l++) {
        if (ij[1] - l > this.minLoopLength) {
            if (RnaUtil.areComplementary(this.sequence[l - 1], this.sequence[ij[1] - 1])) {
                var sigma_prime = JSON.stringify(sigma);
                sigma_prime = JSON.parse(sigma_prime);
                sigma_prime.push([ij[0], l - 1]);
                sigma_prime.push([l + 1, ij[1] - 1]);

                var tmp_P = JSON.stringify(P);
                tmp_P = JSON.parse(tmp_P);
                tmp_P.push([l, ij[1]]);

                var tmp_traces = JSON.stringify(traces);
                tmp_traces = JSON.parse(tmp_traces);

                var NSprime = this.countBasepairs(tmp_P, sigma_prime);

                if (NSprime >= Nmax - delta) {

                    var S_prime = {};
                    S_prime.sigma = sigma_prime;
                    S_prime.P = tmp_P;
                    tmp_traces.unshift([ij, [[ij[0], l - 1], [l + 1, ij[1] - 1]]]);
                    S_prime.traces = tmp_traces;
                    //console.log("ilj:", JSON.stringify(S_prime));
                    // push to the back to keep base pair most prominent to refine
                    R.push(S_prime);
                }

                // check if enough structures found so far
                if (R.length >= maxLengthR) {
                    //console.log("returning R:", JSON.stringify(R));
                    return R;
                }

            }
        }
    }

    // if (i,j) == (i,j-1)
    {
        var sigma_prime = JSON.stringify(sigma);
        sigma_prime = JSON.parse(sigma_prime);
        sigma_prime.unshift([ij[0], ij[1] - 1]);

        var tmp_P = JSON.stringify(P);
        tmp_P = JSON.parse(tmp_P);

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);

        var NSprime = this.countBasepairs(P, sigma_prime);

        if (NSprime >= Nmax - delta) {
            var S_prime = {};
            S_prime.sigma = sigma_prime;
            S_prime.P = tmp_P;
            tmp_traces.unshift([ij, [[ij[0], ij[1] - 1]]]);
            S_prime.traces = tmp_traces;
            //console.log("ij-1:", JSON.stringify(S_prime));
            // push to the front to keep base pair most prominent to refine
            R.unshift(S_prime);
        }

        // check if enough structures found so far
        if (R.length >= maxLengthR) {
            //console.log("returning R:", JSON.stringify(R));
            return R;
        }
    }

    //console.log("returning R:", JSON.stringify(R));
    return R;
}

/**
 * nussinov recursion
 * N(i,j) = max(0, N(i+1,j-1)+1 if bp(i,j), max_{i<=k<j} : N(i,k)+N(k+1,j))
 * @type {DPAlgorithm}
 */
var NussinovDPAlgorithm_Ambiguous2 = Object.create(DPAlgorithm);

NussinovDPAlgorithm_Ambiguous2.Description = "Recursion by Nussinov et al. (1978) with Ambiguous2 decomposition";
NussinovDPAlgorithm_Ambiguous2.Tables = new Array();
NussinovDPAlgorithm_Ambiguous2.Tables.push(Object.create(NussinovMatrix));
NussinovDPAlgorithm_Ambiguous2.Tables[0].latex_representation = "D(i,j) = \\max \\begin{cases} D(i+1,j-1)+1 & S_i,S_j \\text{ compl. base pair} \\\\ \\max_{i\\leq k< j} D(i,k)+D(k+1,j) \\end{cases}";
//NussinovDPAlgorithm_Ambiguous2.Tables[0].latex_representation = "$D(i,j) = \\max(0, D(i+1,j-1)+1 if bp(i,j), \\max_{i<=k<j} : D(i,k)+N(k+1,j))$";
NussinovDPAlgorithm_Ambiguous2.computeMatrix = function (input) {

// resize and initialize matrix
    this.Tables[0].init(input.sequence(), "Ambiguous2");
// store minimal loop length
    var minLL = parseInt(input.loopLength());
    this.Tables[0].minLoopLength = minLL;

// fill matrix by diagonals
// iterate over all substructure spans that can have a base pair
    for (var span = this.Tables[0].minLoopLength; span < this.Tables[0].getDim(); span++) {
        // iterate over all rows
        for (var i = 1; i < this.Tables[0].getDim() - this.Tables[0].minLoopLength; i++) {
            // get column for current span
            var j = i + span;

            // check (i,j) base pair
            if ((j - i > input.loopLength()) && RnaUtil.areComplementary(input.sequence()[i - 1], input.sequence()[j - 1])) {
                // get value for base pair
                this.Tables[0].updateCell(i, j, Object.create(NussinovCellTrace).init([[i + 1, j - 1]], [[i, j]]));
            }
            ;

            // check decomposition into substructures (minLength==2)
            for (var k = i; k < j; k++) {
                // get decomposition value
                this.Tables[0].updateCell(i, j, Object.create(NussinovCellTrace).init([[i, k], [k + 1, j]], []));
            }
            ;
        }
        ;
    }
    ;

    return this.Tables;
};


NussinovDPAlgorithm_Ambiguous2.Tables[0].getSubstructures = function (sigma, P, traces, delta, maxLengthR) {
    var Nmax = this.getCell(1, this.sequence.length).value;
    var R = [];
    var ij = sigma.pop();
    console.log(ij);

    // check for sane interval
    // if i>j dont continue
    if (ij[0] >= ij[1]) {
        //console.log("ij[0] > ij[1]", ij[0], ij[1]);
        var S_prime = {};

        var sigma_prime = JSON.stringify(sigma);
        sigma_prime = JSON.parse(sigma_prime);

        var tmp_P = JSON.stringify(P);
        tmp_P = JSON.parse(tmp_P);

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);

        S_prime.sigma = sigma_prime;
        S_prime.P = tmp_P;
        S_prime.traces = tmp_traces;

        R.push(S_prime);
        //console.log("returning R:", JSON.stringify(R));
        return R;
    }

    // if (i,j) == (i+1,j-1) + bp(ij)
    {
        if (ij[1] - ij[0] > this.minLoopLength) {
            //console.log(this.sequence);
            //console.log(this.sequence[ij[0] - 1], this.sequence[ij[1] - 1]);
            if (RnaUtil.areComplementary(this.sequence[ij[0] - 1], this.sequence[ij[1] - 1])) {
                var sigma_prime = JSON.stringify(sigma);
                sigma_prime = JSON.parse(sigma_prime);
                sigma_prime.push([ij[0] + 1, ij[1] - 1]);

                var tmp_P = JSON.stringify(P);
                tmp_P = JSON.parse(tmp_P);
                tmp_P.push([ij[0], ij[1]]);

                var tmp_traces = JSON.stringify(traces);
                tmp_traces = JSON.parse(tmp_traces);

                var NSprime = this.countBasepairs(tmp_P, sigma_prime);

                if (NSprime >= Nmax - delta) {
                    var S_prime = {};
                    S_prime.sigma = sigma_prime;
                    S_prime.P = tmp_P;
                    tmp_traces.unshift([ij, [[ij[0] + 1, ij[1] - 1]]]);
                    S_prime.traces = tmp_traces;
                    //console.log("i+1,j-1:", JSON.stringify(S_prime));
                    // push to the front to keep base pair most prominent to refine
                    R.unshift(S_prime);
                }
            }
        }

        // check if enough structures found so far
        if (R.length >= maxLengthR) {
            //console.log("returning R:", JSON.stringify(R));
            return R;
        }
    }

    // if (i,j) == (i,l) + (l+1, j)
    for (var l = ij[0]; l < ij[1]; l++) {
        console.log('here');
        var sigma_prime = JSON.stringify(sigma);
        sigma_prime = JSON.parse(sigma_prime);
        sigma_prime.push([ij[0], l]);
        sigma_prime.push([l + 1, ij[1]]);

        var tmp_P = JSON.stringify(P);
        tmp_P = JSON.parse(tmp_P);

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);

        var NSprime = this.countBasepairs(tmp_P, sigma_prime);

        if (NSprime >= Nmax - delta) {

            var S_prime = {};
            S_prime.sigma = sigma_prime;
            S_prime.P = tmp_P;
            tmp_traces.unshift([ij, [[ij[0], l], [l + 1, ij[1]]]]);
            S_prime.traces = tmp_traces;
            //console.log("ilj:", JSON.stringify(S_prime));
            // push to the front to keep base pair most prominent to refine
            R.unshift(S_prime);
        }

        // check if enough structures found so far
        if (R.length >= maxLengthR) {
            //console.log("returning R:", JSON.stringify(R));
            return R;
        }

    }

    console.log("returning R:", JSON.stringify(R));
    return R;
}
;

/**
 * WUCHTY(2nd version) enumerating up to 10 structures
 */
function wuchty_2nd(xmat, delta, formula) {
    console.log("entering wuchty");
    // get maximal number of structures to report
    var maxSOS = 10;
    // call subroutine
    return wuchty_2nd_limited(xmat, delta, formula, maxSOS);
}

/**
 * WUCHTY(2nd version)
 */
function wuchty_2nd_limited(xmat, delta, formula, maxSOS) {
    //if (xmat == undefined)return;
    if (xmat.sequence == undefined)return;
    var seq_length = xmat.sequence.length;
    var Nmax = xmat.getCell(1, seq_length).value;
    console.log("Wuchty beginning\nSequence:", xmat.sequence, "\nNmax:", Nmax, "\nDelta:", delta, "\n");

    var S = {sigma: [[1, seq_length]], P: [], traces: []};
    var R = [S];
    var SOS = [];
    var loop = 0;

    while (R.length != 0) {
        //console.log("\nloop:", ++loop);
        //console.log("R length:" + R.length);

        // Pop R
        var pop_R = R.pop();
        var sigma = pop_R.sigma;
        var P = pop_R.P;
        var t_traces = JSON.stringify(pop_R.traces);
        var traces = JSON.parse(t_traces);
        //console.log("poped R:", JSON.stringify(pop_R));
        //console.log("R_remaining:", JSON.stringify(R));

        var sigma_remaining = 0;
        for (var s in sigma) {//console.log("var s:", sigma[s]);
            if ((sigma[s][0]) <= (sigma[s][1] - xmat.minLoopLength)) sigma_remaining++;
        }

        if (sigma.length == 0 || sigma_remaining == 0) {
            //console.log("no more structs in poped R");
            var temp_sos = {structure: xmat.conv_str(P, seq_length), traces: traces};
            SOS.push(temp_sos);

            //console.log("pushed SOS:", JSON.stringify(SOS));
        }

        else {
            //console.log(formula);
            // compute maximal number of structures still to compute
            var maxLengthR = maxSOS - SOS.length;
            if (maxLengthR < 0) maxLengthR = 0;
            var R_prime = formula.getSubstructures(sigma, P, traces, delta, maxLengthR);
            for (var r in R_prime) {
                R.push(R_prime[r]);
            }
        }
        // check if enough structures found so far
        if (SOS.length >= maxSOS)
            break;
        //console.log("R:", JSON.stringify(R));

    }
    //console.log("SOS:", JSON.stringify(SOS));
    //console.log("\nwuchty end");
    //console.log(SOS);
    return SOS;
}
;


/** Nussinov Structures Count*/

var NussinovDPAlgorithm_structuresCount = Object.create(DPAlgorithm);

NussinovDPAlgorithm_structuresCount.Description = "Nussinov counting";

NussinovDPAlgorithm_structuresCount.Tables = new Array();
NussinovDPAlgorithm_structuresCount.Tables.push(Object.create(NussinovMatrix));

NussinovDPAlgorithm_structuresCount.Tables[0].latex_representation = "C_{i,j} = C_{i,j-1} + \\sum_{i\\leq k <(j-l) \\atop S_k,S_j \\text{ pair}} C_{i,k-1} \\cdot C_{k+1,j-1} \\cdot 1";

// C(i, j) = C(i, j - 1) + sum[k: [i <= k < j - l] && k,j pairs] C(i, k - 1) * C(k + 1, j - 1)
NussinovDPAlgorithm_structuresCount.Tables[0].computeValue = function (i, j) {
    if (i > j + 1 || i < 0 || j < 0 || i >= this.getDim() || j >= this.getDim()) {
        return 0;
    }
    if (i >= j) {
        return 1;
    }
    var res = 0;

    //console.log(i, j);
    // unpaired
    res += this.getValue(i, j - 1);
    for (var k = i; k < j - this.minLoopLength; ++k) {
        if (RnaUtil.areComplementary(this.sequence[k - 1], this.sequence[j - 1])) {
            res += this.getValue(i, k - 1) * this.getValue(k + 1, j - 1) * 1;
        }
    }

    return res;
};


// Invoking getValue for all tables, so that the Dynamic Programming computes all the tables and memoizes them.
NussinovDPAlgorithm_structuresCount.computeMatrix = function (input) {
// resize and initialize matrix
console.log("in struct count");
    this.Tables[0].init(input.sequence(), "structuresCount");
    // store minimal loop length
    var minLL = parseInt(input.loopLength());
    this.Tables[0].minLoopLength = minLL;

    for (var i = 0; i < this.Tables[0].getDim(); i++) {
        for (var j = 0; j < this.Tables[0].getDim(); ++j) {
            // get column for current span
            this.Tables[0].getValue(i, j);
        }
        ;
    }
    ;

    return this.Tables;
};

/****** McCaskill_simple extending NussinovMatrix ************************/

var NussinovDPAlgorithm_McCaskill = Object.create(DPAlgorithm);

NussinovDPAlgorithm_McCaskill.Description = "McCaskill";

NussinovDPAlgorithm_McCaskill.Tables = new Array();
NussinovDPAlgorithm_McCaskill.Tables.push(Object.create(NussinovMatrix)); // Q
NussinovDPAlgorithm_McCaskill.Tables.push(Object.create(NussinovMatrix)); // Qb
NussinovDPAlgorithm_McCaskill.Tables.push(Object.create(NussinovMatrix)); // Pe
NussinovDPAlgorithm_McCaskill.Tables.push(Object.create(NussinovMatrix)); // Pu

NussinovDPAlgorithm_McCaskill.Tables[0].latex_representation = "Q_{i,j} = Q_{i,j-1} + \\sum_{i\\leq k <(j-l)} Q_{i,k-1} \\cdot Q^{b}_{k,j}";
NussinovDPAlgorithm_McCaskill.Tables[1].latex_representation = "Q_{i,j}^{b} = \\begin{cases} Q_{i + 1, j - 1} \\cdot \\exp(-E_{bp}/RT) & \\text{ if }i,j \\text{ can form base pair} \\\\ 0 & \\text{ otherwise}\\end{cases}";
NussinovDPAlgorithm_McCaskill.Tables[2].latex_representation = "P^{bp}_{i, j} = Q^{-1}_{1, n} \\cdot (Q_{1, i - 1} \\cdot Q^{b}_{i, j} \\cdot Q_{j + 1, n})";
NussinovDPAlgorithm_McCaskill.Tables[3].latex_representation = "P^{U}_{i, j} = Q^{-1}_{1, n} \\cdot (Q_{1, i - 1} \\cdot 1 \\cdot Q_{j + 1, n})";

// Q(i,j) = sum[k : [i <= j < j - l] && k,j can pair] Q(i, k - 1) * Qb(k, j)
NussinovDPAlgorithm_McCaskill.Tables[0].computeValue = function (i, j) {
    if (i > j + 1 || i < 0 || j < 0 || i > this.getDim() || j > this.getDim()) {
        return 0;
    }
    if (i >= j) {
        return 1;
    }
    var res = 0;
    res += this.getValue(i, j - 1);
    for (var k = i; k < j - this.minLoopLength; ++k) {
        if (RnaUtil.areComplementary(this.sequence[k - 1], this.sequence[j - 1])) {
            res += this.getValue(i, k - 1) * NussinovDPAlgorithm_McCaskill.Tables[1].getValue(k, j);
        }
    }

    return res;
};

// Qb(i, j) = Q(i + 1, j - 1) * exp(-Eb / RT)
NussinovDPAlgorithm_McCaskill.Tables[1].computeValue = function (i, j) {
    if (i < 0 || j < 0 || i >= this.getDim() || j >= this.getDim() || i >= j - this.minLoopLength) {
        return 0;
    }
    if (RnaUtil.areComplementary(this.sequence[i - 1], this.sequence[j - 1])) {
        return NussinovDPAlgorithm_McCaskill.Tables[0].getValue(i + 1, j - 1) * Math.exp(-this.energy_basepair);
    } else {
        return 0;
    }
};

// Probability that i, j is a base pair.
// Pe(i, j) = Q(1, i - 1) * Qb(i, j) * Q(j + 1, n) / Q(1, n)
NussinovDPAlgorithm_McCaskill.Tables[2].computeValue = function(i, j) {
    if (i < 0 || j < 0) {//} || i >= this.getDim() || j >= this.getDim()) {
        return 0;
    }
    var n = this.getDim() - 1;

    var ret = NussinovDPAlgorithm_McCaskill.Tables[1].getValue(i, j);
    if (i>1) {
        ret *= NussinovDPAlgorithm_McCaskill.Tables[0].getValue(1, i - 1);
    }
    if (j<n) {
        ret *= NussinovDPAlgorithm_McCaskill.Tables[0].getValue(j + 1, n);
    }
    ret /= NussinovDPAlgorithm_McCaskill.Tables[0].getValue(1, n);

    return ret;
};


// Probability that i, j is unpaired.
// Pu(i, j) = Q(1, i - 1) * Q(j + 1, n) / Q(1, n)
NussinovDPAlgorithm_McCaskill.Tables[3].computeValue = function(i, j) {
    if (i < 0 || j < 0 || i > j) {//} || i >= this.getDim() || j >= this.getDim() || i > j) {
        return 0;
    }
    var n = this.getDim() - 1;

    var ret = 1.0;
    if (i>1) {
        ret *= NussinovDPAlgorithm_McCaskill.Tables[0].getValue(1, i - 1);
    }
    if (j<n) {
        ret *= NussinovDPAlgorithm_McCaskill.Tables[0].getValue(j + 1, n);
    }
    ret /= NussinovDPAlgorithm_McCaskill.Tables[0].getValue(1, n);
    return ret;
};

// Invoking getValue for all tables, so that the Dynamic Programming computes all the tables and memoizes them.
NussinovDPAlgorithm_McCaskill.computeMatrix = function (input) {
    this.Tables[0].init(input.sequence(), "McCaskill");
    this.Tables[1].init(input.sequence(), "McCaskill Base");
    this.Tables[2].init(input.sequence(), "McCaskill External Base");
    this.Tables[3].init(input.sequence(), "McCaskill Unpaired");
    // store minimal loop length
    var minLL = parseInt(input.loopLength());
    this.Tables[0].minLoopLength = minLL;

    this.Tables[1].energy_basepair = input.energy();

    for (var i = 0; i <= this.Tables[0].getDim(); i++) {
        for (var j = 0; j <= this.Tables[0].getDim(); ++j) {
            // get column for current span
            this.Tables[0].getValue(i, j);
        }
        ;
    }
    ;

    for (var i = 0; i <= this.Tables[1].getDim(); i++) {
        for (var j = 0; j <= this.Tables[1].getDim(); ++j) {
            // get column for current span
            this.Tables[1].getValue(i, j);
        }
        ;
    }
    ;

    for (var i = 0; i <= this.Tables[2].getDim(); i++) {
        for (var j = 0; j <= this.Tables[2].getDim(); ++j) {
            // get column for current span
            this.Tables[2].getValue(i, j);
        }
        ;
    }
    ;

    for (var i = 0; i <= this.Tables[3].getDim(); i++) {
        for (var j = 0; j <= this.Tables[3].getDim(); ++j) {
            // get column for current span
            this.Tables[3].getValue(i, j);
        }
        ;
    }
    ;

    return this.Tables;
};


var DPAlgorithm_MEA = Object.create(DPAlgorithm);

DPAlgorithm_MEA.Description = "Maximum Expected Accuracy";
DPAlgorithm_MEA.Tables = new Array();
DPAlgorithm_MEA.Tables.push(Object.create(NussinovMatrix));
DPAlgorithm_MEA.Tables.push(Object.create(NussinovMatrix));
DPAlgorithm_MEA.Tables.push(Object.create(NussinovMatrix));
DPAlgorithm_MEA.Tables[0].latex_representation = "M_{i, j} = \\max \\begin{cases} 0 & \\text{i > j} \\\\ M_{i, j - 1} + P^{U}_{j} & \\text{j unpaired} \\\\ M_{i + 1, j - 1} + P^{bp}_{i,j} & \\text{j paired with i} \\end{cases}";
//DPAlgorithm_MEA.Tables[0].latex_representation = "M_{i, j} = \\max \\begin{cases} 0 & i > j \\\\ M_{i, j - 1} + p^{u}_{j} & j unpaired \\\\ M_{i + 1, j - 1} + p^{p}_{i,j} & j paired with i \\\\ \\max_{i \\leq k < j}{M_{i, k} + M_{k + 1, j}} & decomposition \\end{cases}";
//DPAlgorithm_MEA.Tables[0].latex_representation = "D(i,j) = \\max \\begin{cases} D(i+1,j) & S_i \\text{ unpaired} \\\\ D(i,j-1) & S_j \\text{ unpaired} \\\\ D(i+1,j-1)+1 &  S_i,S_j \\text{ compl. base pair and } i+ l< j \\\\ \\max_{i< k< (j-1)} D(i,k)+D(k+1,j) & \\text{decomposition} \\end{cases}";
DPAlgorithm_MEA.Tables[0].updateCell = function (i, j, curVal, curAncestor) {

    var curCell = this.getCell(i, j);
    if (curCell === null || curCell.value <= curVal) {
        // check for new maximal value
        if (curCell === null || curCell.value < curVal) {
            // reset ancestor list
            curCell.traces = [];
            // store new maximum
            curCell.value = curVal;
        }
        ;
        // store this ancestor
        curCell.traces.push(curAncestor);
    }
    ;

}
// "{"parents":[[2,1]],"bps":[[1,2]]}"
DPAlgorithm_MEA.Tables[0].computeValue = function(i, j) {

    if (i > j || i < 0 || j < 0 || i >= this.getDim() || j >= this.getDim()) {
        return 0;
    }
    var unpaired_j = 1.0;
    for (var k = j - 1; k >= 0; --k) {
        unpaired_j -= NussinovDPAlgorithm_McCaskill.Tables[2].getValue(k, j);
    }
    this.updateCell(i, j, this.getValue(i, j - 1) + unpaired_j, Object.create(NussinovCellTrace).init([[i, j - 1]], []));

    if (RnaUtil.areComplementary(this.sequence[i - 1], this.sequence[j - 1])) {
        this.updateCell(i, j, this.getValue(i + 1, j - 1) + NussinovDPAlgorithm_McCaskill.Tables[2].getValue(i, j), Object.create(NussinovCellTrace).init([[i + 1, j - 1]], [i, j]));
    }

    for (var k = i; k < j - this.minLoopLength; ++k) {
        this.updateCell(i, j, this.getValue(i, k) + this.getValue(k + 1, j), Object.create(NussinovCellTrace).init([[i, k], [k + 1, j]], []));
    }

    return this.getCell(i, j).value;
};

DPAlgorithm_MEA.computeMatrix = function(input) {
    
    NussinovDPAlgorithm_McCaskill.computeMatrix(input);

    this.Tables[1] = NussinovDPAlgorithm_McCaskill.Tables[2];
    this.Tables[2] = NussinovDPAlgorithm_McCaskill.Tables[3];

    this.Tables[0].init(input.sequence(), "Maximum Expected Accuracy");
    // store minimal loop length
    var minLL = parseInt(input.loopLength());
    this.Tables[0].minLoopLength = minLL;

    for (var i = 0; i < this.Tables[0].getDim(); i++) {
        for (var j = 0; j < this.Tables[0].getDim(); ++j) {
            // get column for current span
            this.Tables[0].getValue(i, j);
        }
        ;
    }
    ;

    return this.Tables;
};


DPAlgorithm_MEA.Tables[0].getSubstructures = function (sigma, P, traces, delta, maxLengthR) {
    var Nmax = this.getCell(1, this.sequence.length).value;
    var R = [];
    var ij = sigma.pop();
    //console.log(ij);

    // check for sane interval
    // if i>j dont continue
    if (ij[0] > ij[1]) {
        //console.log("ij[0] > ij[1]", ij[0], ij[1]);
        var S_prime = {};

        var sigma_prime = JSON.stringify(sigma);
        sigma_prime = JSON.parse(sigma_prime);

        var tmp_P = JSON.stringify(P);
        tmp_P = JSON.parse(tmp_P);

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);

        S_prime.sigma = sigma_prime;
        S_prime.P = tmp_P;
        S_prime.traces = tmp_traces;

        R.push(S_prime);
        //console.log("returning R:", JSON.stringify(R));
        return R;
    }

    // if (i,j) == (i+1,j-1) + bp(ij)
    {
        if (ij[1] - ij[0] > this.minLoopLength) {
            //console.log(this.sequence);
            //console.log(this.sequence[ij[0] - 1], this.sequence[ij[1] - 1]);
            if (RnaUtil.areComplementary(this.sequence[ij[0] - 1], this.sequence[ij[1] - 1])) {
                var sigma_prime = JSON.stringify(sigma);
                sigma_prime = JSON.parse(sigma_prime);
                sigma_prime.push([ij[0] + 1, ij[1] - 1]);

                var tmp_P = JSON.stringify(P);
                tmp_P = JSON.parse(tmp_P);
                tmp_P.push([ij[0], ij[1]]);

                var tmp_traces = JSON.stringify(traces);
                tmp_traces = JSON.parse(tmp_traces);

                var NSprime = this.countBasepairs(tmp_P, sigma_prime);

                if (NSprime >= Nmax - delta) {
                    var S_prime = {};
                    S_prime.sigma = sigma_prime;
                    S_prime.P = tmp_P;
                    tmp_traces.unshift([ij, [[ij[0] + 1, ij[1] - 1]]]);
                    S_prime.traces = tmp_traces;
                    //console.log("i+1,j-1:", JSON.stringify(S_prime));
                    // push to the front to keep base pair most prominent to refine
                    R.unshift(S_prime);
                }
            }
        }

        // check if enough structures found so far
        if (R.length >= maxLengthR) {
            //console.log("returning R:", JSON.stringify(R));
            return R;
        }
    }

    // if (i,j) == (i,j-1)
    {
        var sigma_prime = JSON.stringify(sigma);
        sigma_prime = JSON.parse(sigma_prime);
        sigma_prime.unshift([ij[0], ij[1] - 1]);

        var tmp_P = JSON.stringify(P);
        tmp_P = JSON.parse(tmp_P);

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);

        var NSprime = this.countBasepairs(P, sigma_prime);

        if (NSprime >= Nmax - delta) {
            var S_prime = {};
            S_prime.sigma = sigma_prime;
            S_prime.P = tmp_P;
            tmp_traces.unshift([ij, [[ij[0], ij[1] - 1]]]);
            S_prime.traces = tmp_traces;
            //console.log("i,j-1:", JSON.stringify(S_prime));
            // push to the front to keep base pair most prominent to refine
            R.unshift(S_prime);
        }

        // check if enough structures found so far
        if (R.length >= maxLengthR) {
            //console.log("returning R:", JSON.stringify(R));
            return R;
        }
    }


    console.log("returning R:", JSON.stringify(R));
    return R;
}
;

var DPAlgorithm_nussiFold = Object.create(DPAlgorithm);

DPAlgorithm_nussiFold.Description = "NussinovFolding";
DPAlgorithm_nussiFold.Tables = new Array();
DPAlgorithm_nussiFold.Tables.push(Object.create(NussinovMatrix));

DPAlgorithm_nussiFold.Tables[0].latex_representation = "D(i,j) = \\max \\begin{cases} D(i,j-1) & S_j \\text{ unpaired} \\\\ \\max_{i\\leq k< (j-l)} D(i,k-1)+D(k+1,j-1)+1 & S_k,S_j \\text{ compl. base pair} \\end{cases}";

DPAlgorithm_nussiFold.computeMatrix = function(input) {

    NussinovDPAlgorithm_Unique.computeMatrix(input);

    this.Tables[0].init(input.sequence(), "NussinovFold");
    this.Tables[0].cells = NussinovDPAlgorithm_Unique.Tables[0].cells;
    return this.Tables;
};

DPAlgorithm_nussiFold.Tables[0].getSubstructures = function (sigma, P, traces, delta, maxLengthR) {
    return NussinovDPAlgorithm_Unique.Tables[0].getSubstructures(sigma, P, traces, delta, maxLengthR);
};

/**
 * global list of available Nussinov algorithm implementations
 */
var availableAlgorithms = {

    /** ambiguous recursion */
    nussinovOriginal: NussinovDPAlgorithm_Ambiguous,//NussinovMatrix_ambiguous,

    /** original unique recursion */
    nussinovUnique: NussinovDPAlgorithm_Unique,//NussinovMatrix_unique,

    /** nussinov neo recursion */
    nussinovAmbiguous2: NussinovDPAlgorithm_Ambiguous2,

    /** McCaskill */
    mcCaskill: NussinovDPAlgorithm_McCaskill,

    /** structure counting */
    nussinovCounting: NussinovDPAlgorithm_structuresCount,

    /** Maximum Expected Accuracy*/
    MaxExpAcc: DPAlgorithm_MEA,

    /** nussinov fold*/
    nussinovFold: DPAlgorithm_nussiFold,

};

