/**
 * @file Main file containing backend algorithms for RNA-algorithms-JS project.
 Main Items contains in theis file are:
 -Matrix class
 -Nussinov algorithms
 -Traceback algorithms
 * @authors "Mostafa Mahmoud Mohamed", "Syed Mohsin Ali", "Martin Mann"
 */
"use strict";

/**
 * Utility class that covers RNA specific functions.
 */

var EPS = 1e-4; // floating point precision error

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
 * data stored within a cell of a nussinov matrix
 */
var NussinovCell = {

// row
    i: -1,

// column
    j: -1,

// value
    value: null,

    logValue: null,

// traces for the current value
    traces: null,

    /**
     * Initialize a cell with the given data and sets traces to an empty list
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
        this.logValue = null;
        // this access for chaining
        return this;
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
 * Nussinov matrix object, stores Sequence and Table. Contains all the utility functionalities for the tables.
 */
var NussinovMatrix = {

    /**
     * Access to the sequence for this matrix
     */
    sequence: null,

    /**
     * Sequence length
     * */
    seq_length: 0,

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
     * The dimensions of the matrix.
     * */
    tablesDimension: 2,

    /**
     * initialize a matrix of dim = (n + 1) x n, where n is the length of the provided sequence
     * @param {string} sequence the RNA sequence (not null or empty)
     * @returns {NussinovMatrix} this
     */
    init: function(sequence, name) {
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
        this.seq_length = sequence.length;

        // create matrix cells
        for (var i = 0; i <= this.seq_length; i++) {
            this.cells[i] = [];
            for (var j = 0; j <= this.seq_length; j++) {
                // create new cell and initialize
                this.cells[i][j] = Object.create(NussinovCell).init(i, j, null);
            }
        }

        return this;
    },

    /**
     * Check if a given tuple is an invalid state or not
     * @param i row #
     * @param j column #
     * @returns {boolean}
     * */
    isInvalidState: function(i, j) {
        if (i < 0 || j < 0 || i > this.seq_length || j > this.seq_length || i > j + 1) {
            return true;
        } else {
            return false;
        }
    },

    /**
     * Compute the cell at a given state in the matrix.
     * It's recommended to make the implementation use the method "updateCell",
     * if it's computing the tracebacks in an optimization problem.
     * TODO: this function has to be overwritten by the instances, before calling computeMatrix/computeAllCells
     * @param i row #
     * @param j column #
     * @returns {NussinovCell}  The computed cell.
     */
    computeCell: function(i, j) {
       // updateCell(i, j);
        return Object.create(NussinovCell).init(i, j, null);
    },

    /**
     * Access a cell at a given state in the matrix. If the cell is null or has no value,
     * then it will be computed using the "computeCell" method.
     * TODO: Implement computeCell
     * @param i row #
     * @param j column #
     * @returns {NussinovCell}  The cell or null if it's an invalid state
     */
    getCell: function (i, j) {
        // check border cases {
        if (this.isInvalidState(i, j)) {
            return null;
        }
        if (this.cells[i][j] === null || this.cells[i][j].value == null) {
            this.cells[i][j] = this.computeCell(i, j);
        }
        return this.cells[i][j];
    },

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
    },

    /**
     * Access the value at a given state in the matrix. If the cell is null or has no value,
     * then the cell will be computed using the "computeCell" method.
     * TODO: Implement computeCell
     * @param i row #
     * @param j column #
     * @returns {float} Cell value or null if invalid cell
     */
    getValue: function (i, j) {
        var cell = this.getCell(i, j);
        if (cell === null) {
            return null;
        }
        return parseFloat(cell.value);
    },

    /**
     * Updates the ancestor list of a given cell if the curVal is higher or
     * equal to the current value within the cell.
     * If the value is equal, curAncestor is added to the list.
     * If the value is smaller than curVal, curAncestor will be set to be the
     * only list entry.
     * TODO: this function can be overwritten by the instances
     * @param curCell The current cell to be updated.
     * @param curAncestor A list of the 4dTraces of the tracebacks at this state.
     */
    updateCell: function (curCell, curAncestor) {
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
            // store this ancestor
            curCell.traces.push(curAncestor);
        }
    },

    /**
     * Compute all the cells of the matrix.
     * TODO: Implement computeCell
     */
    computeAllCells: function() {
        for (var i = 0; i <= this.seq_length; ++i) {
            for (var j = 0; j <= this.seq_length; ++j) {
                this.getCell(i, j);
            }
        }
    },


    /**
     * Fills the matrix according to the recursion.
     * TODO: Implement computeCell
     * TODO: this function has to be overwritten by the instances.
     *
     * @param {input} A Dictionary of the input for the matrix. Should contain all the arguments
     *                needed for initalizing the matrix input properly. Minimally this will be the sequence.
     * @returns {NussinovMatrix} this for call chaining
     */
    computeMatrix: function (input) {
        console.log("WARNING: computeMatrix() not implemented in NussinovMatrix superclass; overwrite in subclass!");

        // resize and init matrix
        this.init(input.sequence(), "Default name");

        // set minimal loop length
        this.minLoopLength = parseInt(input.loopLength());

        this.computeAllCells();

        return this;
    },

    /**
     * Access to the recursion's representation in LaTeX, that is used in by this matrix.
     * @returns {string} latex encoding of the recursion
     */
    getRecursionInLatex: function () {
        return "$$" + this.latex_representation + "$$";
    },

    /**
     * creates a string representation of the matrix
     * @returns {string} matrix as string
     */
    toString: function () {
        var str = this.minLoopLength + "    ";
        for (var i = 0; i < this.seq_length; i++) {
            str += this.sequence[i] + "  ";
        }
        str += "\n";
        for (var i = 1; i <= this.seq_length; i++) {
            // print sequence
            str += this.sequence[i - 1] + " ";
            // print values
            for (var j = 0; j <= this.seq_length; j++) {
                if (j !== 0) {
                    str += ", ";
                }
                str += this.getValue(i, j);
            }
            str += "\n";
        }
        return str;
    },

    /**
     * Compute an Sprime instance that can be used wuchty.
     * Given a trace of a cell @cell and the value of the trace at this cell @NSprime and the optimal value at this
     * cell @Nmax. and returns the SPrime for this trace.
     * It assumes that  NSprime >= Nmax - delta (The error in this trace is less than the remaining potential error)
     * @NSprime: The value of this trace.
     * @NMax: The optimal value at this cell.
     * @cell: The current cell
     * @trace: The current trace being considered
     * @sigma: Old sigma to be cloned
     * @delta: The remaining error potential
     * @P: The basepair structure so far.
     * @traces: A list of pairs of cells and their parents.
     * @returns S_prime for this trace.
     * */
    getSPrime: function(NSprime, Nmax, cell, trace, sigma, delta, P, traces) {
        var sigma_prime = JSON.stringify(sigma);
        sigma_prime = JSON.parse(sigma_prime);
        //console.log(trace);
        for (var i = 0; i < trace.parents.length; ++i) {
            sigma_prime.push(trace.parents[i]);
        }

        var tmp_P = JSON.stringify(P);
        tmp_P = JSON.parse(tmp_P);
        for (var i = 0; i < trace.bps.length; ++i) {
            tmp_P.push(trace.bps[i]);
        }

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);
        tmp_traces.unshift([cell, trace.parents]);


        var S_prime = {};
        S_prime.sigma = sigma_prime;
        S_prime.P = tmp_P;
        S_prime.traces = tmp_traces;
        S_prime.potential = delta - (Nmax - NSprime);
        
        return S_prime;
    },


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

    /**
     * Construct a one string representing the matching base pairs in the sequence.
     * @param x A list of pairs of indices, that represent the matching base pairs.
     * @returns {string} A representation of the matching base pairs in the matching.
     */
    conv_str: function (x, length) {
        var str = "";
        for (var l = 0; l < length; l++) {
            str += ".";
        }
        var linked = this.sequence.indexOf("X");
        if(linked == -1){
            for (var i in x) {
                str = str.substr(0, x[i][0] - 1) + "(" + str.substr(x[i][0], str.length);
                str = str.substr(0, x[i][1] - 1) + ")" + str.substr(x[i][1], str.length);
            }
            return str;
        } else {
            for (var i in x) {
                if(x[i][0]  <= linked && x[i][1] >= linked + 1 + this.minLoopLength) {
                    str = str.substr(0, x[i][0] - 1) + "[" + str.substr(x[i][0], str.length);
                    str = str.substr(0, x[i][1] - 1) + "]" + str.substr(x[i][1], str.length);
                }
                else{
                    str = str.substr(0, x[i][0] - 1) + "(" + str.substr(x[i][0], str.length);
                    str = str.substr(0, x[i][1] - 1) + ")" + str.substr(x[i][1], str.length);

                }
            }
            var st = "";
            for (var l = 0; l < this.minLoopLength + 1; ++l) {
                st += "X";
            }
            str = str.substr(0, linked) + st + str.substr(linked + this.minLoopLength + 1, str.length);
            return str;
        }
    },

};

/**
 * Dynamic programming algorithm.
 *
 * DP Algorithms will work by memoization usually. If so, there's a function getCell/getValue for each of the tables,
 * that computes a value and memoize it if it's not computed(default null cells), and returns the memoized value.
 *
 * The computation is done through the computeCell in each of the tables, this function should be overriden for each
 * table depending on how an entry in the table is computed. The function should include the base case, and should use
 * the other tables by the getValue function and not by accessing the tables directly, in order to ensure the
 * correctness of the memoization.
 *
 * The tables should be usable after invoking the computeMatrix method, this method should be overriden to set the tables
 * parameters and compute all the dynamic programming. Unless there's a special way to do this, it can usually
 * (with memoization) be done by invoking getCell for all the entries of all the tables (or computeAllCells for each table).
 *
 * TODO: How To Use:
 *  * Create the DPAlgorithm instance
 *  * Create new Tables Array, and push the needed tables(NussinovMatrix/NussinovMatrix4d)
 *  * Override latex_representation for each table.
 *  * Override computeCell for each table, and/or isInvalidState, and/or updateCell (In case of storing traceback)
 *      * Remember to use getCell/getValue instead of accessing the cell directly, to preserve the memoization.
 *      * If a given state satisfies isInvalidState, then it should return an invalid default value. Like INF in
 *        minimization algorithms, 0 in counting algorithms, and -INF in maximization algorithms.
 *  * Override getSubstructures (In case of computing tracebacks)
 *  * Override ComputeMatrix, which is the main interface for computing all the tables.
 *      * Can use the instance of other algorithms and clone their tables. (Clone also the needed methods like getCell)
 *      * Calling computeAllCells for each Table should be sufficient most of the time.
 *
 *
 * @type {{Description: Algorithm description,
 *        Tables: Array of tables for all the recursive formulas,
 *        computeMatrix: Compute All the Tables
 *        getRecursionInLatex: Get Latex String describing the recursive equations for all the tables.}}
 */
var DPAlgorithm = {
    /**
     * Algorithm description.
     */
    Description: "Algorithm",

    /**
     * A list of the Tables NussinovMatrix (or NussinovMatrix4d) for each table used by the DP algorithms.
     * Note: Create a new array for each instance.
     */
    Tables: [], // create new array

    /**
     * TODO: Has to be overriden by the instances.
     * @param input A dictionary of the input arguments.
     */
    computeMatrix: function (input) { },

    /**
     * @returns {string} An aligned latex array that contains the latex formula of each table(seperated with empty lines).
     */
    getRecursionInLatex: function () {
        var formula = " \\begin{array} ";
        for (var i = 0; i < this.Tables.length; ++i) {
            formula += " \\\\ \\\\ " + this.Tables[i].latex_representation;
        }
        formula += " \\end{array} ";
        return formula;
    },

};


var NussinovDPAlgorithm_Ambiguous = Object.create(DPAlgorithm);

NussinovDPAlgorithm_Ambiguous.Description = "Ambiguous recursion";
NussinovDPAlgorithm_Ambiguous.Tables = new Array();
NussinovDPAlgorithm_Ambiguous.Tables.push(Object.create(NussinovMatrix));
NussinovDPAlgorithm_Ambiguous.Tables[0].latex_representation = "D(i,j) = \\max \\begin{cases} D(i+1,j) & S_i \\text{ unpaired} \\\\ D(i,j-1) & S_j \\text{ unpaired} \\\\ D(i+1,j-1)+1 & \\text{if } S_i,S_j \\text{ compl. base pair and } i+ l< j \\\\ \\max_{i< k< (j-1)} D(i,k)+D(k+1,j) & \\text{ decomposition} \\end{cases}";

NussinovDPAlgorithm_Ambiguous.Tables[0].computeCell = function(i, j) {

    var curCell = Object.create(NussinovCell).init(i, j, 0);

    if (this.isInvalidState(i, j)) {
        return curCell;
    }
    // i unpaired
    this.updateCell(curCell, Object.create(NussinovCellTrace).init([[i + 1, j]], []));

    // j unpaired
    this.updateCell(curCell, Object.create(NussinovCellTrace).init([[i, j - 1]], []));

    // check (i,j) base pair
    if ((j - i > this.minLoopLength) && RnaUtil.areComplementary(this.sequence[i - 1], this.sequence[j - 1])) {
        // get value for base pair
        this.updateCell(curCell, Object.create(NussinovCellTrace).init([[i + 1, j - 1]], [[i, j]]));
    }

    // check decomposition into substructures (minLength==2)
    for (var k = i + 1; k < (j - 1); k++) {
        // get decomposition value
        this.updateCell(curCell, Object.create(NussinovCellTrace).init([[i, k], [k + 1, j]], []));
    }

    return curCell;
};

NussinovDPAlgorithm_Ambiguous.computeMatrix = function (input) {
// resize and initialize matrix
    this.Tables[0].init(input.sequence(), "ambiguous");
// store minimal loop length
    this.Tables[0].minLoopLength = parseInt(input.loopLength());

    this.Tables[0].computeAllCells();

    return this.Tables;
};

NussinovDPAlgorithm_Ambiguous.Tables[0].getSubstructures = function (sigma, P, traces, delta, maxLengthR) {
    var Nmax = this.getValue(1, this.sequence.length);
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

        var P_prime = JSON.stringify(P);
        P_prime = JSON.parse(P_prime);

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);

        S_prime.sigma = sigma_prime;
        S_prime.P = P_prime;
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

                var P_prime = JSON.stringify(P);
                P_prime = JSON.parse(P_prime);
                P_prime.push([ij[0], ij[1]]);

                var tmp_traces = JSON.stringify(traces);
                tmp_traces = JSON.parse(tmp_traces);

                var NSprime = this.countBasepairs(P_prime, sigma_prime);

                if (NSprime >= Nmax - delta) {
                    var S_prime = {};
                    S_prime.sigma = sigma_prime;
                    S_prime.P = P_prime;
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

        var P_prime = JSON.stringify(P);
        P_prime = JSON.parse(P_prime);

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);

        var NSprime = this.countBasepairs(P_prime, sigma_prime);

        if (NSprime >= Nmax - delta) {
            var S_prime = {};
            S_prime.sigma = sigma_prime;
            S_prime.P = P_prime;
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

        var P_prime = JSON.stringify(P);
        P_prime = JSON.parse(P_prime);

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);

        var NSprime = this.countBasepairs(P_prime, sigma_prime);

        if (NSprime >= Nmax - delta) {
            var S_prime = {};
            S_prime.sigma = sigma_prime;
            S_prime.P = P_prime;
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

        var P_prime = JSON.stringify(P);
        P_prime = JSON.parse(P_prime);

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);

        var NSprime = this.countBasepairs(P_prime, sigma_prime);

        if (NSprime >= Nmax - delta) {

            var S_prime = {};
            S_prime.sigma = sigma_prime;
            S_prime.P = P_prime;
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

    //console.log("returning R:", JSON.stringify(R));
    return R;
}
;


var NussinovDPAlgorithm_Unique = Object.create(DPAlgorithm);

NussinovDPAlgorithm_Unique.Description = "Recursion by Nussinov et al. (1978) with unique decomposition";
NussinovDPAlgorithm_Unique.Tables = new Array();
NussinovDPAlgorithm_Unique.Tables.push(Object.create(NussinovMatrix));
NussinovDPAlgorithm_Unique.Tables[0].latex_representation = "D(i,j) = \\max \\begin{cases} D(i,j-1) & S_j \\text{ unpaired} \\\\ \\max_{i\\leq k< (j-l)} D(i,k-1)+D(k+1,j-1)+1 & \\text{if } S_k,S_j \\text{ compl. base pair} \\end{cases}";

NussinovDPAlgorithm_Unique.Tables[0].computeCell = function(i, j) {

    var curCell = Object.create(NussinovCell).init(i, j, 0);

    if (this.isInvalidState(i, j)) {
        return curCell;
    }
    // j unpaired
    this.updateCell(curCell, Object.create(NussinovCellTrace).init([[i, j - 1]], []));

    // check base pair based decomposition : (k,j) base pair
    for (var k = i; k + this.minLoopLength < j; k++) {
        // check if sequence positions are compatible
        if (RnaUtil.areComplementary(this.sequence[k - 1], this.sequence[j - 1])) {
            this.updateCell(curCell, Object.create(NussinovCellTrace).init([[i, k - 1], [k + 1, j - 1]], [[k, j]]));
        }
    }

    return curCell;
};

NussinovDPAlgorithm_Unique.computeMatrix = function (input) {

// resize and initialize matrix
    this.Tables[0].init(input.sequence(), "unique");
// store minimal loop length
    this.Tables[0].minLoopLength = parseInt(input.loopLength());

    this.Tables[0].computeAllCells();

    return this.Tables;
};


NussinovDPAlgorithm_Unique.Tables[0].getSubstructures = function (sigma, P, traces, delta, maxLengthR) {
    var Nmax = this.getValue(1, this.sequence.length);
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

    // [<[(i, l - 1], (l + 1, j - 1), -1,>, <>, ....]

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

        var NSprime = this.countBasepairs(tmp_P, sigma_prime);

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
NussinovDPAlgorithm_Ambiguous2.Tables[0].latex_representation = "D(i,j) = \\max \\begin{cases} D(i+1,j-1)+1 & \\text{if } S_i,S_j \\text{ compl. base pair} \\\\ \\max_{i\\leq k< j} D(i,k)+D(k+1,j) \\end{cases}";

NussinovDPAlgorithm_Ambiguous2.Tables[0].computeCell = function(i, j) {

    var curCell = Object.create(NussinovCell).init(i, j, 0);

    if (this.isInvalidState(i, j)) {
        return curCell;
    }

    // check (i,j) base pair
    if ((j - i > this.minLoopLength) && RnaUtil.areComplementary(this.sequence[i - 1], this.sequence[j - 1])) {
        // get value for base pair
        this.updateCell(curCell, Object.create(NussinovCellTrace).init([[i + 1, j - 1]], [[i, j]]));
    }

    // check decomposition into substructures (minLength==2)
    for (var k = i; k < j; k++) {
        // get decomposition value
        this.updateCell(curCell, Object.create(NussinovCellTrace).init([[i, k], [k + 1, j]], []));
    }

    return curCell;
};

NussinovDPAlgorithm_Ambiguous2.computeMatrix = function (input) {

// resize and initialize matrix
    this.Tables[0].init(input.sequence(), "Ambiguous2");
// store minimal loop length
    this.Tables[0].minLoopLength = parseInt(input.loopLength());

    this.Tables[0].computeAllCells();

    return this.Tables;
};


NussinovDPAlgorithm_Ambiguous2.Tables[0].getSubstructures = function (sigma, P, traces, delta, maxLengthR) {
    var Nmax = this.getValue(1, this.sequence.length);
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

    // if (i,j) == (i,l) + (l+1, j)
    for (var l = ij[0]; l < ij[1]; l++) {
        //console.log('here');
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
    //console.log("returning R:", JSON.stringify(R));
    return R;
};

var NussinovDPAlgorithm_MostAmbiguous = Object.create(DPAlgorithm);

NussinovDPAlgorithm_MostAmbiguous.Description = "Most ambiguous recursion";
NussinovDPAlgorithm_MostAmbiguous.Tables = new Array();
NussinovDPAlgorithm_MostAmbiguous.Tables.push(Object.create(NussinovMatrix));
NussinovDPAlgorithm_MostAmbiguous.Tables[0].latex_representation = "D(i,j) = \\max \\begin{cases} D(i,j-1) & S_j \\text{ unpaired}  \\\\ \\max_{i \\le (p+l) < r \\le j} D(i,p-1)+1+D(p+1,r-1)+D(r+1,j) & \\text{ decomposition} \\end{cases}";

NussinovDPAlgorithm_MostAmbiguous.Tables[0].computeCell = function(i, j) {

    var curCell = Object.create(NussinovCell).init(i, j, 0);

    if (this.isInvalidState(i, j)) {
        return curCell;
    }
    // j unpaired
    this.updateCell(curCell, Object.create(NussinovCellTrace).init([[i, j - 1]], []));

    // check base pair based decomposition : (p,r) base pair
    for (var p = i; p + this.minLoopLength < j; p++) {
        for(var r = p + this.minLoopLength + 1; r <= j; r++) {
            // check if sequence positions are compatible
            if (RnaUtil.areComplementary(this.sequence[p - 1], this.sequence[r - 1])) {
            this.updateCell(curCell, Object.create(NussinovCellTrace).init([[i, p - 1], [p + 1, r - 1], [r + 1, j]], [[p, r]]));
            }
        }   
    }

    return curCell;
};

NussinovDPAlgorithm_MostAmbiguous.computeMatrix = function (input) {

// resize and initialize matrix
    this.Tables[0].init(input.sequence(), "MostAmbiguous");
// store minimal loop length
    this.Tables[0].minLoopLength = parseInt(input.loopLength());
    
    this.Tables[0].computeAllCells();
    
    return this.Tables;
};

NussinovDPAlgorithm_MostAmbiguous.Tables[0].getSubstructures = function (sigma, P, traces, delta, maxLengthR) {
    var Nmax = this.getValue(1, this.sequence.length);
    var R = [];
    var ij = sigma.pop();

    // check for sane interval
    // if i>j dont continue
    if (ij[0] >= ij[1]) {
        //console.log("ij[0] > ij[1]", ij[0], ij[1]);
        var S_prime = {};

        var sigma_prime = JSON.stringify(sigma);
        sigma_prime = JSON.parse(sigma_prime);

        var p_prime = JSON.stringify(P);
        p_prime = JSON.parse(p_prime);

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);

        S_prime.sigma = sigma_prime;
        S_prime.P = p_prime;
        S_prime.traces = tmp_traces;

        R.push(S_prime);
        //console.log("returning R:", JSON.stringify(R));
        return R;
    }

    // if (i,j) == (i,j-1)
    {
        var sigma_prime = JSON.stringify(sigma);
        sigma_prime = JSON.parse(sigma_prime);
        sigma_prime.unshift([ij[0], ij[1] - 1]);

        var P_prime = JSON.stringify(P);
        P_prime = JSON.parse(P_prime);
        

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);

        var NSprime = this.countBasepairs(P_prime, sigma_prime);

        if (NSprime >= Nmax - delta) {
            var S_prime = {};
            S_prime.sigma = sigma_prime;
            S_prime.P = P_prime;
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

    // trace enclosed base pair (l,m) within i..j range
    for (var l = ij[0]; l < ij[1] - this.minLoopLength; l++) {
       for(var m = l + this.minLoopLength + 1; m <= ij[1]; m++) {
            // check if sequence positions are compatible
            
                if (RnaUtil.areComplementary(this.sequence[l - 1], this.sequence[m - 1])) {
            
                var sigma_prime = JSON.stringify(sigma);
                sigma_prime = JSON.parse(sigma_prime);
                sigma_prime.push([ij[0], l - 1]);
                sigma_prime.push([l + 1, m - 1]);
                sigma_prime.push([m + 1, ij[1]]) // generate sigmaPrime

                var p_prime = JSON.stringify(P);
                p_prime = JSON.parse(p_prime);
                p_prime.push([l, m]);

                var tmp_traces = JSON.stringify(traces);
                tmp_traces = JSON.parse(tmp_traces);

                var NSprime = this.countBasepairs(p_prime, sigma_prime);
                
                // if number of base pairs is within range
                if (NSprime >= Nmax - delta) {

                    var S_prime = {};
                    S_prime.sigma = sigma_prime;
                    S_prime.P = p_prime;
                    tmp_traces.unshift([ij, [[ij[0], l - 1], [l + 1, m - 1], [m + 1, ij[1]]]]);
                    S_prime.traces = tmp_traces;
                    //console.log("ilj:", JSON.stringify(S_prime));
                    // push to the front to keep base pair most prominent to refine
                    R.unshift(S_prime);// add trace to R
                }// end of number of base pairs within range


                // check if enough structures found so far
                if (R.length >= maxLengthR) {
                    //console.log("returning R:", JSON.stringify(R));
                    return R;
                }

            } // end of (l,m) complementary
        }   
    }

    // return all found traces for i..j range
    //console.log("returning R:", JSON.stringify(R));
    return R;

};









/**
 * WUCHTY(2nd version)
 */
function wuchty_unlimited(xmat, delta, formula, maxSOS) {
    //if (xmat == undefined)return;
    if (xmat.sequence == undefined)return;
    var seq_length = xmat.sequence.length;
    var Nmax = xmat.getValue(1, seq_length);
    //console.log("Wuchty beginning\nSequence:", xmat.sequence, "\nNmax:", Nmax, "\nDelta:", delta, "\n");

    var S = {sigma: [[1, seq_length]], P: [], traces: [], potential: delta + EPS};
    var R = [S];
    var SOS = [];
    var Ps = new Set();
    var loop = 0;

    while (R.length != 0) {
        // Pop R
        var pop_R = R.pop();
        var sigma = pop_R.sigma;
        var P = pop_R.P;
        var t_traces = JSON.stringify(pop_R.traces);
        var traces = JSON.parse(t_traces);
        var remaining_potential = pop_R.potential;
        var sigma_remaining = 0;
        for (var s in sigma) {
            // TODO(mostafa): Check the base case condition
            if ((sigma[s][0]) <= (sigma[s][1] - xmat.minLoopLength)) sigma_remaining++;
        }
        if (sigma.length == 0 || sigma_remaining == 0) {
            if (!Ps.has(JSON.stringify(P))) {
                Ps.add(JSON.stringify(P));
                //console.log(Ps);
                var temp_sos = {structure: xmat.conv_str(P, seq_length), traces: traces};
                SOS.push(temp_sos);
            }
        } else {
            //console.log(formula);
            // compute maximal number of structures still to compute
            var maxLengthR = maxSOS - SOS.length;
            if (maxLengthR < 0) maxLengthR = 0;
            var R_prime = formula.getSubstructures(sigma, P, traces, remaining_potential, maxLengthR);
            for (var r in R_prime) {
                R.push(R_prime[r]);
            }
        }
        // check if enough structures found so far
        if (SOS.length >= maxSOS)
            break;

    }
    return SOS;
};


/**
 * WUCHTY(2nd version)
 */
function wuchty_2nd_limited(xmat, delta, formula, maxSOS) {
    //if (xmat == undefined)return;
    if (xmat.sequence == undefined)return;
    var seq_length = xmat.sequence.length;
    var Nmax = xmat.getValue(1, seq_length);
    //console.log("Wuchty beginning\nSequence:", xmat.sequence, "\nNmax:", Nmax, "\nDelta:", delta, "\n");

    var S = {sigma: [[1, seq_length]], P: [], traces: []};
    var R = [S];
    var SOS = [];
    var loop = 0;

    while (R.length != 0) {
        // Pop R
        var pop_R = R.pop();
        var sigma = pop_R.sigma;
        var P = pop_R.P;
        var t_traces = JSON.stringify(pop_R.traces);
        var traces = JSON.parse(t_traces);

        var sigma_remaining = 0;
        for (var s in sigma) {
            // TODO(mostafa): Check the base case condition
            if ((sigma[s][0]) <= (sigma[s][1] - xmat.minLoopLength)) sigma_remaining++;
        }

        if (sigma.length == 0 || sigma_remaining == 0) {
            //console.log("formula used", formula.name);
            if(formula.name == "coFold"){
                var index_link = xmat.sequence.indexOf("X");
                var index_endlink = index_link + xmat.minLoopLength + 1;
                var seq1 = xmat.sequence.substr(0, index_link);
                var seq2 = xmat.sequence.substr(index_endlink, xmat.sequence.length);
                //console.log(index_link, index_endlink, xmat.sequence, seq1, seq2);
                var bps = [];
                if(index_link != -1) {
                    for (var i in P) {
                        if (P[i][0] <= index_link && P[i][1] >= index_link + 1 + xmat.minLoopLength) {
                            //console.log("push", [P[i][0], P[i][1] - index_endlink]);
                            bps.push([P[i][0], P[i][1] - index_endlink]);
                        }
                    }
                }
                //console.log("P", JSON.stringify(P), JSON.stringify(xmat.conv_str(P, seq_length)), "bps:", bps);
                var temp_sos = {
                    structure: xmat.conv_str(P, seq_length),
                    traces: traces,
                    rep4d: repres(visualize4d(seq1, seq2, bps))
                };
            }
            else {
                var temp_sos = {
                    structure: xmat.conv_str(P, seq_length),
                    traces: traces
                };
            }
            SOS.push(temp_sos);
        } else {
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

    }
    return SOS;
};


/** Nussinov Structures Count*/

var NussinovDPAlgorithm_structuresCount = Object.create(DPAlgorithm);

NussinovDPAlgorithm_structuresCount.Description = "Nussinov counting";

NussinovDPAlgorithm_structuresCount.Tables = new Array();
NussinovDPAlgorithm_structuresCount.Tables.push(Object.create(NussinovMatrix));

NussinovDPAlgorithm_structuresCount.Tables[0].latex_representation = "C_{i,j} = C_{i,j-1} + \\sum_{i\\leq k <(j-l) \\atop S_k,S_j \\text{ pair}} C_{i,k-1} \\cdot C_{k+1,j-1}";

// C(i, j) = C(i, j - 1) + sum[k: [i <= k < j - l] && k,j pairs] C(i, k - 1) * C(k + 1, j - 1)
NussinovDPAlgorithm_structuresCount.Tables[0].computeCell = function (i, j) {

    var curCell = Object.create(NussinovCell).init(i, j, 0);

    if (this.isInvalidState(i, j)) {
        return curCell;
    }
    if (i >= j) {
        curCell.value = 1;
        return curCell;
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

    curCell.value = res;
    return curCell;
};


// Invoking getValue for all tables, so that the Dynamic Programming computes all the tables and memoizes them.
NussinovDPAlgorithm_structuresCount.computeMatrix = function (input) {
// resize and initialize matrix
console.log("in struct count");
    this.Tables[0].init(input.sequence(), "structuresCount");
    // store minimal loop length
    this.Tables[0].minLoopLength = parseInt(input.loopLength());

    this.Tables[0].computeAllCells();
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

NussinovDPAlgorithm_McCaskill.Tables[0].latex_representation = "Q_{i,j} = Q_{i,j-1} + \\sum_{i\\leq k <(j-l)} Q_{i,k-1} \\cdot Q^{bp}_{k,j}";
NussinovDPAlgorithm_McCaskill.Tables[1].latex_representation = "Q_{i,j}^{bp} = \\begin{cases} Q_{i + 1, j - 1} \\cdot \\exp(-E_{bp}/RT) & \\text{ if }S_{i},S_{j} \\text{ can form base pair} \\\\ 0 & \\text{ otherwise}\\end{cases}";

NussinovDPAlgorithm_McCaskill.Tables[2].latex_representation = "P^{bp}_{i,j} = \\frac{Q_{1,i-1} \\cdot Q^{bp}_{i,j} \\cdot Q_{j+1,n}}{Q_{1,n}} + \\sum_{p<i,j<q} P^{bp}_{p,q} \\cdot \\frac{\\exp(-E_{bp} / RT) \\cdot Q_{p+1,i-1} \\cdot Q^{bp}_{i,j} \\cdot Q_{j+1,q-1}}{Q^{bp}_{p,q}}";
NussinovDPAlgorithm_McCaskill.Tables[3].latex_representation = "P^{u}_{i,j} = \\frac{Q_{1,i-1} \\cdot 1 \\cdot Q_{j+1,n}}{Q_{1,n}} + \\sum_{p<i,j<q} P^{bp}_{p,q} \\cdot \\frac{\\exp(-E_{bp} / RT) \\cdot Q_{p+1,i-1} \\cdot 1 \\cdot Q_{j+1,q-1}}{Q^{bp}_{p,q}}";

//NussinovDPAlgorithm_McCaskill.Tables[2].latex_representation = "P^{bp}_{i, j} = Q^{-1}_{1, n} \\cdot (Q_{1, i - 1} \\cdot Q^{b}_{i, j} \\cdot Q_{j + 1, n})";
//NussinovDPAlgorithm_McCaskill.Tables[3].latex_representation = "P^{u}_{i, j} = Q^{-1}_{1, n} \\cdot (Q_{1, i - 1} \\cdot 1 \\cdot Q_{j + 1, n})";

// Q(i,j) = sum[k : [i <= j < j - l] && k,j can pair] Q(i, k - 1) * Qb(k, j)

// TODO: Check the base case in algorithms;
NussinovDPAlgorithm_McCaskill.Tables[0].computeCell = function (i, j) {

    var curCell = Object.create(NussinovCell).init(i, j, 0);
    if (i > j + 1 || i < 0 || j < 0 || i > this.seq_length + 1 || j > this.seq_length + 1) {
        return curCell;
    }
    if (i >= j) {
        curCell.value = 1;
        return curCell;
    }
    var res = 0;
    res += this.getValue(i, j - 1);
    for (var k = i; k < j - this.minLoopLength; ++k) {
        if (RnaUtil.areComplementary(this.sequence[k - 1], this.sequence[j - 1])) {
            res += this.getValue(i, k - 1) * NussinovDPAlgorithm_McCaskill.Tables[1].getValue(k, j);
        }
    }
    curCell.value = res;
    return curCell;
};

// Qb(i, j) = Q(i + 1, j - 1) * exp(-Eb / RT)
NussinovDPAlgorithm_McCaskill.Tables[1].computeCell = function (i, j) {

    var curCell = Object.create(NussinovCell).init(i, j, 0);
    if (i < 0 || j < 0 || i > this.seq_length || j > this.seq_length || i >= j - this.minLoopLength) {
        return curCell;
    }
    if (RnaUtil.areComplementary(this.sequence[i - 1], this.sequence[j - 1])) {
        curCell.value = NussinovDPAlgorithm_McCaskill.Tables[0].getValue(i + 1, j - 1) * Math.exp(this.energy_basepair);
    }
    return curCell;
};

// Probability that i, j is a base pair.
// Pe(i, j) = Q(1, i - 1) * Qb(i, j) * Q(j + 1, n) / Q(1, n)
//            + sum[p,q]{Pe(p,q) * Qb(i,j) * Q(p+1,i-1)*Q(j+1,q-1)*exp(-e/RT)/Qb(p,q)}
NussinovDPAlgorithm_McCaskill.Tables[2].computeCell = function(i, j) {

    var curCell = Object.create(NussinovCell).init(i, j, 0);
    if (i < 0 || j < 0) {//} || i > this.seq_length || j > this.seq_length) {
        return curCell;
    }
    if (!RnaUtil.areComplementary(this.sequence[i - 1], this.sequence[j - 1]) || i >= j - this.minLoopLength) {
        return curCell;
    }
    var n = this.seq_length;

    var ret = NussinovDPAlgorithm_McCaskill.Tables[1].getValue(i, j);
    if (i > 1) {
        ret *= NussinovDPAlgorithm_McCaskill.Tables[0].getValue(1, i - 1);
    }
    if (j < n) {
        ret *= NussinovDPAlgorithm_McCaskill.Tables[0].getValue(j + 1, n);
    }
    ret /= NussinovDPAlgorithm_McCaskill.Tables[0].getValue(1, n);

    for (var p = 1; p < i; ++p) {
        for (var q = j + 1; q <= n; ++q) {
            if (RnaUtil.areComplementary(this.sequence[p - 1], this.sequence[q - 1])) {
                var v = this.getValue(p, q) * NussinovDPAlgorithm_McCaskill.Tables[1].getValue(i, j);
                if (p + 1 < i) {
                    v *= NussinovDPAlgorithm_McCaskill.Tables[0].getValue(p + 1, i - 1);
                }
                if (j + 1 < q) {
                    v *= NussinovDPAlgorithm_McCaskill.Tables[0].getValue(j + 1, q - 1);
                }
                // weight of base pair p, q
                v *= Math.exp(this.energy_basepair);
                v /= NussinovDPAlgorithm_McCaskill.Tables[1].getValue(p, q);
                ret += v;
            }
        }
    }
    curCell.value = ret;
    curCell.logValue = -this.energy_normal * Math.log(ret);
    //console.log(curCell);
    return curCell;
};

// Probability that i, j is unpaired.
// Pu(i, j) = Q(1, i - 1) * 1 * Q(j + 1, n) / Q(1, n) + sum[p,q]{Pbn[p,q] * Q[p+1,i-1]*Q[j+1,q-1]*exp(-e)/Qb[p,q]}
NussinovDPAlgorithm_McCaskill.Tables[3].computeCell = function(i, j) {

    var curCell = Object.create(NussinovCell).init(i, j, 0);
    if (i < 0 || j < 0 || i > j) {//} || i > this.seq_length || j > this.seq_length || i > j) {
        return curCell;
    }
    var n = this.seq_length;

    var ret = 1.0;
    if (i > 1) {
        ret *= NussinovDPAlgorithm_McCaskill.Tables[0].getValue(1, i - 1);
    }
    if (j < n) {
        ret *= NussinovDPAlgorithm_McCaskill.Tables[0].getValue(j + 1, n);
    }
    ret /= NussinovDPAlgorithm_McCaskill.Tables[0].getValue(1, n);

    for (var p = 1; p < i; ++p) {
        for (var q = j + 1; q <= n; ++q) {
            if (p < q - this.minLoopLength &&
                RnaUtil.areComplementary(this.sequence[p - 1], this.sequence[q - 1])) {
                var v = NussinovDPAlgorithm_McCaskill.Tables[2].getValue(p, q);
                if (p + 1 < i) {
                    v *= NussinovDPAlgorithm_McCaskill.Tables[0].getValue(p + 1, i - 1);
                }
                if (j + 1 < q) {
                    v *= NussinovDPAlgorithm_McCaskill.Tables[0].getValue(j + 1, q - 1);
                }
                // weight of base pair p,q
                v *= Math.exp(this.energy_basepair);
                v /= NussinovDPAlgorithm_McCaskill.Tables[1].getValue(p, q);

                ret += v;
            }
        }
    }
    curCell.value = ret;
    curCell.logValue = -this.energy_normal * Math.log(ret);
    //console.log(curCell);
    return curCell;
};

NussinovDPAlgorithm_McCaskill.computeMatrix = function (input) {
    this.Tables[0].init(input.sequence(), "McCaskill");
    this.Tables[1].init(input.sequence(), "McCaskill Base");
    this.Tables[2].init(input.sequence(), "McCaskill External Base");
    this.Tables[3].init(input.sequence(), "McCaskill Unpaired");
    // store minimal loop length
    this.Tables[0].minLoopLength = parseInt(input.loopLength());
    this.Tables[1].minLoopLength = parseInt(input.loopLength());
    this.Tables[2].minLoopLength = parseInt(input.loopLength());
    this.Tables[3].minLoopLength = parseInt(input.loopLength());

    //this.Tables[1].energy_basepair = -input.energy() / input.energy_normal();
    //this.Tables[2].energy_basepair = -input.energy() / input.energy_normal();
    //this.Tables[3].energy_basepair = -input.energy() / input.energy_normal();

    for (var i = 1; i < 4; ++i) {
        this.Tables[i].energy_basepair = -input.energy() / input.energy_normal();
        this.Tables[i].energy_normal = input.energy_normal();
    }

    for (var i = 0; i < 4; ++i) {
        this.Tables[i].computeAllCells();
    }
    
    return this.Tables;
};


var DPAlgorithm_MEA = Object.create(DPAlgorithm);

DPAlgorithm_MEA.Description = "Maximum Expected Accuracy";
DPAlgorithm_MEA.Tables = new Array();
DPAlgorithm_MEA.Tables.push(Object.create(NussinovMatrix));
DPAlgorithm_MEA.Tables.push(Object.create(NussinovMatrix));
DPAlgorithm_MEA.Tables.push(Object.create(NussinovMatrix));
DPAlgorithm_MEA.Tables[0].latex_representation = "M_{i, j} = \\max \\begin{cases}M_{i, j - 1} + P^{u}_{j} & S_{j}\\text{ unpaired} \\\\ \\max_{i \\leq k < (j-l)} \\left( M_{i, k - 1} + M_{k + 1, j - 1} + \\gamma \\cdot P^{bp}_{k, j} \\right) & S_{k}\\text{ paired with }S_{j} \\end{cases}";
DPAlgorithm_MEA.Tables[2].latex_representation = "P_{i}^{u} = 1 - \\sum_{k < i}{P^{bp}_{k, i}} - \\sum_{i < j}{P^{bp}_{i, j}}";

//DPAlgorithm_MEA.Tables[0].latex_representation = "M_{i, j} = \\max \\begin{cases} 0 & i > j \\\\ M_{i, j - 1} + p^{u}_{j} & j unpaired \\\\ M_{i + 1, j - 1} + p^{p}_{i,j} & j paired with i \\\\ \\max_{i \\leq k < j}{M_{i, k} + M_{k + 1, j}} & decomposition \\end{cases}";
//DPAlgorithm_MEA.Tables[0].latex_representation = "D(i,j) = \\max \\begin{cases} D(i+1,j) & S_i \\text{ unpaired} \\\\ D(i,j-1) & S_j \\text{ unpaired} \\\\ D(i+1,j-1)+1 & \\text{if } S_i,S_j \\text{ compl. base pair and } i+ l< j \\\\ \\max_{i< k< (j-1)} D(i,k)+D(k+1,j) & \\text{ decomposition} \\end{cases}";
DPAlgorithm_MEA.Tables[0].updateCell = function (curCell, curVal, curAncestor) {

    if (curCell === null || curCell.value <= curVal) {
        // check for new maximal value
        if (curCell === null || curCell.value < curVal) {
            // reset ancestor list
            curCell.traces = [];
            // store new maximum
            curCell.value = curVal;
        }
        // store this ancestor
        curCell.traces.push(curAncestor);
    }
}

DPAlgorithm_MEA.Tables[0].computeCell = function(i, j) {

    var curCell = Object.create(NussinovCell).init(i, j, 0);
    if (this.isInvalidState(i, j) || i > j) {
        return curCell;
    }

    this.updateCell(curCell, this.getValue(i, j - 1) + DPAlgorithm_MEA.Tables[2].getValue(j, j), Object.create(NussinovCellTrace).init([[i, j - 1]], []));
    for (var k = i; k < j - this.minLoopLength; ++k) {
        if (RnaUtil.areComplementary(this.sequence[k - 1], this.sequence[j - 1])) {
            this.updateCell(curCell, this.getValue(i, k - 1) + this.getValue(k + 1, j - 1) + this.gamma * NussinovDPAlgorithm_McCaskill.Tables[2].getValue(k, j), Object.create(NussinovCellTrace).init([[i, k - 1], [k + 1, j - 1]], [k, j]));
        }
    }
    return curCell;
};

DPAlgorithm_MEA.Tables[2].computeCell = function(i, j) {
    var curCell = Object.create(NussinovCell).init(i, j, 0);
    if (i != j || i < 0 || j < 0 || i > this.seq_length || j > this.seq_length) {
        return curCell;
    }
    var ret = 1.0;
    for (var k = 1; k < i; ++k)
        ret -= NussinovDPAlgorithm_McCaskill.Tables[2].getValue(k, i);
    for (var k = j + 1; k <= this.seq_length; ++k)
        ret -= NussinovDPAlgorithm_McCaskill.Tables[2].getValue(j, k);

    curCell.value = ret;
    return curCell;
};

DPAlgorithm_MEA.computeMatrix = function(input) {
    
    NussinovDPAlgorithm_McCaskill.computeMatrix(input);

    this.Tables[1] = NussinovDPAlgorithm_McCaskill.Tables[2];

    this.Tables[0].init(input.sequence(), "Maximum Expected Accuracy");
    this.Tables[2].init(input.sequence(), "Unpaired probabilities");
    // store minimal loop length
    this.Tables[0].minLoopLength = parseInt(input.loopLength());
    this.Tables[2].minLoopLength = parseInt(input.loopLength());
    this.Tables[0].gamma = parseFloat(input.gamma());

    this.Tables[0].computeAllCells();
    for (var i = 0; i <= this.Tables[2].seq_length; ++i) {
        this.Tables[2].getValue(i, i);
    }

    return this.Tables;
};


DPAlgorithm_MEA.Tables[0].getSubstructures = function (sigma, P, traces, delta, maxLengthR) {
    var R = [];
    var ij = sigma.pop();

    if (this.isInvalidState(ij[0], ij[1]) || ij[0] > ij[1]) {//ij[0] > ij[1) {
        var S_prime = {};
        S_prime.sigma = sigma;
        S_prime.P = P;
        S_prime.traces = traces;
        S_prime.potential = delta;
        R.push(S_prime);
        return R;
    }

    var Nmax = this.getValue(ij[0], ij[1]);
    {
        var NSprime = this.getValue(ij[0], ij[1] - 1) + DPAlgorithm_MEA.Tables[2].getValue(ij[1], ij[1]);
        if (NSprime >= Nmax - delta || Math.abs(NSprime - (Nmax - delta)) < EPS) {
            var trace = Object.create(NussinovCellTrace).init([[ij[0], ij[1] - 1]], []);
            R.unshift(this.getSPrime(NSprime, Nmax, ij, trace, sigma, delta, P, traces));
        }
        if (R.length >= maxLengthR) {
            return R;
        }
    }

    // if (i,j) == (i,k -1) + (k+1, j - 1) || (k,j)
    for (var k = ij[0]; k < ij[1] - this.minLoopLength; ++k) {
        if (RnaUtil.areComplementary(this.sequence[k - 1], this.sequence[ij[1] - 1])) {
            var NSprime = this.getValue(ij[0], k - 1) + this.getValue(k + 1, ij[1] - 1) + this.gamma * NussinovDPAlgorithm_McCaskill.Tables[2].getValue(k, ij[1]);//this.countBasepairs(tmp_P, sigma_prime);
            if (NSprime >= Nmax - delta || Math.abs(NSprime - (Nmax - delta)) < EPS) {
                var trace = Object.create(NussinovCellTrace).init([[ij[0], k - 1], [k + 1, ij[1] - 1]], [[k, ij[1]]]);
                R.push(this.getSPrime(NSprime, Nmax, ij, trace, sigma, delta, P, traces));
            }
            if (R.length >= maxLengthR) {
                return R;
            }
        }
    }

    return R;
};


var DPAlgorithm_coFold = Object.create(DPAlgorithm);

DPAlgorithm_coFold.Description = "Co-Fold";
DPAlgorithm_coFold.Tables = new Array();
DPAlgorithm_coFold.Tables.push(Object.create(NussinovMatrix));

DPAlgorithm_coFold.Tables[0].latex_representation = "D(i,j) = \\max \\begin{cases} D(i,j-1) & S_j \\text{ unpaired} \\\\ \\max_{i\\leq k< (j-l)} \\left( D(i,k-1)+D(k+1,j-1)+1 \\right) & \\text{if } S_k,S_j \\text{ compl. base pair} \\end{cases}";

DPAlgorithm_coFold.computeMatrix = function(input) {

    NussinovDPAlgorithm_Unique.computeMatrix(input);
    this.Tables[0].minLoopLength = parseInt(input.loopLength());

    this.Tables[0].init(input.sequence(), "coFold");
    this.Tables[0].cells = NussinovDPAlgorithm_Unique.Tables[0].cells;
    return this.Tables;
};

DPAlgorithm_coFold.Tables[0].getSubstructures = function (sigma, P, traces, delta, maxLengthR) {
    return NussinovDPAlgorithm_Unique.Tables[0].getSubstructures(sigma, P, traces, delta, maxLengthR);
};

/**
 * WUCHTY (generic doesnt give subomtimal structures)
 */
var wuchty = function (xmat) {
    var seq_length = xmat.sequence.length;
    var Nmax = xmat.getValue(1, seq_length);

    var S = {sigma: [[1, seq_length]], P: [], traces: []};
    var R = [S];
    var SOS = [];
    var loop = 0;

    while (R.length != 0) {

        // Pop R
        var pop_R = R.pop();
        var sigma = pop_R.sigma;
        var P = pop_R.P;
        var t_traces = JSON.stringify(pop_R.traces);
        var traces = JSON.parse(t_traces);

        var sigma_remaining = 0;
        for (var s in sigma) {//console.log("var s:", sigma[s]);
            if ((sigma[s][0]) <= (sigma[s][1] - xmat.minLoopLength)) sigma_remaining++;
        }

        if (sigma.length == 0 || sigma_remaining == 0) {
            var temp_sos = {structure: xmat.conv_str(P, seq_length), traces: traces};
            SOS.push(temp_sos);
            if (SOS.length > 10) {
                break;
            }
        }

        else {
            var ij = sigma.pop();

            if (ij[0] > ij[1]) {
                pop_R.traces = traces;
                R.push(pop_R);
                continue;
            }

            var json_ij_traces = JSON.stringify(xmat.getCell(ij[0], ij[1]).traces);
            var ij_traces = JSON.parse(json_ij_traces);

            var S_prime = {};

            for (var t in ij_traces) {
                var S_prime = {};
                var t_trace = JSON.stringify(ij_traces[t]);
                var trace = JSON.parse(t_trace);
                var trace_p = JSON.parse(t_trace);

                // add sigma info in S_prime
                S_prime.sigma = [];
                if (sigma.length > 0) {
                    for (var s in sigma)
                        trace.parents.unshift(sigma[s]);
                    S_prime.sigma = trace.parents;
                }
                else S_prime.sigma = trace.parents;

                // add P(bps) info in S_prime
                S_prime.P = []
                if (P[0] != undefined) {
                    for (var p in P) {
                        S_prime.P.push(P[p]);
                    }
                }
                if (trace.bps[0] != undefined) {
                    S_prime.P.push(trace.bps[0]);
                }

                // add traces info in S_prime
                var temp_trace = [ij];
                temp_trace.push(trace_p.parents);
                if (traces.length == 0) {
                    S_prime.traces = [temp_trace];
                }
                else {
                    var clone_traces = JSON.stringify(traces);
                    var parse_clone_traces = JSON.parse(clone_traces);
                    parse_clone_traces.unshift(temp_trace);
                    S_prime.traces = parse_clone_traces;
                }
                R.push(S_prime);
            }
        }
        ;

    }
    return SOS;
}

///**
// * global list of available Nussinov algorithm implementations
// */
//var availableAlgorithms = {
//
//    /** original unique recursion */
//    nussinovUnique: NussinovDPAlgorithm_Unique,//NussinovMatrix_unique,
//
//    /** ambiguous recursion */
//    nussinovAmbiguous1: NussinovDPAlgorithm_Ambiguous,//NussinovMatrix_ambiguous,
//
//    /** nussinov neo recursion */
//    nussinovAmbiguous2: NussinovDPAlgorithm_Ambiguous2,
//
//    /** McCaskill */
//    mcCaskill: NussinovDPAlgorithm_McCaskill,
//
//    /** structure counting */
//    nussinovCounting: NussinovDPAlgorithm_structuresCount,
//
//    /** Maximum Expected Accuracy*/
//    MaxExpAcc: DPAlgorithm_MEA,
//
//    /** Co-fold*/
//    coFold: DPAlgorithm_coFold,
//
//    hybrid: DPAlgorithm_hybrid,
//
//    rnaup: DPAlgorithm_rnaup
//
//};
