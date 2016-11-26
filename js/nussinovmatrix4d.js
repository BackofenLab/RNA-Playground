/**
 * Created by moni on 02.08.16.
 */
/**
 * Utility class that covers RNA specific functions.
 *
 */
/*
var RnaUtil = {

    /**
     * checks base pair complementary of two nucleotides
     * @param {string} nt1 first nucleotide
     * @param {string} nt2 second nucleotide
     * @returns {boolean} true if complementary; false otherwise
     * /
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
     * /
    isRnaSequence: function (sequence) {
        var isValid =
            // check if sequence given
            (sequence !== null)
                // check RNA alphabet
            && sequence.match("^[ACGU]+$");
        return isValid;
    }


};
*/


var NussinovCell4d = {

// row
    i: -1,

    k: -1,

// column
    j: -1,

    l: -1,

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
    init: function (i, k, j, l, value) {
        // init data
        this.i = i;
        this.k = k;
        this.j = j;
        this.l = l;
        this.value = value;
        this.traces = [];
        // this access for chaining
        return this;
    }

};



/**
 * Ancestor information for a certain traceback
 */
var NussinovCell4dTrace = {

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

var NussinovMatrix4d = {

        /**
         * Access to the sequence for this matrix
         */
        sequence1: null,
        /**
         * Access to the sequence for this matrix
         */
        sequence2: null,
    
        seq1_length: 0,
    
        seq2_length: 0,

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

        tablesDimension: 4,

        /**
         * initialize a matrix of dim = n+1 with indices 0..n, where n is the
         * length of the provided sequence
         * @param {string} sequence the RNA sequence (not null or empty)
         * @returns {NussinovMatrix} this
         */
        init: function (sequence1, sequence2, name) { //initialize matrix

            // reset data
            this.sequence1 = null;
            this.sequence2 = null;
            this.name = null;
            this.cells = [];

            // check input
            if (sequence1 == null || sequence1 === "" || sequence2 == null || sequence2 === ""  || name == null) {
                console.log("Matrix init failed for sequence (", sequence1, sequence2, ")");
                return this;
            }

            // store sequence
            this.sequence1 = sequence1;
            this.sequence2 = sequence2;
            this.name = name;
            this.seq1_length = this.sequence1.length;
            this.seq2_length = this.sequence2.length;

            // create matrix cells

            for (var i = 0; i <= this.seq1_length; i++) {
                this.cells[i] = [];
                for (var k = 0; k <= this.seq1_length; ++k) {
                    this.cells[i][k] = [];
                    for (var j = 0; j <= this.seq2_length; j++) {
                        // create new cell and initialize
                        this.cells[i][k][j] = [];
                        for (var r = 0; r <= this.seq2_length; ++r) {
                            //console.log(i, k, j, r, "init");
                            this.cells[i][k][j][r] = Object.create(NussinovCell4d).init(i, k, j, r, null);
                        }
                    }
                }
                ;
            }
            ;

            return this;
        },

        isBaseCase: function(i, k, j, l) {
            if (i < 0 || j < 0 || k < 0 || l < 0 || k < i || l < j ||
                k > this.sequence1.length || l > this.sequence2.length) {
                return true;
            } else {
                return false;
            }
        },

        computeCell: function(i, k, j, l) {
            // updateCell(i, j);
        },
        /**
         * access whole object at cell location (i,j) in matrix
         * @param {int} i row #.
         * @param {int} j column #.
         * @returns {NussinovCell} cell object or null if not available
         */
        getCell: function (i, k, j, l) {
            // check border cases
            if (this.isBaseCase(i, k, j, l)) {
                return null;
            }
            //console.log(i, k, j, l, this.cells[i][k][j][l]);
            if (this.cells[i][k][j][l] === null || this.cells[i][k][j][l].value === null) {
                this.cells[i][k][j][l] = this.computeCell(i, k, j, l);
            }
            return this.cells[i][k][j][l];
        },

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

        /*
        computeValue: function (i, k, j, l) {
            // Computes the M[i, j] by accessing cells using getValue function for memoization.

            // base case
            return 0;

            // computation
        },
        */
        /**
         * access value of cell at location (i,j) in matrix, and compute it if it's not computed.
         * @param {int} i row #.
         * @param {int} j column #.
         * @returns {int} cell value or 0 if invalid cell
         */
        getValue: function (i, k, j, l) {
            // access cell at location (i,j) in the matrix
            //console.log("get", i, k, j, l);
            var cell = this.getCell(i, k, j, l);
            if (cell === null) {
                return null;
            }

            //console.log("getting", i, k, j, l);
            // check if invalid cell
            ///if (cell.value === null) {
            //    cell.value = this.computeValue(i, k, j, l);
            //}
            // get cell value
            return parseFloat(cell.value);
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
        updateCell: function (curCell, curAncestor) {
            // get cell to update
            // check if something to update
            if (curCell === null) {
                return;
            }

            // init value with number of additional base pairs
            var curVal = curAncestor.bps.length;
            // add scores of ancestor cells
            for (var x = 0; x < curAncestor.parents.length; x++) {
                var i = curAncestor.parents[x][0];
                var k = curAncestor.parents[x][1];
                var j = curAncestor.parents[x][2];
                var l = curAncestor.parents[x][3];
                curVal += this.getValue(i, k, j, l);
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
         */
        computeAllCells: function() {
            for (var i = 0; i <= this.seq1_length; ++i) {
                for (var k = 0; k<= this.seq1_length; ++k) {
                    for (var j = 0; j <= this.seq2_length; ++j) {
                        for (var l = 0; l <= this.seq2_length; ++l) {
                            //console.log(this.cells[i][k][j][l]);
                            //console.log(this.getCell(i, k, j, l));
                        }
                    }
                }
            }
        },

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

            var splitSeq = input.sequence().indexOf('X');
            var sequence1 = input.sequence().substr(0,splitSeq);
            var sequence2 = input.sequence().substr(parseInt(input.loopLength())+splitSeq + 1);
            // resize and init matrix
            this.init(sequence1, sequence2, "Default name");

            // set minimal loop length
            this.minLoopLength = parseInt(input.loopLength());

            this.computeAllCells();

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
            console.log("WARNING: getRecursionInLatex() not implemented in NussinovMatrix4d superclass; overwrite in subclass!");
            return "$$" + this.latex_representation + "$$";

        }
        ,

        conv_str: function(P) {
            var str = "";
            for (var l = 0; l < this.seq1_length + this.seq2_length + 1; l++) {
                str += ".";
            }
            str = str.substr(0, this.seq1_length) + "X" + str.substr(this.seq1_length + 1);
            //str[this.seq1_length] = 'X';
            for (var indx in P) {
                var i = P[indx][0], j = P[indx][1];
                //str[i - 1] = '(';
                //str[str.length - j] = ')';
                str = str.substr(0, i - 1) + "(" + str.substr(i);
                str = str.substr(0, str.length - j) + ")" + str.substr(str.length - j + 1);
            }

            return str;
        },
    
        simpleRepresentation: function() {

            var n1 = this.seq1_length;
            var n2 = this.seq2_length;

            console.log(n1, n2);
            console.log(this.sequence1, this.sequence2);
            var res = "";
            for (var i = 0; i <= n1; ++i)
                for (var k = i; k <= n1; ++k)
                    for (var j = 0; j <= n2; ++j)
                        for (var l = j; l <= n2; ++l)
                            if (this.getValue(i, k, j, l) != null && Math.abs(this.getValue(i, k, j, l)) > 1e-5 ) {
                                res += JSON.stringify(this.getCell(i, k, j, l)) + "\n";
                            }
            
            return res; 

        }

    }
    ;
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




var DPAlgorithm_hybrid = Object.create(DPAlgorithm);

DPAlgorithm_hybrid.Description = "RNA to RNA matching";
DPAlgorithm_hybrid.Tables = new Array();
DPAlgorithm_hybrid.Tables.push(Object.create(NussinovMatrix4d));
DPAlgorithm_hybrid.Tables[0].latex_representation = "D_{i, k}^{j, l} = \\max \\begin{cases} E^{init}(i, j) & \\mathcal{R}^1_i, \\mathcal{R}^2_j  \\text{  pairs}, i = k, j = l \\\\ \\max_{p,q}{ E^{loop}(i, j, p, q) + D_{q, l}^{p, k} } & \\mathcal{R}^1_i, \\mathcal{R}^2_j  \\text{  pairs}, i < k, j < l\\\\ 0 & \\text{otherwise} \\end{cases}";

DPAlgorithm_hybrid.Tables[0].computeCell  = function(i, k, j, l) {

    var curCell = Object.create(NussinovCell4d).init(i, k, j, l, 0);

    // TODO: Check i,j <= 0 or i,j < 0
    if (i <= 0 || j <= 0 || this.isBaseCase(i, k, j, l)) {
    //if (this.isBaseCase(i, k, j, l)) {
        return curCell;
    }
    var ret = 0;
    //console.log(i, k, j, l);
    // TODO: Check the indices of this condition
    if (RnaUtil.areComplementary(this.sequence1[i - 1], this.sequence2[j - 1])) {
        if (i === k && j === l) {
            // Energy init instead of 1
            ret = Math.max(ret, 1);

            this.updateCell(curCell, Object.create(NussinovCell4dTrace).init([], [[i, j]]));
        }

        for (var p = i + 1; p <= k; ++p) {
            for (var q = j + 1; q <= l; ++q) {
                // Energy loop instead of 1
                // if (i < k && j < l)
                if (this.getValue(p, k, q, l) > 0) {
                    ret = Math.max(ret, 1 + this.getValue(p, k, q, l));
                    this.updateCell(curCell, Object.create(NussinovCell4dTrace).init([[p, k, q, l]], [[i, j]]));
                }
            }
        }
    }

    //console.log(i, k, j, l, ret, curCell.value);
    //assert curCell.value = ret;
    return curCell;
};



DPAlgorithm_hybrid.computeMatrix = function(input) {
    var splitSeq = input.sequence().indexOf('X');
    var sequence1 = input.sequence().substr(0,splitSeq);
    var sequence2 = input.sequence().substr(parseInt(input.loopLength())+splitSeq + 1).split("");//.reverse().join("");

    console.log(sequence1, sequence2);
    this.Tables[0].init(sequence1, sequence2, "RNAHybrid");
    //console.log(this.Tables[0].cells);

    this.Tables[0].computeAllCells();

    console.log(this.Tables[0].simpleRepresentation());
    return this.Tables;

};




var DPAlgorithm_rnaup = Object.create(DPAlgorithm);

DPAlgorithm_rnaup.Description = "RNA to RNA matching";
DPAlgorithm_rnaup.Tables = new Array();
DPAlgorithm_rnaup.Tables.push(Object.create(NussinovMatrix4d));
DPAlgorithm_rnaup.Tables.push(Object.create(NussinovMatrix));
DPAlgorithm_rnaup.Tables.push(Object.create(NussinovMatrix));

DPAlgorithm_rnaup.Tables[0].latex_representation = "I_{i, k}^{j, l} = E_{bp} \\cdot D_{i, k}^{j, l} - RT \\cdot (\\log(P^{u_1}_{i,k}) + \\log(P^{u_2}_{j, l}))";//"D_{i, k}^{j, l} = \\max \\begin{cases} E^{init}(i, j) & \\mathcal{R}^1_i, \\mathcal{R}^2_j  pairs, i = k, j = l \\\\ \\max_{p,q}{ E^{loop}(i, j, p, q) + D_{q, l}^{p, k} } & \\mathcal{R}^1_i, \\mathcal{R}^2_j  pairs, i < k, j < l\\\\ 0 & otherwise \\end{cases}";



DPAlgorithm_rnaup.Tables[0].computeCell = function(i, k, j, l) {

    var curCell = Object.create(NussinovCell4d).init(i, k, j, l, 0);

    if (this.isBaseCase(i, k, j, l)) {
        return curCell;
    }
    // I[i1,j1,i2,j2] = Ebp*H[i1,j1,i2,j2] -RT*ln(P^u_1) -RT*ln(P^u_2)
    if (DPAlgorithm_hybrid.Tables[0].getValue(i, k, j, l) == 0) {
        return curCell;
    }
    var logP = Math.log(DPAlgorithm_rnaup.Tables[1].getValue(i, k)) + Math.log(DPAlgorithm_rnaup.Tables[2].getValue(j, l));

    curCell.value = this.energy * DPAlgorithm_hybrid.Tables[0].getValue(i, k, j, l) - this.energy_normal * logP;
    curCell.traces = DPAlgorithm_hybrid.Tables[0].getCell(i, k, j, l).traces;
    return curCell;
};



DPAlgorithm_rnaup.computeMatrix = function(input) {

    DPAlgorithm_hybrid.computeMatrix(input);

    var splitSeq = input.sequence().indexOf('X');
    var sequence1 = input.sequence().substr(0,splitSeq);
    var sequence2 = input.sequence().substr(parseInt(input.loopLength())+splitSeq + 1);

    NussinovDPAlgorithm_McCaskill.computeMatrix({sequence: function(){return sequence1;}, loopLength: input.loopLength, energy: input.energy, energy_normal: input.energy_normal});
    this.Tables[1] = JSON.parse(JSON.stringify(NussinovDPAlgorithm_McCaskill.Tables[3]));
    this.Tables[1].getRecursionInLatex = NussinovDPAlgorithm_McCaskill.Tables[3].getRecursionInLatex;
    this.Tables[1].getDim = NussinovDPAlgorithm_McCaskill.Tables[3].getDim;
    this.Tables[1].getCell = NussinovDPAlgorithm_McCaskill.Tables[3].getCell;
    this.Tables[1].getValue = NussinovDPAlgorithm_McCaskill.Tables[3].getValue;
    this.Tables[1].isBaseCase = NussinovDPAlgorithm_McCaskill.Tables[3].isBaseCase;


    NussinovDPAlgorithm_McCaskill.computeMatrix({sequence: function(){return sequence2;}, loopLength: input.loopLength, energy: input.energy, energy_normal: input.energy_normal});
    this.Tables[2] = JSON.parse(JSON.stringify(NussinovDPAlgorithm_McCaskill.Tables[3]));
    this.Tables[2].getRecursionInLatex = NussinovDPAlgorithm_McCaskill.Tables[3].getRecursionInLatex;
    this.Tables[2].getDim = NussinovDPAlgorithm_McCaskill.Tables[3].getDim;
    this.Tables[2].getCell = NussinovDPAlgorithm_McCaskill.Tables[3].getCell;
    this.Tables[2].getValue = NussinovDPAlgorithm_McCaskill.Tables[3].getValue;
    this.Tables[2].isBaseCase = NussinovDPAlgorithm_McCaskill.Tables[3].isBaseCase;


    console.log(sequence1, sequence2);
    this.Tables[0].init(sequence1, sequence2, "RNAHybrid");
    this.Tables[0].energy = input.energy();
    this.Tables[0].energy_normal = input.energy_normal();

    this.Tables[0].computeAllCells();

    console.log(this.Tables[0].simpleRepresentation());
    return this.Tables;

};


/**
 * WUCHTY (generic doesnt give subomtimal structures)
 */
var wuchty4d = function (xmat) {
    console.log("wuchty4d");
    // TODO: compute NMax
    var sigma_0 = [];
    var NMax = 0;
    for (var i = 1; i <= xmat.seq1_length; ++i) {
        for (var k = i; k <= xmat.seq1_length; ++k) {
            for (var j = 1; j <= xmat.seq2_length; ++j) {
                for (var r = j; r <= xmat.seq2_length; ++r) {
                    console.log("val: ", i, k, j, r, xmat.getValue(i, k, j, r));
                    if (xmat.getValue(i, k, j, r) >= NMax) {
                        console.log("true");
                        if (xmat.getValue(i, k, j, r) > NMax) {
                            sigma_0 = [];
                        }
                        sigma_0.push([[i, k, j, r]]);
                        NMax = xmat.getValue(i, k, j, r);
                    }
                }
            }
        }
    }
    console.log("Sigma0: ", JSON.stringify(sigma_0), NMax);
    var SOS = [];
    var loop = 0;
    for (var sig = 0; sig < sigma_0.length; ++sig) {
        var S = {sigma: sigma_0[sig], P: [], traces: []};
        var R = [S];

        console.log("initial: ", R, S);
        while (R.length != 0) {
            console.log("xx");
            // Pop R
            var pop_R = R.pop();
            var sigma = pop_R.sigma;
            var P = pop_R.P;
            var t_traces = JSON.stringify(pop_R.traces);
            var traces = JSON.parse(t_traces);

            var sigma_remaining = 0;
            for (var s in sigma) {//console.log("var s:", sigma[s]);
                if (!xmat.isBaseCase(sigma[s][0], sigma[s][1], sigma[s][2], sigma[s][3])) sigma_remaining++;
            }
            if (sigma.length == 0 || sigma_remaining == 0) {
                //var temp_sos = {structure: xmat.conv_str(P, seq_length), traces: traces};
                // TODO(mohsin): xmat.conv_str(P), pass it scripts.visualize4d
                console.log('visualize4d', visualize4d(xmat.sequence1, xmat.sequence2, P));
                var temp_sos = {structure: xmat.conv_str(P), traces: traces, rep4d: visualize4d(xmat.sequence1, xmat.sequence2, P)};
                console.log("structures:", P);
                console.log("structures parsed:", xmat.conv_str(P));
                console.log('pushing: ', temp_sos);
                SOS.push(temp_sos);

                console.log('SOS: ', SOS);
                if (SOS.length > 10) {
                    break;
                }
            }

            else {
                var idx = sigma.pop();
                //console.log("idx: ", idx);
                if (xmat.isBaseCase(idx[0], idx[1], idx[2], idx[3])) {
                    pop_R.traces = traces;
                    R.push(pop_R);
                    continue;
                }

                var json_ij_traces = JSON.stringify(xmat.getCell(idx[0], idx[1], idx[2], idx[3]).traces);
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
                    S_prime.P = [];
                    if (P[0] != undefined) {
                        for (var p in P) {
                            S_prime.P.push(P[p]);
                        }
                    }
                    console.log("trace", trace.bps);
                    if (trace.bps[0] != undefined) {
                        S_prime.P.push(trace.bps[0]);
                    }

                    // add traces info in S_prime
                    var temp_trace = [idx];
                    temp_trace.push(trace_p.parents);
                    if (traces.length == 0) {
                        console.log('before2 ', S_prime.traces);
                        S_prime.traces = [temp_trace];
                        console.log('after2 ', S_prime.traces);
                    }
                    else {

                        var clone_traces = JSON.stringify(traces);
                        var parse_clone_traces = JSON.parse(clone_traces);
                        console.log('before', parse_clone_traces);
                        parse_clone_traces.unshift(temp_trace);
                        console.log('after', parse_clone_traces);
                        S_prime.traces = parse_clone_traces;
                    }
                    console.log('my ij', JSON.stringify(ij_traces));
                    R.push(S_prime);
                }
            }
            ;

        }
    }
    console.log('final: ', JSON.stringify(SOS));

    return SOS;
}

var availableAlgorithms = {

    /** ambiguous recursion */
    hybrid: DPAlgorithm_hybrid,//NussinovMatrix_ambiguous,

    rnaup: DPAlgorithm_rnaup,

};

