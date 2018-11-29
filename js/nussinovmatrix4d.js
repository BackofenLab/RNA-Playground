/**
 * @file Main file containing backend algorithms for RNA-algorithms-JS project.
 * @import nussinovmatrix.js
 Main Items contains in this file are:
 -4dMatrix class
 -Nussinov 4d algorithms
 -Traceback algorithms
 * @authors, "Mostafa Mahmoud Mohamed", "Syed Mohsin Ali", "Martin Mann"
 */


/**
 * data stored within a cell of a Nussinov matrix 4d
 */
var NussinovCell4d = {
// Start and End index of string1
    i: -1,
    k: -1,

// Start and End index of string2
    j: -1,
    l: -1,

// value
    value: null,

// traces for the current value
    traces: null,

    /**
     * inits a cell with the given data and sets traces to an empty list
     * @param i the start index of string1
     * @param k the ned index of string1
     * @param j the start index of string2
     * @param l the end index of string2
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
            console.log("ERROR : NussinovCell4dTrace.init : only one value null");
            return this;
        }
        // store sane data
        this.parents = parents;
        this.bps = bps;
        return this;
    }
};

/*
 * Nussinov matrix 4d object, stores 2 Sequences and Table. Contains all the utility functionalities for the tables.
 */
var NussinovMatrix4d = {

    /**
     * Access to the sequence1 for this matrix
     */
    sequence1: null,
    /**
     * Access to the sequence2 for this matrix
     */
    sequence2: null,

    /**
     * The length of sequence 1
     */
    seq1_length: 0,

    /**
     * The length of sequence 2
     */
    seq2_length: 0,

    /**
     * Access name of recursion used
     */
    name: null,

    /**
     * cells of the matrix
     */
    cells: [],

    /**
     * The latex representation of the formula computing the matrix.
     */
    latex_representation: "$$",

    /**
     * The dimensions of the table.
     */
    tablesDimension: 4,

    /**
     * Initialize a 4d matrix of dimensions (L1 + 1) x (L1 + 1) x (L2 + 1) x (L2 + 1)
     * @param {string} sequence1 the RNA sequence1 (not null or empty)
     * @param {string} sequence2 the RNA sequence2 (not null or empty)
     * @param {string} name the name description of the table. (not null or empty)
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

        // Store sequences
        this.sequence1 = sequence1;
        this.sequence2 = sequence2;
        this.name = name;
        this.seq1_length = this.sequence1.length;
        this.seq2_length = this.sequence2.length;

        // Create matrix cells
        for (var i = 0; i <= this.seq1_length; i++) {
            this.cells[i] = [];
            for (var k = 0; k <= this.seq1_length; ++k) {
                this.cells[i][k] = [];
                for (var j = 0; j <= this.seq2_length; j++) {
                    this.cells[i][k][j] = [];
                    for (var r = 0; r <= this.seq2_length; ++r) {
                        // Initialize the 4d Cells, with initial values as null.
                        this.cells[i][k][j][r] = Object.create(NussinovCell4d).init(i, k, j, r, null);
                    }
                }
            }
            ;
        }
        ;

        return this;
    },

    /**
     * Check if a given tuple is an invalid state or not.
     * @param i the start index of string1
     * @param k the ned index of string1
     * @param j the start index of string2
     * @param l the end index of string2
     * @returns {boolean}
     */
    isInvalidState: function(i, k, j, l) {
        if (i < 0 || j < 0 || k < 0 || l < 0 || k < i || l < j ||
            k > this.sequence1.length || l > this.sequence2.length) {
            return true;
        } else {
            return false;
        }
    },

    /**
     * Compute the cell at a given state in the 4d matrix.
     * It's recommended to make the implementation use the method "updateCell",
     * if it's computing the tracebacks in an optimization problem.
     * TODO: this function has to be overwritten by the instances, before calling computeMatrix/computeAllCells
     * @param i
     * @param k
     * @param j
     * @param l
     * @returns {NussinovCell4d}  The computed cell.
     */
    computeCell: function(i, k, j, l) {
        // updateCell(i, k, j, l);
        return Object.create(NussinovCell4d).init(i, k, j, r, null);
    },

    /**
     * Access a cell at a given state in the 4d matrix. If the cell is null or has no value,
     * then it will be computed using the "computeCell" method.
     * TODO: Implement computeCell
     * @param i the start index of string1
     * @param k the ned index of string1
     * @param j the start index of string2
     * @param l the end index of string2
     * @returns {NussinovCell4d}  The cell or null if it's an invalid state
     */
    getCell: function (i, k, j, l) {
        // check border cases
        if (this.isInvalidState(i, k, j, l)) {
            return null;
        }
        if (this.cells[i][k][j][l] === null || this.cells[i][k][j][l].value === null) {
            this.cells[i][k][j][l] = this.computeCell(i, k, j, l);
        }
        return this.cells[i][k][j][l];
    },

    /**
     * Access the value at a given state in the 4d matrix. If the cell is null or has no value,
     * then the cell will be computed using the "computeCell" method.
     * TODO: Implement computeCell
     * @param i the start index of string1
     * @param k the ned index of string1
     * @param j the start index of string2
     * @param l the end index of string2
     * @returns {float} Cell value or null if invalid cell
     */
    getValue: function (i, k, j, l) {
        // access cell at location (i,j) in the matrix
        //console.log("get", i, k, j, l);
        var cell = this.getCell(i, k, j, l);
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
     * TODO: Implement computeCell
     */
    computeAllCells: function() {
        for (var i = 0; i <= this.seq1_length; ++i) {
            for (var k = 0; k <= this.seq1_length; ++k) {
                for (var j = 0; j <= this.seq2_length; ++j) {
                    for (var l = 0; l <= this.seq2_length; ++l) {
                        this.getCell(i, k, j, l);
                    }
                }
            }
        }
    },

    /**
     * Fills the matrix according to the recursion.
     * TODO: Implement computeCell
     * TODO: this function has to be overwritten by the instances.
     *
     * @param {input} A Dictionary of the input for the 4dmatrix. Should contain all the arguments
     *                needed for initalizing the matrix input properly. Minimally this will be the 2 sequences.
     * @returns {NussinovMatrix4d} this for call chaining
     */
    computeMatrix: function (input) {
        console.log("WARNING: computeMatrix() not implemented in NussinovMatrix superclass; overwrite in subclass!");

        var splitSeq = input.sequence().indexOf('X');
        var sequence1 = input.sequence().substr(0,splitSeq);
        var sequence2 = input.sequence().substr(parseInt(input.loopLength())+splitSeq + 1);
        // resize and init matrix
        this.init(sequence1, sequence2, "Default name");

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
     * Construct a one string representing the matching base pairs in the 2 sequences.
     * @param P A list of pairs of indices, that represent the matching base pairs between the 2 sequences.
     * @returns {string} A representation of the matching base pairs in the matching.
     */
    conv_str: function(P) {
        var str = "";
        for (var l = 0; l < this.seq1_length + this.seq2_length + this.minLoopLength + 1; l++) {
            str += ".";
        }
        //str[this.seq1_length] = 'X';
        var linked = "";
        for (var i = 0; i <= this.minLoopLength; ++i) {
            linked += "X";
        }
        str = str.substr(0, this.seq1_length) + linked + str.substr(this.seq1_length + this.minLoopLength + 1);
        for (var indx in P) {
            var i = P[indx][0], j = P[indx][1];
            //str[i - 1] = '(';
            str = str.substr(0, i - 1) + "(" + str.substr(i);
            //str[str.length - j] = ')';
            //str[this.seq1_length + 1 + j - 1] = ')';
            str = str.substr(0, this.seq1_length + this.minLoopLength + j) + ")" + str.substr(this.seq1_length + this.minLoopLength + j + 1);
        }
        return str;
    },

    /**
     * @returns {string} A string representation of non-zero/null cells.
     */
    simpleRepresentation: function() {
        var res = JSON.stringify(this.seq1_length) + " " + JSON.stringify(this.seq2_length) + "\n";
        res = JSON.stringify(this.sequence1) + " " + JSON.stringify(this.sequence2) + "\n";
        for (var i = 0; i <= this.seq1_length; ++i)
            for (var k = i; k <= this.seq1_length; ++k)
                for (var j = 0; j <= this.seq2_length; ++j)
                    for (var l = j; l <= this.seq2_length; ++l)
                        if (this.getValue(i, k, j, l) != null && this.getCell(i, k, j, l).traces.length > 0) {
                            res += JSON.stringify(this.getCell(i, k, j, l)) + "\n";
                        }

        return res;
    }
};


var DPAlgorithm_hybrid = Object.create(DPAlgorithm);

DPAlgorithm_hybrid.Description = "hybrid-only interaction prediction";
DPAlgorithm_hybrid.Tables = new Array();
DPAlgorithm_hybrid.Tables.push(Object.create(NussinovMatrix4d));
DPAlgorithm_hybrid.Tables[0].latex_representation = "D^{i, k}_{j, l} = \\max \\begin{cases} 1 & \\text{if } S^1_i, \\overleftarrow{S^2_j}  \\text{ compl.}, i = k, j = l \\\\ \\max_{\\substack{i<p\\leq k,\\;j<q\\leq l\\\\S^1_p, \\overleftarrow{S^2_q}  \\text{ compl.}}}\\left( 1 + D_{q, l}^{p, k} \\right) & \\text{if } S^1_i, \\overleftarrow{S^2_j}  \\text{ compl.}, i < k, j < l\\\\ 0 & \\text{otherwise} \\end{cases}";
//DPAlgorithm_hybrid.Tables[0].latex_representation = "\\begin{array} \\ D^{i, k}_{j, l} = \\max \\begin{cases} E^{init}(i, j) & S^1_i, S^2_j  \\text{  pairs}, i = k, j = l \\\\ \\max_{p,q}{ E^{loop}(i, j, p, q) + D_{q, l}^{p, k} } & S^1_i, S^2_j  \\text{  pairs}, i < k, j < l\\\\ 0 & \\text{otherwise} \\end{cases} \\\\ \\\\ E^{init} = 1 \\\\ \\\\ E^{loop}_{i, j, p, q} =  \\begin{cases} 1 & \\text{if }S^1_p, S^2_q  \\text{  pairs} \\\\ 0 & \\text{otherwise} \\end{cases} \\end{array}";

DPAlgorithm_hybrid.Tables[0].computeCell = function(i, k, j, l) {
    var curCell = Object.create(NussinovCell4d).init(i, k, j, l, 0);

    if (this.isInvalidState(i, k, j, l)) {
        return curCell;
    }

    if (RnaUtil.areComplementary(this.sequence1[i - 1], this.sequence2[this.seq2_length - j])) {
        if (i === k && j === l) {
            // Energy init instead of 1
            this.updateCell(curCell, Object.create(NussinovCell4dTrace).init([], [[i, this.seq2_length + 1 - j]]));
        }

        for (var p = i + 1; p <= k; ++p) {
            for (var q = j + 1; q <= l; ++q) {
                // Energy loop instead of 1
                if (this.getValue(p, k, q, l) > 0) {
                    // This basepair can be added only if it encloses a structure.
                    this.updateCell(curCell, Object.create(NussinovCell4dTrace).init([[p, k, q, l]], [[i, this.seq2_length + 1 - j]]));
                }
            }
        }
    }

    return curCell;
};

DPAlgorithm_hybrid.computeMatrix = function(input) {
    var splitSeq = input.sequence().indexOf('X');
    var sequence1 = input.sequence().substr(0, splitSeq);
    var sequence2 = input.sequence().substr(parseInt(input.loopLength()) + splitSeq + 1);

    this.Tables[0].init(sequence1, sequence2, "RNAHybrid");
    this.Tables[0].minLoopLength = parseInt(input.loopLength());
    
    this.Tables[0].computeAllCells();

    //console.log(this.Tables[0].simpleRepresentation());
    return this.Tables;

};



var DPAlgorithm_rnaup = Object.create(DPAlgorithm);

DPAlgorithm_rnaup.Description = "accessibility-based RNA-RNA interaction prediction";
DPAlgorithm_rnaup.Tables = new Array();
DPAlgorithm_rnaup.Tables.push(Object.create(NussinovMatrix4d));
DPAlgorithm_rnaup.Tables.push(Object.create(NussinovMatrix));
DPAlgorithm_rnaup.Tables.push(Object.create(NussinovMatrix));

DPAlgorithm_rnaup.Tables[0].latex_representation = "\\begin{array}\\ I^{i, k}_{j, l} &=& \\max \\begin{cases} - (E^{bp}\\cdot D^{i, k}_{j, l} +ED^{1}_{i,k} +ED^{2}_{j, l}) &\\text{if } D^{i, k}_{j, l} > 0 \\\\ -\\infty & \\text{otherwise} \\end{cases} \\\\ \\\\ D^{i, k}_{j, l} &=& \\max \\begin{cases} 1 & \\text{if } S^1_i, \\overleftarrow{S_j^2} \\text{ compl.}, i = k, j = l \\\\ \\underset{\\substack{i<p\\leq k,\\;j<q\\leq l}}{\\max}\\left( 1 + D_{q, l}^{p, k} \\right) & \\text{if } S^1_i, \\overleftarrow{S^2_j} \\text{ compl.}, i < k, j < l\\\\ -\\infty & \\text{otherwise} \\end{cases} \\\\ \\\\ ED^{1}_{i,k} &=& - RT \\cdot \\log(P^{u1}_{i,k}),\\quad\\quad\\quad ED^{2}_{j, l} \\;=\\; - RT \\cdot \\log(P^{u2}_{j, l}) \\end{array}";


DPAlgorithm_rnaup.Tables[0].computeCell = function(i, k, j, l) {

    var curCell = Object.create(NussinovCell4d).init(i, k, j, l, 0);

    if (this.isInvalidState(i, k, j, l)) {
        return curCell;
    }
    if (DPAlgorithm_hybrid.Tables[0].getValue(i, k, j, l) == 0) {
        return curCell;
    }
    var logP = Math.log(DPAlgorithm_rnaup.Tables[1].getValue(i, k)) + Math.log(DPAlgorithm_rnaup.Tables[2].getValue(j, l));

    // I[i,k,j,l] = Ebp*H[i,k,j,l] -RT*ln(P^u_1[i,k]) -RT*ln(reverseP^u_2[j,l])
    // negate value for generic maximization optimization
    curCell.value = - (this.energy * DPAlgorithm_hybrid.Tables[0].getValue(i, k, j, l) - this.energy_normal * logP);
    curCell.traces = DPAlgorithm_hybrid.Tables[0].getCell(i, k, j, l).traces;
    return curCell;
};

DPAlgorithm_rnaup.computeMatrix = function(input) {
    DPAlgorithm_hybrid.computeMatrix(input);

    var splitSeq = input.sequence().indexOf('X');
    var sequence1 = input.sequence().substr(0,splitSeq);
    var sequence2 = input.sequence().substr(parseInt(input.loopLength()) + splitSeq + 1);

    // Clone the matrices of McCaskill algorithm
    NussinovDPAlgorithm_McCaskill.computeMatrix({sequence: function(){return sequence1;}, loopLength: input.loopLength, energy: input.energy, energy_normal: input.energy_normal});
    this.Tables[1] = JSON.parse(JSON.stringify(NussinovDPAlgorithm_McCaskill.Tables[3]));
    this.Tables[1].getRecursionInLatex = NussinovDPAlgorithm_McCaskill.Tables[3].getRecursionInLatex;
    this.Tables[1].getCell = NussinovDPAlgorithm_McCaskill.Tables[3].getCell;
    this.Tables[1].getValue = NussinovDPAlgorithm_McCaskill.Tables[3].getValue;
    this.Tables[1].isInvalidState = NussinovDPAlgorithm_McCaskill.Tables[3].isInvalidState;

    NussinovDPAlgorithm_McCaskill.computeMatrix({sequence: function(){return sequence2;}, loopLength: input.loopLength, energy: input.energy, energy_normal: input.energy_normal});
    this.Tables[2] = JSON.parse(JSON.stringify(NussinovDPAlgorithm_McCaskill.Tables[3]));
    this.Tables[2].getRecursionInLatex = NussinovDPAlgorithm_McCaskill.Tables[3].getRecursionInLatex;
    this.Tables[2].getCell = NussinovDPAlgorithm_McCaskill.Tables[3].getCell;
    this.Tables[2].getValue = NussinovDPAlgorithm_McCaskill.Tables[3].getValue;
    this.Tables[2].isInvalidState = NussinovDPAlgorithm_McCaskill.Tables[3].isInvalidState;
    
    this.Tables[0].init(sequence1, sequence2, "RNAup");
    this.Tables[0].minLoopLength = parseInt(input.loopLength());
    this.Tables[0].energy = input.energy();
    this.Tables[0].energy_normal = input.energy_normal();

    this.Tables[0].computeAllCells();

    //console.log(this.Tables[0].simpleRepresentation());
    return this.Tables;
};


/**
 * WUCHTY (generic doesnt give suboptimal structures)
 * Construct the tracebacks of the optimal solutions. 
 */
var wuchty4d = function (xmat, maxSOS) {
    //console.log("wuchty4d");
    var sigma_0 = [];
    var NMax = 0;
    //console.log(xmat.seq1_length, xmat.seq2_length);
    for (var i = 0; i <= xmat.seq1_length; ++i) {
        for (var k = i; k <= xmat.seq1_length; ++k) {
            for (var j = 0; j <= xmat.seq2_length; ++j) {
                for (var r = j; r <= xmat.seq2_length; ++r) {
                    if (xmat.getValue(i, k, j, r) >= NMax) {
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
    //console.log(NMax);
    //console.log(sigma_0);
    var SOS = [];
    var loop = 0;
    for (var sig = 0; sig < sigma_0.length; ++sig) {
        var S = {sigma: sigma_0[sig], P: [], traces: []};
        var R = [S];
        
        while (R.length != 0) {
            // Pop R
            var pop_R = R.pop();
            var sigma = pop_R.sigma;
            var P = pop_R.P;
            var t_traces = JSON.stringify(pop_R.traces);
            var traces = JSON.parse(t_traces);

            var sigma_remaining = 0;
            for (var s in sigma) {
                if (!xmat.isInvalidState(sigma[s][0], sigma[s][1], sigma[s][2], sigma[s][3])) sigma_remaining++;
            }
            if (sigma.length == 0 || sigma_remaining == 0) {
                var temp_sos = {value: NMax, structure: xmat.conv_str(P), traces: traces, rep4d: repres(visualize4d(xmat.sequence1, xmat.sequence2, P))};
                //console.log("P", JSON.stringify(P), JSON.stringify(xmat.conv_str(P)));
                SOS.push(temp_sos);

                if (SOS.length >= maxSOS) {
                    break;
                }
            } else {
                var idx = sigma.pop();
                if (xmat.isInvalidState(idx[0], idx[1], idx[2], idx[3])) {
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
                    if (trace.bps[0] != undefined) {
                        S_prime.P.push(trace.bps[0]);
                    }

                    // add traces info in S_prime
                    var temp_trace = [idx];
                    temp_trace.push(trace_p.parents);
                    if (traces.length == 0) {
                        S_prime.traces = [temp_trace];
                    } else {
                        var clone_traces = JSON.stringify(traces);
                        var parse_clone_traces = JSON.parse(clone_traces);
                        parse_clone_traces.unshift(temp_trace);
                        S_prime.traces = parse_clone_traces;
                    }
                    R.push(S_prime);
                }
            }
        }

        if (SOS.length >= maxSOS) {
            break;
        }
    }
    
    // check if no interaction stored so far
    if (SOS.length == 0) {
    	// push empty interaction without trace
    	SOS.push( {structure: xmat.conv_str([]), traces: [], rep4d: " "} );
    }
    
    //console.log('final: ', JSON.stringify(SOS));
    return SOS;
}

/**
 * global list of available Nussinov algorithm implementations
 */
var availableAlgorithms = {

    /** original unique recursion */
    nussinovUnique: NussinovDPAlgorithm_Unique,//NussinovMatrix_unique,
    
    /** ambiguous recursion */
    nussinovAmbiguous1: NussinovDPAlgorithm_Ambiguous,//NussinovMatrix_ambiguous,

    /** nussinov neo recursion */
    nussinovAmbiguous2: NussinovDPAlgorithm_Ambiguous2,

    /** Most Ambiguos recursion */
    nussinovMostAmbiguous: NussinovDPAlgorithm_MostAmbiguous, 

    /** McCaskill */
    mcCaskill: NussinovDPAlgorithm_McCaskill,

    /** structure counting */
    nussinovCounting: NussinovDPAlgorithm_structuresCount,

    /** Maximum Expected Accuracy*/
    MaxExpAcc: DPAlgorithm_MEA,

    /** Co-fold*/
    coFold: DPAlgorithm_coFold,

    hybrid: DPAlgorithm_hybrid,

    rnaup: DPAlgorithm_rnaup

};