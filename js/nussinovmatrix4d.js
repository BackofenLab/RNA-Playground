/**
 * Created by moni on 02.08.16.
 */
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
        this.j = j;
        this.k = k;
        this.l = l;
        this.value = value;
        this.traces = [];
        // this access for chaining
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

            // create matrix cells
            var n1 = this.sequence1.length;
            var n2 = this.sequence2.length;

            for (var i = 0; i <= n1; i++) {
                this.cells[i] = [];
                for (var k = 0; k <= n1; ++k) {
                    this.cells[i][k] = [];
                    for (var j = 0; j <= n2; j++) {
                        // create new cell and initialize
                        this.cells[i][k][j] = [];
                        for (var r = 0; r <= n2; ++r) {
                            this.cells[i][k][j][r] = Object.create(NussinovCell4d).init(i, k, j, r, null);
                        }
                    }
                }
                ;
            }
            ;

            return this;
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
            if (i < 0 || j < 0 || k < 0 || l < 0 || k < i || l < j || k >= this.getDim() || l >= this.getDim() || i >= this.getDim() || j >= this.getDim()) {
                return null;
            }
            if (this.cells[i][k][j][l] === null) {
                //this.cells[i][k][j][l] = computeCell(i, k, j, l);
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


        computeValue: function (i, k, j, l) {
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
        getValue: function (i, k, j, l) {
            // access cell at location (i,j) in the matrix
            //console.log(i, k, j, l);
            var cell = this.getCell(i, k, j, l);
            if (cell === null) {
                return null;
            }
            // check if invalid cell
            if (cell.value === null) {
                cell.value = this.computeValue(i, k, j, l);
            }
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
            console.log("WARNING: getRecursionInLatex() not implemented in NussinovMatrix4d superclass; overwrite in subclass!");
            return "$$" + this.latex_representation + "$$";

        }
        ,
    
        simpleRepresentation: function() {

            var n1 = this.sequence1.length;
            var n2 = this.sequence2.length;

            console.log(n1, n2);
            console.log(this.sequence1, this.sequence2);
            var res = "";
            for (var i = 0; i <= n1; ++i)
                for (var k = i; k <= n1; ++k)
                    for (var j = 0; j <= n2; ++j)
                        for (var l = j; l <= n2; ++l)
                            if (this.getValue(i, k, j, l) != null) {
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
DPAlgorithm_hybrid.Tables[0].latex_representation = "D_{i, k}^{j, l} = \\max \\begin{cases} E^{init}(i, j) & \\mathcal{R}^1_i, \\mathcal{R}^2_j  pairs, i = k, j = l \\\\ \\max_{p,q}{ E^{loop}(i, j, p, q) + D_{q, l}^{p, k} } & \\mathcal{R}^1_i, \\mathcal{R}^2_j  pairs, i < k, j < l\\\\ 0 & otherwise \\end{cases}";

DPAlgorithm_hybrid.Tables[0].computeValue = function(i, k, j, l) {
    if (i < 0 || j < 0 || k < 0 || l < 0 || k < i || l < j || i >= this.getDim() || j >= this.getDim() || k >= this.getDim() || l >= this.getDim()) {
        return 0;
    }
    var ret = 0;

    if (RnaUtil.areComplementary(this.sequence1[i - 1], this.sequence2[j - 1])) {
        if (i === k && j === l) {
            // Energy init instead of 1
            ret = Math.max(ret, 1);
        }

        for (var p = i + 1; p <= k; ++p) {
            for (var q = j + 1; q <= l; ++q) {
                // Energy loop instead of 1
                if (i < k && j < l)
                ret = Math.max(ret, 1 + this.getValue(p, k, q, l));
            }
        }
    }

    return ret;
};



DPAlgorithm_hybrid.computeMatrix = function(input) {
    var splitSeq = input.sequence().indexOf('X');
    var sequence1 = input.sequence().substr(0,splitSeq);
    var sequence2 = input.sequence().substr(parseInt(input.loopLength())+splitSeq + 1);

    console.log(sequence1, sequence2);
    this.Tables[0].init(sequence1, sequence2, "RNAHybrid");

    var n1 = sequence1.length;
    var n2 = sequence2.length;
    for (var i = 0; i <= n1; ++i)
        for (var k = i; k <= n1; ++k)
            for (var j = 0; j <= n2; ++j)
                for (var l = j; l <= n2; ++l)
                    this.Tables[0].getValue(i, k, j, l);

    console.log(this.Tables[0].simpleRepresentation());
    return this.Tables;

};




var DPAlgorithm_rnaup = Object.create(DPAlgorithm);

DPAlgorithm_rnaup.Description = "RNA to RNA matching";
DPAlgorithm_rnaup.Tables = new Array();
DPAlgorithm_rnaup.Tables.push(Object.create(NussinovMatrix4d));
DPAlgorithm_rnaup.Tables.push(Object.create(NussinovMatrix));
DPAlgorithm_rnaup.Tables.push(Object.create(NussinovMatrix));

DPAlgorithm_rnaup.Tables[0].latex_representation = "D_{i, k}^{j, l} = \\max \\begin{cases} E^{init}(i, j) & \\mathcal{R}^1_i, \\mathcal{R}^2_j  pairs, i = k, j = l \\\\ \\max_{p,q}{ E^{loop}(i, j, p, q) + D_{q, l}^{p, k} } & \\mathcal{R}^1_i, \\mathcal{R}^2_j  pairs, i < k, j < l\\\\ 0 & otherwise \\end{cases}";



DPAlgorithm_rnaup.Tables[0].computeValue = function(i, k, j, l) {
    if (i < 0 || j < 0 || k < 0 || l < 0 || k < i || l < j || i >= this.getDim() || j >= this.getDim() || k >= this.getDim() || l >= this.getDim()) {
        return 0;
    }
    // I[i1,j1,i2,j2] = Ebp*H[i1,j1,i2,j2] -RT*ln(P^u_1) -RT*ln(P^u_2)

    var logP = Math.log(DPAlgorithm_rnaup.Tables[1].getValue(i, k)) + Math.log(DPAlgorithm_rnaup.Tables[2].getValue(j, l));
    return this.energy * DPAlgorithm_hybrid.Tables[0].getValue(i, k, j, l) - this.energy_normal * logP;
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


    NussinovDPAlgorithm_McCaskill.computeMatrix({sequence: function(){return sequence2;}, loopLength: input.loopLength, energy: input.energy, energy_normal: input.energy_normal});
    this.Tables[2] = JSON.parse(JSON.stringify(NussinovDPAlgorithm_McCaskill.Tables[3]));
    this.Tables[2].getRecursionInLatex = NussinovDPAlgorithm_McCaskill.Tables[3].getRecursionInLatex;
    this.Tables[2].getDim = NussinovDPAlgorithm_McCaskill.Tables[3].getDim;
    this.Tables[2].getCell = NussinovDPAlgorithm_McCaskill.Tables[3].getCell;
    this.Tables[2].getValue = NussinovDPAlgorithm_McCaskill.Tables[3].getValue;


    console.log(sequence1, sequence2);
    this.Tables[0].init(sequence1, sequence2, "RNAHybrid");
    this.Tables[0].energy = input.energy();
    this.Tables[0].energy_normal = input.energy_normal();

    var n1 = sequence1.length;
    var n2 = sequence2.length;
    for (var i = 0; i <= n1; ++i)
        for (var k = i; k <= n1; ++k)
            for (var j = 0; j <= n2; ++j)
                for (var l = j; l <= n2; ++l)
                    this.Tables[0].getValue(i, k, j, l);

    console.log(this.Tables[0].simpleRepresentation());
    return this.Tables;

};


var availableAlgorithms = {

    /** ambiguous recursion */
    hybrid: DPAlgorithm_hybrid,//NussinovMatrix_ambiguous,

    rnaup: DPAlgorithm_rnaup,

};

