/**
 * Returns ONE optimal traceback if the matrix was already computed.
 *
 * @returns a Traceback object or null if the matrix wasnt computed yet
 */

/**
 * Encodes a full/partial traceback for a given cell
 */
var Traceback = {

    /** structure in dot-bracket-notation */
    structure: "",
    /** list of cell traces of the form [source, parent1, parent2, ...] */
    traces: [],
    //partialStructure: [],

    /**
     * set structure value
     */
    setStructure: function (struct) {
        this.structure = struct;
    },

    /**
     * set traces
     */
    setTraces: function (args) {
        if (this.traces == [])
            this.traces = args;
        else this.traces.push(args);
    },

    /**
     * initialize structure with dot string of sequence length
     * @param (int) len length of string to be created.
     */
    init: function (len) {
        if (len <= 0) return "";
        for (var i = 0; i < len; i++) {
            if (this.structure == undefined)
                this.structure = ".";
            else
                this.structure = this.structure + "."
        }
        this.traces = [];
        return this;
    },

    /**
     * Add Brackets to structure for given basepair
     * @param {array} basepair
     */
    addBrackets: function (basepair) {
        if (basepair[1] > this.structure.length || basepair == null) {
            console.log("String too short or empty bp");
            return;
        }

        this.structure = this.structure.substr(0, basepair[0] - 1) + '(' + this.structure.substr(basepair[0], basepair[1] - basepair[0] - 1) + ')' + this.structure.substr(basepair[1]);
        return this.structure;
    }


};

var NussinovMatrix = {
    /**
     * WUCHTY
     */
    wuchty: function (delta) {
        var seq_length = this.sequence.length;
        var Nmax = this.getCell(1, seq_length).value;

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
                if ((sigma[s][0]) <= (sigma[s][1] - this.minLoopLength)) sigma_remaining++;
            }

            if (sigma.length == 0 || sigma_remaining == 0) {
                var temp_sos = {structure: this.conv_str(P, seq_length), traces: traces};
                SOS.push(temp_sos);
            }

            else {
                var ij = sigma.pop();

                if (ij[0] > ij[1]) {
                    pop_R.traces = traces;
                    R.push(pop_R);
                    continue;
                }

                var json_ij_traces = JSON.stringify(this.getCell(ij[0], ij[1]).traces);
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
    },

    getOptimalTraceback: function (i, j, optimalTrace) {
        // check if there was something computed yet
        if (this.sequence === null) {
            return null;
        }
        ;

        // create empty trace to fill
        if (typeof optimalTrace == "undefined") {
            //console.log("optimalTrace object is not defined");
        }
        ;

        if (this.getValue(i, j) === null) {
            return;
        }
        else {
            //console.log(this.getTraces(i,j)[0]);

            for (var p in this.getTraces(i, j)[0].parents) {
                var parent = this.getTraces(i, j)[0].parents[p];
                var basepair = this.getTraces(i, j)[0].bps;

                if (typeof basepair[0] != "undefined") {
                    optimalTrace.addBrackets(basepair[0]);
                }
                ;

                this.getOptimalTraceback(parent[0], parent[1], optimalTrace);

            }
            ;
            optimalTrace.traces.push([[i, j], this.getTraces(i, j)[0].parents]);

        }
        ;

        return optimalTrace;
    }
    ,
}
/****** NussinovMatrix_ambiguous extending NussinovMatrix ************************/

/**
 * Implements an ambiguous recursion
 */
var NussinovMatrix_ambiguous = Object.create(NussinovMatrix);


/**********  OVERWRITING + EXTENSION  ********************/

/**
 * Returns a description for the implemented recursion
 *
 * @returns {string} description of the recursion
 */
NussinovMatrix_ambiguous.getRecursionDescription = function () {
    return "Ambiguous recursion";
};


/**
 * Access to the recursion in LaTeX encoding that is used in the computeMatrix implementation
 *
 * @returns {string} latex encoding of the recursion
 */
NussinovMatrix_ambiguous.getRecursionInLatex = function () {
    return "$D(i,j) = \\max \\begin{cases} D(i+1,j) & S_i \\text{ unpaired} \\\\ D(i,j-1) & S_j \\text{ unpaired} \\\\ D(i+1,j-1)+1 &  S_i,S_j \\text{ compl. base pair and } i+ l< j \\\\ \\max_{i< k< (j-1)} D(i,k)+D(k+1,j) & \\text{decomposition} \\end{cases}$";
};

/**
 * Fills the matrix according to the recursion.
 *
 * @param {string} sequence the RNA sequence to compute the matrix for
 * @param {int} minLoopLength the minimal loop length to be used for computation
 *
 * @returns {NussinovMatrix} this for call chaining
 */
NussinovMatrix_ambiguous.computeMatrix = function (sequence, minLoopLength) {
// resize and initialize matrix
    this.init(sequence, "ambiguous");
// store minimal loop length
    this.minLoopLength = minLoopLength;
    //console.log("computing ambiguos matrix");
// fill matrix by diagonals
// iterate over all substructure spans that can have a base pair
    for (var span = minLoopLength; span < this.getDim(); span++) {
        // iterate over all rows
        for (var i = 1; i < this.getDim() - minLoopLength; i++) {
            // get column for current span
            var j = i + span;

            // i unpaired
            this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i + 1, j]], []));

            // j unpaired
            this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, j - 1]], []));

            // check (i,j) base pair
            if ((j - i > minLoopLength) && RnaUtil.areComplementary(this.sequence[i - 1], this.sequence[j - 1])) {
                // get value for base pair
                this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i + 1, j - 1]], [[i, j]]));
            }
            ;

            // check decomposition into substructures (minLength==2)
            for (var k = i + 1; k < (j - 1); k++) {
                // get decomposition value
                this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, k], [k + 1, j]], []));
            }
            ;
        }
        ;
    }
    ;

    return this;
};

/****** NussinovMatrix_unique extending NussinovMatrix ************************/

/**
 * Implements a non-ambiguous recursion
 */
var NussinovMatrix_unique = Object.create(NussinovMatrix);


/**********  OVERWRITING + EXTENSION  ********************/

/**
 * Returns a description for the implemented recursion
 *
 * @returns {string} description of the recursion
 */
NussinovMatrix_unique.getRecursionDescription = function () {
    return "Recursion by Nussinov et al. (1978) with unique decomposition";
};

/**
 * Access to the recursion in LaTeX encoding that is used in the computeMatrix implementation
 *
 * @returns {string} latex encoding of the recursion
 */
NussinovMatrix_unique.getRecursionInLatex = function () {
    return "$D(i,j) = \\max \\begin{cases} D(i,j-1) & S_j \\text{ unpaired} \\\\ \\max_{i\\leq k< (j-l)} D(i,k-1)+D(k+1,j-1)+1 & S_k,S_j \\text{ compl. base pair} \\end{cases}$";
};

/**
 * Fills the matrix according to the recursion.
 *
 * @param {string} sequence the RNA sequence to compute the matrix for
 * @param {int} minLoopLength the minimal loop length to be used for computation
 *
 * @returns {NussinovMatrix} this for call chaining
 */
NussinovMatrix_unique.computeMatrix = function (sequence, minLoopLength) {

// resize and initialize matrix
    this.init(sequence, "unique");
// store minimal loop length
    this.minLoopLength = minLoopLength;

// fill matrix by diagonals
// iterate over all substructure spans that can have a base pair
    for (var span = minLoopLength; span < this.getDim(); span++) {
        // iterate over all rows
        for (var i = 0; i < this.getDim() - minLoopLength; i++) {
            // get column for current span
            var j = i + span;

            // j unpaired
            this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, j - 1]], []));

            // check base pair based decomposition : (k,j) base pair
            for (var k = i; k + minLoopLength < j; k++) {
                // check if sequence positions are compatible
                if (RnaUtil.areComplementary(this.sequence[k - 1], this.sequence[j - 1])) {
                    this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, k - 1], [k + 1, j - 1]], [[k, j]]));
                }
                ;
            }
            ;
        }
        ;
    }
    ;

    return this;
};

/**
 * getSubstructures(for wuchty)
 */
NussinovMatrix_unique.getSubstructures = function(sigma, P, traces, delta, maxLengthR) {
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
;

/**
 * getSubstructures(for wuchty)
 */
NussinovMatrix_ambiguous.getSubstructures = function (sigma, P, traces, delta, maxLengthR) {
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

    //console.log("returning R:", JSON.stringify(R));
    return R;
}
;



/*****************************************************************************/
/****** NussinovMatrix_structuresCount extending NussinovMatrix ************************/

/**
 * Implements a non-ambiguous recursion
 */
var NussinovMatrix_structuresCount = Object.create(NussinovMatrix);


/**********  OVERWRITING + EXTENSION  ********************/

/**
 * Returns a description for the implemented recursion
 *
 * @returns {string} description of the recursion
 */
NussinovMatrix_structuresCount.getRecursionDescription = function () {
    return "Recursion to count number of total structures.";
};

/**
 * Access to the recursion in LaTeX encoding that is used in the computeMatrix implementation
 *
 * @returns {string} latex encoding of the recursion
 */
NussinovMatrix_structuresCount.getRecursionInLatex = function () {
    return "$C_{i,j} = C_{i,j-1} + \\sum_{i\\leq k <(j-l) \\atop S_k,S_j \\text{ pair}} C_{i,k-1} * C_{k+1,j-1} * 1 $";
};

/**
 * Fills the matrix according to the recursion.
 *
 * @param {string} sequence the RNA sequence to compute the matrix for
 * @param {int} minLoopLength the minimal loop length to be used for computation
 *
 * @returns {NussinovMatrix_structuresCount} this for call chaining
 */
NussinovMatrix_structuresCount.computeMatrix = function (sequence, minLoopLength) {

// resize and initialize matrix
    this.init(sequence, "structuresCount");
    //console.log(sequence, this.sequence.length);
// store minimal loop length
    this.minLoopLength = minLoopLength;

// fill matrix by diagonals
// iterate over all substructure spans that can have a base pair
    for (var span = minLoopLength; span < this.getDim(); span++) {
        // iterate over all rows
        for (var i = 0; i < this.getDim() - minLoopLength; i++) {
            // get column for current span
            var j = i + span;

            // j unpaired
            this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, j - 1]], []));

            // check base pair based decomposition : (k,j) base pair
            for (var k = i; k + minLoopLength < j; k++) {
                // check if sequence positions are compatible
                if (RnaUtil.areComplementary(this.sequence[k - 1], this.sequence[j - 1])) {
                    this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, k - 1], [k + 1, j - 1]], [[k, j]]));
                }
                ;
            }
            ;
        }
        ;
    }
    ;

    return this;
};


 NussinovMatrix_structuresCount.updateCell = function(i, j, curAncestor){
    var val = 1;
    // get cell to update
    var curCell = this.getCell(i, j);
    // check if something to update
    if (curCell === null) {
        return;
    }

    // add scores of ancestor cells
    for (var x = 0; x < curAncestor.parents.length; x++) {
        if(curAncestor.parents[x][0] > 0 && curAncestor.parents[x][1] > 0
            && curAncestor.parents[x][0] < curAncestor.parents[x][1]) {
            //var q = this.getValue(curAncestor.parents[x][0], curAncestor.parents[x][1]) * Math.exp(-1 * (-1));
            var q = this.getValue(curAncestor.parents[x][0], curAncestor.parents[x][1]);
            val *= q;
        }
    }
    // check for case (i,j-1)
    if(curAncestor.bps.length === 0){
        curCell.value = val;
    }
    else {
        curCell.value += val;
    }
    // curCell.traces.push(curAncestor);
}
 ;


//var availableAlgorithms = [NussinovMatrix_ambiguous, NussinovMatrix_unique];

/****** McKaskill_base extending NussinovMatrix ************************/



/**
 * Implements a non-ambiguous recursion
 */
/*
 var McKaskill_base = Object.create(NussinovMatrix);

 /**
 * Returns a description for the implemented recursion
 *
 * @returns {string} description of the recursion
 */
/*
 McKaskill_base.getRecursionDescription = function () {
 return "Recursion to count the energy of all the structures.";
 };
 */
/**
 * Access to the recursion in LaTeX encoding that is used in the computeMatrix implementation
 *
 * @returns {string} latex encoding of the recursion
 */
/*
 McKaskill_base.getRecursionInLatex = function () {
 return "$$Q_{i,j}^{b} = \\begin{cases} Q_{i + 1, j - 1} * \\exp(-E(bp)/RT) & \\text{ if }i,j \\text{ can form base pair} \\\\ 0 & \\text{ otherwise}\\end{cases}$$";
 };
 */

/**********  OVERWRITING + EXTENSION  ********************/
self/*
 McKaskill_base.updateCell = function(i, j){
 // get cell to update
 var curCell = this.getCell(i, j);
 // check if something to update
 if (curCell === null) {
 return;
 }
 curCell.value = McKaskill_simple.getValue(i + 1, j - 1) * Math.exp(1);
 }
 ;
 */


/**
 * Fills the matrix according to the recursion.
 *
 * @param {string} sequence the RNA sequence to compute the matrix for
 * @param {int} minLoopLength the minimal loop length to be used for computation
 *
 * @returns {McKaskill_base} this for call chaining
 */
/*
 McKaskill_base.computeMatrix = function (sequence, minLoopLength) {

 // resize and initialize matrix
 this.init(sequence, "structuresEnergy");
 //console.log(sequence, this.sequence.length);
 // store minimal loop length
 this.minLoopLength = minLoopLength;

 // fill matrix by diagonals
 // iterate over all substructure spans that can have a base pair
 for (var span = minLoopLength; span < this.getDim(); span++) {
 // iterate over all rows
 for (var i = 0; i < this.getDim() - minLoopLength; i++) {
 // get column for current span
 var j = i + span;

 // j unpaired
 this.updateCell(i, j);

 // check base pair based decomposition : (k,j) base pair
 //for (var k = i; k + minLoopLength < j; k++) {
 // check if sequence positions are compatible
 //if (RnaUtil.areComplementary(this.sequence[k - 1], this.sequence[j - 1])) {
 //this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, k - 1], [k + 1, j - 1]], [[k, j]]));
 //}
 ;
 // }
 ;

 }
 ;
 }
 ;

 return this;
 };
 */

/**
 * Implements a non-ambiguous recursion
 */
/*
 var McKaskill_simple = Object.create(NussinovMatrix);
 */
/**
 * Returns a description for the implemented recursion
 *
 * @returns {string} description of the recursion
 */
/*
 McKaskill_simple.getRecursionDescription = function () {
 return "Recursion to count the energy of all the structures.";
 };
 */

/**
 * Access to the recursion in LaTeX encoding that is used in the computeMatrix implementation
 *
 * @returns {string} latex encoding of the recursion
 */
/*
 McKaskill_simple.getRecursionInLatex = function () {
 return "$$Q_{i,j} = Q_{i,j-1} + \\sum_{i\\leq k <(j-l)} Q_{i,k-1} * Q^{b}_{k,j} $$";
 };
 */


/**********  OVERWRITING + EXTENSION  ********************/
/*
 McKaskill_simple.updateCell = function(i, j, curAncestor){
 var val = 1;
 // get cell to update
 var curCell = this.getCell(i, j);
 // check if something to update
 if (curCell === null) {
 return;
 }

 // add scores of ancestor cells
 for (var x = 0; x < curAncestor.parents.length; x++) {
 if(curAncestor.parents[x][0] > 0 && curAncestor.parents[x][1] > 0
 && curAncestor.parents[x][0] < curAncestor.parents[x][1]) {
 val *= this.getValue(curAncestor.parents[x][0], curAncestor.parents[x][1]) * Math.exp(1);
 }
 }
 // check for case (i,j-1)
 if(curAncestor.bps.length === 0){
 curCell.value = val;
 }
 else {
 curCell.value += val;
 }
 // curCell.traces.push(curAncestor);
 }
 ;
 */


/**
 * Fills the matrix according to the recursion.
 *
 * @param {string} sequence the RNA sequence to compute the matrix for
 * @param {int} minLoopLength the minimal loop length to be used for computation
 *
 * @returns {McKaskill_simple} this for call chaining
 */
/*
 McKaskill_simple.computeMatrix = function (sequence, minLoopLength) {

 // resize and initialize matrix
 this.init(sequence, "structuresEnergy");
 //console.log(sequence, this.sequence.length);
 // store minimal loop length
 this.minLoopLength = minLoopLength;

 // fill matrix by diagonals
 // iterate over all substructure spans that can have a base pair
 for (var span = minLoopLength; span < this.getDim(); span++) {
 // iterate over all rows
 for (var i = 0; i < this.getDim() - minLoopLength; i++) {
 // get column for current span
 var j = i + span;

 // j unpaired
 this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, j - 1]], []));

 // check base pair based decomposition : (k,j) base pair
 for (var k = i; k + minLoopLength < j; k++) {
 // check if sequence positions are compatible
 if (RnaUtil.areComplementary(this.sequence[k - 1], this.sequence[j - 1])) {
 this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, k - 1], [k + 1, j - 1]], [[k, j]]));
 }
 ;
 }
 ;
 }
 ;
 }
 ;

 return this;
 };
 */

/***************   TESTING  ********************/
//console.log(NussinovMatrix_structuresCount.getRecursionDescription());
//var m = NussinovMatrix_structuresCount.computeMatrix("CGGGC", 0);
//console.log(m.toString());
//
//{
//    var selectedAlgorithm = "nussinovOriginal";
////console.log("\n####  using", selectedAlgorithm, "  ####\n");
////console.log(availableAlgorithms[selectedAlgorithm].getRecursionDescription() + "\n");
////console.log(availableAlgorithms[selectedAlgorithm].getRecursionInLatex() + "\n");
////console.log(availableAlgorithms[selectedAlgorithm].computeMatrix("CGAGC", 0).toString());
////console.log(availableAlgorithms[selectedAlgorithm].computeMatrix("CGAGC", 1).toString());
////console.log("optimal structure :", availableAlgorithms[selectedAlgorithm].getOptimalTraceback().structure);
//
//    //selectedAlgorithm = "nussinovUnique";
//    console.log("\n####  using", selectedAlgorithm, "  ####\n");
//    //console.log(availableAlgorithms[selectedAlgorithm].getRecursionDescription() + "\n");
//    //console.log(availableAlgorithms[selectedAlgorithm].getRecursionInLatex() + "\n");
//    //console.log(availableAlgorithms[selectedAlgorithm].computeMatrix("GCGCG", 0).toString());
//    //console.log(availableAlgorithms[selectedAlgorithm].computeMatrix("GCGCG", 1).toString());
//    var xmat = availableAlgorithms[selectedAlgorithm].computeMatrix("CGAGCAGC", 0);
//    var opt = Object.create(Traceback);
//    opt.init(xmat.sequence.length);
//
//    var xy = wuchty_2nd(xmat, 0, selectedAlgorithm);
//    console.log("wuchty trace\n");
//
//    for (var i in xy) {
//        console.log(JSON.stringify(xy[i].traces), JSON.stringify(xy[i].structure));
//    }
//
//
//    //var x = xmat.traceallOpt(1, xmat.sequence.length);
//    //var struct = "";
//    //for (var i in x) {
//    //    struct += x[i].structure + ",";
//    //}
//    //console.log();
//    //for (var i in x) {
//        //console.log(JSON.stringify(x[i].traces), JSON.stringify(x[i].structure))
//    //}
//
//}


