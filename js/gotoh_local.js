/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

/**
 * Defines tasks after page-loading.
 */
$(document).ready(function () {
    if (loaded === ALGORITHMS.GOTOH_LOCAL) {  // to avoid self execution on a script import
        gotohLocal.startGotohLocal();
        loaded = ALGORITHMS.NONE;
    }
});

(function () {  // namespace
    // public methods
    namespace("gotohLocal", startGotohLocal, GotohLocal);

    // instances
    var alignmentInstance;
    var gotohLocalInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Function managing objects.
     */
    function startGotohLocal() {
        var subadditiveAlignmentInterface = new interfaces.subadditiveAlignmentInterface.SubadditiveAlignmentInterface();
        subadditiveAlignmentInterface.startSubadditiveAlignmentAlgorithm(GotohLocal, ALGORITHMS.GOTOH_LOCAL);
    }

    /*---- ALGORITHM ----*/
    /**
     * Computes the optimal, local affine alignment.
     * Combination of the Smith-Waterman and Gotoh-Algorithm.
     * @constructor
     * @augments Alignment
     * @see: https://doi.org/10.1016/0022-2836(82)90398-9 and https://doi.org/10.1016/0022-2836(81)90087-5
     *
     * Gotoh, Osamu.
     * "An improved algorithm for matching biological sequences."
     * Journal of molecular biology 162.3 (1982): 705-708.
     *
     * Smith, Temple F., and Michael S. Waterman.
     * "Identification of common molecular subsequences."
     * Journal of molecular biology 147.1 (1981): 195-197.
     */
    function GotohLocal() {
        gotohLocalInstance = this;

        // variables
        this.type = ALGORITHMS.GOTOH_LOCAL;
        this.numberOfTracebacks = 0;

        // instances
        alignmentInstance = new bases.alignment.Alignment(this);

        // public class methods
        this.getInput = getInput;

        this.setInput = setInput;
        this.compute = compute;
        this.getNeighboured = getNeighboured;
        this.getOutput = getOutput;

        this.setIO = setIO;
        this.getSuperclass = getSuperclass;
    }

    /**
     * Returns the input data of the algorithm.
     * @return {Object} - Contains all input data.
     */
    function getInput() {
        return inputData;
    }

    /**
     * Sets the algorithm input for an appropriate algorithm
     * which is using the inputViewmodel properties in its computations.
     * @param inputViewmodel {Object} - The InputViewmodel of an appropriate algorithm.
     */
    function setInput(inputViewmodel) {
        alignmentInstance.setIO(inputData, {});
        alignmentInstance.setSubadditiveAlignmentInput(inputViewmodel);
    }

    /**
     * Starts the computation.
     */
    function compute() {
        initializeMatrices();
        computeMatricesAndScore();
        computeTraceback();
        createAlignments();
        return [inputData, outputData];
    }

    /**
     * Initializes and creates the matrices.
     */
    function initializeMatrices() {
        createMatrices();
        initMatrices();
    }

    /**
     * Creates the matrices without initializing them.
     */
    function createMatrices() {
        outputData.matrix = new Array(inputData.matrixHeight);
        outputData.horizontalGaps = new Array(inputData.matrixHeight);
        outputData.verticalGaps = new Array(inputData.matrixHeight);

        for (var i = 0; i < inputData.matrixHeight; i++) {
            outputData.matrix[i] = new Array(inputData.matrixWidth);
            outputData.horizontalGaps[i] = new Array(inputData.matrixWidth);
            outputData.verticalGaps[i] = new Array(inputData.matrixWidth);
        }
    }

    /**
     * Initializes the default matrix and the gap matrices.
     */
    function initMatrices() {
        outputData.matrix[0][0] = 0;

        for (var i = 1; i < inputData.matrixHeight; i++) {
            outputData.matrix[i][0] = 0;
            outputData.horizontalGaps[i][0] = Number.NEGATIVE_INFINITY;
            outputData.verticalGaps[i][0] = SYMBOLS.GAP;
        }

        for (var j = 1; j < inputData.matrixWidth; j++) {
            outputData.matrix[0][j] = 0;
            outputData.horizontalGaps[0][j] = SYMBOLS.GAP;
            outputData.verticalGaps[0][j] = Number.NEGATIVE_INFINITY;
        }
    }

    /**
     * Computes the matrix by using the recursion function and the score.
     */
    function computeMatricesAndScore() {
        alignmentInstance.setIO(inputData, outputData);
        var maxValue = 0;

        // going through every matrix cell
        for (var i = 1; i < inputData.matrixHeight; i++) {
            var aChar = inputData.sequenceA[i - 1];

            for (var j = 1; j < inputData.matrixWidth; j++) {
                var bChar = inputData.sequenceB[j - 1];

                outputData.matrix[i][j] = alignmentInstance.affineRecursionFunction(aChar, bChar, i, j, Math.max, true);

                // storing maximum
                if (maxValue < outputData.matrix[i][j])
                    maxValue = outputData.matrix[i][j];
            }
        }

        // score is the maximum value
        outputData.score = maxValue;
    }

    /**
     * Initializes the traceback.
     * @override Alignment.computeTraceback()
     */
    function computeTraceback() {
        gotohLocalInstance.numberOfTracebacks = 0;

        // computing all traceback start-positions
        var backtraceStarts = alignmentInstance.getAllMaxPositions(inputData, outputData);

        outputData.tracebackPaths = [];
        outputData.moreTracebacks = false;

        for (var i = 0; i < backtraceStarts.length; i++) {
            var tracebackPaths = alignmentInstance.getLocalTraces([backtraceStarts[i]], inputData, outputData, -1, getNeighboured);
            outputData.tracebackPaths = outputData.tracebackPaths.concat(tracebackPaths);

            if (alignmentInstance.stopTraceback) break;
        }
    }

    /**
     * Returns the neighbours to which you can go from the current cell position used as input.
     * @param position {Vector} - Current cell position in matrix.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     * @param algorithm {Object} - Contains an alignment algorithm.
     * @return {Array} - Contains neighboured positions as Vector-objects.
     * @see Hint: The parameter algorithm is needed!
     * It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function getNeighboured(position, inputData, outputData, algorithm) {
        var neighboured = [];

        if (position.label === MATRICES.VERTICAL)
            return alignmentInstance.getVerticalNeighboured(position, inputData, outputData);
        else if (position.label === MATRICES.HORIZONTAL)
            return alignmentInstance.getHorizontalNeighboured(position, inputData, outputData);

        var left = position.j - 1;
        var up = position.i - 1;

        // retrieve values
        var aChar = left >= 0 ? inputData.sequenceB[left] : SYMBOLS.EMPTY;
        var bChar = up >= 0 ? inputData.sequenceA[up] : SYMBOLS.EMPTY;

        var currentValue = outputData.matrix[position.i][position.j];

        var matchOrMismatch = aChar === bChar ? inputData.match : inputData.mismatch;

        var diagonalValue = left >= 0 && up >= 0 ? outputData.matrix[up][left] : Number.NaN;
        var verticalValue = up >= 0 ? outputData.verticalGaps[position.i][position.j] : Number.NaN;
        var horizontalValue = left >= 0 ? outputData.horizontalGaps[position.i][position.j] : Number.NaN;

        // check
        var isMatchMismatch = currentValue === (diagonalValue + matchOrMismatch);
        var isChangeToP = currentValue === verticalValue;
        var isChangeToQ = currentValue === horizontalValue;

        isMatchMismatch = isMatchMismatch || currentValue === 0 && up >= 0 && left >= 0;
        var isDeletion = currentValue === 0 && up >= 0;  // lower 0 -> cut away
        var isInsertion = currentValue === 0 && left >= 0;  // lower 0 -> cut away

        // add
        if (isMatchMismatch)
            neighboured.push(new bases.alignment.Vector(up, left));

        if (isChangeToP)
            neighboured.push(bases.alignment.create(new bases.alignment.Vector(position.i, position.j), MATRICES.VERTICAL));

        if (isChangeToQ)
            neighboured.push(bases.alignment.create(new bases.alignment.Vector(position.i, position.j), MATRICES.HORIZONTAL));

        if (isInsertion)
            neighboured.push(new bases.alignment.Vector(position.i, left));

        if (isDeletion)
            neighboured.push(new bases.alignment.Vector(up, position.j));

        if (!(isMatchMismatch || isChangeToP || isChangeToQ || isInsertion || isDeletion)
            && (position.i !== 0 || position.j !== 0))
            neighboured.push(new bases.alignment.Vector(0, 0));

        return neighboured;
    }

    /**
     * Creates the alignments.
     * @augments Alignment.createAlignments()
     */
    function createAlignments() {
        alignmentInstance.setIO(inputData, outputData);
        alignmentInstance.createAlignments();
        outputData = alignmentInstance.getOutput();
    }

    /**
     * Returns all algorithm output.
     * @return {Object} - Contains all output of the algorithm.
     */
    function getOutput() {
        return outputData;
    }

    /**
     * Sets the input and output of an algorithm.
     * @param input {Object} - Contains all input data.
     * @param output {Object} - Contains all output data.
     */
    function setIO(input, output) {
        inputData = input;
        outputData = output;
    }

    /**
     * Returns the superclass instance.
     * @return {Object} - Superclass instance.
     */
    function getSuperclass() {
        return alignmentInstance;
    }
}());