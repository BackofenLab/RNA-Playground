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
    if (loaded === ALGORITHMS.SMITH_WATERMAN) {  // to avoid self execution on a script import
        smithWaterman.startSmithWaterman();
        loaded = ALGORITHMS.NONE;
    }
});

(function () {  // namespace
    // public methods
    namespace("smithWaterman", startSmithWaterman, SmithWaterman);

    // instances
    var alignmentInstance;
    var smithWatermanInstance;

    /**
     * Function managing objects.
     */
    function startSmithWaterman() {
        var linearAlignmentInterface = new interfaces.linearAlignmentInterface.LinearAlignmentInterface();
        linearAlignmentInterface.startLinearAlignmentAlgorithm(SmithWaterman, ALGORITHMS.SMITH_WATERMAN);
    }

    /*---- ALGORITHM ----*/
    /**
     * Computes the optimal, local alignment.
     * @constructor
     * @augments Alignment
     * @see https://doi.org/10.1016/0022-2836(81)90087-5
     *
     * Smith, Temple F., and Michael S. Waterman.
     * "Identification of common molecular subsequences."
     * Journal of molecular biology 147.1 (1981): 195-197.
     */
    function SmithWaterman() {
        smithWatermanInstance = this;

        // variables
        this.type = ALGORITHMS.SMITH_WATERMAN;
        this.numberOfTracebacks = 0;

        // inheritance
        alignmentInstance = new bases.alignment.Alignment(this);

        this.setInput = alignmentInstance.setLinearAlignmentInput;
        this.compute = alignmentInstance.compute;
        this.getOutput = alignmentInstance.getOutput;

        this.setIO = alignmentInstance.setIO;

        // public class methods
        this.initializeMatrix = initializeMatrix;
        this.computeMatrixAndScore = computeMatrixAndScore;
        this.recursionFunction = recursionFunction;
        this.computeTraceback = computeTraceback;

        this.getSuperclass = getSuperclass;
    }

    // methods
    /**
     * Initializes the matrix.
     * @augments Alignment.initializeMatrix()
     */
    function initializeMatrix() {
        var inputData = alignmentInstance.getInput();
        var outputData = alignmentInstance.getOutput();

        // initialize
        outputData.matrix[0][0] = 0;

        // initialize left matrix border
        for (var i = 1; i < inputData.matrixHeight; i++)
            outputData.matrix[i][0] = 0;

        // initialize upper matrix border
        for (var j = 1; j < inputData.matrixWidth; j++)
            outputData.matrix[0][j] = 0;
    }

    /**
     * Computes the matrix by using the recursion function and the score.
     * @override Alignment.computeMatrixAndScore()
     */
    function computeMatrixAndScore() {
        var inputData = alignmentInstance.getInput();
        var outputData = alignmentInstance.getOutput();

        var maxValue = 0;

        // going through every matrix cell
        for (var i = 1; i < inputData.matrixHeight; i++) {
            var aChar = inputData.sequenceA[i - 1];

            for (var j = 1; j < inputData.matrixWidth; j++) {
                var bChar = inputData.sequenceB[j - 1];

                outputData.matrix[i][j] = alignmentInstance.recursionFunction(aChar, bChar, i, j);

                // storing maximum
                if (maxValue < outputData.matrix[i][j])
                    maxValue = outputData.matrix[i][j];
            }
        }

        // score is the maximum value
        outputData.score = maxValue;
    }

    /**
     * Computing maximum of the three input values and zero to compute the cell score.
     * @param diagonalValue {number} - First input value.
     * @param upValue {number} - Second input value.
     * @param leftValue {number} - Third input value.
     * @return {number} - Maximum.
     */
    function recursionFunction(diagonalValue, upValue, leftValue) {
        return Math.max(diagonalValue, upValue, leftValue, 0);
    }

    /**
     * Initializes the traceback and starts
     * the traceback from every maximum bigger 0
     * to get all tracebacks.
     * @override Alignment.computeTraceback()
     */
    function computeTraceback() {
        smithWatermanInstance.numberOfTracebacks = 0;

        var inputData = alignmentInstance.getInput();
        var outputData = alignmentInstance.getOutput();

        // computing all traceback start-positions
        var backtraceStarts = alignmentInstance.getAllMaxPositions(inputData, outputData);

        outputData.tracebackPaths = [];
        outputData.moreTracebacks = false;

        for (var i = 0; i < backtraceStarts.length; i++) {
            var tracebackPaths = alignmentInstance.getLocalTraces([backtraceStarts[i]], inputData, outputData, -1, alignmentInstance.getNeighboured);
            outputData.tracebackPaths = outputData.tracebackPaths.concat(tracebackPaths);
        }
    }

    /**
     * Returns the superclass instance.
     * @return {Object} - Superclass instance.
     */
    function getSuperclass() {
        return alignmentInstance;
    }
}());