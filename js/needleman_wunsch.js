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
    if (loaded === ALGORITHMS.NEEDLEMAN_WUNSCH) {  // to avoid self execution on a script import
        needlemanWunsch.startNeedlemanWunsch();
        loaded = ALGORITHMS.NONE;
    }
});

(function () {  // namespace
    // public methods
    namespace("needlemanWunsch", startNeedlemanWunsch, NeedlemanWunsch);

    // instances
    var alignmentInstance;
    var needlemanWunschInstance;

    /**
     * Function managing objects.
     */
    function startNeedlemanWunsch() {
        var linearAlignmentInterface = new interfaces.linearAlignmentInterface.LinearAlignmentInterface();
        linearAlignmentInterface.startLinearAlignmentAlgorithm(NeedlemanWunsch, ALGORITHMS.NEEDLEMAN_WUNSCH);
    }

    /*---- ALGORITHM ----*/
    /**
     * Computes the optimal, global alignment.
     * @constructor
     * @augments Alignment
     * @see https://doi.org/10.1016/0022-2836(70)90057-4
     *
     * Needleman, Saul B., and Christian D. Wunsch.
     * "A general method applicable to the search for similarities in the amino acid sequence of two proteins."
     * Journal of molecular biology 48.3 (1970): 443-453.
     */
    function NeedlemanWunsch() {
        needlemanWunschInstance = this;

        // variables
        this.type = ALGORITHMS.NEEDLEMAN_WUNSCH;
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

        // initialize left upper corner
        outputData.matrix[0][0] = 0;

        // initialize left matrix border
        for (var i = 1; i < inputData.matrixHeight; i++)
            outputData.matrix[i][0] = outputData.matrix[i - 1][0] + inputData.deletion;

        // initialize upper matrix border
        for (var j = 1; j < inputData.matrixWidth; j++)
            outputData.matrix[0][j] = outputData.matrix[0][j - 1] + inputData.insertion;
    }

    /**
     * Computes the matrix by using the recursion function and the score.
     * @override Alignment.computeMatrixAndScore()
     */
    function computeMatrixAndScore() {
        var inputData = alignmentInstance.getInput();
        var outputData = alignmentInstance.getOutput();

        // going through every matrix cell
        for (var i = 1; i < inputData.matrixHeight; i++) {
            var aChar = inputData.sequenceA[i - 1];

            for (var j = 1; j < inputData.matrixWidth; j++) {
                var bChar = inputData.sequenceB[j - 1];

                outputData.matrix[i][j] = alignmentInstance.recursionFunction(aChar, bChar, i, j);
            }
        }

        // score is stored in the right bottom cell
        outputData.score = outputData.matrix[inputData.matrixHeight - 1][inputData.matrixWidth - 1];
    }

    /**
     * Computing maximum or minimum of the three input values to compute the cell score.
     * If the type of calculation is similarity,
     * the maximum will be computed and else the minimum.
     * @param diagonalValue {number} - First input value.
     * @param upValue {number} - Second input value.
     * @param leftValue {number} - Third input value.
     * @return {number} - Maximum or minimum.
     */
    function recursionFunction(diagonalValue, upValue, leftValue) {
        var inputData = alignmentInstance.getInput();

        var value;
        if (inputData.calculationType === ALIGNMENT_TYPES.DISTANCE)
            value = Math.min(diagonalValue, upValue, leftValue);
        else  // inputData.calculationType === ALIGNMENT_TYPES.SIMILARITY
            value = Math.max(diagonalValue, upValue, leftValue);

        return value;
    }

    /**
     * Initializes the traceback.
     * @override Alignment.computeTraceback()
     */
    function computeTraceback() {
        needlemanWunschInstance.numberOfTracebacks = 0;

        var inputData = alignmentInstance.getInput();
        var outputData = alignmentInstance.getOutput();

        var lowerRightCorner = new bases.alignment.Vector(inputData.matrixHeight - 1, inputData.matrixWidth - 1);

        outputData.moreTracebacks = false;
        outputData.tracebackPaths =
            alignmentInstance.getGlobalTraces([lowerRightCorner], inputData, outputData, -1, alignmentInstance.getNeighboured);
    }

    /**
     * Returns the superclass instance.
     * @return {Object} - Superclass instance.
     */
    function getSuperclass() {
        return alignmentInstance;
    }
}());