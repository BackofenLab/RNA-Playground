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
    if (document.title !== UNIT_TEST_WEBTITLE)  // to avoid the execution of the algorithm interfaces during a Unit-Test
        needlemanWunsch.startNeedlemanWunsch();
});

(function () {  // namespace
    // public methods
    namespace("needlemanWunsch", startNeedlemanWunsch, NeedlemanWunsch,
        setIO, compute, initializeMatrix, computeMatrixAndScore, recursionFunction, computeTraceback, getTraces);

    // instances
    var alignmentInstance;
    var needlemanWunschInstance;

    /**
     * Function managing objects.
     */
    function startNeedlemanWunsch() {
        imports();

        var alignmentInterface = new interfaces.alignmentInterface.AlignmentInterface();
        alignmentInterface.startAlignmentAlgorithm(NeedlemanWunsch, ALGORITHMS.NEEDLEMAN_WUNSCH);
    }

    /**
     * Handling imports.
     */
    function imports() {
        $.getScript(PATHS.ALIGNMENT_INTERFACE);
        //$.getScript(PATHS.ALIGNMENT);
    }

    /*---- ALGORITHM ----*/
    /**
     * Computes the optimal, global alignment.
     * @constructor
     * @augments Alignment
     */
    function NeedlemanWunsch() {
        needlemanWunschInstance = this;

        // variables
        this.type = ALGORITHMS.NEEDLEMAN_WUNSCH;
        this.numberOfTracebacks = 0;

        // inheritance
        alignmentInstance = new procedures.bases.alignment.Alignment(this);

        this.setInput = alignmentInstance.setInput;
        this.compute = alignmentInstance.compute;
        this.getOutput = alignmentInstance.getOutput;

        this.setIO = alignmentInstance.setIO;

        // public methods (linking)
        this.initializeMatrix = initializeMatrix;
        this.computeMatrixAndScore = computeMatrixAndScore;
        this.recursionFunction = recursionFunction;
        this.computeTraceback = computeTraceback;

        this.getTraces = getTraces;
    }

    // inheritance
    /**
     * Sets the algorithm input and output for calculation.
     * @param input {Object} - The input structure.
     * @param output {Object} - The output structure.
     * @augments Alignment.setIO(input, output)
     */
    function setIO(input, output) {
        alignmentInstance.setIO(input, output);
    }

    /**
     * Starts computation by starting the superclass computation.
     * @augments Alignment.compute()
     */
    function compute() {
        return alignmentInstance.compute();
    }

    // methods
    /**
     * Initializes the matrix.
     * @override Alignment.initializeMatrix()
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
            var bChar = inputData.sequenceB[i - 1];

            for (var j = 1; j < inputData.matrixWidth; j++) {
                var aChar = inputData.sequenceA[j - 1];

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

        var lowerRightCorner = new procedures.backtracking.Vector(inputData.matrixHeight - 1, inputData.matrixWidth - 1);
        
        outputData.moreTracebacks = false;
        outputData.tracebackPaths = getTraces([lowerRightCorner], inputData, outputData, -1);
    }

    /**
     * Gets tracebacks by starting the traceback procedure
     * with some path containing the first element of the path.
     * @param path {Array} - Array containing the first vector element from which on you want find a path.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     * @param pathLength {number} - Tells after how many edges the procedure should stop.
     * The value -1 indicates arbitrarily long paths.
     * @return {Array} - Array of paths.
     */
    function getTraces(path, inputData, outputData, pathLength) {
        var paths = [];
        var backtracking = new procedures.backtracking.Backtracking();
        traceback(backtracking, paths, path, inputData, outputData, pathLength);
        return paths;
    }

    /**
     * Computing the traceback and stops after it has found a constant number of tracebacks.
     * It sets a flag "moreTracebacks" in the "outputData", if it has stopped before computing all tracebacks.
     * The traceback algorithm executes a recursive,
     * modified deep-first-search (deleting last found path from memory)
     * with special stop criteria on the matrix cells as path-nodes.
     * @param backtracking {Object} - Allows to call up cell neighbours.
     * @param paths {Array} - Array of paths.
     * @param path {Array} - Array containing the first vector element from which on you want find a path.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     * @param pathLength {number} - Tells after how many edges the procedure should stop.
     * @see It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function traceback(backtracking, paths, path, inputData, outputData, pathLength) {
        var currentPosition = path[path.length - 1];
        var neighboured = backtracking.getNeighboured(currentPosition, inputData, outputData, needlemanWunschInstance);

        // going through all successors (initial nodes of possible paths)
        for (var i = 0; i < neighboured.length; i++) {
            if ((neighboured[i].i === 0 && neighboured[i].j === 0)  // stop criteria checks
                || (pathLength !== -1 && path.length >= pathLength)
                || outputData.moreTracebacks) {

                path.push(neighboured[i]);

                // path storage, if MAX_NUMBER_TRACEBACKS is not exceeded
                if (needlemanWunschInstance.numberOfTracebacks < MAX_NUMBER_TRACEBACKS) {
                    paths.push(path.slice());  // creating a shallow copy
                    needlemanWunschInstance.numberOfTracebacks++;
                } else
                    outputData.moreTracebacks = true;

                path.pop();
            } else {
                // executing procedure with a successor
                path.push(neighboured[i]);
                traceback(backtracking, paths, path, inputData, outputData, pathLength);
                path.pop();
            }
        }
    }
}());