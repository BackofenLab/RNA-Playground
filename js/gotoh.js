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
        gotoh.startGotoh();
});

(function () {  // namespace
    // public methods
    namespace("gotoh", startGotoh, Gotoh, getInput, setInput, compute, getTraces, getOutput, setIO);

    // instances
    var alignmentInstance;
    var gotohInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Function managing objects.
     */
    function startGotoh() {
        imports();

        var subadditiveAlignmentInterface = new interfaces.subadditiveAlignmentInterface.SubadditiveAlignmentInterface();
        subadditiveAlignmentInterface.startSubadditiveAlignmentAlgorithm(Gotoh, ALGORITHMS.GOTOH);
    }

    /**
     * Handling imports.
     */
    function imports() {
        $.getScript(PATHS.SUBADDITIVE_ALIGNMENT_INTERFACE);
    }

    /*---- ALGORITHM ----*/
    /**
     * Computes the optimal, global affine alignment.
     * @constructor
     * @augments Alignment
     */
    function Gotoh() {
        gotohInstance = this;

        // variables
        this.type = ALGORITHMS.GOTOH;
        this.numberOfTracebacks = 0;

        // inheritance
        alignmentInstance = new procedures.bases.alignment.Alignment(this);

        this.getInput = getInput;

        this.setInput = setInput;
        this.compute = compute;
        this.getOutput = getOutput;

        this.getTraces = getTraces;

        this.setIO = setIO;
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
        inputData.sequenceA = inputViewmodel.sequence1();
        inputData.sequenceB = inputViewmodel.sequence2();

        inputData.calculationType = inputViewmodel.calculation();

        inputData.baseCosts = inputViewmodel.baseCosts();
        inputData.enlargement = inputViewmodel.enlargement();
        inputData.match = inputViewmodel.match();
        inputData.mismatch = inputViewmodel.mismatch();

        inputData.matrixHeight = inputData.sequenceB.length + 1;
        inputData.matrixWidth = inputData.sequenceA.length + 1;
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
        createComputationMatrix();
        createHorizontalGapCostMatrix();
        createVerticalGapCostMatrix();
    }

    /**
     * Creates the default matrix without initializing it.
     */
    function createComputationMatrix() {
        outputData.matrix = new Array(inputData.matrixHeight);

        for (var i = 0; i < inputData.matrixHeight; i++)
            outputData.matrix[i] = new Array(inputData.matrixWidth);
    }

    /**
     * Creates the matrix for horizontal gaps without initializing it.
     */
    function createHorizontalGapCostMatrix() {
        outputData.horizontalGaps = new Array(inputData.matrixHeight);

        for (var i = 0; i < inputData.matrixHeight; i++)
            outputData.horizontalGaps[i] = new Array(inputData.matrixWidth);
    }

    /**
     * Creates the matrix for vertical gaps without initializing it.
     */
    function createVerticalGapCostMatrix() {
        outputData.verticalGaps = new Array(inputData.matrixHeight);

        for (var i = 0; i < inputData.matrixHeight; i++)
            outputData.verticalGaps[i] = new Array(inputData.matrixWidth);
    }

    /**
     * Initializes the default matrix and the gap matrices.
     */
    function initMatrices() {
        initComputationMatrix();
        initHorizontalGapCostMatrix();
        initVerticalGapCostMatrix();
    }

    /**
     * Initializes the default matrix.
     */
    function initComputationMatrix() {
        outputData.matrix[0][0] = 0;

        for (var i = 1; i < inputData.matrixHeight; i++)
            outputData.matrix[i][0] = inputData.baseCosts + i * inputData.enlargement;

        for (var j = 1; j < inputData.matrixWidth; j++)
            outputData.matrix[0][j] = inputData.baseCosts + j * inputData.enlargement
    }

    /**
     * Initializes the horizontal gap cost matrix.
     */
    function initHorizontalGapCostMatrix() {
        for (var i = 1; i < inputData.matrixHeight; i++)
            outputData.horizontalGaps[i][0] = (ALIGNMENT_TYPES.SIMILARITY === inputData.calculationType)
                ? Number.NEGATIVE_INFINITY : Number.POSITIVE_INFINITY;

        for (var j = 1; j < inputData.matrixWidth; j++)
            outputData.horizontalGaps[0][j] = SYMBOLS.GAP;
    }

    /**
     * Initializes the vertical gap cost matrix.
     */
    function initVerticalGapCostMatrix() {
        for (var i = 1; i < inputData.matrixHeight; i++)
            outputData.verticalGaps[i][0] = SYMBOLS.GAP;

        for (var j = 1; j < inputData.matrixWidth; j++)
            outputData.verticalGaps[0][j] = (ALIGNMENT_TYPES.SIMILARITY === inputData.calculationType)
                ? Number.NEGATIVE_INFINITY : Number.POSITIVE_INFINITY;
    }

    /**
     * Computes the matrix by using the recursion function and the score.
     */
    function computeMatricesAndScore() {
        // going through every matrix cell
        for (var i = 1; i < inputData.matrixHeight; i++) {
            var bChar = inputData.sequenceB[i - 1];

            for (var j = 1; j < inputData.matrixWidth; j++) {
                var aChar = inputData.sequenceA[j - 1];

                if (inputData.calculationType === ALIGNMENT_TYPES.DISTANCE)
                    outputData.matrix[i][j] = recursionFunction(aChar, bChar, i, j, Math.min);
                else  // inputData.calculationType === ALIGNMENT_TYPES.SIMILARITY
                    outputData.matrix[i][j] = recursionFunction(aChar, bChar, i, j, Math.max);
            }
        }

        // score is stored in the right bottom cell
        outputData.score = outputData.matrix[inputData.matrixHeight - 1][inputData.matrixWidth - 1];
    }

    /**
     * Computes the cell score.
     * @param aChar {string} - The current char from the first string.
     * @param bChar {string} - The current char from the second string.
     * @param i {number} - The current vertical position in the matrix.
     * @param j {number} - The current horizontal position in the matrix.
     * @param optimum {Function} - The function which should be used for optimization {Math.min, Math.max}.
     * @return {number} - The value for the cell at position (i,j).
     */
    function recursionFunction(aChar, bChar, i, j, optimum) {
        var matchOrMismatch = aChar === bChar ? inputData.match : inputData.mismatch;

        // gap recursion-functions
        outputData.horizontalGaps[i][j] = horizontalOptimum(optimum, i, j);
        outputData.verticalGaps[i][j] = verticalOptimum(optimum, i, j);

        // default matrix recursion function
        return optimum(
            outputData.horizontalGaps[i][j],
            outputData.matrix[i - 1][j - 1] + matchOrMismatch,
            outputData.verticalGaps[i][j]);
    }

    /**
     * Computes the cell score for the horizontal gap matrix.
     * @param optimum {Function} - The function which should be used for optimization {Math.min, Math.max}.
     * @param i {number} - The current vertical position in the matrix.
     * @param j {number} - The current horizontal position in the matrix.
     * @return {number} - The optimal value.
     */
    function horizontalOptimum(optimum, i, j) {
        return optimum(
            outputData.horizontalGaps[i][j - 1] + inputData.enlargement,
            outputData.matrix[i][j - 1] + inputData.baseCosts + inputData.enlargement);
    }

    /**
     * Computes the cell score for the vertical gap matrix.
     * @param optimum {Function} - The function which should be used for optimization {Math.min, Math.max}.
     * @param i {number} - The current vertical position in the matrix.
     * @param j {number} - The current horizontal position in the matrix.
     * @return {number} - The optimal value.
     */
    function verticalOptimum(optimum, i, j) {
        return optimum(
            outputData.verticalGaps[i - 1][j] + inputData.enlargement,
            outputData.matrix[i - 1][j] + inputData.baseCosts + inputData.enlargement);
    }

    /**
     * Initializes the traceback.
     * @override Alignment.computeTraceback()
     */
    function computeTraceback() {
        gotohInstance.numberOfTracebacks = 0;

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
        var neighboured = backtracking.getMultiNeighboured(currentPosition, inputData, outputData);

        // going through all successors (initial nodes of possible paths)
        for (var i = 0; i < neighboured.length; i++) {
            if ((neighboured[i].i === 0 && neighboured[i].j === 0)
                || (pathLength !== -1 && path.length >= pathLength)
                || outputData.moreTracebacks) {
                path.push(neighboured[i]);

                // path storage, if MAX_NUMBER_TRACEBACKS is not exceeded
                if (gotohInstance.numberOfTracebacks < MAX_NUMBER_TRACEBACKS) {
                    paths.push(path.slice());  // creating a shallow copy
                    gotohInstance.numberOfTracebacks++;
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
}());