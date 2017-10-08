/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("procedures.bases.alignment",
        Alignment, getInput, setInput, compute, recursionFunction, createAlignments, getOutput, setIO);

    // instances
    var childInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Contains functions to compute optimal alignments.
     * It is used by algorithms global and local alignment algorithms as superclass.
     * @param child - The child algorithm which inherits from this class.
     * @constructor
     */
    function Alignment(child) {
        childInstance = child;

        // public methods (linking)
        this.getInput = getInput;

        this.setInput = setInput;
        this.compute = compute;
        this.recursionFunction = recursionFunction;
        this.createAlignments = createAlignments;
        this.getOutput = getOutput;

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
     * @param inputViewmodel {InputViewmodel}
     * - The InputViewmodel of an appropriate algorithm (Needleman-Wunsch, Smith-Waterman).
     */
    function setInput(inputViewmodel) {
        inputData.sequenceA = inputViewmodel.sequence1();
        inputData.sequenceB = inputViewmodel.sequence2();

        inputData.calculationType = inputViewmodel.calculation();

        inputData.deletion = inputViewmodel.deletion();
        inputData.insertion = inputViewmodel.insertion();
        inputData.match = inputViewmodel.match();
        inputData.mismatch = inputViewmodel.mismatch();

        inputData.matrixHeight = inputData.sequenceB.length + 1;
        inputData.matrixWidth = inputData.sequenceA.length + 1;
    }

    /**
     * Starts the computation.
     */
    function compute() {
        initializeMatrix();
        computeMatrixAndScore();
        computeTraceback();
        createAlignments();
        return [inputData, outputData];
    }

    /**
     * Initializes and creates the matrix.
     */
    function initializeMatrix() {
        createMatrix();
        childInstance.initializeMatrix();
    }

    /**
     * Creates the matrix without initializing them.
     */
    function createMatrix() {
        outputData.matrix = new Array(inputData.matrixHeight);

        for (var i = 0; i < inputData.matrixHeight; i++)
            outputData.matrix[i] = new Array(inputData.matrixWidth);
    }

    /**
     * Computes the matrix by using the recursion function and the score.
     * @abstract
     */
    function computeMatrixAndScore() {
        childInstance.computeMatrixAndScore();
    }

    /**
     * Computes the cell score.
     * @param aChar {string} - The current char from the first string.
     * @param bChar {string} - The current char from the second string.
     * @param i {number} - The current vertical position in the matrix.
     * @param j {number} - The current horizontal position in the matrix.
     * @return {number} - The value for the cell at position (i,j).
     */
    function recursionFunction(aChar, bChar, i, j) {
        var matchOrMismatch = aChar === bChar ? inputData.match : inputData.mismatch;

        var diagonalValue = outputData.matrix[i - 1][j - 1] + matchOrMismatch;
        var upValue = outputData.matrix[i - 1][j] + inputData.deletion;
        var leftValue = outputData.matrix[i][j - 1] + inputData.insertion;

        return childInstance.recursionFunction(diagonalValue, upValue, leftValue);
    }

    /**
     * Initializes the traceback.
     * @abstract
     */
    function computeTraceback() {
        childInstance.computeTraceback();
    }

    /**
     * Creates the alignments.
     */
    function createAlignments() {
        outputData.alignments = [];
        var numTracebacks = outputData.tracebackPaths.length;

        for (var i = 0; i < numTracebacks; i++) {
            var alignment = createAlignment(outputData.tracebackPaths[i]);
            outputData.alignments.push(alignment);
        }
    }

    /**
     * Creates an alignment by going through the path array of vectors.
     * @param path {Array} - Path of vectors which is used to create the alignment.
     * @return {[alignedSequenceA, matchOrMismatchString, alignedSequenceB]} - The pair of strings which have to be displayed.
     * @see: It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function createAlignment(path) {
        path.reverse();  // allows more intuitive calculations from left to right

        var alignedSequenceA = SYMBOLS.EMPTY;
        var matchOrMismatchString = SYMBOLS.EMPTY;
        var alignedSequenceB = SYMBOLS.EMPTY;

        var currentPositionA = path[0].j;
        var currentPositionB = path[0].i;

        // going through each element of the path and look on the differences between vectors
        // to find out the type of difference vector (arrow)
        for (var i = 1; i < path.length; i++) {
            var aChar = inputData.sequenceA[currentPositionA];
            var bChar = inputData.sequenceB[currentPositionB];

            if (path[i].i - path[i - 1].i > 0 && path[i].j - path[i - 1].j > 0) {  // diagonal case
                alignedSequenceA += aChar;
                matchOrMismatchString += aChar === bChar ? SYMBOLS.STAR : SYMBOLS.VERTICAL_BAR;
                alignedSequenceB += bChar;

                currentPositionA++;
                currentPositionB++;
            } else if (path[i].j - path[i - 1].j > 0) {  // horizontal case
                alignedSequenceA += aChar;
                matchOrMismatchString += SYMBOLS.SPACE;
                alignedSequenceB += SYMBOLS.GAP;

                currentPositionA++;
            } else if (path[i].i - path[i-1].i > 0) {  // vertical case
                // Hint: for Gotoh really "else if" is needed because you can switch between matrices
                alignedSequenceA += SYMBOLS.GAP;
                matchOrMismatchString += SYMBOLS.SPACE;
                alignedSequenceB += bChar;

                currentPositionB++;
            }
        }

        return [alignedSequenceA, matchOrMismatchString, alignedSequenceB];
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
