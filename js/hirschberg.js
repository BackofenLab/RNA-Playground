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
    if (loaded === ALGORITHMS.HIRSCHBERG) {  // to avoid self execution on a script import
        hirschberg.startHirschberg();
        loaded = ALGORITHMS.NONE;
    }
});

(function () {  // namespace
    // public methods
    namespace("hirschberg", startHirschberg, Hirschberg, getInput, setInput, compute, getOutput, setIO, getSuperclass);

    // instances
    var alignmentInstance;
    var needlemanWunschInstance;

    var hirschbergInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Function managing objects.
     */
    function startHirschberg() {
        var linearAlignmentInterface = new interfaces.linearAlignmentInterface.LinearAlignmentInterface();
        linearAlignmentInterface.startLinearAlignmentAlgorithm(Hirschberg, ALGORITHMS.HIRSCHBERG);
    }

    /*---- ALGORITHM ----*/
    /**
     * Simulates the Hirschberg algorithm
     * with the Needleman-Wunsch algorithm.
     * @constructor
     * @augments Alignment
     * @see https://dl.acm.org/citation.cfm?doid=360825.360861
     *
     * Hirschberg, Daniel S.
     * "A linear space algorithm for computing maximal common subsequences."
     * Communications of the ACM 18.6 (1975): 341-343.
     */
    function Hirschberg() {
        hirschbergInstance = this;

        // variables
        this.type = ALGORITHMS.HIRSCHBERG;
        this.numberOfTracebacks = 0;

        // inheritance
        alignmentInstance = new bases.alignment.Alignment(this);
        needlemanWunschInstance = new needlemanWunsch.NeedlemanWunsch();

        this.setIO = alignmentInstance.setIO;

        // public class methods
        this.getInput = getInput;

        this.setInput = setInput;
        this.compute = compute;
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
        alignmentInstance.setLinearAlignmentInput(inputViewmodel);
    }

    /**
     * Starts the computation.
     */
    function compute() {
        initializeStructs();

        // initialize
        var input = {};
        initializeInput(input);

        hirschbergInstance.numberOfIterations = 0;
        outputData.maxNumberIterations = false;

        computeAllRecursionData(input, [1]);

        return [inputData, outputData];
    }

    /**
     * Creates the input structs used to visualize algorithm.
     */
    function initializeStructs() {
        outputData.firstSequences = [];
        outputData.secondSequences = [];
        outputData.definedVerticalMatrixPositions = [];
        outputData.definedHorizontalMatrixPositions = [];

        outputData.forwardRows = [];
        outputData.mirroredBackwardRows = [];
        outputData.addedRows = [];
        outputData.minimum = [];  // stores [i=ceil(matrix.Height/2) + c , j + d]
        outputData.recursionNumbersContainer = [];  // stores the position within the computation tree
    }

    /**
     * Initializes the given input with the read in inputData.
     * @param input {Object} - The input which has to be initialized.
     */
    function initializeInput(input) {
        input.sequenceA = inputData.sequenceA;
        input.sequenceB = inputData.sequenceB;

        input.sequenceAPositions = createCharacterPositions(input.sequenceA);
        input.sequenceBPositions = createCharacterPositions(input.sequenceB);

        input.calculationType = inputData.calculationType;

        input.deletion = inputData.deletion;
        input.insertion = inputData.deletion;
        input.match = inputData.match;
        input.mismatch = inputData.mismatch;

        input.matrixHeight = inputData.sequenceA.length + 1;
        input.matrixWidth = inputData.sequenceB.length + 1;
    }

    /**
     * Creates for each character a position in an array.
     * Hint: Gives more information about the position of a submatrix in the initial matrix.
     * @param sequence {string} - The sequence in which for every character a position is created.
     */
    function createCharacterPositions(sequence) {
        var positions = [];

        for (var i = 0; i < sequence.length; i++) {
            positions.push(i);
        }

        return positions;
    }

    /**
     * Executes a recursive
     * deep-first-search (deleting last found path from memory)
     * to find all possible submatrices.
     * @param input {Object} - The initialized Needleman-Wunsch input structure.
     * @param recursionNumbers {Array} - Contains information about recursion (i.e. [1, 1, 2] is upper, upper, lower).
     * @see: Restricted to one path for better runtime! So, first founded minimum chosen for splitting.
     */
    function computeAllRecursionData(input, recursionNumbers) {
        debugger;
        if (input.matrixHeight === 1)
            return;

        // [1] find trace-cell
        var forwardMatrix = computeForwardSequenceMatrix(input);
        var backwardMatrix = computeBackwardSequenceMatrix(deepCopy(input));  // deep copy, because else reversed strings are saved

        var minimumRowPosI = Math.ceil(input.matrixHeight / 2);

        var forwardRow = forwardMatrix[minimumRowPosI];
        var backwardRow = backwardMatrix[minimumRowPosI];
        var mirroredBackwardRow = backwardRow.slice().reverse();  // create a new mirrored row

        var sumRow = addRows(forwardRow, mirroredBackwardRow);

        var minimumColumnPosJ = findMinimum(sumRow);

        createDataCopy(input, forwardRow, mirroredBackwardRow, sumRow, minimumRowPosI, minimumColumnPosJ, recursionNumbers);

        // [2] divide and conquer
        recursionNumbers.push(1);
        computeAllRecursionData(initializedUpperMatrixInput(deepCopy(input), minimumRowPosI, minimumColumnPosJ), recursionNumbers);
        recursionNumbers.pop();

        recursionNumbers.push(2);
        computeAllRecursionData(initializedLowerMatrixInput(deepCopy(input), minimumRowPosI, minimumColumnPosJ), recursionNumbers);
        recursionNumbers.pop();
    }

    /**
     * Computes the Needleman-Wunsch matrix with the input sequences in usual order.
     * @param input {Object} - The initialized Needleman-Wunsch input structure.
     * @return {Object} - Output data of Needleman-Wunsch.
     */
    function computeForwardSequenceMatrix(input) {
        needlemanWunschInstance.setIO(input, {});
        return needlemanWunschInstance.compute()[1].matrix;
    }

    /**
     * Creates a deep copy of the struct used as input.
     * @param input {Object} - The object which should be copied.
     * @return {Object} - Copy of the Object.
     */
    function deepCopy(input) {
        return jQuery.extend(true, {}, input);
    }

    /**
     * Computes the Needleman-Wunsch matrix with the input sequences in reversed order.
     * @param input {Object} - The initialized Needleman-Wunsch input structure.
     * @return {Object} - Output data of Needleman-Wunsch.
     */
    function computeBackwardSequenceMatrix(input) {
        input.sequenceA = input.sequenceA.split(SYMBOLS.EMPTY).reverse().join(SYMBOLS.EMPTY);
        input.sequenceB = input.sequenceB.split(SYMBOLS.EMPTY).reverse().join(SYMBOLS.EMPTY);

        needlemanWunschInstance.setIO(input, {});
        return needlemanWunschInstance.compute()[1].matrix;
    }

    /**
     * Adds two rows (arrays) of same length like vectors and returns a new row (array).
     * @param row1 {Array} - The array of numbers.
     * @param row2 {Array} - The array of numbers.
     * @return {Array} - The sum of two rows.
     */
    function addRows(row1, row2) {
        var sumRow = [];

        for (var i = 0; i < row1.length; i++) {
            var row1Value = row1[i];
            var row2Value = row2[i];

            sumRow.push(row1Value + row2Value);
        }

        return sumRow;
    }

    /**
     * Returns the position of the first minimum.
     * @param row {Array} - The array in which it is searched for the minimum.
     * @return {number} - The first minimum.
     */
    function findMinimum(row) {
        var minimumValue = Number.POSITIVE_INFINITY;
        var minimumPosition = -1;

        for (var i = 0; i < row.length; i++) {
            var currentValue = row[i];

            if (currentValue < minimumValue) {
                minimumValue = currentValue;
                minimumPosition = i;
            }
        }

        return minimumPosition;
    }

    /**
     * Creates a copy of the given recursion data to display it later on.
     */
    function createDataCopy(input, forwardRow, mirroredBackwardRow, sumRow, minimumRowPosI, minimumColumnPosJ, recursionNumbers) {
        outputData.firstSequences.push(input.sequenceA);
        outputData.secondSequences.push(input.sequenceB);
        outputData.definedVerticalMatrixPositions.push(input.sequenceAPositions);
        outputData.definedHorizontalMatrixPositions.push(input.sequenceBPositions);

        outputData.forwardRows.push(forwardRow);
        outputData.mirroredBackwardRows.push(mirroredBackwardRow);
        outputData.addedRows.push(sumRow);
        outputData.minimum.push([minimumRowPosI, minimumColumnPosJ]);
        outputData.recursionNumbersContainer.push(recursionNumbers);
    }

    /**
     * Initializes the input of the next recursion for the upper created matrix.
     * @param input {Object} - The input used to execute the algorithm.
     * @param minimumRowPosI {number} - The i-position on which the algorithm does a split.
     * @param minimumColumnPosJ {number} - The j-position on which the algorithm does a split.
     */
    function initializedUpperMatrixInput(input, minimumRowPosI, minimumColumnPosJ) {
        /**
         * _ 0 1 2 ... j
         * 0            |
         * 1            |
         * .            |
         * .____________|
         * i
         */
        input.sequenceA = input.sequenceA.split(SYMBOLS.EMPTY).slice(0, minimumRowPosI).join(SYMBOLS.EMPTY);  // without i
        input.sequenceB = input.sequenceB.split(SYMBOLS.EMPTY).slice(0, minimumColumnPosJ + 1).join(SYMBOLS.EMPTY);  // with j

        input.sequenceAPositions = input.sequenceAPositions.slice(0, minimumRowPosI);  // without i
        input.sequenceBPositions = input.sequenceBPositions.slice(0, minimumColumnPosJ + 1);  // with j

        input.matrixHeight = input.sequenceA.length + 1;
        input.matrixWidth = input.sequenceB.length + 1;

        return input;
    }

    /**
     * Initializes the input of the next recursion for the lower created matrix.
     * @param input {Object} - The input used to execute the algorithm.
     * @param minimumRowPosI {number} - The i-position on which the algorithm does a split.
     * @param minimumColumnPosJ {number} - The j-position on which the algorithm does a split.
     */
    function initializedLowerMatrixInput(input, minimumRowPosI, minimumColumnPosJ) {
        /**
         * _ 0 1 2 ... j ... m
         * 0
         * 1
         * .
         * .
         * i          ________
         * i+1       |        |
         * .         |        |
         * .         |        |
         * n         |        |
         */

         input.sequenceA = input.sequenceA.split(SYMBOLS.EMPTY).slice(minimumRowPosI + 1).join(SYMBOLS.EMPTY);  // without i
         input.sequenceB = input.sequenceB.split(SYMBOLS.EMPTY).slice(minimumColumnPosJ).join(SYMBOLS.EMPTY);  // with j

         input.sequenceAPositions = input.sequenceAPositions.slice(minimumRowPosI + 1);  // without i
         input.sequenceBPositions = input.sequenceBPositions.slice(minimumColumnPosJ);  // with j

         input.matrixHeight = inputData.sequenceA.length + 1;
         input.matrixWidth = inputData.sequenceB.length + 1;

         return input;
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
        return multiSequenceAlignmentInstance;
    }
}());