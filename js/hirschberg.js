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

        // to compute the global row
        var currentGlobalRow = Math.ceil(input.sequenceA.length / 2);
        computeAllRecursionData(input, [1], [currentGlobalRow]);

        debugger;
        return [inputData, outputData];
    }

    /**
     * Creates the input structs used to visualize algorithm.
     */
    function initializeStructs() {
        outputData.firstSequences = [];
        outputData.secondSequences = [];
        outputData.firstSequencePositions = [];
        outputData.secondSequencePositions = [];

        outputData.forwardRows = [];
        outputData.mirroredBackwardRows = [];
        outputData.addedRows = [];
        outputData.minimum = [];  // stores [i=ceil(matrix.Height/2), j]
        outputData.currentGlobalRow = [];  // the
        outputData.recursionNumbersContainer = [];  // stores the position within the computation tree
    }

    /**
     * Initializes the given input with the read in inputData.
     * @param input {Object} - The input which has to be initialized.
     */
    function initializeInput(input) {
        input.sequenceA = inputData.sequenceA;
        input.sequenceB = inputData.sequenceB;

        input.matrixHeight = inputData.sequenceA.length + 1;
        input.matrixWidth = inputData.sequenceB.length + 1;

        input.sequenceAPositions = createMatrixPositions(input.matrixHeight, 1);
        input.sequenceBPositions = createMatrixPositions(input.matrixWidth, 1);

        input.calculationType = inputData.calculationType;

        input.deletion = inputData.deletion;
        input.insertion = inputData.deletion;
        input.match = inputData.match;
        input.mismatch = inputData.mismatch;
    }

    /**
     * Creates an array of positions.
     * Hint: Gives more information about the position of a submatrix in the initial matrix.
     * @param numberPositions {number} - The number of positions to create.
     */
    function createMatrixPositions(numberPositions, start) {
        var positions = [];

        for (var i = start; i < numberPositions; i++) {
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
     * @param rowTree {Array} - Binary tree storing current row at its last position.
     * @see: Restricted to one path for better runtime! So, first founded minimum chosen for splitting.
     */
    function computeAllRecursionData(input, recursionNumbers, rowTree) {
        debugger;
        var isRightNode = recursionNumbers[recursionNumbers.length-1] === HIRSCHBERG_RIGHT_NODE;

        // [1] find trace-cell
        var forwardMatrix = computeForwardSequenceMatrix(input);
        var backwardMatrix = computeBackwardSequenceMatrix(shallowCopy(input));  // shallow copy, because else reversed strings are saved

        var minimumRowPosI = Math.ceil(input.sequenceAPositions.length / 2);

        var forwardRow = forwardMatrix[minimumRowPosI];
        var backwardRow = backwardMatrix[minimumRowPosI];
        var mirroredBackwardRow = backwardRow.slice().reverse();  // create a new mirrored row

        var sumRow = addRows(forwardRow, mirroredBackwardRow);

        var minimumColumnPosJ = findMinimum(input, sumRow, isRightNode);

        createDataCopy(input, forwardRow, mirroredBackwardRow, sumRow, minimumRowPosI, minimumColumnPosJ, recursionNumbers, rowTree);

        if (input.sequenceAPositions.length <= 1)
            return;

        // [2] divide and conquer
        recursionNumbers.push(1);
        rowTree.push(getNextRow(rowTree, false));
        computeAllRecursionData(initializedUpperMatrixInput(shallowCopy(input), minimumRowPosI, minimumColumnPosJ), recursionNumbers, rowTree);
        recursionNumbers.pop();
        rowTree.pop();

        recursionNumbers.push(2);
        rowTree.push(getNextRow(rowTree, true));
        computeAllRecursionData(initializedLowerMatrixInput(shallowCopy(input), minimumRowPosI, minimumColumnPosJ), recursionNumbers, rowTree);
        recursionNumbers.pop();
        rowTree.pop()
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
    function shallowCopy(input) {
        return jQuery.extend(false, {}, input);
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
     * Returns the first minimum position in the columns.
     * @param row {Array} - The array in which it is searched for the minimum.
     * @param right {boolean} - Tells if we have to search the minimum from right to left.
     * @return {number} - The first minimum.
     */
    function findMinimum(input, row, right) {
        var minimumValue = Number.POSITIVE_INFINITY;
        var minimumPosition = -1;

        if (right) {  // search minimum from right side
            for (var i = row.length - 1; i >= 0; i--) {
                var currentValue = row[i];

                if (currentValue < minimumValue) {
                    minimumValue = currentValue;
                    minimumPosition = i;
                }
            }
        } else {  // else from left side
            for (var i = 0; i < row.length; i++) {
                var currentValue = row[i];

                if (currentValue < minimumValue) {
                    minimumValue = currentValue;
                    minimumPosition = i;
                }
            }
        }

        return minimumPosition;
    }

    /**
     * Creates a copy of the given recursion data to display it later on.
     */
    function createDataCopy(input, forwardRow, mirroredBackwardRow, sumRow, minimumRowPosI, minimumColumnPosJ, recursionNumbers, rowTree) {
        outputData.firstSequences.push(input.sequenceA);
        outputData.secondSequences.push(input.sequenceB);
        outputData.firstSequencePositions.push(input.sequenceAPositions);
        outputData.secondSequencePositions.push(input.sequenceBPositions);

        outputData.forwardRows.push(forwardRow);
        outputData.mirroredBackwardRows.push(mirroredBackwardRow);
        outputData.addedRows.push(sumRow);

        outputData.minimum.push([minimumRowPosI, minimumColumnPosJ]);
        outputData.currentGlobalRow.push(rowTree[rowTree.length-1]);

        outputData.recursionNumbersContainer.push(recursionNumbers.slice());
    }

    /**
     * Returns the next row.
     * @param rowTree {Array} - Binary tree storing current row at its last position.
     * @param right {boolean} - Tells if we have to compute the row for a right or a left node.
     */
    function getNextRow(rowTree, right) {
        var lastRow;
        var nextRow;

        if (right)
            nextRow = rowTree[rowTree.length - 1] + 1
        else {
            lastRow = rowTree[rowTree.length - 1] - 1;
            nextRow = Math.ceil(lastRow / 2);
        }

        return nextRow;
    }

    /**
     * Initializes the input of the next recursion for the upper created matrix.
     * @param input {Object} - The input used to execute the algorithm.
     * @param minimumRowPosI {number} - The i-position on which the algorithm does a split.
     * @param minimumColumnPosJ {number} - The j-position on which the algorithm does a split.
     */
    function initializedUpperMatrixInput(input, minimumRowPosI, minimumColumnPosJ) {
        input.sequenceA = input.sequenceA.split(SYMBOLS.EMPTY).slice(0, minimumRowPosI - 1).join(SYMBOLS.EMPTY);
        input.sequenceB = input.sequenceB.split(SYMBOLS.EMPTY).slice(0, minimumColumnPosJ).join(SYMBOLS.EMPTY);

        input.sequenceAPositions = input.sequenceAPositions.slice(0, minimumRowPosI - 1);
        input.sequenceBPositions = input.sequenceBPositions.slice(0, minimumColumnPosJ);

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
        // Hint: There is one row and one column more than affiliated sequence size.
        input.sequenceA = input.sequenceA.split(SYMBOLS.EMPTY).slice(minimumRowPosI).join(SYMBOLS.EMPTY);
        input.sequenceB = input.sequenceB.split(SYMBOLS.EMPTY).slice(minimumColumnPosJ - 1).join(SYMBOLS.EMPTY);

        input.sequenceAPositions = input.sequenceAPositions.slice(minimumRowPosI);
        input.sequenceBPositions = input.sequenceBPositions.slice(minimumColumnPosJ - 1);

        input.matrixHeight = input.sequenceA.length + 1;
        input.matrixWidth = input.sequenceB.length + 1;

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