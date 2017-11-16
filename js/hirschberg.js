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
        var recursionData = [];
        var currentData = [];

        // initialize
        var input = {};
        initializeInput(input);
        hirschbergInstance.numberOfIterations = 0;
        outputData.maxNumberIterations = false;
        computeAllRecursionData(input, currentData, recursionData);

        return [inputData, outputData];
    }

    /**
     * Initializes the given input with the read in inputData.
     * @param input {Object} - The input which has to be initialized.
     */
    function initializeInput(input) {
        input.sequenceA =  inputData.sequenceA;
        input.sequenceB = inputData.sequenceB;

        input.calculationType = inputData.calculationType;

        input.matrixHeight = inputData.sequenceA.length + 1;
        input.matrixWidth = inputData.sequenceB.length + 1;
    }

    /**
     * Executes a recursive
     * deep-first-search (deleting last found path from memory)
     * to find all possible sub-matrices.
     * @param input {Object} - The initialized Needleman-Wunsch input structure.
     * @param currentData {Object} - Stores the data from the current recursion. At the beginning it is empty.
     * @param recursionData {Object} - Stores the data from all iterations.
     * @see: Restricted to one path for better runtime! So, first founded minimum chosen for splitting.
     */
    function computeAllRecursionData(input, currentData, recursionData) {
        // [1] find trace-cell
        var forwardMatrix = computeForwardSequenceMatrix(input);
        var backwardMatrix = computeBackwardSequenceMatrix(input);

        var forwardRow = forwardMatrix[Math.ceil(input.matrixHeight / 2)];
        var backwardRow = forwardMatrix[Math.ceil(input.matrixHeight / 2)];
        var mirroredBackwardRow = backwardRow.splice().reverse();  // create a new mirrored row

        var addedRows = addRows(forwardRow, backwardRow);
        var minimumColumnPosJ = findMinimum(addedRows);

        // [2] divide and conquer
    }

    /**
     * Computes the Needleman-Wunsch matrix with the input sequences in usual order.
     * @param input {Object} - The initialized Needleman-Wunsch input structure.
     * @return {Object} - Output data of Needleman-Wunsch.
     */
    function computeForwardSequenceMatrix(input) {
    }

    /**
     * Computes the Needleman-Wunsch matrix with the input sequences in reversed order.
     * @param input {Object} - The initialized Needleman-Wunsch input structure.
     * @return {Object} - Output data of Needleman-Wunsch.
     */
    function computeBackwardSequenceMatrix(input) {
    }

    /**
     * Adds two rows (arrays) like vectors and returns a new row (array).
     * @param row1 {Array} - The array of numbers.
     * @param row2 {Array} - The array of numbers.
     * @return {Array} - The sum of two rows.
     */
    function addRows(row1, row2) {
        return [];
    }

    /**
     * Returns the position of the first minimum.
     * @param row {Array} - The array in which it is searched for the minimum.
     * @return {number} - The first minimum.
     */
    function findMinimum(row) {
        return 0;
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