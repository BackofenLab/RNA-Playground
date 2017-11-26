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

        debugger;
        if (input.sequenceAPositions.length > 1 && input.sequenceBPositions.length > 1) {
            computeAllRecursionData(input, [HIRSCHBERG_UPPER_NODE]);
            processDiscoveredTracecells();
            createAlignments();
        }

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

        outputData.forwardMatrices = [];
        outputData.backwardMatrices = [];

        outputData.forwardRows = [];
        outputData.mirroredBackwardRows = [];
        outputData.addedRows = [];
        outputData.relativeSplittingPoint = [];  // stores local [i=floor(matrix.Height/2), j]
        outputData.allMinimaPosJ = [];
        outputData.tracecellLines = [];
        outputData.globalPositionsI = [];

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

        input.calculationType = ALIGNMENT_DEFAULTS.CALCULATION_HIRSCHBERG;

        input.deletion = inputData.deletion;
        input.insertion = inputData.insertion;
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
     * @see: Restricted to one path for better runtime! So, first founded minimum chosen for splitting.
     */
    function computeAllRecursionData(input, recursionNumbers) {
        if (input.sequenceAPositions.length <= 1)
            return;

        // [1] find trace-cell
        var forwardMatrix = computeForwardSequenceMatrix(input);
        var backwardMatrix = computeBackwardSequenceMatrix(shallowCopy(input));  // shallow copy, because else reversed strings are saved

        var splittingPosI = Math.floor(input.sequenceAPositions.length / 2);

        var forwardRow = forwardMatrix[splittingPosI];
        var backwardRow = backwardMatrix[(backwardMatrix.length - 1) - splittingPosI];
        var mirroredBackwardRow = backwardRow.slice().reverse();  // create a new mirrored row

        var sumRow = addRows(forwardRow, mirroredBackwardRow);

        var allMinima = findAllMinima(sumRow);
        var splittingPosJ = allMinima[0];

        createDataCopy(input, forwardMatrix, backwardMatrix, forwardRow, mirroredBackwardRow, sumRow,
            splittingPosI, splittingPosJ, allMinima, recursionNumbers);

        // [2] divide and conquer
        recursionNumbers.push(HIRSCHBERG_UPPER_NODE);
        computeAllRecursionData(initializedUpperMatrixInput(shallowCopy(input), splittingPosI, splittingPosJ), recursionNumbers);
        recursionNumbers.pop();

        recursionNumbers.push(HIRSCHBERG_LOWER_NODE);
        computeAllRecursionData(initializedLowerMatrixInput(shallowCopy(input), splittingPosI, splittingPosJ), recursionNumbers);
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
     * Returns all minimum positions of the given row.
     * @param row {Array} - The array in which it is searched for the minimum.
     * @return {Array} - The minimum positions.
     */
    function findAllMinima(row) {
        var minimumValue = findMinima(row);
        var minimumPositions = [];

        for (var j = 0; j < row.length; j++) {
            if (row[j] === minimumValue)
                minimumPositions.push(j);
        }

        return minimumPositions;
    }

    /**
     * Returns the minimum value of the given row.
     * @param row {Array} - The array in which it is searched for the minimum.
     * @return {number} - The first minimum.
     */
    function findMinima(row) {
        var minimumValue = Number.POSITIVE_INFINITY;

        var currentValue;

        for (var i = 0; i < row.length; i++) {
            currentValue = row[i];

            if (currentValue < minimumValue)
                minimumValue = currentValue;
        }

        return minimumValue;
    }

    /**
     * Creates a copy of the given recursion data to display it later on.
     * @param input {Object} - The input with which Needleman-Wunsch was executed.
     * @param forwardMatrix {Array} - The matrix computed with Needleman-Wunsch and ordinary order of sequences.
     * @param backwardMatrix {Array} - The matrix computed with Needleman-Wunsch and reversed order of sequences.
     * @param forwardRow {Array} - The row computed with Needleman-Wunsch and ordinary order of sequences.
     * @param mirroredBackwardRow {Array} - The reversed row computed with Needleman-Wunsch and reversed order of sequences.
     * @param sumRow {Array} - The sum array of forwardRow and mirroredBackwardRow.
     * @param splittingPosI {number} - The i-position on which the algorithm does a split.
     * @param splittingPosJ {number} - The j-position on which the algorithm does a split.
     * @param allMinima {Array} - The relative j-positions.
     * @param recursionNumbers {Array} - Contains information about recursion (i.e. [1, 1, 2] is upper, upper, lower).
     */
    function createDataCopy(input, forwardMatrix, backwardMatrix, forwardRow, mirroredBackwardRow, sumRow,
                            splittingPosI, splittingPosJ, allMinima, recursionNumbers) {
        outputData.firstSequences.push(input.sequenceA);
        outputData.secondSequences.push(input.sequenceB);
        outputData.firstSequencePositions.push(input.sequenceAPositions);
        outputData.secondSequencePositions.push(input.sequenceBPositions);

        outputData.forwardMatrices.push(forwardMatrix);
        outputData.backwardMatrices.push(backwardMatrix);

        outputData.forwardRows.push(forwardRow);
        outputData.mirroredBackwardRows.push(mirroredBackwardRow);
        outputData.addedRows.push(sumRow);

        outputData.relativeSplittingPoint.push([splittingPosI, splittingPosJ]);
        outputData.allMinimaPosJ.push(allMinima);

        outputData.recursionNumbersContainer.push(recursionNumbers.slice());

        var currentRound = outputData.recursionNumbersContainer.length;
        var globalPosI = outputData.firstSequencePositions[currentRound-1][outputData.relativeSplittingPoint[currentRound-1][0]-1];
        var globalMinimaPositions = getGlobalRowTracecells(allMinima, currentRound, globalPosI);

        outputData.tracecellLines.push(globalMinimaPositions);
        outputData.globalPositionsI.push(globalMinimaPositions[0].i);
    }

    /**
     * Returns the global minima positions from a row.
     * @param allMinima {Array} - The relative j-positions.
     * @param currentRound {number} - The current recursion round.
     * @param globalPosI {number} - The global i-positions.
     * @param {Array} - The global positions with minima from a line.
     */
    function getGlobalRowTracecells(allMinima, currentRound, globalPosI) {
        var lineTracecells = [];

        for (var j = 0; j < allMinima.length; j++) {
            var globalPosJ = outputData.secondSequencePositions[currentRound-1][allMinima[j]-1];

            if (globalPosJ === undefined) {
                var firstDefinedPosition = outputData.secondSequencePositions[currentRound-1][0];
                globalPosJ = firstDefinedPosition-1;
            }

            lineTracecells.push(new bases.alignment.Vector(globalPosI, globalPosJ));
        }

        return lineTracecells;
    }

    /**
     * Initializes the input of the next recursion for the upper created matrix.
     * @param input {Object} - The input used to execute the algorithm.
     * @param minimumRowPosI {number} - The i-position on which the algorithm does a split.
     * @param minimumColumnPosJ {number} - The j-position on which the algorithm does a split.
     */
    function initializedUpperMatrixInput(input, minimumRowPosI, minimumColumnPosJ) {
        input.sequenceA = input.sequenceA.split(SYMBOLS.EMPTY).slice(0, minimumRowPosI).join(SYMBOLS.EMPTY);
        input.sequenceB = input.sequenceB.split(SYMBOLS.EMPTY).slice(0, minimumColumnPosJ).join(SYMBOLS.EMPTY);

        input.sequenceAPositions = input.sequenceAPositions.slice(0, minimumRowPosI);
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
        input.sequenceB = input.sequenceB.split(SYMBOLS.EMPTY).slice(minimumColumnPosJ).join(SYMBOLS.EMPTY);

        input.sequenceAPositions = input.sequenceAPositions.slice(minimumRowPosI);
        input.sequenceBPositions = input.sequenceBPositions.slice(minimumColumnPosJ);

        input.matrixHeight = input.sequenceA.length + 1;
        input.matrixWidth = input.sequenceB.length + 1;

        return input;
    }

    /**
     * Processes the tracecell-lines to get back one traceback.
     */
    function processDiscoveredTracecells() {
        addEndings();
        reorderTraceCellLines();  // for faster access during path creation

        var lowerRightCorner = new bases.alignment.Vector(inputData.matrixHeight - 1, inputData.matrixWidth - 1);
        outputData.tracebackPaths = computeTraceback([lowerRightCorner]);
    }

    /**
     * The start-row and the end-row are not contained in the tracecell-lines
     * and have to be added to create a path.
     */
    function addEndings() {
        var matrix = outputData.forwardMatrices[0];
        var reversedStringsMatrix = outputData.backwardMatrices[0];

        addLine(matrix, reversedStringsMatrix, 0);  // first line
        addLine(matrix, reversedStringsMatrix, inputData.sequenceA.length);  // last line
    }

    /**
     * Adds a line to tracecell-lines.
     * @param matrix {Array} - The matrix for the strings in right order.
     * @param reversedStringsMatrix {Array} - The matrix for the reversed strings.
     * @param row {number} - The row which ahould be added.
     */
    function addLine(matrix, reversedStringsMatrix, row) {
        var allMinima = getMinimaPositions(matrix, reversedStringsMatrix, row);

        outputData.tracecellLines.push(allMinima);
        outputData.globalPositionsI.push(row);
    }

    /**
     * Returns all horizontal positions with a minima from a given row.
     * @param matrix {Array} - The matrix for the strings in right order.
     * @param reversedStringsMatrix {Array} - The matrix for the reversed strings.
     * @param posI {number} - The vertical position from which you want all minima.
     * @return {Array} - All minima.
     */
    function getMinimaPositions(matrix, reversedStringsMatrix, posI) {
        var row = matrix[posI];
        var backwardRow = reversedStringsMatrix[(reversedStringsMatrix.length - 1) - posI];
        var reversedBackwardRow = backwardRow.slice().reverse();

        var sumRow = addRows(row, reversedBackwardRow);

        if (posI === 0)
            return getGlobalRowTracecells(findAllMinima(sumRow), 1, posI);
        return getGlobalRowTracecells(findAllMinima(sumRow), 1, posI);
    }

    /**
     * The tracecell-lines are not in the right order and have to be reordered.
     */
    function reorderTraceCellLines() {
        var switches = [];

        outputData.globalPositionsI.sort(function (a, b) {
            switches.push(a - b);
            return switches[switches.length-1];
        });

        var i = 0;

        outputData.tracecellLines.sort(function (a,b) {
            return switches[i++];
        });
    }

    /**
     * Afterwards we have to go from bottom bottom right to top left
     * of an imaginary matrix and select always one neighbour cell
     * to get back a traceback and such way find an alignment.
     */
    function computeTraceback(path) {
        var paths = [];
        globalTraceback(paths, path, getNeighboured);
        return paths;
    }

    /**
     * Returns one neighbour to which you can go from the current cell position used as input.
     * @param position {Object} - Current cell position in matrix.
     * @return {Array} - Contains neighboured positions as Vector-objects.
     */
    function getNeighboured(position) {
        var neighboured = [];

        // retrieve neighbours (order is important)
        var horizontalTraceCell = getNextHorizontalTraceCell(position);
        var matchTraceCell = getMatchCell(position);
        var verticalTraceCell = getNextVerticalTraceCell(position);

        // add
        if (horizontalTraceCell !== undefined)
            neighboured.push(horizontalTraceCell);
        else if (matchTraceCell !== undefined)
            neighboured.push(matchTraceCell);
        else
            neighboured.push(verticalTraceCell);

        return neighboured;
    }

    /**
     * Return the most left trace cell.
     * @param position {Object} - The current vector position.
     * @return {Object} - The most left trace cell up from the given position.
     */
    function getNextHorizontalTraceCell(position) {
        var tracecellLine = outputData.tracecellLines[position.i];

        if (tracecellLine !== undefined) {
            var nextPosition = tracecellLine[0];

            return nextPosition.j !== position.j ? nextPosition : undefined;
        }

        return undefined;
    }

    /**
     * Looks in the next line above for a match.
     * @param position {Object} - The current vector position.
     * @return {Object} - The most left trace cell up from the given position.
     */
    function getMatchCell(position) {
        var left = position.j - 1;
        var up = position.i - 1;

        var tracecellLine = outputData.tracecellLines[up];

        if (tracecellLine !== undefined ) {
            var horizontalPos = -1;  // in tracecell-lines, not in the imaginary matrix

            for (var j = 0; j < tracecellLine.length; j++) {
                if (tracecellLine[j].j === left) {  // test for match or cell above
                    horizontalPos = j;
                    break;
                }
            }

            return horizontalPos !== -1 ? tracecellLine[horizontalPos] : undefined;
        }

        return undefined;
    }

    /**
     * Return the most top trace cell.
     * @param position {Object} - The current vector position.
     * @return {Object} - The most left trace cell up from the given position.
     */
    function getNextVerticalTraceCell(position) {
        debugger;
        var tracecelLines = outputData.tracecellLines;

        var verticalPos = -1;  // in tracecell-lines, not in the imaginary matrix
        var horizontalPos = -1;  // in tracecell-lines, not in the imaginary matrix

        for (var i = 0; i < tracecelLines.length; i++) {
            var tracecellLine = tracecelLines[i];
            var mostRightCell = tracecellLine[tracecellLine.length-1];

            if (mostRightCell.j === position.j) {
                verticalPos = i;
                horizontalPos = tracecellLine.length-1;  // in tracecell-lines, not in the imaginary matrix
                break;
            }
        }

        return verticalPos !== -1 && verticalPos !== position.i ? tracecelLines[verticalPos][horizontalPos] : undefined;
    }

    /**
     * Computing one global traceback.
     * with special stop criteria on the matrix cells as path-nodes.
     * @param paths {Array} - Array of paths.
     * @param path {Array} - Array containing the first vector element from which on you want find a path.
     * @param neighbourFunction {Function} - The function which have to be used to retrieve neighbours.
     * @see It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function globalTraceback(paths, path, neighbourFunction) {
        var currentPosition = path[path.length - 1];
        var neighboured = neighbourFunction(currentPosition);

        if ((neighboured[0].i === 0 && neighboured[0].j === 0)) {  // stop criteria checks
            path.push(neighboured[0]);
            paths.push(path.slice());  // creating a shallow copy
            path.pop();
        } else {  // executing procedure with a successor
            path.push(neighboured[0]);
            globalTraceback(paths, path, neighbourFunction);
            path.pop();
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

    /**
     * Returns the superclass instance.
     * @return {Object} - Superclass instance.
     */
    function getSuperclass() {
        return multiSequenceAlignmentInstance;
    }
}());