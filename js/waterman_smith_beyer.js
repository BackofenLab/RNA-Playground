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
    if (loaded === ALGORITHMS.WATERMAN_SMITH_BEYER) {  // to avoid self execution on a script import
        watermanSmithBeyer.startWatermanSmithBeyer();
        loaded = ALGORITHMS.NONE;
    }
});

(function () {  // namespace
    // public methods
    namespace("watermanSmithBeyer", startWatermanSmithBeyer, WatermanSmithBeyer);

    // instances
    var alignmentInstance;
    var watermanSmithBeyerInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Function managing objects.
     */
    function startWatermanSmithBeyer() {
        var subadditiveAlignmentInterface = new interfaces.subadditiveAlignmentInterface.SubadditiveAlignmentInterface();
        subadditiveAlignmentInterface.startSubadditiveAlignmentAlgorithm(WatermanSmithBeyer, ALGORITHMS.WATERMAN_SMITH_BEYER);
    }

    /*---- ALGORITHM ----*/
    /**
     * Computes the optimal, global subadditive alignment.
     * @constructor
     * @see https://doi.org/10.1016/0001-8708(76)90202-4
     *
     * Waterman, Michael S., Temple F. Smith, and William A. Beyer.
     * "Some biological sequence metrics."
     * Advances in Mathematics 20.3 (1976): 367-387.
     */
    function WatermanSmithBeyer() {
        watermanSmithBeyerInstance = this;

        // variables
        this.type = ALGORITHMS.WATERMAN_SMITH_BEYER;
        this.numberOfTracebacks = 0;

        // instances
        alignmentInstance = new bases.alignment.Alignment(this);

        // public class methods
        this.getInput = getInput;

        this.setInput = setInput;
        this.compute = compute;
        this.gapFunction = gapFunction;
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
        inputData.subadditiveFunction = inputViewmodel.subadditiveFunction();
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
     * Initializes and creates the matrices.
     */
    function initializeMatrix() {
        createMatrix();
        initMatrix();
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
     * Initializes the matrix by distinguishing between three possible subadditive functions.
     */
    function initMatrix() {
        outputData.matrix[0][0] = 0;

        for (var i = 1; i < inputData.matrixHeight; i++)
            outputData.matrix[i][0] = gapFunction(i);

        for (var j = 1; j < inputData.matrixWidth; j++)
            outputData.matrix[0][j] = gapFunction(j);
    }

    /**
     * The gap function used to compute values.
     * @param k {number} - The integer for which the gap function-value should be computed.
     * @return {number} - The gap costs.
     */
    function gapFunction(k) {
        switch (inputData.subadditiveFunction) {
            case SUBADDITIVE_FUNCTIONS.AFFINE:
                return inputData.baseCosts + inputData.enlargement * k;
            case SUBADDITIVE_FUNCTIONS.LOGARITHMIC:
                return inputData.baseCosts + inputData.enlargement * Math.log(k);
            case SUBADDITIVE_FUNCTIONS.QUADRATIC:
                return inputData.baseCosts + inputData.enlargement * Math.pow(k, 2);
        }
    }

    /**
     * Computes the matrix by using the recursion function and the score.
     */
    function computeMatrixAndScore() {
        // going through every matrix cell
        for (var i = 1; i < inputData.matrixHeight; i++) {
            var aChar = inputData.sequenceA[i - 1];

            for (var j = 1; j < inputData.matrixWidth; j++) {
                var bChar = inputData.sequenceB[j - 1];

                if (inputData.calculationType === ALIGNMENT_TYPES.DISTANCE)
                    outputData.matrix[i][j] = recursionFunction(aChar, bChar, i, j, Math.min);
                else
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

        // recursion function
        return optimum(
            horizontalOptimum(optimum, i, j),
            outputData.matrix[i - 1][j - 1] + matchOrMismatch,
            verticalOptimum(optimum, i, j));
    }

    /**
     * Computes horizontal gap score.
     * @param optimum {Function} - The function which should be used for optimization {Math.min, Math.max}.
     * @param i {number} - The current vertical position in the matrix.
     * @param j {number} - The current horizontal position in the matrix.
     * @return {number} - The optimal value.
     */
    function horizontalOptimum(optimum, i, j) {
        var optimumValue;
        var value;

        if (optimum === Math.min) {
            optimumValue = Number.POSITIVE_INFINITY;

            for (var k = 1; k <= j; k++) {
                value = outputData.matrix[i][j - k] + gapFunction(k);

                if (value < optimumValue)
                    optimumValue = value;
            }
        }
        else {
            optimumValue = Number.NEGATIVE_INFINITY;

            for (var k = 1; k <= j; k++) {
                value = outputData.matrix[i][j - k] + gapFunction(k);

                if (value > optimumValue)
                    optimumValue = value;
            }
        }

        return optimumValue;
    }

    /**
     * Computes vertical gap score.
     * @param optimum {Function} - The function which should be used for optimization {Math.min, Math.max}.
     * @param i {number} - The current vertical position in the matrix.
     * @param j {number} - The current horizontal position in the matrix.
     * @return {number} - The optimal value.
     */
    function verticalOptimum(optimum, i, j) {
        var optimumValue;
        var value;

        if (optimum === Math.min) {
            optimumValue = Number.POSITIVE_INFINITY;

            for (var k = 1; k <= i; k++) {
                value = outputData.matrix[i - k][j] + gapFunction(k);

                if (value < optimumValue)
                    optimumValue = value;
            }
        }
        else {
            optimumValue = Number.NEGATIVE_INFINITY;

            for (var k = 1; k <= i; k++) {
                value = outputData.matrix[i - k][j] + gapFunction(k);

                if (value > optimumValue)
                    optimumValue = value;
            }
        }

        return optimumValue;
    }

    /**
     * Initializes the traceback.
     * @override Alignment.computeTraceback()
     */
    function computeTraceback() {
        watermanSmithBeyerInstance.numberOfTracebacks = 0;

        var lowerRightCorner = new bases.alignment.Vector(inputData.matrixHeight - 1, inputData.matrixWidth - 1);

        outputData.moreTracebacks = false;
        outputData.tracebackPaths =
            alignmentInstance.getGlobalTraces([lowerRightCorner], inputData, outputData, -1, getNeighboured);
    }

    /**
     * Returns the neighbours to which you can go from the current cell position used as input.
     * @param position {Object} - Current cell position in matrix.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     * @param algorithm {Object} - Contains an alignment algorithm.
     * @return {Array} - Contains neighboured positions as Vector-objects.
     */
    function getNeighboured(position, inputData, outputData, algorithm) {
        var neighboured = [];

        var left = position.j - 1;
        var up = position.i - 1;

        // retrieve values
        var aChar = left >= 0 ? inputData.sequenceB[left] : SYMBOLS.EMPTY;
        var bChar = up >= 0 ? inputData.sequenceA[up] : SYMBOLS.EMPTY;

        var currentValue = outputData.matrix[position.i][position.j];

        var matchOrMismatch = aChar === bChar ? inputData.match : inputData.mismatch;
        var horizontalK = searchHorizontalMatchPosition(algorithm, currentValue, position, outputData);
        var verticalK = searchVerticalMatchPosition(algorithm, currentValue, position, outputData);

        var diagonalValue = left >= 0 && up >= 0 ? outputData.matrix[up][left] : Number.NaN;
        var upValue = up >= 0 && position.j === 0 ? outputData.matrix[up][position.j] : Number.NaN;
        var leftValue = left >= 0 && position.i === 0 ? outputData.matrix[position.i][left] : Number.NaN;

        // check
        var isMatchMismatch = alignmentInstance.differenceLowerEpsilon(currentValue, diagonalValue + matchOrMismatch, EPSILON);
        var isHorizontal = !isNaN(horizontalK);  // if a position exists to which we can horizontally jump
        var isVertical = !isNaN(verticalK);  // if a position exists to which we can vertically jump

        var isDeletion = currentValue === upValue + inputData.enlargement;
        var isInsertion = currentValue === leftValue + inputData.enlargement;

        // add
        if (isMatchMismatch)
            neighboured.push(new bases.alignment.Vector(up, left));

        if (isHorizontal)
            neighboured.push(new bases.alignment.Vector(position.i, horizontalK));

        if (isVertical)
            neighboured.push(new bases.alignment.Vector(verticalK, position.j));

        if (isInsertion)
            neighboured.push(new bases.alignment.Vector(position.i, left));

        if (isDeletion)
            neighboured.push(new bases.alignment.Vector(up, position.j));

        if (!(isMatchMismatch || isHorizontal || isVertical || isInsertion || isDeletion)
            && (position.i !== 0 || position.j !== 0))
            neighboured.push(new bases.alignment.Vector(0, 0));

        return neighboured;
    }

    /**
     * Computes the vertical position from which you get to the currentValue.
     * @param algorithm {Object} - Contains an alignment algorithm.
     * @param currentValue {number} - The value from the current cell.
     * @param position {Object} - Current cell position in matrix.
     * @param outputData {Object} - Contains all output data.
     * @return {number} - The matching position. You get back NaN if such position does not exists.
     */
    function searchVerticalMatchPosition(algorithm, currentValue, position, outputData) {
        if (position.j > 0) {
            for (var k = 1; k <= position.i; k++) {
                if (alignmentInstance.differenceLowerEpsilon(outputData.matrix[position.i - k][position.j] + algorithm.gapFunction(k), currentValue, EPSILON))
                    return position.i - k;
            }
        }

        return Number.NaN;
    }

    /**
     * Computes the horizontal position from which you get to the currentValue.
     * @param algorithm {Object} - Contains an alignment algorithm.
     * @param currentValue {number} - The value from the current cell.
     * @param position {Object} - Current cell position in matrix.
     * @param outputData {Object} - Contains all output data.
     * @return {number} - The matching position. You get back NaN if such position does not exists.
     */
    function searchHorizontalMatchPosition(algorithm, currentValue, position, outputData) {
        if (position.i > 0) {
            for (var k = 1; k <= position.j; k++) {
                if (alignmentInstance.differenceLowerEpsilon(outputData.matrix[position.i][position.j - k] + algorithm.gapFunction(k), currentValue, EPSILON))
                    return position.j - k;
            }
        }

        return Number.NaN;
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