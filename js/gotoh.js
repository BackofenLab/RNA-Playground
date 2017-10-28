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
    if (loaded === ALGORITHMS.GOTOH) {  // to avoid self execution on a script import
        gotoh.startGotoh();
        loaded = ALGORITHMS.NONE;
    }
});

(function () {  // namespace
    // public methods
    namespace("gotoh", startGotoh, Gotoh, getInput, setInput, compute,
        recursionFunction, getNeighboured, getVerticalNeighboured, getHorizontalNeighboured, getOutput, setIO, getSuperclass);

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
        var subadditiveAlignmentInterface = new interfaces.subadditiveAlignmentInterface.SubadditiveAlignmentInterface();
        subadditiveAlignmentInterface.startSubadditiveAlignmentAlgorithm(Gotoh, ALGORITHMS.GOTOH);
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
        alignmentInstance = new bases.alignment.Alignment(this);

        // public class methods
        this.getInput = getInput;

        this.setInput = setInput;
        this.compute = compute;
        this.recursionFunction = recursionFunction;
        this.getNeighboured = getNeighboured;
        this.getVerticalNeighboured = getVerticalNeighboured;
        this.getHorizontalNeighboured = getHorizontalNeighboured;
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
        inputData.computeOneAlignment = false;  // extension to speed up Feng-Doolittle, default value is false
        alignmentInstance.setIO(inputData, {});
        alignmentInstance.setSubadditiveAlignmentInput(inputViewmodel);
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
        outputData.matrix = new Array(inputData.matrixHeight);
        outputData.horizontalGaps = new Array(inputData.matrixHeight);
        outputData.verticalGaps = new Array(inputData.matrixHeight);

        for (var i = 0; i < inputData.matrixHeight; i++) {
            outputData.matrix[i] = new Array(inputData.matrixWidth);
            outputData.horizontalGaps[i] = new Array(inputData.matrixWidth);
            outputData.verticalGaps[i] = new Array(inputData.matrixWidth);
        }
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
                    outputData.matrix[i][j] = recursionFunction(aChar, bChar, i, j, Math.min, false);
                else  // inputData.calculationType === ALIGNMENT_TYPES.SIMILARITY
                    outputData.matrix[i][j] = recursionFunction(aChar, bChar, i, j, Math.max, false);
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
     * @param local {boolean} - Tells if the local recursion function should be used.
     * @return {number} - The value for the cell at position (i,j).
     */
    function recursionFunction(aChar, bChar, i, j, optimum, local) {
        var matchOrMismatch = aChar === bChar ? inputData.match : inputData.mismatch;
        if (aChar === SYMBOLS.NONE || bChar === SYMBOLS.NONE) matchOrMismatch = 0;  // extension for Feng-Doolittle

        // gap recursion-functions
        outputData.horizontalGaps[i][j] = horizontalOptimum(optimum, i, j);
        outputData.verticalGaps[i][j] = verticalOptimum(optimum, i, j);

        // default matrix recursion function
        if (local)
            return optimum(
                outputData.horizontalGaps[i][j],
                outputData.matrix[i - 1][j - 1] + matchOrMismatch,
                outputData.verticalGaps[i][j],
                0);

        // else global
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

        var lowerRightCorner = new bases.alignment.Vector(inputData.matrixHeight - 1, inputData.matrixWidth - 1);

        outputData.moreTracebacks = false;
        outputData.tracebackPaths = alignmentInstance.getGlobalTraces([lowerRightCorner], inputData, outputData, -1, getNeighboured);
    }

    /**
     * Returns the neighbours to which you can go from the current cell position used as input.
     * @param position {Object} - Current cell position in matrix.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     * @param algorithm {Object} - Contains an alignment algorithm.
     * @return {Array} - Contains neighboured positions as Vector-objects.
     * @see Hint: The parameter algorithm is needed!
     * It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function getNeighboured(position, inputData, outputData, algorithm) {
        var neighboured = [];

        if (position.label === MATRICES.VERTICAL)
            return getVerticalNeighboured(position, inputData, outputData);
        else if (position.label === MATRICES.HORIZONTAL)
            return getHorizontalNeighboured(position, inputData, outputData);

        var left = position.j - 1;
        var up = position.i - 1;

        // retrieve values
        var aChar = left >= 0 ? inputData.sequenceA[left] : SYMBOLS.EMPTY;
        var bChar = up >= 0 ? inputData.sequenceB[up] : SYMBOLS.EMPTY;

        var currentValue = outputData.matrix[position.i][position.j];

        var matchOrMismatch = aChar === bChar ? inputData.match : inputData.mismatch;
        if (aChar === SYMBOLS.NONE || bChar === SYMBOLS.NONE) matchOrMismatch = 0;  // extension for Feng-Doolittle

        var diagonalValue = left >= 0 && up >= 0 ? outputData.matrix[up][left] : Number.NaN;
        var verticalValue = up >= 0 ? outputData.verticalGaps[position.i][position.j] : Number.NaN;
        var horizontalValue = left >= 0 ? outputData.horizontalGaps[position.i][position.j] : Number.NaN;

        var upValue = up >= 0 && position.j === 0 ? outputData.matrix[up][position.j] : Number.NaN;
        var leftValue = left >= 0 && position.i === 0 ? outputData.matrix[position.i][left] : Number.NaN;

        // check
        var isMatchMismatch = currentValue === (diagonalValue + matchOrMismatch);
        var isChangeToP = currentValue === verticalValue;
        var isChangeToQ = currentValue === horizontalValue;

        var isDeletion = currentValue === upValue + inputData.enlargement;
        var isInsertion = currentValue === leftValue + inputData.enlargement;

        // add
        if (isMatchMismatch)
            neighboured.push(new bases.alignment.Vector(up, left));

        if (isChangeToP)
            neighboured.push(bases.alignment.create(new bases.alignment.Vector(position.i, position.j), MATRICES.VERTICAL));

        if (isChangeToQ)
            neighboured.push(bases.alignment.create(new bases.alignment.Vector(position.i, position.j), MATRICES.HORIZONTAL));

        if (isInsertion)
            neighboured.push(new bases.alignment.Vector(position.i, left));

        if (isDeletion)
            neighboured.push(new bases.alignment.Vector(up, position.j));

        if (!(isMatchMismatch || isChangeToP || isChangeToQ || isInsertion || isDeletion)
            && (position.i !== 0 || position.j !== 0))
            neighboured.push(new bases.alignment.Vector(0, 0));

        return neighboured;
    }

    /**
     * Returns the neighbours to which you can go from the current cell position in the matrix for vertical gap costs.
     * @param position {Object} - Current cell position in matrix.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     * @return {Array} - Contains neighboured positions as Vector-objects.
     * @see: It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function getVerticalNeighboured(position, inputData, outputData) {
        var neighboured = [];

        var up = position.i - 1;

        // retrieve values
        var currentValue = outputData.verticalGaps[position.i][position.j];

        var pUpValue = Number.NaN;
        var xUpValue = Number.NaN;

        if (position.i >= 0 && up >= 0) {
            pUpValue = outputData.verticalGaps[up][position.j];
            xUpValue = outputData.matrix[up][position.j];
        }

        // check
        var isUpInP = currentValue === pUpValue + inputData.enlargement;
        var isUpInX = currentValue === xUpValue + inputData.baseCosts + inputData.enlargement;

        // add
        if (isUpInP)
            neighboured.push(bases.alignment.create(new bases.alignment.Vector(up, position.j), MATRICES.VERTICAL));

        if (isUpInX)
            neighboured.push(new bases.alignment.Vector(up, position.j));

        return neighboured;
    }

    /**
     * Returns the neighbours to which you can go from the current cell position in the matrix for horizontal gap costs.
     * @param position {Object} - Current cell position in matrix.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     * @return {Array} - Contains neighboured positions as Vector-objects.
     * @see: It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function getHorizontalNeighboured(position, inputData, outputData) {
        var neighboured = [];

        var left = position.j - 1;

        // retrieve values
        var currentValue = outputData.horizontalGaps[position.i][position.j];

        var qLeftValue = Number.NaN;
        var xLeftValue = Number.NaN;

        if (position.i >= 0 && left >= 0) {
            qLeftValue = outputData.horizontalGaps[position.i][left];
            xLeftValue = outputData.matrix[position.i][left];
        }

        // check
        var isLeftInQ = currentValue === qLeftValue + inputData.enlargement;
        var isLeftInX = currentValue === xLeftValue + inputData.baseCosts + inputData.enlargement;

        // add
        if (isLeftInQ)
            neighboured.push(bases.alignment.create(new bases.alignment.Vector(position.i, left), MATRICES.HORIZONTAL));

        if (isLeftInX)
            neighboured.push(new bases.alignment.Vector(position.i, left));

        return neighboured;
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