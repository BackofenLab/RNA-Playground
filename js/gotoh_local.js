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
    if (loaded === ALGORITHMS.GOTOH_LOCAL) {  // to avoid self execution on a script import
        gotohLocal.startGotohLocal();
        loaded = ALGORITHMS.NONE;
    }
});

(function () {  // namespace
    // public methods
    namespace("gotohLocal", startGotohLocal, GotohLocal, getInput, setInput, compute, getNeighboured, getOutput, setIO, getSuperclass);

    // instances
    var alignmentInstance;
    var gotohInstance;
    var gotohLocalInstance;
    var smithWatermanInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Function managing objects.
     */
    function startGotohLocal() {
        imports();

        var subadditiveAlignmentInterface = new interfaces.subadditiveAlignmentInterface.SubadditiveAlignmentInterface();
        subadditiveAlignmentInterface.startSubadditiveAlignmentAlgorithm(GotohLocal, ALGORITHMS.GOTOH_LOCAL);
    }

    /**
     * Handling imports.
     */
    function imports() {
        $.getScript(PATHS.SUBADDITIVE_ALIGNMENT_INTERFACE);
    }

    /*---- ALGORITHM ----*/
    /**
     * Computes the optimal, local affine alignment.
     * @constructor
     * @augments Alignment
     */
    function GotohLocal() {
        gotohLocalInstance = this;

        // variables
        this.type = ALGORITHMS.GOTOH;
        this.numberOfTracebacks = 0;

        // inheritance
        alignmentInstance = new bases.alignment.Alignment(this);
        gotohInstance = new gotoh.Gotoh();
        smithWatermanInstance = new smithWaterman.SmithWaterman();

        // public class methods
        this.getInput = getInput;

        this.setInput = setInput;
        this.compute = compute;
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
        inputData.sequenceA = inputViewmodel.sequence1();
        inputData.sequenceB = inputViewmodel.sequence2();

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
        outputData.matrix[0][0] = 0;

        for (var i = 1; i < inputData.matrixHeight; i++) {
            outputData.matrix[i][0] = 0;
            outputData.horizontalGaps[i][0] = Number.NEGATIVE_INFINITY;
            outputData.verticalGaps[i][0] = SYMBOLS.GAP;
        }

        for (var j = 1; j < inputData.matrixWidth; j++) {
            outputData.matrix[0][j] = 0;
            outputData.horizontalGaps[0][j] = SYMBOLS.GAP;
            outputData.verticalGaps[0][j] = Number.NEGATIVE_INFINITY;
        }
    }

    /**
     * Computes the matrix by using the recursion function and the score.
     */
    function computeMatricesAndScore() {
        var maxValue = 0;

        // going through every matrix cell
        for (var i = 1; i < inputData.matrixHeight; i++) {
            var bChar = inputData.sequenceB[i - 1];

            for (var j = 1; j < inputData.matrixWidth; j++) {
                var aChar = inputData.sequenceA[j - 1];

                outputData.matrix[i][j] = gotohInstance.recursionFunction(aChar, bChar, i, j, Math.max, local);

                // storing maximum
                if (maxValue < outputData.matrix[i][j])
                    maxValue = outputData.matrix[i][j];
            }
        }

        // score is the maximum value
        outputData.score = maxValue;
    }

    /**
     * Initializes the traceback.
     * @override Alignment.computeTraceback()
     */
    function computeTraceback() {
        gotohLocalInstance.numberOfTracebacks = 0;

        // computing all traceback start-positions
        var backtraceStarts = smithWatermanInstance.getAllMaxPositions(inputData, outputData);

        outputData.tracebackPaths = [];
        outputData.moreTracebacks = false;

        for (var i = 0; i < backtraceStarts.length; i++) {
            var tracebackPaths = alignmentInstance.getLocalTraces([backtraceStarts[i]], inputData, outputData, -1, getNeighboured);
            outputData.tracebackPaths = outputData.tracebackPaths.concat(tracebackPaths);
        }
    }

    /**
     * Returns the neighbours to which you can go from the current cell position used as input.
     * @param position {Vector} - Current cell position in matrix.
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
            return gotohInstance.getVerticalNeighboured(position, inputData, outputData);
        else if (position.label === MATRICES.HORIZONTAL)
            return gotohInstance.getHorizontalNeighboured(position, inputData, outputData);

        var left = position.j - 1;
        var up = position.i - 1;

        // retrieve values
        var aChar = left >= 0 ? inputData.sequenceA[left] : SYMBOLS.EMPTY;
        var bChar = up >= 0 ? inputData.sequenceB[up] : SYMBOLS.EMPTY;

        var currentValue = outputData.matrix[position.i][position.j];

        var matchOrMismatch = aChar === bChar ? inputData.match : inputData.mismatch;

        var diagonalValue = left >= 0 && up >= 0 ? outputData.matrix[up][left] : Number.NaN;
        var verticalValue = up >= 0 ? outputData.verticalGaps[position.i][position.j] : Number.NaN;
        var horizontalValue = left >= 0 ? outputData.horizontalGaps[position.i][position.j] : Number.NaN;

        var upValue = up >= 0 && position.j === 0 ? outputData.matrix[up][position.j] : Number.NaN;
        var leftValue = left >= 0 && position.i === 0 ? outputData.matrix[position.i][left] : Number.NaN;

        // check
        var isMatchMismatch = currentValue === (diagonalValue + matchOrMismatch);
        var isChangeToP = currentValue === verticalValue;
        var isChangeToQ = currentValue === horizontalValue;

        isMatchMismatch = isMatchMismatch || currentValue === 0 && up >= 0 && left >= 0;
        var isDeletion = currentValue === 0 && up >= 0;  // lower 0 -> cut away
        var isInsertion = currentValue === 0 && left >= 0;  // lower 0 -> cut away

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