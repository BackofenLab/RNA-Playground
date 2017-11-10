/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("bases.alignment", Vector, create,
        Alignment, getInput, setLinearAlignmentInput, setSubadditiveAlignmentInput, compute, recursionFunction,
        affineRecursionFunction, getGlobalTraces, getLocalTraces, getAllMaxPositions, createAlignments, getOutput, setIO, getLastChild,
        getNeighboured, getVerticalNeighboured, getHorizontalNeighboured, differenceLowerEpsilon);

    // instances
    var alignmentInstance;
    var childInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Creates a vector with a label, which is by the default the default matrix label.
     * @param i {number} - First matrix axis.
     * @param j {number} - Second matrix axis.
     * @constructor
     */
    function Vector(i, j) {
        this.i = i;
        this.j = j;
        this.label = MATRICES.DEFAULT;
    }

    /**
     * Changes the label of a vector. It is used to create a vector with a specific label.
     * @example
     * create(new Vector(i,j), MATRICES.HORIZONTAL)
     * @param Vector - The vector of which the label should be changed.
     * @param label {string} - The string label from "defaults.js".
     * @return {Vector} - The vector with the changed label.
     */
    function create(Vector, label) {
        Vector.label = label;
        return Vector;
    }

    /**
     * Contains functions to compute optimal alignments.
     * It is used by global and local alignment algorithms as superclass
     * to avoid code duplicates.
     * @param child - The child algorithm which inherits from this class.
     * @constructor
     */
    function Alignment(child) {
        alignmentInstance = this;
        childInstance = child;

        // public methods
        this.getInput = getInput;

        this.setLinearAlignmentInput = setLinearAlignmentInput;
        this.setSubadditiveAlignmentInput = setSubadditiveAlignmentInput;
        this.compute = compute;
        this.recursionFunction = recursionFunction;
        this.affineRecursionFunction = affineRecursionFunction;
        this.getGlobalTraces = getGlobalTraces;
        this.getLocalTraces = getLocalTraces;
        this.getAllMaxPositions = getAllMaxPositions;
        this.getNeighboured = getNeighboured;
        this.getVerticalNeighboured = getVerticalNeighboured;
        this.getHorizontalNeighboured = getHorizontalNeighboured;
        this.createAlignments = createAlignments;
        this.getOutput = getOutput;

        this.setIO = setIO;
        this.getLastChild = getLastChild;

        this.differenceLowerEpsilon = differenceLowerEpsilon;
    }

    /**
     * Returns the input data of the algorithm.
     * @return {Object} - Contains all input data.
     */
    function getInput() {
        return inputData;
    }

    /**
     * Sets the algorithm input for an appropriate linear alignment algorithm
     * which is using the inputViewmodel properties in its computations.
     * @param inputViewmodel {Object} - The InputViewmodel of an appropriate algorithm (NW, SW, AEP).
     */
    function setLinearAlignmentInput(inputViewmodel) {
        inputData.sequenceA = inputViewmodel.sequence1();
        inputData.sequenceB = inputViewmodel.sequence2();

        inputData.calculationType = inputViewmodel.calculation();

        inputData.deletion = inputViewmodel.gap();
        inputData.insertion = inputViewmodel.gap();
        inputData.match = inputViewmodel.match();
        inputData.mismatch = inputViewmodel.mismatch();

        inputData.matrixHeight = inputData.sequenceB.length + 1;
        inputData.matrixWidth = inputData.sequenceA.length + 1;
    }

    /**
     * Sets the algorithm input for an appropriate subadditive alignment algorithm
     * which is using the inputViewmodel properties in its computations.
     * @param inputViewmodel {Object} - The InputViewmodel of an appropriate algorithm (G, WSB).
     */
    function setSubadditiveAlignmentInput(inputViewmodel) {
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
     * Computes the cell score.
     * @param aChar {string} - The current char from the first string.
     * @param bChar {string} - The current char from the second string.
     * @param i {number} - The current vertical position in the matrix.
     * @param j {number} - The current horizontal position in the matrix.
     * @param optimum {Function} - The function which should be used for optimization {Math.min, Math.max}.
     * @param local {boolean} - Tells if the local recursion function should be used.
     * @return {number} - The value for the cell at position (i,j).
     */
    function affineRecursionFunction(aChar, bChar, i, j, optimum, local) {
        var matchOrMismatch = aChar === bChar ? inputData.match : inputData.mismatch;

        if (aChar === SYMBOLS.NONE || bChar === SYMBOLS.NONE) matchOrMismatch = 0;  // extension for Feng-Doolittle

        if (inputData.substitutionFunction !== undefined)  // extension for T-Coffee
            matchOrMismatch = inputData.substitutionFunction(i,j);

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
     * @abstract
     */
    function computeTraceback() {
        childInstance.computeTraceback();
    }

    /**
     * Gets the global tracebacks by starting the traceback procedure
     * with some path containing the first element of the path.
     * @param path {Array} - Array containing the first vector element from which on you want find a path.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     * @param pathLength {number} - Tells after how many edges the procedure should stop.
     * The value -1 indicates arbitrarily long paths.
     * @param neighbourFunction {Function} - The function which have to be used to retrieve neighbours.
     * @return {Array} - Array of paths.
     */
    function getGlobalTraces(path, inputData, outputData, pathLength, neighbourFunction) {
        alignmentInstance.stopTraceback = false;

        if (childInstance.getTraces !== undefined)
            return childInstance.getTraces(path, inputData, outputData, pathLength, neighbourFunction);

        var paths = [];
        globalTraceback(paths, path, inputData, outputData, pathLength, neighbourFunction);
        return paths;
    }

    /**
     * Computing the global traceback and stops after it has found a constant number of tracebacks.
     * It sets a flag "moreTracebacks" in the "outputData", if it has stopped before computing all tracebacks.
     * The traceback algorithm executes a recursive,
     * modified deep-first-search (deleting last found path from memory)
     * with special stop criteria on the matrix cells as path-nodes.
     * @param paths {Array} - Array of paths.
     * @param path {Array} - Array containing the first vector element from which on you want find a path.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     * @param pathLength {number} - Tells after how many edges the procedure should stop.
     * @param neighbourFunction {Function} - The function which have to be used to retrieve neighbours.
     * @see It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function globalTraceback(paths, path, inputData, outputData, pathLength, neighbourFunction) {
        var currentPosition = path[path.length - 1];
        var neighboured = neighbourFunction(currentPosition, inputData, outputData, childInstance);

        // going through all successors (initial nodes of possible paths)
        for (var i = 0; i < neighboured.length; i++) {
            if ((neighboured[i].i === 0 && neighboured[i].j === 0)  // stop criteria checks
                || (pathLength !== -1 && path.length >= pathLength)
                || outputData.moreTracebacks
                || alignmentInstance.stopTraceback) {

                if (alignmentInstance.stopTraceback)
                    break;

                if (inputData.computeOneAlignment)  // extension to speed up Feng-Doolittle
                    alignmentInstance.stopTraceback = true;

                path.push(neighboured[i]);

                // path storage, if MAX_NUMBER_TRACEBACKS is not exceeded
                if (childInstance.numberOfTracebacks < MAX_NUMBER_TRACEBACKS) {
                    paths.push(path.slice());  // creating a shallow copy
                    childInstance.numberOfTracebacks++;
                } else
                    outputData.moreTracebacks = true;

                path.pop();
            } else {
                // executing procedure with a successor
                path.push(neighboured[i]);
                globalTraceback(paths, path, inputData, outputData, pathLength, neighbourFunction);
                path.pop();
            }
        }
    }

    /**
     * Gets tracebacks by starting the traceback procedure
     * with some path containing the first element of the path.
     * @param path {Array} - Array containing the first vector element from which on you want find a path.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     * @param pathLength {number} - Tells after how many edges the procedure should stop.
     * The value -1 indicates arbitrarily long paths.
     * @param neighbourFunction {Function} - The function which have to be used to retrieve neighbours.
     * @return {Array} - Array of paths.
     */
    function getLocalTraces(path, inputData, outputData, pathLength, neighbourFunction) {
        alignmentInstance.stopTraceback = false;

        var paths = [];
        localTraceback(paths, path, inputData, outputData, pathLength, neighbourFunction);
        return paths;
    }

    /**
     * Computing the local traceback and stops after it has found a constant number of tracebacks.
     * It sets a flag "moreTracebacks" in the "outputData", if it has stopped before computing all tracebacks.
     * The traceback algorithm executes a recursive,
     * modified deep-first-search (deleting last found path from memory)
     * with special stop criteria on the matrix cells as path-nodes.
     * @param paths {Array} - Array of paths.
     * @param path {Array} - Array containing the first vector element from which on you want find a path.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     * @param pathLength {number} - Tells after how many edges the procedure should stop.
     * @param neighbourFunction {Function} - The function which have to be used to retrieve neighbours.
     * @see It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function localTraceback(paths, path, inputData, outputData, pathLength, neighbourFunction) {
        var currentPosition = path[path.length - 1];
        var neighboured = neighbourFunction(currentPosition, inputData, outputData, childInstance);

        // going through all successors (initial nodes of possible paths)
        for (var i = 0; i < neighboured.length; i++) {
            if (outputData.matrix[neighboured[i].i][neighboured[i].j] === 0  // stop criteria checks
                || (pathLength !== -1 && path.length >= pathLength)
                || outputData.moreTracebacks) {

                if (alignmentInstance.stopTraceback)
                    break;

                if (inputData.computeOneAlignment)  // extension to speed up T-Coffee's local alignment
                    alignmentInstance.stopTraceback = true;

                path.push(neighboured[i]);

                // path storage, if MAX_NUMBER_TRACEBACKS is not exceeded
                if (childInstance.numberOfTracebacks < MAX_NUMBER_TRACEBACKS) {
                    paths.push(path.slice());  // creating a shallow copy
                    childInstance.numberOfTracebacks++;
                } else
                    outputData.moreTracebacks = true;

                path.pop();
            } else {
                // executing procedure with a successor
                path.push(neighboured[i]);
                localTraceback(paths, path, inputData, outputData, pathLength, neighbourFunction);
                path.pop();
            }
        }
    }

    /**
     * Returning all maximums of the computed matrix.
     * @param inputData {Object} - Containing information about the output matrix.
     * @param outputData {Object} - Containing the output matrix.
     * @return {Array} - Array of vectors (max-positions).
     */
    function getAllMaxPositions(inputData, outputData) {
        var maxPositions = [];

        if (outputData.score > 0) {  // only positions bigger 0 can be start positions (because local alignments never lower 0)
            for (var i = 0; i < inputData.matrixHeight; i++) {
                for (var j = 0; j < inputData.matrixWidth; j++) {
                    if (outputData.matrix[i][j] === outputData.score) {
                        maxPositions.push(new bases.alignment.Vector(i, j));
                    }
                }
            }
        }

        return maxPositions;
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
     * @return {[alignedSequenceA, matchOrMismatchString, alignedSequenceB]} - The triple of strings which have to be displayed.
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
            var horizontalDifference = path[i].j - path[i - 1].j;
            var verticalDifference = path[i].i - path[i - 1].i;

            if (verticalDifference === 1 && horizontalDifference === 1) {  // diagonal case
                var aChar = inputData.sequenceA[currentPositionA];
                var bChar = inputData.sequenceB[currentPositionB];

                alignedSequenceA += aChar;
                matchOrMismatchString += aChar === bChar ? SYMBOLS.STAR : SYMBOLS.VERTICAL_BAR;
                alignedSequenceB += bChar;

                currentPositionA++;
                currentPositionB++;
            } else if (horizontalDifference > 0) {  // horizontal case
                for (var k = 0; k < horizontalDifference; k++) {
                    alignedSequenceA += inputData.sequenceA[currentPositionA];
                    matchOrMismatchString += SYMBOLS.SPACE;
                    alignedSequenceB += SYMBOLS.GAP;

                    currentPositionA++;
                }
            } else if (verticalDifference > 0) {  // vertical case
                // Hint: for Gotoh really "else if" is needed because you can switch between matrices
                for (var k = 0; k < verticalDifference; k++) {
                    alignedSequenceA += SYMBOLS.GAP;
                    matchOrMismatchString += SYMBOLS.SPACE;
                    alignedSequenceB += inputData.sequenceB[currentPositionB];

                    currentPositionB++;
                }
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

    /**
     * Returns the child which is has currently worked with that class.
     * @return {Object} - The child object.
     */
    function getLastChild() {
        return childInstance;
    }

    /**
     * Returns the neighbours to which you can go from the current cell position used as input.
     * This function is used by Smith-Waterman and Needleman-Wunsch.
     * Other algorithms specify their neighbour-function within the algorithm.
     * @param position {Vector} - Current cell position in matrix.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     * @param algorithm {Object} - Contains an alignment algorithm.
     * @return {Array} - Contains neighboured positions as Vector-objects.
     * @see: It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function getNeighboured(position, inputData, outputData, algorithm) {
        var neighboured = [];

        var left = position.j - 1;
        var up = position.i - 1;

        // retrieve values
        var aChar = left >= 0 ? inputData.sequenceA[left] : SYMBOLS.EMPTY;
        var bChar = up >= 0 ? inputData.sequenceB[up] : SYMBOLS.EMPTY;

        var currentValue = outputData.matrix[position.i][position.j];

        var matchOrMismatch = aChar === bChar ? inputData.match : inputData.mismatch;

        var diagonalValue = left >= 0 && up >= 0 ? outputData.matrix[up][left] : NaN;
        var upValue = up >= 0 ? outputData.matrix[up][position.j] : NaN;
        var leftValue = left >= 0 ? outputData.matrix[position.i][left] : NaN;

        // check
        var isMatchMismatch = differenceLowerEpsilon(currentValue, diagonalValue + matchOrMismatch, EPSILON);
        var isDeletion = differenceLowerEpsilon(currentValue, upValue + inputData.deletion, EPSILON);
        var isInsertion = differenceLowerEpsilon(currentValue, leftValue + inputData.insertion, EPSILON);

        if (algorithm.type === ALGORITHMS.SMITH_WATERMAN) {
            isMatchMismatch = isMatchMismatch || currentValue === 0 && up >= 0 && left >= 0;
            isDeletion = isDeletion || currentValue === 0 && up >= 0;  // lower 0 -> cut away
            isInsertion = isInsertion || currentValue === 0 && left >= 0;  // lower 0 -> cut away
        }

        // add
        if (isMatchMismatch)
            neighboured.push(new Vector(up, left));

        if (isDeletion)
            neighboured.push(new Vector(up, position.j));

        if (isInsertion)
            neighboured.push(new Vector(position.i, left));

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
     * Tests if the difference between to values is lower some parameter Epsilon.
     * Hint: It would maybe work without this function. So, this function is only for security reasons.
     * @param value1 - The first value.
     * @param value2 - The second value.
     * @param epsilon - The small number you test against.
     * @return {boolean} - The
     */
    function differenceLowerEpsilon(value1, value2, epsilon) {
        var difference = value1 - value2;

        return Math.abs(difference) < epsilon;
    }
}());
