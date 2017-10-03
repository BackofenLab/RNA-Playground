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

    function getInput() {
        return inputData;
    }

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

    function compute() {
        debugger;
        initializeMatrix();
        computeMatrixAndScore();
        computeTraceback();
        createAlignments();
        return [inputData, outputData];
    }

    function initializeMatrix() {
        createMatrix();
        childInstance.initializeMatrix();
    }

    function createMatrix() {
        outputData.matrix = new Array(inputData.matrixHeight);

        for (var i = 0; i < inputData.matrixHeight; i++)
            outputData.matrix[i] = new Array(inputData.matrixWidth);
    }

    function computeMatrixAndScore() {
        childInstance.computeMatrixAndScore();
    }

    function recursionFunction(aChar, bChar, i, j) {
        var matchOrMismatch = aChar === bChar ? inputData.match : inputData.mismatch;

        var diagonalValue = outputData.matrix[i - 1][j - 1] + matchOrMismatch;
        var upValue = outputData.matrix[i - 1][j] + inputData.deletion;
        var leftValue = outputData.matrix[i][j - 1] + inputData.insertion;

        return childInstance.recursionFunction(diagonalValue, upValue, leftValue);
    }

    function computeTraceback() {
        childInstance.computeTraceback();
    }

    function createAlignments() {
        outputData.alignments = [];
        var numTracebacks = outputData.tracebackPaths.length;

        for (var i = 0; i < numTracebacks; i++) {
            var alignment = createAlignment(outputData.tracebackPaths[i]);
            outputData.alignments.push(alignment);
        }
    }

    /**
     * Hint: It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     * @param path - Path which is used to create the alignment.
     * @return {[alignedSequenceA, matchOrMismatchString, alignedSequenceB]}
     */
    function createAlignment(path) {
        path.reverse();  // allows more intuitive calculations from left to right

        var alignedSequenceA = SYMBOLS.EMPTY;
        var matchOrMismatchString = SYMBOLS.EMPTY;
        var alignedSequenceB = SYMBOLS.EMPTY;

        var currentPositionA = path[0].j;
        var currentPositionB = path[0].i;

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

        return [alignedSequenceA.toUpperCase(), matchOrMismatchString, alignedSequenceB.toUpperCase()];
    }

    function getOutput() {
        return outputData;
    }

    function setIO(input, output) {
        inputData = input;
        outputData = output;
    }
}());
