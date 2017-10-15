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
    debugger;
    if (document.title !== UNIT_TEST_WEBTITLE)  // to avoid the execution of the algorithm interfaces during a Unit-Test
        arslanEgeciougluPevzner.startArslanEgeciougluPevzner();
});

(function () {  // namespace
    // public methods
    namespace("arslanEgeciougluPevzner", startArslanEgeciougluPevzner, ArslanEgeciougluPevzner,
        getInput, setInput, compute, getNeighboured, getOutput, setIO, getSuperclass);

    // instances
    var alignmentInstance;
    var smithWatermanInstance;
    var arslanEgeciougluPevznerInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Function managing objects.
     */
    function startArslanEgeciougluPevzner() {
        imports();

        var alignmentInterface = new interfaces.alignmentInterface.AlignmentInterface();
        alignmentInterface.startAlignmentAlgorithm(ArslanEgeciougluPevzner, ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER);
    }

    /**
     * Handling imports.
     */
    function imports() {
        $.getScript(PATHS.ALIGNMENT_INTERFACE);
    }

    /*---- ALGORITHM ----*/
    /**
     * Computes the optimal, normalized local alignment.
     * @constructor
     */
    function ArslanEgeciougluPevzner() {
        arslanEgeciougluPevznerInstance = this;

        // variables
        this.type = ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER;
        this.numberOfIterations = 0;
        this.lambda = 0;
        this.lastLambda = Number.POSITIVE_INFINITY;

        // inheritance
        alignmentInstance = new bases.alignment.Alignment(this);
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
     * @param inputViewmodel {InputViewmodel} - The InputViewmodel of an appropriate algorithm.
     */
    function setInput(inputViewmodel) {
        inputData.sequenceA = inputViewmodel.sequence1();
        inputData.sequenceB = inputViewmodel.sequence2();

        inputData.calculationType = inputViewmodel.calculation();

        inputData.deletion = inputViewmodel.deletion();
        inputData.insertion = inputViewmodel.insertion();
        inputData.match = inputViewmodel.match();
        inputData.mismatch = inputViewmodel.mismatch();
        inputData.length = inputViewmodel.length();

        inputData.matrixHeight = inputData.sequenceB.length + 1;
        inputData.matrixWidth = inputData.sequenceA.length + 1;
    }

    /**
     * Starts the computation by initializing
     * the Smith-Waterman input with the inputData.
     */
    function compute() {
        var iterationData = [];
        var currentData = [];

        // initialize
        var input = {};
        initializeInput(input);
        arslanEgeciougluPevznerInstance.numberOfIterations = 0;
        computeAllIterationData(input, currentData, iterationData);

        // storage of output
        outputData.iterationData = iterationData;
        return [inputData, outputData];
    }

    /**
     * Initializes the given input with the read in inputData.
     * @param input {Object} - The input which has to be initialized.
     */
    function initializeInput(input) {
        input.sequenceA = inputData.sequenceA;
        input.sequenceB =  inputData.sequenceB;
        input.calculationType = inputData.calculationType;
        input.matrixHeight = inputData.sequenceB.length + 1;
        input.matrixWidth = inputData.sequenceA.length + 1;
    }

    /**
     * Executes a recursive
     * deep-first-search (deleting last found path from memory)
     * to find all possible iteration paths.
     * So, there are multiple alignments which are found by Smith-Waterman
     * and every alignment has to be tried out
     * or you maybe won't found all possible alignments with final Lambda.
     * @param input {Object} - The initialized Waterman-Smith input structure.
     * @param currentData {Object} - Stores the data from the current iteration. At the beginning it is empty.
     * @param iterationData {Object} - Stores the data from all iterations.
     */
    function computeAllIterationData(input, currentData, iterationData) {
        // [1,4] computes Smith-Waterman with the given Lambda
        var ioData = computeSmithWaterman(input, arslanEgeciougluPevznerInstance.lambda);

        var alignments = ioData[1].alignments;

        // going through all alignments (initial nodes of possible paths)
        for (var i = 0; i < alignments.length; i++) {
            if (arslanEgeciougluPevznerInstance.lambda === arslanEgeciougluPevznerInstance.lastLambda
                || arslanEgeciougluPevznerInstance.numberOfIterations >= MAX_NUMBER_ITERATIONS) {  // stop criteria checks

                iterationData.push(currentData.slice());  // shallow copy
                arslanEgeciougluPevznerInstance.lambda = 0;
                arslanEgeciougluPevznerInstance.lastLambda = Number.POSITIVE_INFINITY;
                arslanEgeciougluPevznerInstance.numberOfIterations++;
            } else { // executing procedure with an alignment
                // [2] compute score and alignment length
                var saData = computeScoreAndLength(alignments, i);

                var score = saData[0];
                var alignmentLength = saData[1];

                // [3] compute scale-factor lambda
                arslanEgeciougluPevznerInstance.lastLambda = arslanEgeciougluPevznerInstance.lambda;
                arslanEgeciougluPevznerInstance.lambda = score / (alignmentLength + inputData.length);

                currentData.push(
                    getDataCopy(score, alignmentLength, arslanEgeciougluPevznerInstance.lambda, alignments, ioData[1].matrix, ioData[1].tracebackPaths, i));
                computeAllIterationData(input, currentData, iterationData);
                currentData.pop();
            }
        }
    }

    /**
     * Computes Smith-Waterman output with the given lambda.
     * @param input {Object} - The initialized Waterman-Smith input structure.
     * @param lambda {number} - The last computed normalized score lambda.
     * @return {Object} - Output data of Smith-Waterman with the given parameter lambda.
     */
    function computeSmithWaterman(input, lambda) {
        input.deletion = inputData.deletion - lambda;
        input.insertion = inputData.deletion - lambda;
        input.match = inputData.match - 2 * lambda;
        input.mismatch = inputData.mismatch - 2 * lambda;

        smithWatermanInstance.setIO(input, {});

        return smithWatermanInstance.compute();
    }

    /**
     * Computes score and length.
     * @param alignments {Object} - The alignments computed with Smith-Waterman.
     * @param i {number} - The selected alignments index.
     * @return {[number,number]} - Score and length.
     */
    function computeScoreAndLength(alignments, i) {
        return [getScore(alignments[i]), getAlignmentLength(alignments[i])];
    }

    /**
     * Computes the score.
     * @param alignment {Object} - The alignment computed with Smith-Waterman.
     * @return {number} - Score.
     */
    function getScore(alignment) {
        var alignedSequenceA = alignment[0];
        var alignedSequenceB = alignment[2];

        var alignmentLength = alignedSequenceA.length;

        var score = 0;
        for (var i = 0; i < alignmentLength; i++) {
            if (alignedSequenceA[i] === SYMBOLS.GAP)
                score += inputData.insertion;
            else if (alignedSequenceB[i] === SYMBOLS.GAP)
                score += inputData.deletion;
            else if (alignedSequenceA[i] === alignedSequenceB[i])
                score += inputData.match;
            else if (alignedSequenceA[i] !== alignedSequenceB[i])
                score += inputData.mismatch;
        }

        return score;
    }

    /**
     * Computes the alignment length.
     * @param alignment {Object} - The alignment computed with Smith-Waterman.
     * @return {number} - Score.
     */
    function getAlignmentLength(alignment) {
        var alignedSequenceA = alignment[0];
        var alignedSequenceB = alignment[2];

        var numOfCharactersA = countCharacters(alignedSequenceA);
        var numOfCharactersB = countCharacters(alignedSequenceB);

        return numOfCharactersA + numOfCharactersB;
    }

    /**
     * Returns the number of non-gaps in the given sequence.
     * @param sequence {string} - Sequence in which the characters should be counted.
     * @return {number} - Number of non-gaps.
     */
    function countCharacters(sequence) {
        var numCharacters = 0;

        for (var i = 0; i < sequence.length; i++) {
            if (sequence[i] !== SYMBOLS.GAP)
                numCharacters++;
        }

        return numCharacters;
    }

    /**
     * Creates a copy of the given iteration data to display it later on.
     * @param score - The Smith-Waterman score you want store.
     * @param alignmentLength - The alignment length you want store.
     * @param lambda - The normalized score you want store.
     * @param alignments - The alignments you want store.
     * @param alignments - The matrix you want store.
     * @return {Object} - [score, lambda, deletion, insertion, match, mismatch, alignments matrix]
     */
    function getDataCopy(score, alignmentLength, lambda, alignments, matrix, tracebackPaths, alignmentNumber) {
        return [
            score,
            alignmentLength,
            lambda,
            inputData.deletion - lambda,
            inputData.deletion - lambda,
            inputData.match - 2 * lambda,
            inputData.mismatch - 2 * lambda,
            alignments.slice(),
            matrix.slice(),
            tracebackPaths.slice(),
            alignmentNumber
        ];
    }

    /**
     * Returns the neighbours to which you can go from the current cell position used as input.
     * @param position {Vector} - Current cell position in matrix.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     * @param algorithm {Object} - Contains an alignment algorithm.
     * @return {Array} - Contains neighboured positions as Vector-objects.
     * @see: It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function getNeighboured(position, inputData, outputData, algorithm) {
        return smithWatermanInstance.getNeighboured(position, inputData, outputData, algorithm);
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