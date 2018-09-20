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
    if (loaded === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER) {  // to avoid self execution on a script import
        arslanEgeciougluPevzner.startArslanEgeciougluPevzner();
        loaded = ALGORITHMS.NONE;
    }
});

(function () {  // namespace
    // public methods
    namespace("arslanEgeciougluPevzner", startArslanEgeciougluPevzner, ArslanEgeciougluPevzner);

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
        var linearAlignmentInterface = new interfaces.linearAlignmentInterface.LinearAlignmentInterface();
        linearAlignmentInterface.startLinearAlignmentAlgorithm(ArslanEgeciougluPevzner, ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER);
    }

    /*---- ALGORITHM ----*/
    /**
     * Computes the optimal, normalized local alignment.
     * @constructor
     * @see https://doi.org/10.1093/bioinformatics/17.4.327
     *
     * Arslan, Abdullah N., Omer Egecioglu, and Pavel A. Pevzner.
     * "A new approach to sequence comparison: normalized sequence alignment."
     * Bioinformatics 17.4 (2001): 327-337.
     */
    function ArslanEgeciougluPevzner() {
        arslanEgeciougluPevznerInstance = this;

        // variables
        this.type = ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER;
        this.numberOfIterations = 0;
        this.lambda = 0;
        this.lastLambda = Number.POSITIVE_INFINITY;

        // instances (do not change order)
        alignmentInstance = new bases.alignment.Alignment(this);
        smithWatermanInstance = new smithWaterman.SmithWaterman();

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
     * @augments Alignment.setLinearAlignmentInput()
     */
    function setInput(inputViewmodel) {
        alignmentInstance.setIO(inputData, {});
        alignmentInstance.setLinearAlignmentInput(inputViewmodel);
        inputData.length = inputViewmodel.length();
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
        arslanEgeciougluPevznerInstance.lambda = 0;
        arslanEgeciougluPevznerInstance.lastLambda = Number.POSITIVE_INFINITY;
        arslanEgeciougluPevznerInstance.numberOfIterations = 0; // because last one is counted
        outputData.maxNumberIterations = false;
        computeAllIterationData(input, currentData, iterationData);

        // storage of output
        outputData.iterationData = iterationData;
        outputData.maxNumberIterations = arslanEgeciougluPevznerInstance.numberOfIterations > MAX_NUMBER_ITERATIONS;
        return [inputData, outputData];
    }

    /**
     * Initializes the given input with the read in inputData.
     * @param input {Object} - The input which has to be initialized.
     */
    function initializeInput(input) {
        input.sequenceB = inputData.sequenceB;
        input.sequenceA = inputData.sequenceA;
        input.calculationType = inputData.calculationType;
        input.matrixHeight = inputData.sequenceA.length + 1;
        input.matrixWidth = inputData.sequenceB.length + 1;
    }

    /**
     * Executes a recursive
     * deep-first-search (deleting last found path from memory)
     * to find potentially all possible iteration paths.
     * Hint: For scores bigger 0 the algorithm have always to converge (Dinkelbach).
     * @param input {Object} - The initialized Waterman-Smith input structure.
     * @param currentData {Object} - Stores the data from the current iteration. At the beginning it is empty.
     * @param iterationData {Object} - Stores the data from all iterations.
     * @see: Restricted to one path for better runtime!
     */
    function computeAllIterationData(input, currentData, iterationData) {
        // [1,4] computes Smith-Waterman with the given Lambda
        var ioData = computeSmithWaterman(input, arslanEgeciougluPevznerInstance.lambda);
        var alignments = ioData[1].alignments;

        if (alignments.length > 0) {
            // going through all alignments (initial nodes of possible paths)
            for (var i = 0; i < 1; i++) { // RESTRICTION: change "1" to "alignments.length" in for-loop to get all paths and remove MAX_NUMBER_ITERATIONS
                if (arslanEgeciougluPevznerInstance.lambda === arslanEgeciougluPevznerInstance.lastLambda
                    || arslanEgeciougluPevznerInstance.numberOfIterations >= MAX_NUMBER_ITERATIONS) {  // stop criteria checks

                    iterationData.push(currentData.slice());  // shallow copy
                    arslanEgeciougluPevznerInstance.lambda = 0;
                    arslanEgeciougluPevznerInstance.lastLambda = Number.POSITIVE_INFINITY;
                    arslanEgeciougluPevznerInstance.numberOfIterations--;  // because else break up is counted twice
                } else { // executing procedure with an alignment
                    // [2] compute score and alignment length
                    var saData = computeScoreAndLength(alignments, i);

                    var score = saData[0];
                    var alignmentLength = saData[1];

                    // [3] compute scale-factor lambda
                    arslanEgeciougluPevznerInstance.lastLambda = arslanEgeciougluPevznerInstance.lambda;
                    arslanEgeciougluPevznerInstance.lambda = score / (alignmentLength + inputData.length);

                    currentData.push(
                        getDataCopy(score, alignmentLength, arslanEgeciougluPevznerInstance.lambda, alignments,
                            ioData[1].matrix, ioData[1].tracebackPaths, i, ioData[1].moreTracebacks));
                    computeAllIterationData(input, currentData, iterationData);
                    currentData.pop();
                }
            }
            arslanEgeciougluPevznerInstance.numberOfIterations++;
        } else { // case: DELETION: -2,INSERTION: -2, MATCH: 3, MISMATCH: -1, LENGTH: 0
            var score = 0;
            var alignmentLength = 0;
            arslanEgeciougluPevznerInstance.lambda = score / (alignmentLength + inputData.length);  // 0/0 = 0 would also make sense
            currentData.push(
                getDataCopy(score, alignmentLength, arslanEgeciougluPevznerInstance.lambda, alignments,
                    ioData[1].matrix, ioData[1].tracebackPaths, 0, ioData[1].moreTracebacks));
            iterationData.push(currentData.slice());  // shallow copy
            arslanEgeciougluPevznerInstance.numberOfIterations++;
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
     * @return {[number, number]} - Score and length.
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
     * @param score {number} - The Smith-Waterman score you want store.
     * @param alignmentLength {number} - The alignment length you want store.
     * @param lambda {number} - The normalized score you want store.
     * @param alignments {Array} - The alignments you want store.
     * @param matrix {Array} - The matrix you want store.
     * @param tracebackPaths {Array} - The tracebackPaths you want store.
     * @param alignmentNumber {number} - The number of the alignment to which parameters like score and lambda are stored.
     * @param moreTracebacks {boolean} - Tells if the algorithm has aborted before all alignments were calculated.
     * @return {Object} - [score, lambda, deletion, insertion, match, mismatch, alignments matrix]
     */
    function getDataCopy(score, alignmentLength, lambda, alignments, matrix, tracebackPaths, alignmentNumber, moreTracebacks) {
        return [
            score,
            alignmentLength,
            lambda,
            inputData.deletion - lambda,
            inputData.insertion - lambda,
            inputData.match - 2 * lambda,
            inputData.mismatch - 2 * lambda,
            alignments.slice(),
            matrix.slice(),
            tracebackPaths.slice(),
            alignmentNumber,
            moreTracebacks
        ];
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