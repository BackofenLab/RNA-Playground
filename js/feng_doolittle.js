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
    if (loaded === ALGORITHMS.FENG_DOOLITTLE) {  // to avoid self execution on a script import
        fengDoolittle.startFengDoolittle();
        loaded = ALGORITHMS.NONE;
    }
});

(function () {  // namespace
    // public methods
    namespace("fengDoolittle", startFengDoolittle, FengDoolittle, getInput, setInput, compute, getOutput, setIO, getSuperclass);

    // instances
    var alignmentInstance;
    var fengDoolittleInstance;
    var gotohInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Function managing objects.
     */
    function startFengDoolittle() {
        var subadditiveAlignmentInterface = new interfaces.subadditiveAlignmentInterface.SubadditiveAlignmentInterface();
        subadditiveAlignmentInterface.startSubadditiveAlignmentAlgorithm(GotohLocal, ALGORITHMS.GOTOH_LOCAL);
    }

    /*---- ALGORITHM ----*/
    /**
     * Computes the optimal, local affine alignment.
     * @constructor
     * @augments Alignment
     * @see: The superclass "alignmentInstance" have to be created as last instance
     * or the childInstance in the superclass will be probably wrong!
     */
    function FengDoolittle() {
        fengDoolittleInstance = this;

        // variables
        this.type = ALGORITHMS.FENG_DOOLITTLE;
        this.numberOfTracebacks = 0;

        // inheritance
        gotohInstance = new gotoh.Gotoh();
        alignmentInstance = new bases.alignment.Alignment(this);

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
        inputData.sequences = inputViewmodel.sequences();

        inputData.baseCosts = inputViewmodel.baseCosts();
        inputData.enlargement = inputViewmodel.enlargement();
        inputData.match = inputViewmodel.match();
        inputData.mismatch = inputViewmodel.mismatch();

        inputData.numOfStartClusters = inputData.sequences.length;
    }

    /**
     * Starts the computation.
     */
    function compute() {
        computePairwiseData();
        computeDistancesFromSimilarities();
        createDistanceMatrix();
        createPhylogeneticTree();
        createProgressiveAlignments();
        return [inputData, outputData];
    }

    /**
     * Computes scores (similarities), the number of gaps and the alignment lengths between all sequences.
     */
    function computePairwiseData() {
        var gotohInput = {};
        initializeInput(gotohInput);

        outputData.sequencePairs = [];
        outputData.similarities = [];
        outputData.alignmentLengths = [];
        outputData.gapNumbers = [];

        for (var i = 1; i < inputData.sequences.length; i++) {
            for (var j = 0; j < i; j++) {
                var sequenceA = inputData.sequences[j];
                var sequenceB = inputData.sequences[i];

                var ioData = computeGotoh(gotohInput, sequenceA, sequenceB);
                outputData.sequencePairs.push([sequenceA, sequenceB]);
                outputData.similarities.push(ioData[1].score);
                outputData.alignmentLengths.push(getAlignmentLength(sequenceA, sequenceB));
                outputData.gapNumbers.push(getNumberOfGaps(ioData[1].alignments[0]));
            }
        }
    }

    /**
     * Initializes the given input with the read in inputData.
     * @param input {Object} - The input which has to be initialized.
     */
    function initializeInput(input) {
        input.baseCosts = inputData.baseCosts;
        input.enlargement = inputData.enlargement;
        input.match = inputData.match;
        input.mismatch = inputData.mismatch;

        input.calculationType = inputData.calculationType;

        input.matrixHeight = inputData.sequenceB.length + 1;
        input.matrixWidth = inputData.sequenceA.length + 1;
    }

    /**
     * Computes Gotoh with the given input sequences and the function.
     * @param input {Object} - The initialized Gotoh input structure.
     * @return {Object} - Output data of Gotoh with the given sequences in the input.
     */
    function computeGotoh(input, sequenceA, sequenceB) {
        input.sequenceA = sequenceA;
        input.sequenceB = sequenceB;

        input.matrixHeight = inputData.sequenceB.length + 1;
        input.matrixWidth = inputData.sequenceA.length + 1;

        gotohInstance.setIO(input, {});

        return gotohInstance.compute();
    }

    /**
     * Computes the length of the alignment (length of both sequences together without gaps).
     * Hint: Usually you have to look into the alignment, but in the case of global sequences this is not needed.
     * @param sequenceA - The first (not aligned) sequence.
     * @param sequenceB - The second (not aligned) sequence.
     * @return {number} - The length of the alignment.
     */
    function getAlignmentLength(sequenceA, sequenceB) {
        return sequenceA.length + sequenceB;
    }

    /**
     * Returns the number of gaps in the alignment.
     * @param alignment {[alignedSequenceA, matchOrMismatchString, alignedSequenceB]}
     * - The triple of strings for which the number of gaps has to be computed.
     * @return {number} - The number of gaps in the alignment.
     */
    function getNumberOfGaps(alignment) {
        return countNumberOfGaps(alignment[0]) + countNumberOfGaps(alignment[2]);
    }

    /**
     * Counts the number of gaps in the given sequence.
     * @param sequence {string} - The sequence for which the number of gaps has to be counted.
     * @return {number} - The number of gaps.
     */
    function countNumberOfGaps(sequence) {
        var numGaps = 0;

        for (var i = 0; i < sequence.length; i++) {
            if (sequence[i] === SYMBOLS.GAP)
                numGaps++;
        }

        return numGaps;
    }

    /**
     * Converting the pairwise similarities into distances
     * by using the Feng-Doolittle formulas.
     * Hint: The factor 100 is really
     * inside the logarithm to scale S^eff(a,b) between (0,100].
     * Hint 2: The formula is proven with the help of several sources!
     * @example:
     * D(a,b) = -ln(S^eff(a,b) * 100)
     * where
     * S^eff(a,b) = [S(a,b) - S^rand(a,b)] / [S^max(a,b) - S^rand(a,b)]
     * @see: https://doi.org/10.1007/PL00006155
     * Feng, Da-Fei, and Russell F. Doolittle.
     * "Converting amino acid alignment scores into measures of evolutionary time:
     * a simulation study of various relationships." Journal of molecular evolution 44.4 (1997): 361-370.
     *
     * Hint: S^eff(a,b) should scale between [0,1]
     * where 1 means identical and 0 means totally non-identical,
     * but in reality S^rand can be bigger as S(a,b) and we can get negative values of S^eff(a,b).
     * To avoid negative values the approach from the paper linked above is used.
     *
     * The idea:
     * if [S(a,b) - S^rand(a,b)] < 0
     * then you set [S(a,b) - S^rand(a,b)] = 0.001
     *
     * Hint: for [S(a,b) - S^rand(a,b)] == 0 it was not defined,
     * but for simplicity we also use [S(a,b) - S^rand(a,b)] = 0.001
     */
    function computeDistancesFromSimilarities() {
        outputData.distances = [];

        // going over all similarities
        for (var i = 0; i < outputData.similarities.length; i++) {
            var score = outputData.similarities[i];

            // retrieve parameters
            var alignmentLength = outputData.alignmentLengths[i];
            var sequences = outputData.sequencePairs[i];
            var numOfGaps = outputData.gapNumbers[i];

            var scoreRandom = getExpectedScore(alignmentLength, sequences[0], sequences[1], numOfGaps);
            var scoreMax = getAverageMaximumScore(sequences[0], sequences[1]);

            var scoreEffective = 0;
            if (score - scoreRandom > 0)
                scoreEffective = (score - scoreRandom) / (scoreMax - scoreRandom);
            else
                scoreEffective = FENG_DOOLITTLE_CONSTANT / (scoreMax - scoreRandom);

            outputData.distances.push(-Math.log(scoreEffective));  // natural logarithm
        }
    }

    /**
     * Computes the expected score for aligning two sequences
     * by using the approximative formula from 1996 of Feng and Doolittle.
     * Hint: Usually random shuffling is used (1987),
     * but then the algorithm would become non-deterministic and it couldn't be tested.
     * Hint 2: The formula could be computed more efficient with character frequency tables,
     * but for their computation usually sets are needed which are not fully implemented in [the here used
     * version of] Javascript.
     * @param alignmentLength - The length of the alignment (different definitions possible).
     * @param sequenceA - The first (not aligned) sequence.
     * @param sequenceB - The second (not aligned) sequence.
     * @param numOfGaps - The number of gaps in the alignment of sequence a and b.
     * @see: https://doi.org/10.1016/S0076-6879(96)66023-6
     * Feng, Da-Fei and Doolittle, Russell F. «[21] Progressive alignment of amino
     * acid sequences and construction of phylogenetic trees from them». In: Methods in
     * enzymology 266 (1996), pp. 368–382
     * @example:
     * S^rand(a,b) = [1/L(a,b)] * [\sum_{i in a} \sum_{j in b} s(i,j) N_a(i) N_b(j)] - N_{a,b}("_")
     * Hint:
     * Usually a parameter d is multiplied with N_{a,b}("_") to avoid very high S^rand which
     * can lead to a negative S^eff(a,b) = [S(a,b) - S^rand(a,b)] / ...
     * So, the correct formula would be
     * S^rand(a,b) = [1/L(a,b)] * [\sum_{i in a} \sum_{j in b} s(i,j) N_a(i) N_b(j)] - N_{a,b}("_") * d
     * but for simplicity this parameter is removed and another approach used.
     * @return {number} - The expected score.
     */
    function getExpectedScore(alignmentLength, sequenceA, sequenceB, numOfGaps) {
        var doubleSum = 0;

        for (var i = 0; i < sequenceA.length; i++) {
            for (var j = 0; j < sequenceB.length; j++) {
                var aChar = sequenceA[i];
                var bChar = sequenceB[j];

                var similarity = aChar === bChar ? inputData.match : inputData.mismatch;
                var frequencyInA = getFrequency(aChar, sequenceA);
                var frequencyInB = getFrequency(bChar, sequenceB);
                doubleSum += similarity * frequencyInA * frequencyInB;
            }
        }

        return (1/alignmentLength) * doubleSum - numOfGaps;  // expected - parameter -> expected score with a shift
    }

    /**
     * Returns the absolute frequency of a char in a sequence.
     * @param char {string} - The char of which the frequency is computed.
     * @param sequence {string} - The sequence in which is counted.
     * @return {number} - The absolute frequency.
     */
    function getFrequency(char, sequence) {
        var charFrequency = 0;

        for (var i = 0; i < sequence.length; i++) {
            if (sequence[i] === char)
                charFrequency++;
        }

        return charFrequency;
    }

    /**
     * Computes the average score of aligning both sequences with themselves.
     * @param sequenceA {string} - The first sequence which was aligned.
     * @param sequenceB {string} - The second sequence which was aligned.
     * @example:
     * S^max(a,b) = [S(a,a) + S(b,b)] / 2
     */
    function getAverageMaximumScore(sequenceA, sequenceB) {
        return (getMaximumScore(sequenceA) + getMaximumScore(sequenceB)) / 2;
    }

    /**
     * Returns the maximum score for a sequence.
     * @param sequence {string} - The sequence for which the maximum score is computed.
     * @return {number} - The maximum score.
     * @example: S(a,a)
     */
    function getMaximumScore(sequence) {
        var score = 0;

        for (var i = 0; i < sequence.length; i++) {
            score += inputData.match;
        }

        return score;
    }

    /**
     * Creates dependency between cluster names and distances.
     * So, a function dist(a,b) which giving you
     * for two cluster-names a and b the distance
     * (needed for the clustering algorithm).
     * Hint: The matrix for visualization purposes is created
     * in the corresponding HTML and the corresponding interface files.
     * Hint 2: Diagonal is not needed because it is full of zeros and not needed.
     * @see: It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function createDistanceMatrix() {
        var clusterNames = getClusterNames();
        var nameBySequence = {};

        for (var i = 0; i < clusterNames.length; i++) {
            nameBySequence[outputData.sequences[i]] = clusterNames[i];
        }

        outputData.distanceMatrix = {};

        // right half upper diagonal
        for (var i = 0; i < outputData.distances.length; i++) {
            var sequencePair = outputData.sequencePairs[i];

            var firstClusterName = nameBySequence[sequencePair[0]];
            var secondClusterName = nameBySequence[sequencePair[1]];

            outputData.distanceMatrix[[firstClusterName, secondClusterName]] = outputData.distances[i];
        }
    }

    /**
     * Returns names for clusters associated with the distance data.
     * Hint: After all characters are depleted,
     * a number is concatenated to the character
     * to make this function generic.
     * @example:
     * CLUSTER NAMES:
     * a, b, c, ..., z,         FIRST EPISODE
     * a2, b2, c2, ..., z2,     SECOND EPISODE
     * a3, b3, ...              THIRD ...
     */
    function getClusterNames() {
        var clusterNames = [];
        var currentEpisode = 1;

        // for every pairwise distance we need a symbol
        for (var i = 0; i < outputData.distances.length; i++) {
            if (i < CLUSTER_NAMES.length)
                clusterNames.push(CLUSTER_NAMES[i]);  // add a, b, c, ..., z

            if (i >= CLUSTER_NAMES.length && i % CLUSTER_NAMES.length === 0)  // out of characters
                currentEpisode++;  // new episode

            if (i >= CLUSTER_NAMES.length)  // out of characters -> a2, b2, c2, ..., z2, a3, b3, ...
                clusterNames.push(CLUSTER_NAMES[i % CLUSTER_NAMES.length] + SYMBOLS.EMPTY + currentEpisode);
        }

        return clusterNames;
    }

    /**
     * Using a clustering algorithm like UPGMA (Group Average)
     * the algorithm computes a visualization of the distance matrix
     * in form of a binary guide tree.
     */
    function createPhylogeneticTree() {
    }

    /**
     * By going through a guide tree,
     * the algorithm generates progressive alignments.
     */
    function createProgressiveAlignments() {
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