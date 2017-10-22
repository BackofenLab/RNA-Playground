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
    }

    /**
     * Starts the computation.
     */
    function compute() {
        computePairwiseData();
        convertSimilaritiesToDistances();
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
                outputData.alignmentLengths.push(getAlignmentLength(ioData[1].alignments[0]));
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
     * Computes the length of the alignment.
     * @see: The alignment length can be defined differently.
     * It does not make sense to define it like in the AEP algorithm
     * as the number of characters of both aligned sequences (without gaps),
     * because here global alignments are used (GOTOH) and so the alignment length would be
     * always the length of both sequences (what does not make sense to call it then alignment length).
     * So, here the MSA-length definition used (because Feng-Doolittle is a MSA algorithm):
     * length is the number of columns in the global alignment of the sequences.
     * @param alignment {[alignedSequenceA, matchOrMismatchString, alignedSequenceB]}
     * - The triple of strings for which the length has to be computed.
     * @return {number} - The length of the alignment.
     */
    function getAlignmentLength(alignment) {
        return alignment[0].length;
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
     * @param sequence - The sequence for which the number of gaps has to be counted.
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
     * @example:
     * D(a,b) = -ln(S^eff(a,b))
     * where
     * S^eff(a,b) = [S(a,b) - S^rand(a,b)] / [S^max(a,b) - S^rand(a,b)]
     */
    function convertSimilaritiesToDistances() {
        outputData.distances = [];
    }

    /**
     * Computes the expected score for aligning two sequences
     * by using the approximative formula from 1996 of Feng and Doolittle.
     * Hint: Usually random shuffling is used (1987),
     * but then the algorithm would become non-deterministic and it couldn't be tested.
     * @param alignmentLength - The length of the alignment (different definitions possible).
     * @param sequenceA - The first sequence which was aligned.
     * @param sequenceB - The second sequence which was aligned.
     * @param charFreqTableA - The frequency table for characters in sequence a.
     * @param charFreqTableB - The frequency table for characters in sequence b.
     * @param numOfGaps - The number of gaps in the alignment of sequence a and b.
     * @see:
     * Feng, Da-Fei and Doolittle, Russell F. «[21] Progressive alignment of amino
     * acid sequences and construction of phylogenetic trees from them». In: Methods in
     * enzymology 266 (1996), pp. 368–382
     * @example:
     * S^rand(a,b) = [1/L(a,b)] * [\sum_{i in a} \sum_{j in b} s(i,j) N_a(i) N_b(j)] - N_{a,b}("_")
     */
    function getExpectedScore(alignmentLength, sequenceA, sequenceB, charFreqTableA, charFreqTableB, numOfGaps) {
    }

    /**
     * Computes the average score of aligning both sequences with themselves.
     * @example:
     * S^max(a,b) = [S(a,a) + S(b,b)] / 2
     */
    function getMaximumScore(sequenceA, sequenceB) {
    }

    /**
     * Creates a distance matrix for computed pairwise scores.
     */
    function createDistanceMatrix() {
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