/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("bases.multiSequenceAlignment", MultiSequenceAlignment, getInput, setInput,
        computePairwiseData, initializeInput, computeWithAlgorithm, getOutput, setIO, getLastChild);

    // instances
    var alignmentInstance;
    var childInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Contains functions to compute multi-sequence alignments.
     * It is used multi-sequence alignment algorithms as superclass
     * to avoid code duplicates.
     * @param child - The child algorithm which inherits from this class.
     * @constructor
     */
    function MultiSequenceAlignment(child) {
        alignmentInstance = this;
        childInstance = child;

        // public methods
        this.getInput = getInput;
        this.setInput = setInput;
        this.computePairwiseData = computePairwiseData;
        this.initializeInput = initializeInput;
        this.computeWithAlgorithm = computeWithAlgorithm;
        this.getOutput = getOutput;

        this.setIO = setIO;
        this.getLastChild = getLastChild;
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
    function setInput(inputViewmodel) {
        inputData.sequences = inputViewmodel.sequences();
        inputData.arrayPositionsOfRemovedSequences = getDuplicatePositions(inputData.sequences);

        inputData.calculationType = inputViewmodel.calculation();

        inputData.baseCosts = inputViewmodel.baseCosts();
        inputData.enlargement = inputViewmodel.enlargement();
        inputData.match = inputViewmodel.match();
        inputData.mismatch = inputViewmodel.mismatch();
    }

    /**
     * Returns the non-first position of sequences which are multiple times in the sequences.
     * @param sequences - The sequences which are checked for duplicates.
     */
    function getDuplicatePositions(sequences) {
        var positions = {};
        var duplicatePositions = [];
        var uniqueSequences = getUniqueElements(sequences);

        for (var i = 0; i < uniqueSequences.length; i++) {
            positions[uniqueSequences[i]] = getAllPositions(uniqueSequences[i], sequences);

            if (positions[uniqueSequences[i]].length > 1) {  // if a sequence contained multiple times in sequences
                // removing first position, because it's not needed (not removed)
                // and copy the other positions into duplicatePositions (to remove the associated sequences later)
                var nonFirstPositions = positions[uniqueSequences[i]].slice(1);  // first position removed and rest is copied
                duplicatePositions = duplicatePositions.concat(nonFirstPositions)
            }
        }

        return duplicatePositions;
    }

    /**
     * Returns unique elements.
     * @param array {Array} - The sequence in which is counted.
     * @return {Array} - The input array without duplicates.
     */
    function getUniqueElements(array) {
        var elements = [];

        for (var i = 0; i < array.length; i++) {
            if (elements.indexOf(array[i]) === -1)
                elements.push(array[i]);
        }

        return elements;
    }

    /**
     * Returns all positions of a given sequence in the sequences.
     * @param sequence {string} - The sequence to search for.
     * @param sequences {Array} - The array of sequences in which is searched.
     * @return {Array} - The array of all positions of the sequence in the sequence array.
     */
    function getAllPositions(sequence, sequences) {
        var positions = [];

        for (var i = 0; i < sequences.length; i++) {
            if (sequence === sequences[i])
                positions.push(i);
        }

        return positions;
    }

    /**
     * Computes scores (similarities),
     * the number of gaps, the alignment lengths
     * and so on between all sequences.
     * @param {Object} - The algorithm with which the alignment data is computed.
     */
    function computePairwiseData(algorithm) {
        var gotohInput = {};
        initializeInput(gotohInput);

        outputData.alignmentLengths = [];
        outputData.alignmentsAndScores = [];
        outputData.gapNumbers = [];
        outputData.gapStarts = [];
        outputData.sequencePairs = [];
        outputData.similarities = [];

        for (var i = 1; i < inputData.sequences.length; i++) {
            if (inputData.arrayPositionsOfRemovedSequences.indexOf(i) === -1) {  // only if the sequence is not a duplicate
                for (var j = 0; j < i; j++) {
                    var sequenceA = inputData.sequences[j];
                    var sequenceB = inputData.sequences[i];

                    var ioData = computeWithAlgorithm(algorithm, gotohInput, sequenceA, sequenceB);
                    var alignment = getAlignment(ioData);

                    outputData.alignmentsAndScores[[sequenceA, sequenceB]] = [alignment, ioData[1].score];  // for faster access

                    outputData.alignmentLengths.push(getAlignmentLength(alignment));
                    outputData.gapNumbers.push(getNumberOfGaps(alignment));
                    outputData.gapStarts.push(getNumberOfGapsStarts(alignment));
                    outputData.sequencePairs.push([sequenceA, sequenceB]);
                    outputData.similarities.push(ioData[1].score);
                }
            }
        }
    }

    /**
     * Returns the alignment given the input/output data object.
     * @param ioData {Object} - The input/output object
     * @return {Object} - The alignment.
     */
    function getAlignment(ioData) {
        return ioData[1].alignments[0] !== undefined ? ioData[1].alignments[0] : EMPTY_ALIGNMENT;
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
    }

    /**
     * Computes Gotoh with the given input sequences and the function.
     * @param {Object} - The algorithm with which the alignment data is computed.
     * @param input {Object} - The initialized Gotoh input structure.
     * @param sequenceA {string} - The first sequence.
     * @param sequenceA {string} - The second sequence.
     * @return {Object} - Output data of Gotoh with the given sequences in the input.
     */
    function computeWithAlgorithm(algorithm, input, sequenceA, sequenceB) {
        input.sequenceA = sequenceA;
        input.sequenceB = sequenceB;

        input.matrixHeight = input.sequenceB.length + 1;
        input.matrixWidth = input.sequenceA.length + 1;

        input.computeOneAlignment = true;  // speed up for Feng-Doolittle

        algorithm.setIO(input, {});

        return algorithm.compute();
    }

    /**
     * Computes the length of the alignment.
     * Hint: The alignment length can be defined differently.
     * It does not make sense to define it like in the AEP algorithm
     * as the number of characters of both aligned sequences (without gaps),
     * because here global alignments are used (GOTOH) and so the alignment length would be
     * always the length of both sequences (what does not make sense to call it then alignment length).
     * So, here the MSA-length definition used (Feng-Doolittle is a MSA algorithm):
     * length is the number of columns in the global alignment of the sequences.
     * Hint 2: Verified with original paper.
     * @param alignment {[alignedSequenceA, matchOrMismatchString, alignedSequenceB]}
     * - The triple of strings for which the length has to be computed.
     * @see: https://doi.org/10.1016/S0076-6879(96)66023-6
     * Feng, Da-Fei, and Russell F. Doolittle.
     * "[21] Progressive alignment of amino acid sequences and construction of phylogenetic trees from them."
     * Methods in enzymology 266 (1996): 368-382.
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
     * Returns the number of gap starts in the alignment.
     * @param alignment {[alignedSequenceA, matchOrMismatchString, alignedSequenceB]}
     * - The triple of strings for which the number of gaps has to be computed.
     * @return {number} - The number of gaps in the alignment.
     */
    function getNumberOfGapsStarts(alignment) {
        return countNumberOfGapsStarts(alignment[0]) + countNumberOfGapsStarts(alignment[2]);
    }

    /**
     * Counts the number of gap starts in the given sequence.
     * @param sequence {string} - The sequence for which the number of gaps has to be counted.
     * @return {number} - The number of gaps.
     */
    function countNumberOfGapsStarts(sequence) {
        var numGapsStarts = 0;

        for (var i = 0; i < sequence.length; i++) {
            if (sequence[i] === SYMBOLS.GAP) {
                numGapsStarts++;

                while (sequence[i] === SYMBOLS.GAP && i < sequence.length)  // iterate to the gap end
                    i++;
            }
        }

        return numGapsStarts;
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
}());
