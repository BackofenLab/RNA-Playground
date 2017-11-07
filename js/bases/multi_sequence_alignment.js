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
        computePairwiseData, initializeInput, computeWithAlgorithm, computeDistancesFromSimilarities,
        createDistanceMatrix, getOutput, setIO, getLastChild);

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
        this.computeDistancesFromSimilarities = computeDistancesFromSimilarities;
        this.createDistanceMatrix = createDistanceMatrix;
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
     * Converting the pairwise similarities into distances
     * by using the Feng-Doolittle formulas.
     * Hint: The factor 100 is really
     * inside the logarithm to scale S^eff(a,b) between (0,100].
     * Hint 2: Verified with original paper and several sources.
     * @example:
     * D(a,b) = -ln(S^eff(a,b) * 100)
     * where
     * S^eff(a,b) = [S(a,b) - S^rand(a,b)] / [S^max(a,b) - S^rand(a,b)]
     *
     * Hint: The factor 100 was omitted in the computation (to avoid negative distances),
     * because it works only with very big S^max(a,b) (very long sequences).
     * Also, the real Feng-Doolittle uses substitution matrices instead of a simple scoring-function.
     *
     * @see: https://doi.org/10.1007/PL00006155
     * Feng, Da-Fei, and Russell F. Doolittle.
     * "Converting amino acid alignment scores into measures of evolutionary time:
     * a simulation study of various relationships." Journal of molecular evolution 44.4 (1997): 361-370.
     *
     * S^eff(a,b) should scale between [0,1]
     * where 1 means identical and 0 means totally non-identical,
     * but in reality S^rand can be bigger as S(a,b) and we can get negative values of S^eff(a,b).
     * To avoid negative values the approach from the paper linked above is used.
     *
     * The idea:
     * if [S(a,b) - S^rand(a,b)] < 0
     * then you set [S(a,b) - S^rand(a,b)] = 0.001 (other values are allowed)
     *
     * Hint: for [S(a,b) - S^rand(a,b)] == 0, it was not defined,
     * but for simplicity it is used [S(a,b) - S^rand(a,b)] = FENG_DOOLITTLE_CONSTANT
     *
     * Hint 3: Because of floating point errors,
     * it is also maybe possible, that [S^max(a,b) - S^rand(a,b)] <= 0.
     * In this case [S^max(a,b) - S^rand(a,b)] is set to FENG_DOOLITTLE_CONSTANT
     * ([S^max(a,b) - S^rand(a,b)] <= 0 is mathematically not possible,
     * because duplicate sequences are removed)
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
            var numOfGapStarts = outputData.gapStarts[i];

            var scoreRandom = getApproximatedRandomScore(alignmentLength, sequences[0], sequences[1], numOfGaps, numOfGapStarts);
            var scoreMax = getAverageMaximumScore(sequences[0], sequences[1]);

            var dividend = score - scoreRandom;
            var divisor = scoreMax - scoreRandom;

            var scoreEffective = 0;
            if (dividend <= 0) dividend = FENG_DOOLITTLE_CONSTANT;
            if (divisor <= 0) divisor = FENG_DOOLITTLE_CONSTANT;  // mathematically not possible, but because of FLOATING POINT ERRORS

            scoreEffective = dividend / divisor;

            outputData.distances.push(-Math.log(scoreEffective));  // natural logarithm
        }
    }

    /**
     * Computes an approximation of the random score for aligning two sequences
     * by using the approximative formula from 1996 of Feng and Doolittle.
     * Hint: Usually random shuffling is used (1987),
     * but then the algorithm would become non-deterministic and it couldn't be tested so easy.
     * @param alignmentLength - The length of the alignment (number of columns).
     * @param sequenceA - The first (not aligned) sequence.
     * @param sequenceB - The second (not aligned) sequence.
     * @param numOfGaps - The number of gaps in the alignment of sequence a and b.
     * @param numOfGapStarts - The number of gap starts in the alignment of sequence a and b.
     * @see: https://doi.org/10.1016/S0076-6879(96)66023-6
     * Feng, Da-Fei and Doolittle, Russell F. «[21] Progressive alignment of amino
     * acid sequences and construction of phylogenetic trees from them».
     * In: Methods in enzymology 266 (1996), pp. 368–382
     *
     * Hint 2: The penalty is a negated value in the paper!
     * So, for example a PAM250-matrix has a positive penalty like d=8.
     * This is the reason for writing a "+" instead of "-" in the formula:
     * ... + N_{a,b}(gaps) * \beta
     *
     * Hint 3: Mismatches s(i,j) have been omitted in the original formula! Probably this was a mistake.
     *
     * @example:
     * S^rand(a,b)
     * = [1/L(a,b)] * [\sum_{i in A(a,b)} \sum_{j in A(b)} s(i,i) N_a(i) N_b(j)]
     * + N_{a,b}(gaps) * beta
     * + N_{a,b}(gap-starts) * alpha
     * @return {number} - The expected score.
     */
    function getApproximatedRandomScore(alignmentLength, sequenceA, sequenceB, numOfGaps, numOfGapStarts) {
        var doubleSum = 0;

        var aChars = getUniqueChars(sequenceA);
        var bChars = getUniqueChars(sequenceB);

        for (var i = 0; i < aChars.length; i++) {
            var aChar = aChars[i];

            for (var j = 0; j < bChars.length; j++) {
                var bChar = bChars[j];

                var similarity = aChar === bChar ? inputData.match : inputData.mismatch;
                var frequencyInA = getFrequency(aChar, sequenceA);
                var frequencyInB = getFrequency(bChar, sequenceB);
                doubleSum += similarity * frequencyInA * frequencyInB;
            }
        }

        return (1/alignmentLength) * doubleSum + numOfGaps * inputData.enlargement + numOfGapStarts * inputData.baseCosts;
    }

    /**
     * Returns unique characters of both sequences.
     * @param sequence {string} - The first sequence in which is searched.
     * @return {Array} - Array of characters without duplicates.
     */
    function getUniqueChars(sequence) {
        var chars = [];

        for (var i = 0; i < sequence.length; i++) {
            if (chars.indexOf(sequence[i]) === -1)
                chars.push(sequence[i]);
        }

        return chars;
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
        return inputData.match * sequence.length;
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
        outputData.clusterNames = getClusterNames();

        var names = {};

        var currentNamePos = 0;
        for (var i = 0; i < inputData.sequences.length; i++) {
            if (inputData.arrayPositionsOfRemovedSequences.indexOf(i) === -1) {  // if sequence is a used sequence
                names[inputData.sequences[i]] = outputData.clusterNames[currentNamePos];
                currentNamePos++;
            }
        }

        outputData.distanceMatrix = {};
        outputData.sequencePairNames = [];

        // right half upper diagonal
        for (var i = 0; i < outputData.distances.length; i++) {
            var sequencePair = outputData.sequencePairs[i];

            var firstClusterName = names[sequencePair[0]];
            var secondClusterName = names[sequencePair[1]];

            outputData.distanceMatrix[[firstClusterName, secondClusterName]] = outputData.distances[i];
            outputData.sequencePairNames.push([firstClusterName, secondClusterName]);
        }

        outputData.distanceMatrixLength = inputData.sequences.length - inputData.arrayPositionsOfRemovedSequences.length;
    }

    /**
     * Returns names for clusters associated with the distance data.
     * Hint: After all characters are depleted,
     * a number is concatenated to the character
     * to make this function generic.
     * Hint: Duplicates do not get a name
     * by checking the array of removed positions (if not removed, then create name).
     * @example:
     * CLUSTER NAMES:
     * a, b, c, ..., z,         FIRST EPISODE
     * a2, b2, c2, ..., z2,     SECOND EPISODE
     * a3, b3, ...              THIRD ...
     * @return {Array} - The cluster names.
     */
    function getClusterNames() {
        var clusterNames = [];
        var currentEpisode = 1;

        // for every pairwise distance we need a symbol
        for (var i = 0; i < inputData.sequences.length; i++) {
            if (i < CLUSTER_NAMES.length && inputData.arrayPositionsOfRemovedSequences.indexOf(i) === -1)
                clusterNames.push(CLUSTER_NAMES[i]);  // add a, b, c, ..., z

            if (i >= CLUSTER_NAMES.length && i % CLUSTER_NAMES.length === 0)  // out of characters
                currentEpisode++;  // new episode

            // out of characters -> a2, b2, c2, ..., z2, a3, b3, ...
            if (i >= CLUSTER_NAMES.length && inputData.arrayPositionsOfRemovedSequences.indexOf(i) === -1)
                clusterNames.push(CLUSTER_NAMES[i % CLUSTER_NAMES.length] + SYMBOLS.EMPTY + currentEpisode);
        }

        return clusterNames;
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
