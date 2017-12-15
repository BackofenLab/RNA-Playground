/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace ("getIndividualSequenceNames" is set static because creation of a full instance to get just the sequence names is too inefficient)
    namespace("bases.multiSequenceAlignment", MultiSequenceAlignment, getIndividualSequenceNames);

    // instances
    var multiSequenceAlignmentInstance;
    var alignmentInstance;
    var childInstance;
    var gotohInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Contains functions to compute multi-sequence alignments.
     * It is used by multi-sequence alignment algorithms as superclass
     * to avoid code duplicates.
     * @param child {Object} - The child algorithm which inherits from this class.
     * @constructor
     */
    function MultiSequenceAlignment(child) {
        multiSequenceAlignmentInstance = this;
        childInstance = child;
        alignmentInstance = new bases.alignment.Alignment(this);
        gotohInstance = new gotoh.Gotoh();

        // public class methods
        this.getInput = getInput;
        this.setInput = setInput;
        this.computePairwiseData = computePairwiseData;
        this.initializeInput = initializeInput;
        this.computeWithAlgorithm = computeWithAlgorithm;
        this.computeDistancesFromSimilarities = computeDistancesFromSimilarities;
        this.createDistanceMatrix = createDistanceMatrix;
        this.getPhylogeneticTree = getPhylogeneticTree;
        this.createProgressiveAlignment = createProgressiveAlignment;
        this.getBestAlignment = getBestAlignment;
        this.createGroup = createGroup;
        this.replaceGapsWithPlaceHolder = replaceGapsWithPlaceHolder;
        this.replacePlaceholdersWithGaps = replacePlaceholdersWithGaps;
        this.getAffineSumOfPairsScore = getAffineSumOfPairsScore;
        this.getIndividualSequenceNames = getIndividualSequenceNames;
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
        inputData.initialNamingIndex = inputData.sequences.length;

        inputData.calculationType = inputViewmodel.calculation();

        inputData.baseCosts = inputViewmodel.baseCosts();
        inputData.enlargement = inputViewmodel.enlargement();
        inputData.match = inputViewmodel.match();
        inputData.mismatch = inputViewmodel.mismatch();
    }

    /**
     * Returns the non-first position of sequences which are multiple times in the sequences.
     * @param sequences {Array} - The sequences which are checked for duplicates.
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
                duplicatePositions = duplicatePositions.concat(nonFirstPositions);
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
     * and so on between all sequences with the given algorithm.
     * Hint: If Gotoh (Local) computes the pairwise data, then the output data has a "Local" at the end of every parameter
     * and only alignments and scores are stored.
     * @param algorithm {Object} - The algorithm with which the pairwise data should be computed [Gotoh or Gotoh (Local)].
     */
    function computePairwiseData(algorithm) {
        var global = algorithm.type === ALGORITHMS.GOTOH;
        alignmentInstance.setLastChild(algorithm);  // or the Alignment will work with wrong child algorithm

        var algorithmInput = {};
        initializeInput(algorithmInput, global);

        if (global) {
            outputData.alignmentLengths = [];
            outputData.alignmentsAndScores = [];
            outputData.gapNumbers = [];
            outputData.gapStarts = [];
            outputData.sequencePairs = [];
            outputData.similarities = [];
        } else {  // local
            outputData.alignmentsAndScoresLocal = [];
            outputData.tracebacks = [];
            outputData.maxNumberOfAlignments = 0;
        }

        for (var j = 1; j < inputData.sequences.length; j++) {
            if (inputData.arrayPositionsOfRemovedSequences.indexOf(j) === -1) {  // only if the sequence is not a duplicate
                for (var i = 0; i < j; i++) {
                    var sequenceA = inputData.sequences[i];
                    var sequenceB = inputData.sequences[j];

                    var ioData = computeWithAlgorithm(algorithm, algorithmInput, {}, sequenceA, sequenceB);
                    var alignment = getAlignment(ioData);

                    if (global) {
                        // for faster access: hash-table
                        outputData.alignmentsAndScores[[sequenceA, sequenceB]]
                            = [alignment, ioData[1].score, ioData[1].alignments];  // extension for T-Coffee: all alignments stored

                        outputData.alignmentLengths.push(getAlignmentLength(alignment));
                        outputData.gapNumbers.push(getNumberOfGaps(alignment));
                        outputData.gapStarts.push(getNumberOfGapsStarts(alignment));
                        outputData.sequencePairs.push([sequenceA, sequenceB]);
                        outputData.similarities.push(ioData[1].score);
                    } else {  // local
                        outputData.alignmentsAndScoresLocal[[sequenceA, sequenceB]]
                            = [alignment, ioData[1].score, ioData[1].alignments];  // extension for T-Coffee: all alignments stored

                        outputData.tracebacks[[sequenceA, sequenceB]] = ioData[1].tracebackPaths;
                        outputData.maxNumberOfAlignments++;
                    }
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
     * @param global {boolean} - Tells if global data should be computed or local.
     */
    function initializeInput(input, global) {
        if (global) {
            input.baseCosts = inputData.baseCosts;
            input.enlargement = inputData.enlargement;
            input.match = inputData.match;
            input.mismatch = inputData.mismatch;
        } else {
            input.baseCosts = inputData.baseCostsLocal;
            input.enlargement = inputData.enlargementLocal;
            input.match = inputData.matchLocal;
            input.mismatch = inputData.mismatchLocal;
        }

        input.calculationType = inputData.calculationType;
        input.substitutionFunction = inputData.substitutionFunction;  // extension for T-Coffee
    }

    /**
     * Computes Gotoh or Gotoh (Local) with the given input sequences and the function.
     * @param algorithm {Object} - The algorithm with which the alignment data is computed.
     * @param input {Object} - The initialized Gotoh input structure.
     * @param output {Object} - The initialized Gotoh output structure.
     * @param sequenceA {string} - The first sequence.
     * @param sequenceB {string} - The second sequence.
     * @return {Object} - Output data of Gotoh with the given sequences in the input.
     */
    function computeWithAlgorithm(algorithm, input, output, sequenceA, sequenceB) {
        input.sequenceA = sequenceA;
        input.sequenceB = sequenceB;

        input.matrixHeight = input.sequenceA.length + 1;
        input.matrixWidth = input.sequenceB.length + 1;

        if (algorithm.type === ALGORITHMS.FENG_DOOLITTLE)
            input.computeOneAlignment = true;  // speed up for Feng-Doolittle and maybe later other multi-alignment algorithms
        else {
            input.computeOneAlignment = false;
            input.numberOfAlignmentsPerSequencePair = inputData.globalAlignmentsPerSequencePair;
            input.numberOfLocalAlignmentsPerSequencePair = inputData.localAlignmentsPerSequencePair;
        }
        input.recomputeTraceback = output.tracebackPaths === undefined;
        algorithm.setIO(input, output);

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
     * @param alignment {Object}
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
     * @param alignment {Object}
     * - The triple of strings for which the number of gaps has to be computed.
     * @return {number} - The number of gaps in the alignment.
     */
    function getNumberOfGaps(alignment) {
        return countNumberOfGaps(alignment[0]) + countNumberOfGaps(alignment[2]);
    }

    /**
     * Counts the number of gaps in the given sequence.
     * @param sequence {alignedSequenceA|matchOrMismatchString|alignedSequenceB} - The sequence for which the number of gaps has to be counted.
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
     * @param alignment {Object}
     * - The triple of strings for which the number of gaps has to be computed.
     * @return {number} - The number of gaps in the alignment.
     */
    function getNumberOfGapsStarts(alignment) {
        return countNumberOfGapsStarts(alignment[0]) + countNumberOfGapsStarts(alignment[2]);
    }

    /**
     * Counts the number of gap starts in the given sequence.
     * @param sequence {alignedSequenceA|matchOrMismatchString|alignedSequenceB} - The sequence for which the number of gaps has to be counted.
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
     * @param alignmentLength {number} - The length of the alignment (number of columns).
     * @param sequenceA {string} - The first (not aligned) sequence.
     * @param sequenceB {string} - The second (not aligned) sequence.
     * @param numOfGaps {number} - The number of gaps in the alignment of sequence a and b.
     * @param numOfGapStarts {number} - The number of gap starts in the alignment of sequence a and b.
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

        return (1 / alignmentLength) * doubleSum + numOfGaps * inputData.enlargement + numOfGapStarts * inputData.baseCosts;
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

        outputData.nameOfSequence = {};

        var currentNamePos = 0;
        for (var i = 0; i < inputData.sequences.length; i++) {
            if (inputData.arrayPositionsOfRemovedSequences.indexOf(i) === -1) {  // if sequence is a used sequence
                outputData.nameOfSequence[inputData.sequences[i]] = outputData.clusterNames[currentNamePos];
                currentNamePos++;
            }
        }

        outputData.distanceMatrix = {};
        outputData.sequencePairNames = [];

        // right half upper diagonal
        for (var i = 0; i < outputData.distances.length; i++) {
            var sequencePair = outputData.sequencePairs[i];

            var firstClusterName = outputData.nameOfSequence[sequencePair[0]];
            var secondClusterName = outputData.nameOfSequence[sequencePair[1]];

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
     * Using a clustering algorithm like UPGMA (Group Average)
     * the algorithm returns the binary guide tree branches in creation order.
     * @return {Object} - The tree branches.
     */
    function getPhylogeneticTree() {
        inputData.numOfStartClusters = inputData.sequences.length - inputData.arrayPositionsOfRemovedSequences.length;

        var clustering = new agglomerativeClustering.AgglomerativeClustering();
        inputData.clusteringSubalgorithm = CLUSTERING_ALGORITHMS.UPGMA;
        clustering.setIO(inputData, outputData);
        var ioData = clustering.compute();
        return ioData[1].treeBranches;
    }

    /**
     * By going through the guide tree branches (in correct merging order),
     * the algorithm generates progressive alignments.
     * Hint: A repeated post order traversal of the tree would be less efficient.
     * This is why just an iteration through the branches is done.
     * @param treeBranches {Object} - The tree branches which are defining the order for the merging process.
     * @see: It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function createProgressiveAlignment(treeBranches) {
        initializeGroups();

        for (var i = 0; i < treeBranches.length; i++) {
            var treeBranch = treeBranches[i];

            var leftChildName = treeBranch.rightChild.name;
            var rightChildName = treeBranch.leftChild.name;
            var groupName = treeBranch.name;

            alignGroups(leftChildName, rightChildName, groupName);
        }

        var alignmentWithPlaceHolders = outputData.groups[outputData.allClusterNames[outputData.allClusterNames.length - 1]];
        outputData.progressiveAlignment = replacePlaceholdersWithGaps(alignmentWithPlaceHolders);
    }

    /**
     * Initializes the groups-structure,
     * which is storing the created groups during progressive alignment
     * with the start sequences.
     */
    function initializeGroups() {
        outputData.groups = {};

        var currentNamePos = 0;
        for (var i = 0; i < inputData.sequences.length; i++) {
            if (inputData.arrayPositionsOfRemovedSequences.indexOf(i) === -1) {
                outputData.groups[outputData.clusterNames[currentNamePos]] = [inputData.sequences[i]];
                currentNamePos++;
            }
        }

        // needed for visualization
        outputData.firstGroups = [];
        outputData.firstGroupsNames = [];
        outputData.secondGroups = [];
        outputData.secondGroupsNames = [];
        outputData.guideAlignments = [];
        outputData.guideAlignmentsNames = [];
        outputData.joinedGroups = [];
        outputData.joinedGroupNames = [];

        // T-Coffee extension: substitution values are dependant on the groups -> sequences read by substitution function
        outputData.currentFirstGroup = [];  // needed to recompute substitution values for column
        outputData.currentSecondGroup = [];  // needed to recompute substitution values for column
        outputData.currentGroupName = [];  // extension for T-Coffee: to avoid recomputation of matrices
        outputData.groupMatrices = {};  // extension for T-Coffee: store previously computed matrices
        outputData.groupTracebackPaths = {};  // extension for T-Coffee: store previously computed traceback paths
    }

    /**
     * Aligns two groups or sequences to each other.
     * Groups are aligned by aligning every member of the one group
     * with very member of the other group.
     * Afterwards the best alignment is chosen and used as a guide alignment
     * to which the other group members are added,
     * such that characters which were aligned before in the group block,
     * are still together.
     * After the two groups were merged together,
     * all gaps SYMBOLS.GAP are replaced with a placeholder element SYMBOLS.NONE,
     * to preserve the gap (avoid gap movement).
     * @param leftChildName {string} - The name of the left child group or sequence.
     * @param rightChildName {string} - The name of the right child group or sequence
     * @param groupName {string} - The name of the new group.
     * Hint: For the sequence alignment the previously computed sequence-scores are used.
     */
    function alignGroups(leftChildName, rightChildName, groupName) {
        var group1Sequences = outputData.groups[leftChildName];
        var group2Sequences = outputData.groups[rightChildName];
        outputData.currentGroupName = groupName;  // extension for T-Coffee: look for #optimization

        var bestAlignment = getBestAlignment(group1Sequences, group2Sequences);

        if (childInstance.type === ALGORITHMS.NOTREDAME_HIGGINS_HERINGA)  // extension for T-Coffee
            outputData.groupMatrices[groupName] = outputData.currentMatrix;  // store previously calculated matrices (for a Unit-Test)

        outputData.groups[groupName] = createGroup(group1Sequences, group2Sequences, bestAlignment);

        // for visualization of steps
        outputData.firstGroups.push(group1Sequences);
        outputData.firstGroupsNames.push(leftChildName);
        outputData.secondGroups.push(group2Sequences);
        outputData.secondGroupsNames.push(rightChildName);
        outputData.guideAlignments.push(bestAlignment);
        outputData.joinedGroups.push(outputData.groups[groupName]);
        outputData.joinedGroupNames.push(groupName);
    }

    /**
     * Returns the best alignment (with the highest score) after
     * a pair-wise alignment of the group sequences of both groups
     * or after a look-up in the previously computed alignments for the start-sequences.
     * @param group1Sequences {Array} - The sequences of the first group.
     * @param group2Sequences {Array} - The sequences of the second group.
     * @return {[alignedSequenceA, matchOrMismatchString, alignedSequenceB]} - A triple of strings.
     */
    function getBestAlignment(group1Sequences, group2Sequences) {
        var maxScore = Number.NEGATIVE_INFINITY;
        var maxAlignment = [];
        var maxSequence1 = SYMBOLS.EMPTY;
        var maxSequence2 = SYMBOLS.EMPTY;

        // T-Coffee extension: substitution values are dependant on the groups -> sequences read by substitution function
        outputData.currentFirstGroup = group1Sequences;
        outputData.currentSecondGroup = group2Sequences;

        // iterate through all sequences and search for an alignment with maximum similarity
        for (var i = 0; i < group1Sequences.length; i++) {
            var sequence1 = group1Sequences[i];

            for (var j = 0; j < group2Sequences.length; j++) {
                var sequence2 = group2Sequences[j];
                var asData = getAlignmentAndScore(sequence1, sequence2);

                if (asData[1] > maxScore) {
                    maxScore = asData[1];
                    maxAlignment = asData[0];

                    maxSequence1 = sequence1;
                    maxSequence2 = sequence2;
                }
            }
        }

        if (outputData.nameOfSequence !== undefined) {  // not all multi-sequence-alignment algorithms need this
            var firstSequence
                = outputData.nameOfSequence[maxSequence1.replace(MULTI_SYMBOLS.NONE, SYMBOLS.EMPTY).replace(MULTI_SYMBOLS.GAP, SYMBOLS.EMPTY)];
            var secondSequence
                = outputData.nameOfSequence[maxSequence2.replace(MULTI_SYMBOLS.NONE, SYMBOLS.EMPTY).replace(MULTI_SYMBOLS.GAP, SYMBOLS.EMPTY)];

            outputData.guideAlignmentsNames.push(firstSequence + SYMBOLS.ALIGN + secondSequence);
        }

        return maxAlignment;
    }

    /**
     * Returns the alignment and score.
     * @param sequence1 {string} - The first sequence.
     * @param sequence2 {string} - The second sequence.
     * @return {[string, number]} - The alignment and score.
     */
    function getAlignmentAndScore(sequence1, sequence2) {
        if (outputData.alignmentsAndScores !== undefined) {  // optimization for Feng-Doolittle
            var alignmentAndScore = outputData.alignmentsAndScores[[sequence1, sequence2]];  // constant time!

            if (alignmentAndScore !== undefined && childInstance.type !== ALGORITHMS.NOTREDAME_HIGGINS_HERINGA)  // T-Coffee extension
                return alignmentAndScore;
        }

        var input = {};
        var output = {};
        initializeInput(input, true);

        if (childInstance.type === ALGORITHMS.NOTREDAME_HIGGINS_HERINGA)  // T-Coffee extension
            initializeOutput(output);

        var ioData = computeWithAlgorithm(gotohInstance, input, output, sequence1, sequence2);

        if (childInstance.type === ALGORITHMS.NOTREDAME_HIGGINS_HERINGA) {  // optimization for T-Coffee
            // extension for T-Coffee: needed for a Unit-Test
            // if we have a recalculation (look at #optimization), then the matrix is not computed and undefined, we use the last one
            outputData.currentMatrix = ioData[1].matrix === undefined ? outputData.currentMatrix : ioData[1].matrix;

            // #optimization:
            // traceback paths are computed once for a group and then reused, because they do not change anymore
            // look into Unit-Test: Notredame-Higgins-Heringa.pdf -> matrix ab~c on page 6
            outputData.groupTracebackPaths[outputData.currentGroupName] = ioData[1].tracebackPaths;
        }

        return [ioData[1].alignments[0], ioData[1].score];
    }

    /**
     * If traceback path and matrix have been already computed,
     * then the traceback path is not recomputed.
     * @param output {Object} - Output data which have to be initialized.
     */
    function initializeOutput(output) {
        var parentGroupName = outputData.currentGroupName;
        output.tracebackPaths = outputData.groupTracebackPaths[parentGroupName];

        if (output.tracebackPaths !== undefined && output.tracebackPaths[0] !== undefined)
            output.tracebackPaths[0].reverse();  // because the createAlignment-function in Alignment has reversed the array for the later visualization
    }

    /**
     * Creates a group out of the sequences with the help of the guide-alignment
     * and replaces all gaps SYMBOLS.GAP with the placeholder symbol SYMBOLS.NONE.
     * @param group1Sequences {Array} - The sequences of the first group.
     * @param group2Sequences {Array} - The sequences of the second group.
     * @param guideAlignment {Object} - The alignment which is used to align the other group members.
     * @return {Object} - The MSA group (MSA alignment).
     * @see: It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function createGroup(group1Sequences, group2Sequences, guideAlignment) {
        var guideSequence1 = guideAlignment[0];
        var guideSequence2 = guideAlignment[2];

        if (group1Sequences.length === 1 && group2Sequences.length === 1)  // if already aligned, because only two sequences
            return replaceGapsWithPlaceHolder([guideSequence1, guideSequence2]);

        // first all group members of group1Sequences are added to the guideSequence1
        // and then all group members of group2Sequences are added to the guideSequence2
        // and afterwards both groups are joined
        var firstGroup = getGuidedGroup(group1Sequences, guideSequence1);
        var secondGroup = getGuidedGroup(group2Sequences, guideSequence2);

        return replaceGapsWithPlaceHolder(firstGroup.concat(secondGroup));
    }

    /**
     * Adds the group sequences appropriately to the guide sequence and returns the group.
     * @param groupSequences {Array} - The sequences which should be added to the guide sequence.
     * @param guideSequence {Array} - The sequence which is used to align the sequences of the group.
     * @return {Array} - The group in which group sequences have been added to the guide sequence.
     * @see: It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function getGuidedGroup(groupSequences, guideSequence) {
        if (groupSequences.length === 1)  // if only one element in group
            return [guideSequence];

        var currentPosition = 0;
        var alignedSequence = SYMBOLS.EMPTY;
        var guidedGroup = [];

        // going through each member of the group
        // and add new gaps into the group-sequence
        // accordingly to the "new" gaps in the guide sequence
        for (var i = 0; i < groupSequences.length; i++) {
            var currentSequence = groupSequences[i];

            // going through the guide-sequence and adding step by step a new symbol
            for (var j = 0; j < guideSequence.length; j++) {
                var symbol = guideSequence[j];

                if (symbol !== SYMBOLS.GAP) {
                    alignedSequence += currentSequence[currentPosition];
                    currentPosition++;
                } else {
                    alignedSequence += SYMBOLS.GAP;
                }
            }

            guidedGroup.push(alignedSequence);  // adding aligned groupSequences[i] to the group
            currentPosition = 0;
            alignedSequence = SYMBOLS.EMPTY;
        }

        return guidedGroup;
    }

    /**
     * Returns an array of sequences in which the gaps SYMBOLS.GAP of sequences are replaced with placeholders SYMBOLS.NONE.
     * @param group {Array} - The group in which gaps are replaced with placeholders.
     * @return {Array} - The array of aligned sequences in which the gaps are replaced by placeholders.
     */
    function replaceGapsWithPlaceHolder(group) {
        for (var i = 0; i < group.length; i++) {
            group[i] = group[i].replace(MULTI_SYMBOLS.GAP, SYMBOLS.NONE);
        }

        return group;
    }

    /**
     * Returns an array of sequences in which the placeholders SYMBOLS.NONE of sequences are replaced with gaps SYMBOLS.GAP.
     * @param group {Array} - The group in which placeholders are replaced with gaps.
     * @return {Array} - The array of aligned sequences in which the gaps are replaced by placeholders.
     */
    function replacePlaceholdersWithGaps(group) {
        for (var i = 0; i < group.length; i++) {
            group[i] = group[i].replace(MULTI_SYMBOLS.NONE, SYMBOLS.GAP);
        }

        return group;
    }

    /**
     * Computes the affine sum-of-pairs score for a given alignment.
     * The Gotoh-score between every alignment is computed.
     * These pairwise scores are summed up.
     * Hint: The score for overlying gaps is zero.
     * @param inputData {Object} - Contains all input data.
     * @param alignment {Object} - The alignment for which the sum-of-pairs score is computed.
     * @return {number} - The sum-of-pairs score.
     */
    function getAffineSumOfPairsScore(inputData, alignment) {
        var score = 0;

        // every alignment with every other alignment
        for (var i = 1; i < alignment.length; i++) {
            for (var j = 0; j < i; j++) {

                var sequences = removeOverlyingGaps(alignment[j], alignment[i]);

                var sequenceA = sequences[0];
                var sequenceB = sequences[1];
                var gapSize = 0;

                for (var k = 0; k < sequenceB.length; k++) {
                    if (sequenceB[k] === SYMBOLS.GAP) {
                        gapSize = getGapSize(k, sequenceB);
                        score += inputData.baseCosts + gapSize * inputData.enlargement;
                        k += gapSize - 1;
                    } else if (sequenceA[k] === SYMBOLS.GAP) {
                        gapSize = getGapSize(k, sequenceA);
                        score += inputData.baseCosts + gapSize * inputData.enlargement;
                        k += gapSize - 1;
                    } else
                        score += sequenceB[k] === sequenceA[k] ? inputData.match : inputData.mismatch;
                }
            }
        }

        return score;
    }

    /**
     * Removes overlying gaps in sequences, because they are not counted.
     * @param sequence1 {string} - The first sequence.
     * @param sequence2 {string} - The second sequence.
     * @return {[string, string]} - The alignment without overlying gaps.
     */
    function removeOverlyingGaps(sequence1, sequence2) {
        var newSequence1 = SYMBOLS.EMPTY;
        var newSequence2 = SYMBOLS.EMPTY;

        for (var i = 0; i < sequence1.length; i++) {
            if (sequence1[i] === sequence2[i] && sequence1[i] === SYMBOLS.GAP)  // if overlying gaps
                continue;

            newSequence1 += sequence1[i];
            newSequence2 += sequence2[i];
        }

        return [newSequence1, newSequence2];
    }

    /**
     * Returns the size of gap.
     * @param gapStartPosition {number} - The start position.
     * @param sequence {string} - The sequence with the gap.
     * @return {number} - The gap size.
     */
    function getGapSize(gapStartPosition, sequence) {
        var gapSize = 0;

        var i = gapStartPosition;

        while (i < sequence.length && sequence[i] === SYMBOLS.GAP) {
            gapSize++;
            i++;
        }

        return gapSize;
    }

    /**
     * Returns the names of the individual group members.
     * @param groupName {string} - The group name which encoding the group members.
     * @param commaSeparated {boolean} - Tells if between a name-string and a number should be a comma or not.
     * @return {Array} - The array with the group members.
     */
    function getIndividualSequenceNames(groupName, commaSeparated) {
        var names = [];

        if (groupName !== undefined) {
            for (var i = 0; i < groupName.length; i++) {
                var character = groupName[i];
                var number = SYMBOLS.EMPTY;

                while (i + 1 < groupName.length && groupName[i + 1].match(CHARACTER.NUMBER)) {
                    number += groupName[i + 1];
                    i++;
                }

                if (commaSeparated)
                    names.push(number.length > 0 ? character + SYMBOLS.COMMA + number : character);
                else
                    names.push(number.length > 0 ? character + number : character);
            }
        }

        return names;
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
