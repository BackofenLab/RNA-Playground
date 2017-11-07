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
    namespace("fengDoolittle", startFengDoolittle, FengDoolittle, getInput, setInput, compute, computePairwiseData,
        getOutput, setIO, getSuperclass);

    // instances
    var multiSequenceAlignmentInstance;

    var fengDoolittleInstance;
    var gotohInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Function managing objects.
     */
    function startFengDoolittle() {
        var multiSequenceInterface = new interfaces.multiSequenceInterface.MultiSequenceInterface();
        multiSequenceInterface.startMultiSequenceInterface(FengDoolittle, ALGORITHMS.FENG_DOOLITTLE);
    }

    /*---- ALGORITHM ----*/
    /**
     * Computes global affine multi-alignments (non-optimal approach).
     * Extended (= affine) version of the original Feng-Doolittle approaches
     * between 1986 and 1997.
     * @constructor
     * @augments MultiSequenceAlignment
     */
    function FengDoolittle() {
        fengDoolittleInstance = this;

        // variables
        this.type = ALGORITHMS.FENG_DOOLITTLE;

        // instances (do not change order)
        multiSequenceAlignmentInstance = new bases.multiSequenceAlignment.MultiSequenceAlignment(this);
        gotohInstance = new gotoh.Gotoh();

        // public class methods
        this.getInput = getInput;

        this.setInput = setInput;
        this.compute = compute;
        this.computePairwiseData = computePairwiseData;
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
        multiSequenceAlignmentInstance.setIO(inputData, {});
        multiSequenceAlignmentInstance.setInput(inputViewmodel);
    }

    /**
     * Starts the computation.
     */
    function compute() {
        multiSequenceAlignmentInstance.setIO(inputData, outputData);

        computePairwiseData();
        computeDistancesFromSimilarities();
        createDistanceMatrix();
        createProgressiveAlignment(getPhylogeneticTree());

        outputData.score = formats.scoringFunctions.getAffineSumOfPairsScore(inputData , outputData.progressiveAlignment);
        return [inputData, outputData];
    }

    /**
     * Computes scores (similarities),
     * the number of gaps, the alignment lengths
     * and so on between all sequences.
     */
    function computePairwiseData() {
        multiSequenceAlignmentInstance.computePairwiseData(gotohInstance);
    }

    /**
     * Converting the pairwise similarities into distances
     * by using the Feng-Doolittle formulas.
     */
    function computeDistancesFromSimilarities() {
        multiSequenceAlignmentInstance.computeDistancesFromSimilarities();
    }

    /**
     * Creates dependency between cluster names and distances.
     * So, a function dist(a,b) which giving you
     * for two cluster-names a and b the distance
     * (needed for the clustering algorithm).
     */
    function createDistanceMatrix() {
        multiSequenceAlignmentInstance.createDistanceMatrix();
    }

    /**
     * Using a clustering algorithm like UPGMA (Group Average)
     * the algorithm returns the binary guide tree branches in creation order.
     * @return {Object} - The tree branches.
     */
    function getPhylogeneticTree() {
        inputData.numOfStartClusters = inputData.sequences.length - inputData.arrayPositionsOfRemovedSequences.length;

        var clustering = new upgma.Upgma();
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

            var leftChildName = treeBranch.leftChild.name;
            var rightChildName = treeBranch.rightChild.name;
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
        outputData.secondGroups = [];
        outputData.guideAlignments = [];
        outputData.joinedGroups = [];
        outputData.joinedGroupNames = [];
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
        var group1Sequences = getGroupSequences(leftChildName);
        var group2Sequences = getGroupSequences(rightChildName);

        var bestAlignment = getBestAlignment(group1Sequences, group2Sequences);

        outputData.groups[groupName] = createGroup(group1Sequences, group2Sequences, bestAlignment);

        // for visualization of steps
        outputData.guideAlignments.push(bestAlignment);
        outputData.firstGroups.push(group1Sequences);
        outputData.secondGroups.push(group2Sequences);
        outputData.joinedGroups.push(outputData.groups[groupName]);
        outputData.joinedGroupNames.push(groupName);
    }

    /**
     * Returns the sequences of a group.
     * @param groupName {Array} - The name of the group.
     */
    function getGroupSequences(groupName) {
        return outputData.groups[groupName];
    }

    /**
     * Returns the best alignment (with the highest score) after
     * a pair-wise alignment of the group sequences of both groups
     * or after a look-up in the previously computed alignments for the start-sequences.
     * @param group1Sequences {Array} - The sequences of the first group.
     * @param group2Sequences {Array} - The sequences of the second group.
     */
    function getBestAlignment(group1Sequences, group2Sequences) {
        var maxScore = Number.NEGATIVE_INFINITY;
        var maxAlignment = [];

        // iterate through all sequences and search for search for alignment with maximum similarity
        for (var i = 0; i < group1Sequences.length; i++) {
            var sequence1 = group1Sequences[i];

            for (var j = 0; j < group2Sequences.length; j++) {
                var sequence2 = group2Sequences[j];
                var asData = getAlignmentAndScore(sequence1, sequence2);

                if (asData[1] > maxScore) {
                    maxScore = asData[1];
                    maxAlignment = asData[0];
                }
            }
        }

        return maxAlignment;
    }

    /**
     * Returns the initial previously computed alignment of the both sequences.
     * It avoids a time-intensive recomputation.
     * @param sequence1 {string} - The first sequence.
     * @param sequence2 {string} - The second sequence.
     * @return {string} - The initial alignment.
     */
    function getInitialAlignment(sequence1, sequence2) {
        return outputData.alignmentsAndScores[[sequence1, sequence2]];
    }

    /**
     * Returns the alignment and score.
     * @param sequence1 {string} - The first sequence.
     * @param sequence2 {string} - The second sequence.
     * @return {[string, number]} - The alignment and score.
     */
    function getAlignmentAndScore(sequence1, sequence2) {
        var alignmentAndScore = outputData.alignmentsAndScores[[sequence1, sequence2]];  // constant time!

        if (alignmentAndScore !== undefined)
            return alignmentAndScore;

        var input = {};
        multiSequenceAlignmentInstance.initializeInput(input);
        var ioData = multiSequenceAlignmentInstance.computeWithAlgorithm(gotohInstance, input, sequence1, sequence2);

        return [ioData[1].alignments[0], ioData[1].score];
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
        return multiSequenceAlignmentInstance;
    }
}());