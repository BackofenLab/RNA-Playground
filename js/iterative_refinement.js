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
    if (loaded === ALGORITHMS.ITERATIVE_REFINMENT) {  // to avoid self execution on a script import
        iterativeRefinement.startIterativeRefinement();
        loaded = ALGORITHMS.NONE;
    }
});

(function () {  // namespace
    // public methods
    namespace("iterativeRefinement", startIterativeRefinement, IterativeRefinement);

    // instances
    var multiSequenceAlignmentInstance;

    var iterativeRefinementInstance;
    var fengDoolittleInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Function managing objects.
     */
    function startIterativeRefinement() {
        var multiSequenceInterface = new interfaces.multiSequenceInterface.MultiSequenceInterface();
        multiSequenceInterface.startMultiSequenceInterface(IterativeRefinement, ALGORITHMS.ITERATIVE_REFINMENT);
    }

    /*---- ALGORITHM ----*/
    /**
     * Computes an affine multi-sequence-alignment (non-optimal approach) with Feng-Doolittle
     * and afterwards post-processes the multi-sequence-alignment with an iterative refinement.
     * @constructor
     * @augments MultiSequenceAlignment
     * @see http://www.bioinf.uni-freiburg.de/Lehre/Courses/2017_SS/V_BioinfoII/ (just used to get some ideas)
     */
    function IterativeRefinement() {
        iterativeRefinementInstance = this;

        // variables
        this.type = ALGORITHMS.ITERATIVE_REFINMENT;

        // instances (do not change order)
        multiSequenceAlignmentInstance = new bases.multiSequenceAlignment.MultiSequenceAlignment(this);
        fengDoolittleInstance = new fengDoolittle.FengDoolittle();

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
     */
    function setInput(inputViewmodel) {
        multiSequenceAlignmentInstance.setIO(inputData, {});
        multiSequenceAlignmentInstance.setInput(inputViewmodel);

        inputData.iterativeRefinementSubalgorithm = inputViewmodel.selectedApproach()[0];  // because it's array from which we choose
    }

    /**
     * Starts the computation.
     */
    function compute() {
        initializeStructs();

        var ioData = computeFengDoolittle();
        var startMsa = ioData.progressiveAlignment.slice();
        var names = getNamesInOrderAddedToMSA(ioData.treeBranches[ioData.treeBranches.length - 1]);

        for (var i = 0; i < names.length; i++) {
            // [0] create copy of the alignment before a realignment
            var msa = startMsa;

            // [1] remove sequence from MSA
            var mrData  // MSA and removed sequence
                = getMsaAndRemovedSequence(names, i, msa, ioData.joinedGroupNames[ioData.joinedGroupNames.length - 1]);

            // [2] realign the removed sequence
            var msaRefined = getRealignment(mrData[0], mrData[1]);  // ([MSA], [removed sequence])

            // [3] compute score of the MSA and refined MSA (replace startMsa with refinedMsa if [refinedMsa score] > [msa score])
            startMsa = getBetterMultiSequenceAlignment(msa, msaRefined);
        }

        storeAlignmentData(ioData[1].progressiveAlignment, startMsa, ioData[1].score,
            multiSequenceAlignmentInstance.getAffineSumOfPairsScore(inputData , outputData.refinedProgressiveAlignment));
        return [inputData, outputData];
    }

    /**
     * Initializes structs used in the algorithm.
     * Structs store information only if realignment is accepted (score increased).
     */
    function initializeStructs() {
        outputData.removedSequences = [];
        outputData.removedSequencesNames = [];

        outputData.remainingAlignment = [];
        outputData.remainingAlignmentName = [];

        outputData.guideAlignments = [];
        outputData.guideAlignmentsNames = [];

        outputData.acceptedRealignments = [];
        outputData.acceptedRealignmentsNames = [];
        outputData.acceptedRealignmentsScores = [];
    }

    /**
     * Computes Feng-Doolittle output.
     */
    function computeFengDoolittle(input) {
        fengDoolittleInstance.setIO(inputData, {});
        return fengDoolittleInstance.compute();
    }

    /**
     * Returns the names of sequences in order added to the MSA by doing an post-order-traversal.
     * @param tree - The phylogenetic tree from which you want get the representation.
     * @return {Array} - The names of the sequences in order they have been added to the MSA.
     */
    function getNamesInOrderAddedToMSA(tree) {
        var names = [];
        postOrder(tree, names);
        return names;
    }

    /**
     * Does a post-order traversal to get the sequence names.
     * @param node {Object} - The node from which on traversal takes place. At the beginning it is the root node.
     * @param names {Array} - The array which is filled with names.
     */
    function postOrder(node, names) {
        if (node === undefined)
            return;

        postOrder(node.leftChild, names);
        postOrder(node.rightChild, names);

        var isLeaf = node.leftChild === undefined && node.rightChild === undefined;

        if (isLeaf) names.push(node.name);
    }

    /**
     * Removes a sequence by its name from the MSA and then
     * @param names {Array} - The names of the sequences.
     * @param namePos {number} - The position which within the names array which is used to get a sequence to remove.
     * @param progressiveAlignment {Array} - The progressive alignment from which is removed.
     * @param lastGroupName {string} - The name of the progressive alignment.
     * @return {[removedSequence, remainingAlignment]}
     */
    function getMsaAndRemovedSequence(names, namePos, progressiveAlignment, lastGroupName) {
        var removedSequenceName = names[namePos];

        var removedSequence = getRemovedSequence(progressiveAlignment, lastGroupName, removedSequenceName);
        var remainingAlignment = removeSequenceFromMSA(progressiveAlignment, lastGroupName, removedSequenceName);

        outputData.removedSequences.push(removedSequence);
        outputData.removedSequencesNames.push(removedSequenceName);

        outputData.remainingAlignment.push(remainingAlignment[0]);
        outputData.remainingAlignment.push(remainingAlignment[1]);  // push the name

        return [removedSequence, remainingAlignment];
    }

    /**
     * Returns the sequence removed from given multi-sequence-alignment.
     * @param multiSequenceAlignment - The multi-sequence-alignment from which a sequence should be removed.
     * @param msaSequenceNames - The names of sequences in the multiSequenceAlignment.
     * @param sequenceName - The name of sequence which should be removed from the multi-sequence-alignment.
     * @return {string} - Returns the removed sequence.
     */
    function getRemovedSequence(multiSequenceAlignment, msaSequenceNames, sequenceName) {
        return undefined;
    }

    /**
     * Removes a sequence from the MSA by sequence-name.
     * @param multiSequenceAlignment - The multi-sequence-alignment from which a sequence should be removed.
     * @param msaSequenceNames - The names of sequences in the multiSequenceAlignment.
     * @param sequenceName - The name of sequence which should be removed from the multi-sequence-alignment.
     * @return [multiSequenceAlignment, msaSequencesNames] - The multi-sequence-alignment in which the sequence with the given name is removed.
     */
    function removeSequenceFromMSA(multiSequenceAlignment, msaSequenceNames, sequenceName) {
        // retrieve separate names from msaSequenceNames
        // remove sequenceName from msaSequenceNames and the associated sequence from the multiSequenceAlignment
        return undefined;
    }

    /**
     * Returns the realignment of the removed sequence with the input MSA.
     * @param msa {Array} - The multi-sequence alignment to which the removed sequence should be realigned.
     * @param removedSequence - The removed sequence which should be realigned.
     * @return {Array} - The MSA in which the removed sequence is realigned.
     */
    function getRealignment(msa, removedSequence) {
        var realignment = [];

        if (inputData.iterativeRefinementSubalgorithm === ITERATIVE_REFINEMENT_STRATEGIES.MIN_DISTANCE_PAIR)
            realignment = getMinDistanceRealignment(msa, removedSequence);
        else if (inputData.iterativeRefinementSubalgorithm === ITERATIVE_REFINEMENT_STRATEGIES.ONE_VS_ALL)
            realignment = getOneVsAllRealignment(msa, removedSequence);

        return realignment;
    }

    /**
     * Returns the realignment of the removed sequence
     * with the input MSA using the minimum distance approach (described in the algorithm HTML-file).
     * @param msa {Array} - The multi-sequence alignment to which the removed sequence should be realigned.
     * @param removedSequence - The removed sequence which should be realigned.
     * @return {Array} - The MSA in which the removed sequence is realigned.
     */
    function getMinDistanceRealignment(msa, removedSequence) {
    }

    /**
     * Returns the realignment of the removed sequence
     * with the input MSA using the minimum distance approach (described in the algorithm HTML-file).
     * @param msa {Array} - The multi-sequence alignment to which the removed sequence should be realigned.
     * @param removedSequence - The removed sequence which should be realigned.
     * @return {Array} - The MSA in which the removed sequence is realigned.
     */
    function getOneVsAllRealignment(msa, removedSequence) {
    }

    /**
     * Computes the Sum-Of-Pairs score of both input MSA and returns the one with the higher score.
     * @param msa {Array} - The original multi-sequence alignment.
     * @param refinedMsa {Array} - The refined multi-sequence alignment.
     * @return {Array} - The better MSA.
     */
    function getBetterMultiSequenceAlignment(msa, refinedMsa) {
        // compute affine scores of both MSA

        // check which is higher and return the better MSA
        return undefined;
    }

    /**
     * Stores the alignment data from input.
     * @param msa {Array} - The original multi-sequence alignment.
     * @param refinedMsa {Array} - The refined multi-sequence alignment.
     * @param score {number} - The original multi-sequence alignment score.
     * @param refinedScore {number} - The refined multi-sequence alignment score.
     */
    function storeAlignmentData(msa, refinedMsa, score, refinedScore) {
        outputData.progressiveAlignment = progressiveAlignment;
        outputData.refinedProgressiveAlignment = refinedAlignment;

        outputData.score = score;
        outputData.refinedScore = refinedScore;
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