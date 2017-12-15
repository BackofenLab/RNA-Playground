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
    var gotohInstance;

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
     * @see http://www.bioinf.uni-freiburg.de/Lehre/Courses/2017_SS/V_BioinfoII/ (source for algorithm)
     * and https://doi.org/10.1016/0022-2836(87)90316-0 (not implemented, but original idea)
     *
     * Barton, Geoffrey J., and Michael J.E. Sternberg.
     * "A strategy for the rapid multiple alignment of protein sequences: confidence levels from tertiary structure comparisons."
     * Journal of molecular biology 198.2 (1987): 327-337.
     */
    function IterativeRefinement() {
        iterativeRefinementInstance = this;

        // variables
        this.type = ALGORITHMS.ITERATIVE_REFINMENT;

        // instances (do not change order)
        multiSequenceAlignmentInstance = new bases.multiSequenceAlignment.MultiSequenceAlignment(this);
        gotohInstance = new gotoh.Gotoh();
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
        inputData.roundRobinOrder = inputViewmodel.selectedOrder()[0];
    }

    /**
     * Starts the computation.
     */
    function compute() {
        initializeStructs();

        var ioData = computeFengDoolittle();
        var startMsa = ioData[1].progressiveAlignment.slice();
        var startMsaName = ioData[1].joinedGroupNames[ioData[1].joinedGroupNames.length - 1];
        var startMsaSequenceNames = multiSequenceAlignmentInstance.getIndividualSequenceNames(startMsaName, false);

        var names = getNamesInOrder(ioData[1].treeBranches[ioData[1].treeBranches.length - 1], ioData[1].clusterNames);

        if (names.length > 1) {
            for (var i = 0; i < names.length; i++) {
                var removedSequenceName = names[i];

                // create a copy of the alignment before a realignment
                var msa = startMsa.slice();  // shallow copy, because the elements itself are not modified
                var msaSequenceNames = startMsaSequenceNames.slice();

                // [1] remove sequence from MSA
                var mrData = getMsaAndRemovedSequence(removedSequenceName, msa, msaSequenceNames); // ([MSA], [removed sequence])

                // [2] realign the removed sequence
                var msaRefinedWithName = getRealignment(mrData[1][0], mrData[1][1], mrData[0], removedSequenceName, ioData[1].distanceMatrix, ioData[1].nameOfSequence);

                // [3] compute score of the MSA and refined MSA (replace startMsa with refinedMsa if [refinedMsa score] > [startMsa score])
                var msaWithName = getBetterMultiSequenceAlignment([startMsa, startMsaName], msaRefinedWithName);
                startMsa = msaWithName[0];
                startMsaName = msaWithName[1];
                startMsaSequenceNames = multiSequenceAlignmentInstance.getIndividualSequenceNames(startMsaName, false);
            }
        }

        // storing data from Feng-Doolittle and refined alignment
        storeAlignmentData(ioData, startMsa, multiSequenceAlignmentInstance.getAffineSumOfPairsScore(inputData, startMsa), startMsaName);

        return [inputData, outputData];
    }

    /**
     * Initializes structs used in the algorithm.
     * Structs store information only if realignment is accepted (score increased).
     */
    function initializeStructs() {
        outputData.firstGroups = [];  // removed sequences
        outputData.firstGroupsNames = [];

        outputData.secondGroups = [];  // remaining alignments
        outputData.secondGroupsNames = [];

        outputData.guideAlignments = [];   // guide alignments
        outputData.guideAlignmentsNames = [];

        outputData.joinedGroups = [];  // realignments
        outputData.joinedGroupNames = [];
        outputData.realignmentsScores = [];

        outputData.accepted = [];
    }

    /**
     * Computes Feng-Doolittle output.
     */
    function computeFengDoolittle() {
        fengDoolittleInstance.setIO(inputData, {});
        return fengDoolittleInstance.compute();
    }

    /**
     * Returns the names of sequences in order added to the MSA by doing an post-order-traversal.
     * @param tree {Object} - The phylogenetic tree from which you want get the representation.
     * @param sequenceNames {Array} - The names of the sequences.
     * @return {Array} - The names of the sequences in order they have been added to the MSA.
     */
    function getNamesInOrder(tree, sequenceNames) {
        var names = [];

        if (inputData.roundRobinOrder === ITERATIVE_REFINEMENT_ORDERS.INPUT_ORDER)
            names = sequenceNames;
        if (inputData.roundRobinOrder === ITERATIVE_REFINEMENT_ORDERS.LEFT_FIRST)
            postOrder(tree, names, false);
        else if (inputData.roundRobinOrder === ITERATIVE_REFINEMENT_ORDERS.RIGHT_FIRST)
            postOrder(tree, names, true);

        return names;  // needed because else you get reverse order
    }

    /**
     * Does a post-order traversal to get the sequence names.
     * @param node {Object} - The node from which on traversal takes place. At the beginning it is the root node.
     * @param names {Array} - The array which is filled with names.
     * @param rightFirst {boolean} - Tells if we first look on the right node or not.
     */
    function postOrder(node, names, rightFirst) {
        if (node === undefined)
            return;

        if (rightFirst) {
            postOrder(node.leftChild, names, rightFirst);
            postOrder(node.rightChild, names, rightFirst);
        } else {
            postOrder(node.rightChild, names, rightFirst);
            postOrder(node.leftChild, names, rightFirst);
        }

        var isLeaf = node.rightChild === undefined && node.leftChild === undefined;

        if (isLeaf) names.push(node.name);
    }

    /**
     * Removes a sequence by its name from the MSA and then
     * @param removedSequenceName {string} - The name of the sequence which should be removed.
     * @param progressiveAlignment {Array} - The progressive alignment from which is removed.
     * @param msaSequenceNames {Array} - The names of sequences in the multi-sequence-alignment.
     * @return {[removedSequence, remainingAlignment]}
     */
    function getMsaAndRemovedSequence(removedSequenceName, progressiveAlignment, msaSequenceNames) {
        var removedSequence = getSequenceByName(progressiveAlignment, msaSequenceNames, removedSequenceName);
        var remainingAlignment = removeSequenceFromMSA(progressiveAlignment, msaSequenceNames, removedSequenceName);

        outputData.firstGroups.push([removedSequence.replace(MULTI_SYMBOLS.GAP, SYMBOLS.EMPTY).replace(MULTI_SYMBOLS.NONE, SYMBOLS.EMPTY)]);
        outputData.firstGroupsNames.push(removedSequenceName);

        outputData.secondGroups.push(remainingAlignment[0]);
        outputData.secondGroupsNames.push((remainingAlignment[1].toString()).replace(MULTI_SYMBOLS.COMMA, SYMBOLS.EMPTY));  // push the name

        return [removedSequence, remainingAlignment];
    }

    /**
     * Returns the sequence from given multi-sequence-alignment.
     * @param multiSequenceAlignment {Array} - The multi-sequence-alignment from which a sequence should be returned.
     * @param msaSequenceNames {Array} - The names of sequences in the multi-sequence-alignment.
     * @param sequenceName {string} - The name of the sequence which should be returned from the multi-sequence-alignment.
     * @return {string} - Returns the sequence.
     */
    function getSequenceByName(multiSequenceAlignment, msaSequenceNames, sequenceName) {
        var indexInSequenceNames = msaSequenceNames.indexOf(sequenceName);
        return multiSequenceAlignment[indexInSequenceNames];
    }

    /**
     * Removes a sequence from the MSA by sequence-name.
     * @param multiSequenceAlignment {Array} - The multi-sequence-alignment from which a sequence should be removed.
     * @param msaSequenceNames {Array} - The names of sequences in the multiSequenceAlignment.
     * @param sequenceName {string} - The name of the sequence which should be removed from the multi-sequence-alignment.
     * @return [multiSequenceAlignment, msaSequencesNames] - The multi-sequence-alignment in which the sequence with the given name is removed.
     */
    function removeSequenceFromMSA(multiSequenceAlignment, msaSequenceNames, sequenceName) {
        var indexInSequenceNames = msaSequenceNames.indexOf(sequenceName);

        // remove sequenceName from msaSequenceNames and the associated sequence from the multiSequenceAlignment
        msaSequenceNames.splice(indexInSequenceNames, 1);
        multiSequenceAlignment.splice(indexInSequenceNames, 1);
        multiSequenceAlignment = removeGapOnlyColumns(multiSequenceAlignment);

        return [multiSequenceAlignment, msaSequenceNames];
    }

    /**
     * Removes gap-only columns.
     * Hint: Is a part of the original algorithm by Geoffrey J.Barton and Michael J. E. Sternberg.
     * @param msa {Array} - The msa from which gap only columns should be removed.
     * @return {Array} - The msa without gap-only columns.
     */
    function removeGapOnlyColumns(msa) {
        var storeColumn = [];

        // find out which gap columns have to be stored
        for (var j = 0; j < msa[0].length; j++) {  // iteration over each column
            var isGapOnly = true;

            for (var i = 0; i < msa.length; i++) { // iteration over each sequence
                if (msa[i][j] !== SYMBOLS.GAP && msa[i][j] !== SYMBOLS.NONE) {  // if not only gaps in column
                    isGapOnly = false;
                    break;
                }
            }

            storeColumn.push(!isGapOnly);
        }

        // create new msa without gap-only columns
        var decentMsa = [];  // without gap-only columns
        for (var i = 0; i < msa.length; i++)
            decentMsa[i] = SYMBOLS.EMPTY;

        for (var j = 0; j < storeColumn.length; j++) {

            if (storeColumn[j]) {
                for (var i = 0; i < msa.length; i++) { // iteration over each sequence
                    decentMsa[i] += msa[i][j];
                }
            }
        }

        return decentMsa;
    }

    /**
     * Returns the realignment of the removed sequence with the input MSA.
     * @param msa {Array} - The multi-sequence alignment to which the removed sequence should be realigned.
     * @param msaSequenceNames {Array} - The names of sequences in the multi-sequence-alignment.
     * @param removedSequence {string} - The removed sequence which should be realigned.
     * @param removedSequenceName {string} - The name of the sequence which is removed.
     * @param distanceMatrix {Object} - The distances which are to get the element with lowest distance to the removed sequence.
     * @param sequenceIdentificator {Object} - An object which identifies a sequence and returns its name.
     * @return {Array} - The MSA in which the removed sequence is realigned and possibly the MSA-score contained.
     */
    function getRealignment(msa, msaSequenceNames, removedSequence, removedSequenceName, distanceMatrix, sequenceIdentificator) {
        var realignment = [];

        if (inputData.iterativeRefinementSubalgorithm === ITERATIVE_REFINEMENT_STRATEGIES.MIN_DISTANCE_PAIR)
            realignment = getMinDistanceRealignment(msa, msaSequenceNames, removedSequence, removedSequenceName, distanceMatrix);
        else if (inputData.iterativeRefinementSubalgorithm === ITERATIVE_REFINEMENT_STRATEGIES.ONE_VS_ALL)
            realignment = getOneVsAllRealignment(msa, msaSequenceNames, removedSequence, removedSequenceName);
        else if (inputData.iterativeRefinementSubalgorithm === ITERATIVE_REFINEMENT_STRATEGIES.PAIRWISE_BEST_PAIR)
            realignment = getPairwiseBestRealignment(msa, msaSequenceNames, removedSequence, removedSequenceName, sequenceIdentificator);

        return realignment;
    }

    /**
     * Returns the realignment of the removed sequence
     * with the input MSA using the minimum distance approach (described in the algorithm HTML-file).
     * @param msa {Array} - The multi-sequence alignment to which the removed sequence should be realigned.
     * @param msaSequenceNames {Array} - The names of sequences in the multi-sequence-alignment.
     * @param removedSequence {string} - The removed sequence which should be realigned.
     * @param removedSequenceName {string} - The name of the sequence which is removed.
     * @param distanceMatrix {Object} - The distances which are to get the element with lowest distance to the removed sequence.
     * @return {Array} - The MSA in which the removed sequence is realigned.
     */
    function getMinDistanceRealignment(msa, msaSequenceNames, removedSequence, removedSequenceName, distanceMatrix) {
        // get nearest sequence
        var nearestElementName = getNearestElement(removedSequenceName, distanceMatrix);
        var nearestSequence = getSequenceByName(msa, msaSequenceNames, nearestElementName);

        // realign removed sequence with nearest sequence
        var cleanSequence = removedSequence.replace(MULTI_SYMBOLS.GAP, SYMBOLS.EMPTY).replace(MULTI_SYMBOLS.NONE, SYMBOLS.EMPTY);
        var realignment = getAffineRealignment(cleanSequence, multiSequenceAlignmentInstance.replaceGapsWithPlaceHolder([nearestSequence])[0]);

        var msaName = removedSequenceName + (msaSequenceNames.toString()).replace(MULTI_SYMBOLS.COMMA, SYMBOLS.EMPTY);

        // only for visualization
        outputData.guideAlignments.push(realignment);
        outputData.guideAlignmentsNames.push(removedSequenceName + SYMBOLS.ALIGN + nearestElementName);

        // add realignment to MSA and possibly fill out with new gaps
        return [multiSequenceAlignmentInstance.createGroup([cleanSequence], msa, realignment), msaName];
    }

    /**
     * Given a name of an element, it is searched for nearest partner of this element.
     * @param name {string} - The name for which the nearest element should be found.
     * @param distanceMatrix {Object} - The distances which are to get the element with lowest distance to the removed sequence.
     * @return {string} - The nearest element from distance matrix.
     */
    function getNearestElement(name, distanceMatrix) {
        var keys = Object.keys(distanceMatrix);
        var names = getNamesFromKeys(keys);

        var minDistance = Number.POSITIVE_INFINITY;
        var minElement = SYMBOLS.EMPTY;

        for (var i = 0; i < names.length; i++) {
            if (names[i] !== name) {  // no diagonals
                var distance = distanceMatrix[[name, names[i]]];

                if (distance === undefined)
                    distance = distanceMatrix[[names[i], name]];

                if (distance < minDistance) {
                    minDistance = distance;
                    minElement = names[i];
                }
            }
        }

        return minElement;
    }

    /**
     * Returns the individual arguments of all given keys.
     * @example: A key has this form: "[a,b]" and the arguments are "a" and "b" are stored and later on returned with all other arguments.
     * @param keys {Object} - The keys from which the individual arguments are returned.
     * @return {Array} - The individual arguments of the keys.
     */
    function getNamesFromKeys(keys) {
        var names = [];

        for (var i = 0; i < keys.length; i++) {
            var key = keys[i].split(SYMBOLS.COMMA);

            if (names.indexOf(key[0]) === -1)
                names.push(key[0]);

            if (names.indexOf(key[1]) === -1)
                names.push(key[1]);
        }

        return names;
    }

    /**
     * Returns an affine realignment.
     * @param sequence1 {string} - The first string which should be aligned.
     * @param sequence2 {string} - The second string which should be aligned.
     * @return {[alignedSequenceB, matchOrMismatchString, alignedSequenceA]} - The alignment.
     */
    function getAffineRealignment(sequence1, sequence2) {
        inputData.sequenceA = sequence1;
        inputData.sequenceB = sequence2;

        inputData.computeOneAlignment = true;  // extension to speed up computation
        inputData.recomputeTraceback = true;  // an alignment is not stored, we have to recompute

        inputData.matrixHeight = inputData.sequenceA.length + 1;
        inputData.matrixWidth = inputData.sequenceB.length + 1;

        gotohInstance.setIO(inputData, {});
        return gotohInstance.compute()[1].alignments[0];
    }

    /**
     * Returns the realignment of the removed sequence
     * with the input MSA using the minimum distance approach (described in the algorithm HTML-file).
     * @param msa {Array} - The multi-sequence alignment to which the removed sequence should be realigned.
     * @param msaSequenceNames {Array} - The names of sequences in the multi-sequence-alignment.
     * @param removedSequence {string} - The removed sequence which should be realigned.
     * @param removedSequenceName {string} - The name of the sequence which is removed.
     * @return {Array} - The MSA in which the removed sequence is realigned and the score contained.
     */
    function getOneVsAllRealignment(msa, msaSequenceNames, removedSequence, removedSequenceName) {
        multiSequenceAlignmentInstance.setIO(inputData, outputData);

        var bestScore = Number.NEGATIVE_INFINITY;
        var bestPairwiseAlignment = SYMBOLS.EMPTY;
        var bestPairwisePartnerName = SYMBOLS.EMPTY;
        var bestMsa = [];

        // realign removed sequence once with every sequence
        var cleanSequence = removedSequence.replace(MULTI_SYMBOLS.GAP, SYMBOLS.EMPTY).replace(MULTI_SYMBOLS.NONE, SYMBOLS.GAP);

        for (var i = 0; i < msaSequenceNames.length; i++) {  // remaining sequence names without the sequence name of the removed sequence
            var sequence = getSequenceByName(msa, msaSequenceNames, msaSequenceNames[i]);
            var pairwiseRealignment = getAffineRealignment(cleanSequence, multiSequenceAlignmentInstance.replaceGapsWithPlaceHolder([sequence])[0]);
            var refinedMsa = multiSequenceAlignmentInstance.createGroup([cleanSequence], msa, pairwiseRealignment);
            var refinedMsaScore
                = multiSequenceAlignmentInstance.getAffineSumOfPairsScore(inputData, multiSequenceAlignmentInstance.replacePlaceholdersWithGaps(refinedMsa));

            if (refinedMsaScore > bestScore) {
                bestScore = refinedMsaScore;
                bestPairwiseAlignment = pairwiseRealignment;
                bestPairwisePartnerName = msaSequenceNames[i];
                bestMsa = refinedMsa;
            }
        }

        var msaName = removedSequenceName + (msaSequenceNames.toString()).replace(MULTI_SYMBOLS.COMMA, SYMBOLS.EMPTY);

        // only for visualization
        outputData.guideAlignments.push(bestPairwiseAlignment);
        outputData.guideAlignmentsNames.push(removedSequenceName + SYMBOLS.ALIGN + bestPairwisePartnerName);

        return [bestMsa, msaName, bestScore];
    }

    /**
     * Returns the realignment of the removed sequence
     * with the input MSA using the minimum distance approach (described in the algorithm HTML-file).
     * @param msa {Array} - The multi-sequence alignment to which the removed sequence should be realigned.
     * @param msaSequenceNames {Array} - The names of sequences in the multi-sequence-alignment.
     * @param removedSequence {string} - The removed sequence which should be realigned.
     * @param removedSequenceName {string} - The name of the sequence which is removed.
     * @param sequenceIdentificator {Object} - An object which identifies a sequence and returns its name.
     * @return {Array} - The MSA in which the removed sequence is realigned.
     */
    function getPairwiseBestRealignment(msa, msaSequenceNames, removedSequence, removedSequenceName, sequenceIdentificator) {
        multiSequenceAlignmentInstance.setIO(inputData, outputData);

        // realign removed sequence with best sequence
        var cleanSequence = removedSequence.replace(MULTI_SYMBOLS.GAP, SYMBOLS.EMPTY).replace(MULTI_SYMBOLS.NONE, SYMBOLS.GAP);
        var bestAlignment
            = multiSequenceAlignmentInstance.getBestAlignment([cleanSequence], multiSequenceAlignmentInstance.replaceGapsWithPlaceHolder(msa));

        var bestElementName
            = sequenceIdentificator[bestAlignment[2].replace(MULTI_SYMBOLS.NONE, SYMBOLS.EMPTY).replace(MULTI_SYMBOLS.GAP, SYMBOLS.EMPTY)];

        var msaName = removedSequenceName + (msaSequenceNames.toString()).replace(MULTI_SYMBOLS.COMMA, SYMBOLS.EMPTY);

        // only for visualization
        outputData.guideAlignments.push(bestAlignment);
        outputData.guideAlignmentsNames.push(removedSequenceName + SYMBOLS.ALIGN + bestElementName);

        return [multiSequenceAlignmentInstance.createGroup([cleanSequence], msa, bestAlignment), msaName];
    }

    /**
     * Computes the Sum-Of-Pairs score of both input MSA and returns the one with the higher score.
     * @param msaWithName {Array} - The original multi-sequence alignment with its name.
     * @param refinedMsaWithName {Array} - The refined multi-sequence alignment with its name and possibly its score.
     * @return {Array} - The better MSA.
     */
    function getBetterMultiSequenceAlignment(msaWithName, refinedMsaWithName) {
        var msa = multiSequenceAlignmentInstance.replacePlaceholdersWithGaps(msaWithName[0]);
        var refinedMsa = multiSequenceAlignmentInstance.replacePlaceholdersWithGaps(refinedMsaWithName[0]);

        // compute affine scores of both MSA
        var msaScore = multiSequenceAlignmentInstance.getAffineSumOfPairsScore(inputData, msa);
        var refinedMsaScore = refinedMsaWithName.length !== 3 ?  // score not contained
            multiSequenceAlignmentInstance.getAffineSumOfPairsScore(inputData, refinedMsa) : refinedMsaWithName[2];

        // only for visualization
        outputData.joinedGroups.push(refinedMsa);
        outputData.joinedGroupNames.push(refinedMsaWithName[1]);
        outputData.realignmentsScores.push(refinedMsaScore);

        // check which is higher and return the better MSA
        if (msaScore >= refinedMsaScore) {
            outputData.accepted.push(false);
            return msaWithName;
        }

        // only for visualization
        outputData.accepted.push(true);

        return refinedMsaWithName;
    }

    /**
     * Stores the alignment data from input.
     * @param fengDoolittleData {[input, output]} - The input and output data of Feng-Doolittle.
     * @param refinedMsa {Array} - The refined multi-sequence alignment.
     * @param refinedScore {number} - The refined multi-sequence alignment score.
     * @param refinedMsaName {string} - The refined multi-sequence alignment name.
     */
    function storeAlignmentData(fengDoolittleData, refinedMsa, refinedScore, refinedMsaName) {
        outputData.progressiveAlignment = fengDoolittleData[1].progressiveAlignment;
        outputData.score = fengDoolittleData[1].score;

        outputData.newickString = fengDoolittleData[1].newickString;
        outputData.progressiveAlignmentName = fengDoolittleData[1].joinedGroupNames[fengDoolittleData[1].joinedGroupNames.length - 1];

        outputData.distanceMatrix = fengDoolittleData[1].distanceMatrix;
        outputData.distanceMatrixLength = fengDoolittleData[1].distanceMatrixLength;
        outputData.remainingClusters = fengDoolittleData[1].remainingClusters;

        outputData.refinedProgressiveAlignment = refinedMsa;
        outputData.refinedScore = refinedScore;
        outputData.refinedProgressiveAlignmentName = refinedMsaName;
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