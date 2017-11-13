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
    if (loaded === ALGORITHMS.NOTREDAME_HIGGINS_HERINGA) {  // to avoid self execution on a script import
        notredameHigginsHeringa.startNotredameHigginsHeringa();
        loaded = ALGORITHMS.NONE;
    }
});

(function () {  // namespace
    // public methods
    namespace("notredameHigginsHeringa", startNotredameHigginsHeringa, NotredameHigginsHeringa,
        getInput, setInput, compute, getOutput, setIO, getSuperclass);

    // instances
    var multiSequenceAlignmentInstance;

    var gotohInstance;
    var gotohLocalInstance;
    var notredameHigginsHeringaInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Function managing objects.
     */
    function startNotredameHigginsHeringa() {
        var multiSequenceInterface = new interfaces.multiSequenceInterface.MultiSequenceInterface();
        multiSequenceInterface.startMultiSequenceInterface(NotredameHigginsHeringa, ALGORITHMS.NOTREDAME_HIGGINS_HERINGA);
    }

    /*---- ALGORITHM ----*/
    /**
     * Computes affine multi-alignments (non-optimal approach) with T-Coffee.
     * @constructor
     * @augments MultiSequenceAlignment
     */
    function NotredameHigginsHeringa() {
        notredameHigginsHeringaInstance = this;

        // variables
        this.type = ALGORITHMS.NOTREDAME_HIGGINS_HERINGA;

        // instances (do not change order)
        multiSequenceAlignmentInstance = new bases.multiSequenceAlignment.MultiSequenceAlignment(this);

        gotohInstance = new gotoh.Gotoh();
        gotohLocalInstance = new gotohLocal.GotohLocal();

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
        reinitializeInputOutput();
        multiSequenceAlignmentInstance.setIO(inputData, {});
        multiSequenceAlignmentInstance.setInput(inputViewmodel);

        inputData.useLocalLibrary = inputViewmodel.useLocalLibrary();

        inputData.baseCostsLocal = inputViewmodel.baseCostsLocal();
        inputData.enlargementLocal = inputViewmodel.enlargementLocal();
        inputData.matchLocal = inputViewmodel.matchLocal();
        inputData.mismatchLocal = inputViewmodel.mismatchLocal();

        //inputData.maxNumberOptimalAlignments = inputViewmodel.maxNumberOptimalAlignments();
        //inputData.maxNumberOptimalAlignmentsLocal = inputViewmodel.maxNumberOptimalAlignmentsLocal();
    }

    /**
     * Reinitializes the input and the output before a recomputation with the algorithm.
     * It is needed, because else previously computed data can disturb newly computed data.
     */
    function reinitializeInputOutput() {
        inputData = {};
        outputData = {};
    }

    /**
     * Starts the computation.
     */
    function compute() {
        debugger;
        computePrimaryLibraries();
        computeCombinedWeightPrimaryLibrary();
        computeExtendedWeightPrimaryLibrary();
        startProgressiveAlignment();

        outputData.score = formats.scoringFunctions.getAffineSumOfPairsScore(inputData , outputData.progressiveAlignment);
        return [inputData, outputData];
    }

    /**
     * Computes the set of pairwise alignments for local and global alignments.
     */
    function computePrimaryLibraries() {
        computePairwiseGlobalAlignmentData();
        if (inputData.useLocalLibrary) computePairwiseLocalAlignmentData();
    }

    /**
     * Computes scores (similarities),
     * the number of gaps, the alignment lengths
     * and so on between all sequences.
     */
    function computePairwiseGlobalAlignmentData() {
        multiSequenceAlignmentInstance.setIO(inputData, outputData);
        multiSequenceAlignmentInstance.computePairwiseData(gotohInstance);
    }

    function computePairwiseLocalAlignmentData() {
        multiSequenceAlignmentInstance.setIO(inputData, outputData);
        multiSequenceAlignmentInstance.computePairwiseData(gotohLocalInstance);
    }

    /**
     * Computes the weights for both libraries
     * and combines both libraries to one big library (signal addition).
     */
    function computeCombinedWeightPrimaryLibrary() {
        // Hint: conversion of double sequences into alignment edges is directly done during the computations of pairwise weights
        outputData.primaryGlobalWeightLib = computePairwiseWeights(outputData, true);

        if (inputData.useLocalLibrary) {
            outputData.primaryLocalWeightLib = computePairwiseWeights(outputData, false);
            addSignals();
        } else  // only global library
            outputData.primaryWeightLib = outputData.primaryGlobalWeightLib;
    }

    /**
     * Computes the sequence identity.
     * So, how much is identical between two sequences
     * with respect to the smaller sequence.
     * @param output {Object} - The output which is used to compute a primary library.
     * @param global {boolean} - Tells if local or global data should be used for computation.
     */
    function computePairwiseWeights(output, global) {
        var primaryWeightLib = {};
        outputData.numberOfPairs = 0;  // for visualization

        // iterate over each sequence a and sequence b to compute structure primLib^{a,b}(i,j) = {L_{1,3}, L_{2,4}, ..., L_{5,7}}
        // Hint: a and b are the aligned sequences
        for (var j = 1; j < inputData.sequences.length; j++) {
            if (inputData.arrayPositionsOfRemovedSequences.indexOf(j) === -1) {  // only if the sequence is not a duplicate
                for (var i = 0; i < j; i++) {
                    var sequenceA = inputData.sequences[i];
                    var sequenceB = inputData.sequences[j];

                    var alignment = [];
                    var traceback = [];
                    var sequenceIdentities = {};

                    if (global) {
                        alignment = output.alignmentsAndScores[[sequenceA, sequenceB]][0];
                        sequenceIdentities = getSequenceIdentities(alignment, undefined);  // alignment = [alignedSequenceA, matchMismatchString, alignedSequenceB]
                    }
                    else {
                        debugger;
                        alignment = output.alignmentsAndScoresLocal[[sequenceA, sequenceB]][0];
                        traceback = output.tracebacks[[sequenceA, sequenceB]];
                        sequenceIdentities = getSequenceIdentities(alignment, traceback);  // alignment = [alignedSequenceA, matchMismatchString, alignedSequenceB]
                    }
                    primaryWeightLib[[sequenceA, sequenceB]] = sequenceIdentities;
                    outputData.numberOfPairs++;
                }
            }
        }

        return primaryWeightLib;
    }

    /**
     * Computes a structure L with elements
     * containing the sequence identity.
     * @param alignment - The alignment for which you want compute a structure that returns position sequence identities.
     * @param traceback {Array} - The traceback of the alignment which tells the defined positions.
     * @return {Object}
     */
    function getSequenceIdentities(alignment, traceback) {
        var global = traceback === undefined;  // global alignments do not need the traceback to get defined positions

        var sequenceA = alignment[0];
        var sequenceB = alignment[2];

        var sequenceLength = sequenceA.length;  // OR: sequenceB.length
        var numCharactersInA = 0;
        var numCharactersInB = 0;
        var numMatches = 0;
        var numMatchesOrMismatches = 0;
        var sequenceIdentity = 0;

        var L = {};

        // iterate over each position to create keys L_{i,j}
        for (var k = 0; k < sequenceLength; k++) {
            var currentPosition = global ? undefined : traceback[k+1];  // "+1" to jump over local end-position (with value 0)

            if (sequenceA[k] === SYMBOLS.GAP) {  // means there is no gap in sequence b
                numCharactersInB++;
            } else if (sequenceB[k] === SYMBOLS.GAP) {  // means there is no gap in sequence a
                numCharactersInA++;
            } else {  // match or mismatch
                numCharactersInA++;
                numCharactersInB++;
                numMatchesOrMismatches++;

                numMatches += sequenceA[k] === sequenceB[k] ? 1 : 0;  // if match, then increment

                sequenceIdentity = (100 * numMatches) / numMatchesOrMismatches;

                if (global)
                    L[[numCharactersInA, numCharactersInB]] = sequenceIdentity;  // creating key with non-final value
                else  // local
                    L[[currentPosition.i, currentPosition.j]] = sequenceIdentity;
            }
        }

        var definedKeys = Object.keys(L);

        // set final weight: L^{ii'}(jj') = L^{ii'}(jj') + weight(A)
        // where initial L^{ii'}(jj') = 0 and weight(A) = seqIdentity(A)
        for (var i = 0; i < definedKeys.length; i++)
            L[definedKeys[i]] = sequenceIdentity;  // overwriting wrong value

        return L;
    }

    /**
     * The elements of both weight libraries added into one primary weight library.
     * If the key (alignment) [a, b] of an element primLib^{a,b} is present in both libraries,
     * then the weights for this key [a, b] retrieved in both libraries and added together.
     * @see: https://doi.org/10.1006/jmbi.2000.4042
     * Notredame, Cédric, Desmond G. Higgins, and Jaap Heringa.
     * "T-Coffee: A novel method for fast and accurate multiple sequence alignment."
     * Journal of molecular biology 302.1 (2000): 205-217.
     *
     * Hint: Described on p.207: Combination of the libraries
     */
    function addSignals() {
        debugger;

        var globalAlignmentKeys = Object.keys(outputData.primaryGlobalWeightLib);  // global alignments {{a,b}, {a,d}, ... {d,f}}
        var localAlignmentKeys = Object.keys(outputData.primaryLocalWeightLib);  // local alignments {{a,b}, {a,g}, ... {d,f}}
        var commonAlignmentKeys = getCommonArguments(globalAlignmentKeys, localAlignmentKeys);  // common alignments {{a,b}, {d,f}}

        outputData.primaryWeightLib = {};

        // append global, non common elements into primary weight library -> primLib^{a,b} =  {L_{1,3}, L_{2,4}, ..., L_{5,7}}
        outputData.primaryWeightLib
            = appendToLibrary(outputData.primaryWeightLib, outputData.primaryGlobalWeightLib, globalAlignmentKeys, commonAlignmentKeys);

        // append local, non common elements into primary weight library -> primLib^{a,b} = {L_{2,3}, L_{2,5}, ..., L_{5,5}}
        outputData.primaryWeightLib
            = appendToLibrary(outputData.primaryWeightLib, outputData.primaryLocalWeightLib, localAlignmentKeys, commonAlignmentKeys);

        // append common elements into primary weight library
        for (var i = 0; i < commonAlignmentKeys.length; i++) {
            var globalWeights = outputData.primaryGlobalWeightLib[commonAlignmentKeys[i]];  // {L_{1,3}, L_{2,4}, ..., L_{5,7}}
            var localWeights = outputData.primaryLocalWeightLib[commonAlignmentKeys[i]];  // {L_{2,3}, L_{2,5}, ..., L_{5,5}}

            var globalWeightKeys = Object.keys(globalWeights);  // {{1,3}, {2,4}, ..., {5,7}}
            var localWeightKeys = Object.keys(localWeights);  // {{2,3}, {2,5}, ..., {5,5}}
            var commonWeightKeys = getCommonArguments(globalWeightKeys, localWeightKeys);  // {{3,3}, {3,5}, ..., {5,4}}

            var weights = {};

            // append global, non common weights
            weights = appendToLibrary(weights, globalWeights, globalWeightKeys, commonWeightKeys);

            // append local, non common weights
            weights = appendToLibrary(weights, localWeights, localWeightKeys, commonWeightKeys);

            // add common weights of both libraries together and append
            for (var j = 0; j < commonWeightKeys.length; j++)
                weights[[commonWeightKeys[j]]] = globalWeights[commonWeightKeys[j]] + localWeights[commonWeightKeys[j]];

            outputData.primaryWeightLib[commonAlignmentKeys[i]] = weights;
        }
    }

    /**
     * Returns the common keys from two arrays.
     * @param keys1 {Array} - The keys from the first structure.
     * @param keys2 {Array} - The keys from the second structure.
     * @return {Array} - The common keys.
     */
    function getCommonArguments(keys1, keys2) {
        var commonKeys = [];

        var shorterSequenceLength = 0;
        var shorterKey = undefined;
        var longerKey = undefined;

        if (keys1.length < keys2.length) {
            shorterSequenceLength = keys1.length;
            shorterKey = keys1;
            longerKey = keys2;
        } else {
            shorterSequenceLength = keys2.length;
            shorterKey = keys2;
            longerKey = keys1;
        }


        for (var i = 0; i < shorterSequenceLength; i++) {
            if (longerKey.indexOf(shorterKey[i]) >= 0)  // if (key from keys1 contained in keys2)
                commonKeys.push(shorterKey[i]);
        }

        return commonKeys;
    }

    /**
     * Append elements from one structure into a second structure by using the given keys.
     * @param libraryWriteTo {Object} - The structure in which elements are stored.
     * @param libraryReadFrom {Object} - The structure from which elements are read.
     * @param keys {Array} - Keys, which are read from the "readFrom"-structure.
     * @param nonKeys {Array} - Keys of the elements, that are not written in the "WriteTo"-structure (read, but not written).
     * @return {Object} - The structure in which elements were appended.
     */
    function appendToLibrary(libraryWriteTo, libraryReadFrom, keys, nonKeys) {
        // append non common elements into libraryWriteTo
        for (var i = 0; i < keys.length; i++) {
            if (nonKeys.indexOf(keys[i]) === -1)  // if (not contained in nonKeys)
                libraryWriteTo[keys[i]] = libraryReadFrom[keys[i]];
        }

        return libraryWriteTo;
    }

    /**
     * The weights in the primary weight library are recomputed
     * to add consistency-information.
     * @example:
     * \forall(ab) \forall(ij):
     * EL^{ab}(ij) = L^{ab}(ij) + \sum_{x \in S_aligned} \sum_{k in pos(x)} min(L^{ax}(ik), L^{xb}(kj))
     * is computed.
     */
    function computeExtendedWeightPrimaryLibrary() {
        outputData.extendedWeightLib = {};  // extLib^{a,b}

        var outerPrimLibKeys = Object.keys(outputData.primaryWeightLib);  // [(a,b), (a,c), ..., (d,f)] where [character] = [ALIGNED SEQUENCE]
        var alignmentSequenceNames = getIndividualArguments(outerPrimLibKeys);  // [a, b, c, ..., f]

        // iterate over each element in primary library (so over all ii')
        for (var i = 0; i < outerPrimLibKeys.length; i++) {
            var alignmentKey = outerPrimLibKeys[i].split(SYMBOLS.COMMA);  // [a,b]
            var leftAlignmentKeyArgument = alignmentKey[0];  // a -> [character] = [ALIGNED SEQUENCE]
            var rightAlignmentKeyArgument = alignmentKey[1];  // b

            var innerPrimLibKeys = Object.keys(outputData.primaryWeightLib[outerPrimLibKeys[i]]);  // [(2,3), (2,5), ..., (5,7)]

            outputData.extendedWeightLib[outerPrimLibKeys[i]] = {};  // extLib^{a,b}(i,j)

            for (var j = 0; j < innerPrimLibKeys.length; j++) {  // iteration over each weight
                var weightKey = innerPrimLibKeys[j].split(SYMBOLS.COMMA);  // [i,j]

                var leftWeightKeyArgument = weightKey[0];  // i
                var rightWeightKeyArgument = weightKey[1];  // j

                var sum = computeExtensionSum(
                    alignmentSequenceNames, leftAlignmentKeyArgument, rightAlignmentKeyArgument,
                    leftWeightKeyArgument, rightWeightKeyArgument);

                outputData.extendedWeightLib[outerPrimLibKeys[i]][innerPrimLibKeys[j]]
                    = outputData.primaryWeightLib[outerPrimLibKeys[i]][innerPrimLibKeys[j]] + sum;
            }
        }
    }

    /**
     * Returns the individual arguments of a set of key-pairs.
     * @param keyPairs - The key-pairs from which the arguments are retrieved.
     * @return {Array} - The set of individual arguments.
     */
    function getIndividualArguments(keyPairs) {
        var args = [];  // "arguments" is Javascript keyword

        for (var i = 0; i < keyPairs.length; i++) {
            var keyPair = keyPairs[i].split(SYMBOLS.COMMA);
            var leftArg = keyPair[0];
            var rightArg = keyPair[1];

            if (args.indexOf(leftArg) === -1)  // if (not contained)
                args.push(leftArg);
            if (args.indexOf(rightArg) === -1)  // if (not contained)
                args.push(rightArg);
        }

        return args;
    }

    /**
     * Computes the extension.
     * @example:
     * \sum_{x \in S_aligned} \sum_{k in pos(x)} min(L^{ax}(ik), L^{xb}(kj))
     * @param alignmentSequenceNames {Array} - The different aligned sequence-strings: [a, b, c, ..., f]
     * @param leftAlignmentKeyArgument {string} - The aligned left sequence which is used in a key-pair (a,b).
     * @param rightAlignmentKeyArgument {string} - The aligned sequence which is used in a key-pair (a,b).
     * @param leftWeightKeyArgument {number} - The number which is used in a key-pair (i,j).
     * @param rightWeightKeyArgument {number} - The number which is used in a key-pair (i,j).
     * @return {number} - The extension sum which is added to a library weight to add consistency information.
     */
    function computeExtensionSum(alignmentSequenceNames, leftAlignmentKeyArgument, rightAlignmentKeyArgument,
                                 leftWeightKeyArgument, rightWeightKeyArgument) {
        var sum = 0;

        // iterate overall aligned sequence x
        for (var x = 0; x < alignmentSequenceNames.length; x++) {

            if (alignmentSequenceNames[x] !== leftAlignmentKeyArgument
                && alignmentSequenceNames[x] !== rightAlignmentKeyArgument) {  // just an optimization (x in S\{a,b})

                var positionSequenceNames = getPositions(alignmentSequenceNames[x]);

                // iterate overall positions in aligned sequence x
                for (var k = 0; k < positionSequenceNames.length; k++) {

                    sum += getMinimumWeight(
                        alignmentSequenceNames, leftAlignmentKeyArgument, rightAlignmentKeyArgument,
                        positionSequenceNames, leftWeightKeyArgument, rightWeightKeyArgument, x, k);
                }
            }
        }

        return sum;
    }

    /**
     * Return the possible positions in an aligned sequence without counting gap positions.
     * @example A__T returns [1,2]
     * @param sequence {string} - The sequence in which the non-gap positions are searched.
     * @return {Array} - The non-gap positions of the input sequence.
     */
    function getPositions(sequence) {
        var positions = [];

        var position = 1;   // "+1" because counting starts in this algorithm with 1

        for (var i = 0; i < sequence.length; i++) {
            if (sequence[i] !== SYMBOLS.GAP) {
                positions.push(position);
                position++;
            }
        }

        return positions;
    }

    /**
     * Computes the minimum weight.
     * @example:
     * min(L^{ax}(ik), L^{xb}(kj))
     * where L^{ab}(ij) = L^{ba}(ji)
     * @param alignmentSequenceNames {Array} - The different aligned sequence-strings: [a, b, c, ..., f]
     * @param leftAlignmentKeyArgument {string} - The aligned left sequence which is used in a key-pair (a,b).
     * @param rightAlignmentKeyArgument {string} - The aligned sequence which is used in a key-pair (a,b).
     * @param positionSequenceNames {Array} - The different positions within a aligned sequence: [1, 2, ..., n]
     * @param leftWeightKeyArgument {number} - The number which is used in a key-pair (i,j).
     * @param rightWeightKeyArgument {number} - The number which is used in a key-pair (i,j).
     * @param x {number} - The position within alignmentSequenceNames which should be called.
     * @param k {number} - The position within positionSequenceNames which should be called.
     * @return {number} - The weight which is added to a library weight to add consistency information.
     */
    function getMinimumWeight(alignmentSequenceNames, leftAlignmentKeyArgument, rightAlignmentKeyArgument,
                              positionSequenceNames, leftWeightKeyArgument, rightWeightKeyArgument, x, k) {
        var weightStruct1 = outputData.primaryWeightLib[[leftAlignmentKeyArgument, alignmentSequenceNames[x]]];  // {ax}
        var weightStruct2 = outputData.primaryWeightLib[[alignmentSequenceNames[x], rightAlignmentKeyArgument]];  // {xb}

        var symmetry1 = false;  // it holds: L^{ab}(ij) = L^{ba}(ji)
        var symmetry2 = false;

        if (weightStruct1 === undefined) {
            weightStruct1 = outputData.primaryWeightLib[[alignmentSequenceNames[x], leftAlignmentKeyArgument]];  // {xa}
            symmetry1 = true;
        }

        if (weightStruct2 === undefined) {
            weightStruct2 = outputData.primaryWeightLib[[rightAlignmentKeyArgument, alignmentSequenceNames[x]]];  // {bx}
            symmetry2 = true;
        }

        var weight1 = getWeight(weightStruct1, symmetry1, positionSequenceNames, leftWeightKeyArgument, false, k);
        var weight2 = getWeight(weightStruct2, symmetry2, positionSequenceNames, rightWeightKeyArgument, true, k);

        return Math.min(weight1, weight2);
    }

    /**
     * Returns the weight.
     * @example
     * L^{ax}(ik) or L^{xb}(kj)
     * @param weightStruct {Object} - The structure L which is used to output a weight.
     * @param symmetry {boolean} - Tells if arguments have to be swapped in their order (e.g. ik to ki).
     * @param positionSequenceNames {Array} - The different positions within a aligned sequence: [1, 2, ..., n]
     * @param weightKeyArgument {number} - The number which is used in a key-pair (i,j).
     * @param right {boolean} - Tells if weightKeyArgument is a right or left weight key argument.
     * @param k {number} - The position within positionSequenceNames which should be called.
     * @return {number} - The weight of the weight structure at a specific position.
     */
    function getWeight(weightStruct, symmetry, positionSequenceNames, weightKeyArgument, right, k) {
        var weight = 0;

        symmetry = right ? !symmetry : symmetry;

        if (weightStruct !== undefined) {
            if (symmetry)
                weight = weightStruct[[positionSequenceNames[k], weightKeyArgument]];  // (ki) or (kj)
            else
                weight = weightStruct[[weightKeyArgument, positionSequenceNames[k]]];  // (ik) or (jk)

            weight = weight !== undefined ? weight : 0;
        }

        return weight;
    }

    /**
     * Starts the creation of a progressive alignment under the computed, extended weight library as a scoring-function.
     */
    function startProgressiveAlignment() {
        // use all pairwise scores, gaps, ... from alignments in step 1 to compute distances
        // problem: we have local and global alignment data
        // ideas:
        // (1)  take average scores, gaps, ... and compute then distances for a distance matrix
        // (2)  just use global alignment data for the distance matrix
        // (3)  compute distance matrix first with global alignment data,
        //      then with local alignment data and take then average distance matrix = (D_local + D_global) / 2
        //
        // but: second possibility (2) taken, because
        // [1] in other sources it is said that local alignments less important
        // http://www.bioinfbook.org/ : Bioinformatics and Functional Genomics 2nd Edition p.193
        // (not all, only 10 local alignments for library computation used)
        // [2] we want create "global" progressive alignments

        multiSequenceAlignmentInstance.setIO(inputData, outputData);
        multiSequenceAlignmentInstance.computeDistancesFromSimilarities();
        multiSequenceAlignmentInstance.createDistanceMatrix();
        createProgressiveAlignment(multiSequenceAlignmentInstance.getPhylogeneticTree());
    }

    /**
     * By going through the guide tree branches (in correct merging order),
     * the algorithm generates progressive alignments.
     * Hint: A repeated post order traversal of the tree would be less efficient.
     * This is why just an iteration through the branches is done.
     * @param treeBranches {Object} - The tree branches which are defining the order for the merging process.
     * @see: https://doi.org/10.1006/jmbi.2000.4042
     * Notredame, Cédric, Desmond G. Higgins, and Jaap Heringa.
     * "T-Coffee: A novel method for fast and accurate multiple sequence alignment."
     * Journal of molecular biology 302.1 (2000): 205-217.
     */
    function createProgressiveAlignment(treeBranches) {
        // store
        var enlargement = inputData.enlargement;
        var baseCosts = inputData.baseCosts;

        inputData.substitutionFunction = substitutionFunction;
        inputData.enlargement = 0;  // p.210 top-left
        inputData.baseCosts = 0;

        multiSequenceAlignmentInstance.setIO(inputData, outputData);
        multiSequenceAlignmentInstance.createProgressiveAlignment(treeBranches);

        // restore
        inputData.enlargement = enlargement;
        inputData.baseCosts = baseCosts;
    }

    /**
     * Returns the position specific scoring.
     * @param i {number} - The position in the first sequence.
     * @param j {number} - The position in the second sequence.
     * @returns {number} - The position specific score.
     * @see: https://doi.org/10.1006/jmbi.2000.4042
     * Notredame, Cédric, Desmond G. Higgins, and Jaap Heringa.
     * "T-Coffee: A novel method for fast and accurate multiple sequence alignment."
     * Journal of molecular biology 302.1 (2000): 205-217.
     *
     * Hint: Described on p.209: Progressive Alignment Strategy
     * When aligning group sequences with other group sequences,
     * the average library score of the group columns is used.
     * Else usual extended library score is used.
     */
    function substitutionFunction(i, j) {
        var group1 = outputData.currentFirstGroup;
        var group2 = outputData.currentSecondGroup;

        return getPairwiseColumnWeight(group1, group2, i, j) / (group1.length * group2.length);
    }

    /**
     * Returns the average weight of two groups at a specific position.
     * @param group1 {Array} - The group in which the score for a column is computed.
     * @param group2 {Array} - The group in which the score for a column is computed.
     * @param i {number} - The row in which is looked up.
     * @param j {number} - The column in which is looked up.
     * @return {number} - The weight.
     */
    function getPairwiseColumnWeight(group1, group2, i, j) {
        var weight = 0;

        for (var k = 0; k < group1.length; k++) {
            for (var l = 0; l < group2.length; l++) {
                var preSequenceA = group1[k];  // contains neutral symbol SYMBOLS.NONE
                var preSequenceB = group2[l];

                var sequenceA = preSequenceA.replace(MULTI_SYMBOLS.NONE, SYMBOLS.EMPTY);  // without neutral symbol SYMBOLS.NONE
                var sequenceB = preSequenceB.replace(MULTI_SYMBOLS.NONE, SYMBOLS.EMPTY);

                var L = outputData.extendedWeightLib[[sequenceA, sequenceB]];

                var switchArguments = false;

                if (L === undefined) {
                    L = outputData.extendedWeightLib[[sequenceB, sequenceA]];
                    switchArguments = true;  // L^{a,b}(i,j) = L^{b,a}(j,i)
                }

                if (L !== undefined) {
                    var argI = i - getNumberOfNeutrals(preSequenceA, i);
                    var argJ = j - getNumberOfNeutrals(preSequenceB, j);

                    if (switchArguments)
                        weight += L[[argJ,argI]] !== undefined ? L[[argJ,argI]] : 0;
                    else
                        weight += L[[argI,argJ]] !== undefined ? L[[argI,argJ]] : 0;  // not defined positions have always score 0
                }
            }
        }

        return weight;
    }

    /**
     * Returns the number of neutrals in the sequence up to a certain position.
     * @param sequence {string} - The sequence in which neutral symbols are counted.
     * @param position {number} - The position up to which neutral symbols are inclusively counted.
     * @return {number} - The number of neutral symbols in the alignment.
     */
    function getNumberOfNeutrals(sequence, position) {
        var numOfNeutrals = 0;

        for (var i = 0; i <= position; i++) {
            if (sequence[i] === SYMBOLS.NONE)
                numOfNeutrals++;
        }

        return numOfNeutrals;
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