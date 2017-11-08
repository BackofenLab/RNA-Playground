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
    var globalOutputData = {};
    var localOutputData = {};

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
        multiSequenceAlignmentInstance.setIO(inputData, {});
        multiSequenceAlignmentInstance.setInput(inputViewmodel);

        inputData.maxNumberOptimalAlignments = inputViewmodel.maxNumberOptimalAlignments();
        inputData.maxNumberOptimalAlignmentsLocal = inputViewmodel.maxNumberOptimalAlignmentsLocal();
    }

    /**
     * Starts the computation.
     */
    function compute() {
        computePrimaryLibraries();
        computeCombinedWeightPrimaryLibrary();
        computeExtendedWeightPrimaryLibrary();
        startProgressiveAlignment();
        return [inputData, outputData];
    }

    /**
     * Computes the set of pairwise alignments for local and global alignments.
     */
    function computePrimaryLibraries() {
        computePairwiseGlobalAlignmentData();
        //computePairwiseLocalAlignmentData();
    }

    /**
     * Computes scores (similarities),
     * the number of gaps, the alignment lengths
     * and so on between all sequences.
     */
    function computePairwiseGlobalAlignmentData() {
        multiSequenceAlignmentInstance.setIO(inputData, globalOutputData);
        multiSequenceAlignmentInstance.computePairwiseData(gotohInstance);
    }

    function computePairwiseLocalAlignmentData() {
        /*
        // later
        inputData.baseCosts = inputViewmodel.baseCostsLocal();
        inputData.enlargement = inputViewmodel.enlargementLocal();
        inputData.match = inputViewmodel.matchLocal();
        inputData.mismatch = inputViewmodel.mismatchLocal();
        */
        multiSequenceAlignmentInstance.setIO(inputData, localOutputData);
        multiSequenceAlignmentInstance.computePairwiseData(gotohLocalInstance);
    }

    /**
     * Computes the weights for both libraries
     * and combines both libraries to one big library (signal addition).
     */
    function computeCombinedWeightPrimaryLibrary() {
        // Hint: conversion of double sequences into alignment edges is directly done during the computations of pairwise weights
        outputData.primaryGlobalWeightLib = computePairwiseWeights(globalOutputData);
        // later
        // outputData.primaryLocalWeightLib = computePairwiseWeights(localOutputData);
        // addSignals();
    }

    /**
     * Computes the sequence identity.
     * So, how much is identical between two sequences
     * with respect to the smaller sequence.
     * @param {Object} - The output which is used to compute a primary library.
     */
    function computePairwiseWeights(output) {
        var primaryWeightLib = {};

        // iterate over each sequence a and sequence b to compute structure primLib^{a,b}(i,j) = {L_{1,3}, L_{2,4}, ..., L_{5,7}}
        // Hint: a and b are the aligned sequences
        for (var i = 0; i < inputData.sequences.length; i++) {
            for (var j = 0; j < i; j++) {
                var sequenceA = inputData.sequences[j];
                var sequenceB = inputData.sequences[i];

                var alignment = output.alignmentsAndScores[[sequenceA, sequenceB]][0];
                var sequenceIdentities = getSequenceIdentities(alignment);  // alignment = [alignedSequenceA, matchMismatchString, alignedSequenceB]
                primaryWeightLib[[alignment[0], alignment[2]]] = sequenceIdentities;
            }
        }

        return primaryWeightLib;
    }

    /**
     * Computes a structure L with elements
     * containing the sequence identity.
     * @param alignment - The alignment for which you want compute
     * a structure that returns position sequence identities.
     * @return {Object}
     */
    function getSequenceIdentities(alignment) {
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
            if (sequenceA[k] === SYMBOLS.GAP) {  // means there is a gap in sequence b
                numCharactersInB++;
            } else if (sequenceB[k] === SYMBOLS.GAP) {  // means there is a gap in sequence a
                numCharactersInA++;
            } else {  // match or mismatch
                numCharactersInA++;
                numCharactersInB++;
                numMatchesOrMismatches++;

                numMatches += sequenceA[k] === sequenceB[k] ? 1 : 0;  // if match, then increment

                sequenceIdentity = (100 * numMatches) / numMatchesOrMismatches;
                L[[numCharactersInA, numCharactersInB]] = sequenceIdentity;  // creating key with non-final value
            }
        }

        var definedKeys = Object.keys(L);  //

        // set final weight: L^{ii'}(jj') = L^{ii'}(jj') + weight(A)
        // where initial L^{ii'}(jj') = 0 and weight(A) = seqIdentity(A)
        for (var i = 0; i < definedKeys.length; i++) {
            L[definedKeys[i]] = sequenceIdentity;  // overwriting wrong value
        }

        return L;
    }

    /**
     * The elements of both weight libraries added into one primary weight library.
     * If the key (alignment) [a, b] of an element primLib^{a,b} is present in both libraries,
     * then the weights for this key [a, b] retrieved in both libraries and added together.
     */
    function addSignals() {
        var globalAlignmentKeys = Object.keys(outputData.primaryGlobalWeightLib);  // global alignments {{a,b}, {a,d}, ... {d,f}}
        var localAlignmentKeys = Object.keys(outputData.primaryLocalWeightLib);  // local alignments {{a,b}, {a,g}, ... {d,f}}
        var commonAlignmentKeys = getCommonArguments(globalAlignmentKeys, localAlignmentKeys);  // common alignments {{a,b}, {d,f}}

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
            weights = appendToLibrary(weights, globalWeights, globalWeightKeys, commonWeightKeys);

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

        var shorterSequenceLength = keys1.length < keys2.length ? keys1.length : keys2.length;

        for (var i = 0; i < shorterSequenceLength; i++) {
            if (keys2.indexOf(keys1[i]) >= 0)  // if (key from keys1 contained in keys2)
                commonKeys.push(keys1[i]);
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
     */
    function computeExtendedWeightPrimaryLibrary() {
        outputData.extendedWeightLib = {};  // extLib^{a,b}
        var outerPrimLibKeys = Object.keys(outputData.primaryWeightLib);  // {{a,b}, {a,c}, ..., {d,f}}
        var alignmentSequenceNames = getIndividualArguments(outerPrimLibKeys);  // [a, b, c, ..., f]

        // iterate over each element in primary library (so over all ii')
        for (var i = 0; i < outerPrimLibKeys.length; i++) {
            var alignmentKey = outerPrimLibKeys[i].split(SYMBOLS.COMMA);  // [a,b]
            var leftAlignmentKeyArgument = alignmentKey[0];  // a
            var rightAlignmentKeyArgument = alignmentKey[1];  // b

            var innerPrimLibKeys = Object.keys(outputData.primaryWeightLib[outerPrimLibKeys[i]]);  // {{2,3}, {2,5}, ..., {5,7}}

            var positionSequenceNames = getIndividualArguments(innerPrimLibKeys);

            outputData.extendedWeightLib[outerPrimLibKeys[i]] = {};  // extLib^{a,b}(i,j)

            // computation of following formula
            // EL^{ab}(ij) = L^{ab}(ij) + \sum_{x \in S_aligned} \sum_{k in pos(x)} min(L^{ax}(ik), L^{xb}(kj))
            for (var j = 0; j < innerPrimLibKeys.length; j++) {  // iteration over each weight
                var weightKey = innerPrimLibKeys[j].split(SYMBOLS.COMMA);  // [i,j]
                var leftWeightKeyArgument = weightKey[0];  // i
                var rightWeightKeyArgument = weightKey[1];  // j

                var sum = 0;

                // iterate overall aligned sequence x
                for (var x = 0; x < alignmentSequenceNames.length; x++) {

                    // iterate overall positions in aligned sequence x
                    for (var k = 0; k < positionSequenceNames.length; k++) {
                        var alignment1 = outputData.primaryWeightLib[leftAlignmentKeyArgument, alignmentSequenceNames[x]];  // {ax}
                        var alignment2 = outputData.primaryWeightLib[alignmentSequenceNames[x], rightAlignmentKeyArgument];  // {xb}

                        var weight1 = 0;
                        var weight2 = 0;

                        if (alignment1 !== undefined) {
                            weight1 = alignment1[leftWeightKeyArgument, positionSequenceNames[k]];  // (ik)
                            weight1 = weight1 !== undefined ? weight1 : 0;
                        }

                        if (alignment2 !== undefined) {
                            weight2 = alignment2[positionSequenceNames[k], rightWeightKeyArgument];  // (kj)
                            weight2 = weight2 !== undefined ? weight2 : 0;
                        }

                        sum += Math.min(weight1, weight2);
                    }
                }

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
            var leftArg = keyPairs[i][0];
            var rightArg = keyPairs[i][1];

            if (args.indexOf(leftArg) === -1)  // if (not contained)
                args.push(leftArg);
            if (args.indexOf(rightArg) === -1)  // if (not contained)
                args.push(rightArg);
        }

        return args;
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
        // [2] we want create global progressive alignments

        multiSequenceAlignmentInstance.setIO(inputData, globalOutputData);
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
     * @see: It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function createProgressiveAlignment() {
        inputData.substitutionFunction = substitutionFunction;
        inputData.enlargement = 0;
        inputData.gapCosts = 0;

        multiSequenceAlignmentInstance.setIO(inputData, globalOutputData);
    }

    /**
     * Returns the position specific scoring.
     * @param i {number} - The position in the first sequence.
     * @param j {number} - The position in the second sequence.
     * @returns {number} - The position specific score.
     */
    function substitutionFunction(i, j) {
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