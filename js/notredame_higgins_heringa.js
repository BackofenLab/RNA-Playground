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

    var fengDoolittleInstance;
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

        fengDoolittleInstance = new fengDoolittle.FengDoolittle();
        gotohInstance = new gotoh.Gotoh();
        gotohLocalInstance = new gotohLocalInstance.GotohLocal();

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
    }

    /**
     * Starts the computation.
     */
    function compute() {
        computePrimaryLibraries();
        computeCombinedWeightPrimaryLibrary();
        computeExtendedWeightPrimaryLibrary();
        createProgressiveAlignment();
        return [inputData, outputData];
    }

    /**
     * Computes the set of pairwise alignments for local and global alignments.
     */
    function computePrimaryLibraries() {
        computePairwiseGlobalAlignmentData();
        computePairwiseLocalAlignmentData();
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
        outputData.primaryLocalWeightLib = computePairwiseWeights(localOutputData);
        addSignals();
    }

    /**
     * Computes the sequence identity.
     * So, how much is identical between two sequences
     * with respect to the smaller sequence.
     * Hint: It could be computed during computation of alignment data,
     * but for better understanding and less code complexity
     * nearly everything computed in order defined in the original paper.
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
     * Computes a structure L dependant on an argument of the form [i,j]
     * which returns the position specific sequence identity.
     * @param alignment - The alignment for which you want compute
     * a structure that returns position specific sequence identities.
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

        var L = {};

        // iterate over each position to compute structures L_{i,j}
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

                L[[numCharactersInA, numCharactersInB]] = (100 * numMatches) / numMatchesOrMismatches;
            }
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
    }

    /**
     * Creates a progressive alignment with Feng-Doolittle
     * and the computed, extended weight library.
     */
    function createProgressiveAlignment() {
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