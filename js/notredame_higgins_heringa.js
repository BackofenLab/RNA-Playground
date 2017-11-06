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
     * If the alignment [a, b] of an element L^{a,b} is contained in both libraries,
     * then the weights of this element [a, b] retrieved in both libraries and added together.
     */
    function addSignals() {

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