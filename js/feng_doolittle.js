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
    namespace("fengDoolittle", startFengDoolittle, FengDoolittle);

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
     * @see https://doi.org/10.1016/S0076-6879(96)66023-6
     *
     * Feng, Da-Fei, and Russell F. Doolittle.
     * "[21] Progressive alignment of amino acid sequences and construction of phylogenetic trees from them."
     * Methods in enzymology 266 (1996): 368-382.
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

        multiSequenceAlignmentInstance.computePairwiseData(gotohInstance);  // computes scores, number of gaps, alignment lengths, ...
        multiSequenceAlignmentInstance.computeDistancesFromSimilarities(); // converting the pairwise similarities into distances
        multiSequenceAlignmentInstance.createDistanceMatrix();  // creates a function dist(a,b) between cluster names and distances
        multiSequenceAlignmentInstance.createProgressiveAlignment(multiSequenceAlignmentInstance.getPhylogeneticTree());

        outputData.score = multiSequenceAlignmentInstance.getAffineSumOfPairsScore(inputData, outputData.progressiveAlignment);
        return [inputData, outputData];
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