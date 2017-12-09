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
        //outputData.score = multiSequenceAlignmentInstance.getAffineSumOfPairsScore(inputData , outputData.progressiveAlignment);
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