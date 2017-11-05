/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("bases.multiSequenceAlignment", MultiSequenceAlignment, getInput, setInput, getOutput, setIO, getLastChild);

    // instances
    var alignmentInstance;
    var childInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Contains functions to compute multi-sequence alignments.
     * It is used multi-sequence alignment algorithms as superclass
     * to avoid code duplicates.
     * @param child - The child algorithm which inherits from this class.
     * @constructor
     */
    function MultiSequenceAlignment(child) {
        alignmentInstance = this;
        childInstance = child;

        // public methods
        this.getInput = getInput;
        this.setInput = setInput;
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

        inputData.calculationType = inputViewmodel.calculation();

        inputData.baseCosts = inputViewmodel.baseCosts();
        inputData.enlargement = inputViewmodel.enlargement();
        inputData.match = inputViewmodel.match();
        inputData.mismatch = inputViewmodel.mismatch();
    }

    /**
     * Returns the non-first position of sequences which are multiple times in the sequences.
     * @param sequences - The sequences which are checked for duplicates.
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
                duplicatePositions = duplicatePositions.concat(nonFirstPositions)
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
