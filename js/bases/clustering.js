/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("bases.clustering", Clustering, setInput, compute, getOutput, setIO);

    // instances
    var childInstance;
    var clusteringInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Contains functions to compute (hierarchical) agglomerative clustering.
     * It is used by agglomerative clustering algorithms as superclass
     * to avoid code duplicates.
     * @constructor
     */
    function Clustering(child) {
        clusteringInstance = this;

        // variables
        this.cardinalities = {};  // needed for computations (for example in UPGMA)
        this.nameIndex = 0;
        this.treeParts = [];  // parts of the hierarchical tree

        // inheritance
        childInstance = child;

        // public methods (linking)
        this.setInput = setInput;

        this.compute = compute;
        this.getOutput = getOutput;

        this.setIO = setIO;
    }

    /**
     * Sets the input data of the algorithm.
     * @return {Object} - Contains all input data.
     */
    function setInput(input) {
        inputData = input;
    }

    /**
     * Starts the computation.
     * Hint: Because the distance matrix is changing during the following procedure
     * and it has to be visualized later on, a copy is made which is written back.
     *
     */
    function compute() {
        var distanceMatrixCopy = jQuery.extend(true, {}, outputData.distanceMatrix);  // deep copy

        var numOfIterations = inputData.sequences.length - 1;  // always lower by one in hierarchical clustering algorithms

        initializeCardinalities(numOfIterations);

        for (var i = 0; i < numOfIterations; i++) {
            var minimum = determineMatrixMinimum();  // minimum = [key, value]
            var newClusterName = mergeClusters(minimum.clusterName1, minimum.clusterName2);
            clusteringInstance
                .treeParts.push([minimum.clusterName1, minimum.clusterName2, newClusterName, minimum.distance / 2]);
            recomputeDistances(newClusterName);
        }

        outputData.distanceMatrix = distanceMatrixCopy;  // write-back
        return [inputData, outputData];
    }

    /**
     * Initializes the size parameters of the clusters.
     * @param numOfIterations - The number of iterations the algorithm will do.
     */
    function initializeCardinalities(numOfIterations) {
        clusteringInstance.nameIndex = numOfIterations + 1;  // current cluster-name set-size + 1

        for (var i = 0; i < clusteringInstance.nameIndex; i++)
            clusteringInstance.cardinalities[CLUSTER_NAMES[i]] = 1;
    }

    /**
     * Determines the entry with the lowest distance
     * in the remaining entries.
     * @return {Object} - The matrix-entry which contains the minimum value.
     * @see: It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function determineMatrixMinimum() {
        var minKey = SYMBOLS.EMPTY;
        var minValue = Number.POSITIVE_INFINITY;

        var keys = Object.keys(outputData.distanceMatrix);

        // searching for minimum by going over all keys
        for (var i = 0; i < keys.length; i++) {
            var key = keys[i].split(SYMBOLS.COMMA);  // string to array: Javascript stores keys as strings
            var value = outputData.distanceMatrix[key];

            if (value < minValue) {  // hint: diagonals were not computed because they are zero
                minKey = key;  // name: [a,b] with a,b cluster names
                minValue = value;  // distance
            }
        }

        // create structure for better understandable access
        var minimum = {};
        minimum.clusterName1 = minKey[0];
        minimum.clusterName2 = minKey[1];
        minimum.distance = minValue;
        return minimum;
    }

    /**
     * Computes the union of the two clusters which form the minimum.
     * @param clusterName1 - The name of the first cluster.
     * @param clusterName2 - The name of the second cluster.
     * @return {string} - The name of the new cluster.
     */
    function mergeClusters(clusterName1, clusterName2) {
        var newClusterName = createNewCluster(clusterName1, clusterName2);
        delete outputData.distanceMatrix[[clusterName1, clusterName2]];  // removing key

        return newClusterName;
    }

    /**
     * Creates the union of the two given clusters.
     * @param clusterName1 - The name of the first cluster.
     * @param clusterName2 - The name of the second cluster.
     * @return {string} - The name of the new cluster.
     * @see: It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function createNewCluster(clusterName1, clusterName2) {
        var newClusterName = getNextClusterName();

        var firstClusterCardinality = clusteringInstance.cardinalities[clusterName1];
        var SecondClusterCardinality = clusteringInstance.cardinalities[clusterName2];
        clusteringInstance.cardinalities[newClusterName] = firstClusterCardinality + SecondClusterCardinality;

        return newClusterName;
    }

    /**
     * Returns the next name of a cluster.
     * Hint: After all characters are depleted,
     * a number is concatenated to the character
     * to make this function generic.
     * @example: (with 26 characters)
     * CLUSTER NAMES:
     * a, b, c, ..., z,         FIRST EPISODE   (0 <= index < 26)
     * a2, b2, c2, ..., z2,     SECOND EPISODE  (26 <= index < 52)
     * a3, b3, ...              THIRD ...       (52 <= index < 78)
     * @return {string} - Cluster name.
     */
    function getNextClusterName() {
        var clusterName = SYMBOLS.EMPTY;

        if (clusteringInstance.nameIndex >= CLUSTER_NAMES.length) {
            var quotient = Math.floor(clusteringInstance.nameIndex / CLUSTER_NAMES.length);
            var remainder = clusteringInstance.nameIndex % CLUSTER_NAMES.length;

            clusterName = CLUSTER_NAMES[remainder] + SYMBOLS.EMPTY + quotient;
        } else
            clusterName = CLUSTER_NAMES[clusteringInstance.nameIndex];

        clusteringInstance.nameIndex++;  // to get next time the next cluster name

        return clusterName;
    }

    /**
     * Recomputes the distance to the other clusters with the approach of the child class.
     * @param clusterName {string} - The name of the new cluster.
     */
    function recomputeDistances(clusterName) {
        childInstance.recomputeDistances(clusterName);
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
}());
