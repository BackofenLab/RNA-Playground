/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("bases.clustering", Clustering, setInput, compute, getMatrixKeys, getOutput, setIO);

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
        this.nameIndex = 0;

        // inheritance
        childInstance = child;

        // public methods
        this.setInput = setInput;
        this.compute = compute;
        this.getMatrixKeys = getMatrixKeys;
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
     */
    function compute() {
        var distanceMatrixCopy = jQuery.extend(true, {}, outputData.distanceMatrix);  // deep copy

        var numOfIterations = inputData.numOfStartClusters - 1;  // always lower by one in hierarchical clustering algorithms

        initializeStructs();
        initializeCardinalities(numOfIterations);

        for (var i = 0; i < numOfIterations; i++) {
            var minimum = determineMatrixMinimum();
            var newClusterName = mergeClusters(minimum.cluster1Name, minimum.cluster2Name);
            var subtree = appendToTree(minimum.cluster1Name, minimum.cluster2Name, newClusterName, minimum.distance / 2);
            childInstance.computeDistances(subtree);
        }

        outputData.distanceMatrix = distanceMatrixCopy;  // write-back
        outputData.newickString = formats.newickFormat.getEncoding(outputData.tree);
        return [inputData, outputData];
    }

    /**
     * Initializes structs used in the algorithm.
     */
    function initializeStructs() {
        outputData.cardinalities = {};  // needed for distance computations (for example in UPGMA)
        outputData.remainingClusterNames = ["a", "b", "c", "d", "e"];  // have to be fixed
        outputData.removedKeys = [];
        outputData.tree = [];  // hierarchical tree
    }

    /**
     * Initializes the size parameters of the clusters.
     * @param numOfIterations - The number of iterations the algorithm will do.
     */
    function initializeCardinalities(numOfIterations) {
        clusteringInstance.nameIndex = numOfIterations + 1;  // current cluster-name set-size + 1

        for (var i = 0; i < clusteringInstance.nameIndex; i++)
            outputData.cardinalities[CLUSTER_NAMES[i]] = 1;
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

        var keys = getMatrixKeys(outputData.distanceMatrix);

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
        minimum.cluster1Name = minKey[0];
        minimum.cluster2Name = minKey[1];
        minimum.distance = minValue;
        return minimum;
    }

    /**
     * Returns the remaining acceptable keys from the matrix.
     * @param distanceMatrix - The associative array
     * from which the cluster names should be determined.
     * @return {Array} - Array of cluster names.
     */
    function getMatrixKeys(distanceMatrix) {
        var keys = Object.keys(distanceMatrix);
        var remainingKeys = [];

        for (var i = 0; i < keys.length; i++) {
            if (outputData.removedKeys.indexOf(keys[i]) === -1)  // if (not contained)
                remainingKeys.push(keys[i]);
        }

        return remainingKeys;
    }

    /**
     * Computes the union of the two clusters which form the minimum.
     * @param cluster1Name {string} - The name of the first cluster.
     * @param cluster2Name {string} - The name of the second cluster.
     * @return {string} - The name of the new cluster.
     */
    function mergeClusters(cluster1Name, cluster2Name) {
        var newClusterName = createNewCluster(cluster1Name, cluster2Name)
        removeEntriesWith(cluster1Name, cluster2Name);
        return newClusterName;
    }

    /**
     * Creates the union of the two given clusters.
     * @param cluster1Name {string} - The name of the first cluster.
     * @param cluster2Name {string} - The name of the second cluster.
     * @return {string} - The name of the new cluster.
     * @see: It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function createNewCluster(cluster1Name, cluster2Name) {
        var newClusterName = getNextClusterName();

        var firstClusterCardinality = outputData.cardinalities[cluster1Name];
        var SecondClusterCardinality = outputData.cardinalities[cluster2Name];
        outputData.cardinalities[newClusterName] = firstClusterCardinality + SecondClusterCardinality;

        return newClusterName;
    }

    /**
     * Storing distance matrix entries
     * which are not allowed to use anymore for a minimum value.
     * Hint: It is not allowed to remove entries
     * from the distance matrix because it is possible that the entries
     * will be needed for later calculations (with for example UPGMA).
     * Because of this, "removed" entries are just stored
     * and no more used for calculations.
     * @param cluster1Name {string} - The name of the first cluster.
     * @param cluster2Name {string} - The name of the second cluster.
     */
    function removeEntriesWith(cluster1Name, cluster2Name) {
        var keys = Object.keys(outputData.distanceMatrix);

        for (var i = 0; i < keys.length; i++) {
            var key = keys[i].split(SYMBOLS.COMMA);

            var firstKeyPart = key[0];
            var secondKeyPart = key[1];

            if (firstKeyPart === cluster1Name || firstKeyPart === cluster2Name
                || secondKeyPart === cluster1Name || secondKeyPart === cluster2Name)
                outputData.removedKeys.push(keys[i]);
        }

        removeFromRemaining(cluster1Name);
        removeFromRemaining(cluster2Name);
    }

    /**
     * Removes cluster name from remaining names.
     * @param clusterName - The name which should be removed.
     */
    function removeFromRemaining(clusterName) {
        var index = outputData.remainingClusterNames.indexOf(clusterName);

        if (index >= 0)
            outputData.remainingClusterNames.splice(index, 1);
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
     * Appends a node with the given parameters to the hierarchical tree.
     * @param cluster1Name {string} - The name of the first cluster.
     * @param cluster2Name {string} - The name of the second cluster.
     * @param newClusterName {string} - The name of the new cluster.
     * @param distance {number} - The distance between cluster 1 and cluster 2.
     * @return
     */
    function appendToTree(cluster1Name, cluster2Name, newClusterName, distance) {
        // create node
        var node = {};
        node.leftChild = getNode(cluster1Name, distance);
        node.rightChild = getNode(cluster2Name, distance);
        node.name = newClusterName;
        node.value = distance;  // non-final evolutionary distance for edge above this node

        outputData.tree.push(node);
        return outputData.tree[outputData.tree.length-1];
    }

    /**
     * Creates a new node with the given name or uses an existing if it exists.
     * @param name {string} - The name of the node which have to be created.
     * @param value {number} - The computed distance value for the node.
     * @return {Object} - The node with the given parameters.
     */
    function getNode(name, value) {
        for (var i = 0; i < outputData.tree.length; i++) {
            if (outputData.tree[i].name === name) {
                var node = outputData.tree.splice(i, 1)[0];  // removes and returns the removed element
                node.value = value - node.value;  // computing final evolutionary distance for edge above the node
                return node;
            }
        }

        // create node
        var node = {};
        node.leftChild = undefined;
        node.rightChild = undefined;
        node.name = name;
        node.value = value;  // the value of a edge above a node

        return node;
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
