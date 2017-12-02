/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("bases.clustering", Clustering);

    // instances
    var childInstance;
    var clusteringInstance;
    var newickEncoderInstance;
    var parserInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Contains functions to compute a hierarchical clustering (at the moment only hierarchical agglomerative clustering).
     * It is used by agglomerative clustering algorithms as superclass
     * to avoid code duplicates.
     * @constructor
     */
    function Clustering(child) {
        clusteringInstance = this;
        newickEncoderInstance = new formats.newickEncoder.NewickEncoder();
        parserInstance = new formats.csvParser.CsvParser();

        // variables
        this.nameIndex = 0;  // only really needed, if getNextClusterName() function is used
        this.remainingClusterNames = [];
        this.removedKeys = [];
        this.treeParts = [];  // storing the current incurred tree parts

        // inheritance
        childInstance = child;

        // public class methods
        this.setInput = setInput;
        this.getInput = getInput;
        this.compute = compute;
        this.getMatrixKeys = getMatrixKeys;
        this.getOutput = getOutput;

        this.setIO = setIO;
    }

    /**
     * Sets the algorithm input for an appropriate algorithm
     * which is using the inputViewmodel properties in its computations.
     * @param inputViewmodel {Object} - The InputViewmodel of an appropriate algorithm.
     */
    function setInput(inputViewmodel) {
        var noError = inputViewmodel.errorInput() === SYMBOLS.EMPTY;

        inputData.clusteringSubalgorithm = inputViewmodel.selectedApproach()[0];  // because it's array from which we choose

        outputData.distanceMatrix = getDistanceMatrix(noError, inputViewmodel.csvTable());
        inputData.numOfStartClusters = getNumOfStartClusters(noError, inputViewmodel.csvTable());

        inputData.initialNamingIndex = inputData.numOfStartClusters;
        outputData.distanceMatrixLength = inputData.numOfStartClusters;
        outputData.clusterNames = getColumnNames(noError, outputData.distanceMatrix);
    }

    /**
     * Returns the distance matrix object created by a CSV-parser.
     * @param isCorrectData {boolean} - If false, it returns an empty object and else the distance matrix object.
     * @param csvData {string} - The csv string which has to be converted into a cluster algorithm distance matrix.
     * @return {Object} - The distance matrix.
     */
    function getDistanceMatrix(isCorrectData, csvData) {
        var distanceMatrix = {};

        if (isCorrectData)
            distanceMatrix = parserInstance.getMatrix(csvData, getClusterNames);

        return distanceMatrix;
    }

    /**
     * Returns names for clusters associated with the data.
     * Hint: After all characters are depleted,
     * a number is concatenated to the character
     * to make this function generic.
     * @param number {number} - The number of names you want create.
     * @example:
     * CLUSTER NAMES:
     * a, b, c, ..., z,         FIRST EPISODE
     * a2, b2, c2, ..., z2,     SECOND EPISODE
     * a3, b3, ...              THIRD ...
     * @return {Array} - The cluster names.
     */
    function getClusterNames(number) {
        var clusterNames = [];
        var currentEpisode = 1;

        // for every pairwise distance we need a symbol
        for (var i = 0; i < number; i++) {
            if (i < CLUSTER_NAMES.length)
                clusterNames.push(CLUSTER_NAMES[i]);  // add a, b, c, ..., z

            if (i >= CLUSTER_NAMES.length && i % CLUSTER_NAMES.length === 0)  // out of characters
                currentEpisode++;  // new episode

            // out of characters -> a2, b2, c2, ..., z2, a3, b3, ...
            if (i >= CLUSTER_NAMES.length)
                clusterNames.push(CLUSTER_NAMES[i % CLUSTER_NAMES.length] + SYMBOLS.EMPTY + currentEpisode);
        }

        return clusterNames;
    }

    /**
     * Returns the number of start clusters given the CSV data.
     * @param isCorrectData {boolean} - If false, it returns an empty object and else the distance matrix object.
     * @param csvData {string} - The csv string which has to be converted into a cluster algorithm distance matrix.
     * @return {Object} - The distance matrix.
     */
    function getNumOfStartClusters(isCorrectData, csvData) {
        var numOfStartClusters = 0;

        if (isCorrectData)
            numOfStartClusters = parserInstance.getNumOfTableLines(csvData);

        return numOfStartClusters;
    }

    /**
     * Returns the names from a distance matrix.
     * @param isCorrectData {boolean} - If false, it returns an empty object and else the distance matrix object.
     * @param distanceMatrix {Object} - The distance matrix object from which the column names are returned.
     * @return {Array} - The name of the columns.
     */
    function getColumnNames(isCorrectData, distanceMatrix) {
        var names = [CLUSTER_NAMES[0]];

        if (isCorrectData) {
            var namePairs = Object.keys(distanceMatrix);

            for (var i = 0; i < namePairs.length; i++) {
                var namesOfPair = namePairs[i].split(SYMBOLS.COMMA);
                var name1 = namesOfPair[0];
                var name2 = namesOfPair[1];

                if (names.indexOf(name1) === -1)
                    names.push(name1);

                if (names.indexOf(name2) === -1)
                    names.push(name2);
            }
        }

        return names;
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
        initializeCardinalities();

        for (var i = 0; i < numOfIterations; i++) {
            var minimum = determineMatrixMinimum();
            var newClusterName = mergeClusters(minimum.cluster1Name, minimum.cluster2Name);
            var subtree = appendToTree(minimum.cluster1Name, minimum.cluster2Name, newClusterName, minimum.distance / 2);
            computeDistances(subtree, i, numOfIterations);
        }

        getMatrixKeys(outputData.distanceMatrix);  // only for visualization called again, to store also the last matrix
        outputData.distanceMatrix = distanceMatrixCopy;  // write-back
        outputData.newickString = newickEncoderInstance.getEncoding(outputData.treeBranches[outputData.treeBranches.length-1]);
        return [inputData, outputData];
    }

    /**
     * Initializes structs used in the algorithm.
     */
    function initializeStructs() {
        debugger;
        clusteringInstance.remainingClusterNames = outputData.clusterNames.slice();  // shallow copy (because they won't be changed)
        clusteringInstance.removedKeys = [];
        clusteringInstance.treeParts = [];

        outputData.cardinalities = {};  // needed for distance computations (for example in UPGMA)

        // stores the created tree branches in the order
        // in which they were created to avoid a repeated tree traversal during the progressive alignment
        // and for possible step by step visualizations of tree-growthment with libraries
        outputData.treeBranches = [];

        outputData.allClusterNames = outputData.clusterNames.slice();
        outputData.remainingClusters = [jQuery.extend(true, [], clusteringInstance.remainingClusterNames)];
        outputData.distanceMatrices = [];
        outputData.keys = [];
        outputData.minimums = [];
    }

    /**
     * Initializes the size-parameters of the clusters.
     */
    function initializeCardinalities() {
        clusteringInstance.nameIndex = inputData.initialNamingIndex;  // do not change that!

        for (var i = 0; i < clusteringInstance.nameIndex; i++)
            outputData.cardinalities[outputData.clusterNames[i]] = 1;
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

        outputData.minimums.push(Math.round(minValue*10)/10);  // only for visualization

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
            if (clusteringInstance.removedKeys.indexOf(keys[i]) === -1)  // if (not contained)
                remainingKeys.push(keys[i]);
        }

        outputData.distanceMatrices.push(jQuery.extend(true, {}, outputData.distanceMatrix));  // only for visualization the matrix of each round is stored
        outputData.keys.push(remainingKeys);

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
        outputData.allClusterNames.push(newClusterName);
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
        var newClusterName = cluster1Name + cluster2Name;  // getNextClusterName();  // alternative name generation

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
                clusteringInstance.removedKeys.push(keys[i]);
        }

        removeFromRemaining(cluster1Name);
        removeFromRemaining(cluster2Name);
    }

    /**
     * Removes cluster name from remaining names.
     * @param clusterName - The name which should be removed.
     */
    function removeFromRemaining(clusterName) {
        var index = clusteringInstance.remainingClusterNames.indexOf(clusterName);

        if (index >= 0)
            clusteringInstance.remainingClusterNames.splice(index, 1);
    }

    /** ALTERNATIVE from lecture WS 2016/2017.
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
    /*
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
    */

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

        outputData.treeBranches.push(node);
        clusteringInstance.treeParts.push(node);
        return clusteringInstance.treeParts[clusteringInstance.treeParts.length-1];
    }

    /**
     * Creates a new node with the given name or uses an existing if it exists.
     * @param name {string} - The name of the node which have to be created.
     * @param value {number} - The computed distance value for the node.
     * @return {Object} - The node with the given parameters.
     */
    function getNode(name, value) {
        for (var i = 0; i < clusteringInstance.treeParts.length; i++) {
            if (clusteringInstance.treeParts[i].name === name) {
                var node = clusteringInstance.treeParts.splice(i, 1)[0];  // removes and returns the removed element
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
     * Computes the distance of the new cluster to the other clusters.
     * Hint: It is really UPGMA and not WPGMA!
     * Calculated distances in UPGMA are unweighted
     * with respect to the cluster-sizes. From this the "unweighted"-term results.
     * @example:
     * dist(c, k = i union j) = [|i|* dist(c, i) + |j|* dist(c, j)] / [|i|+|j|]
     * @param subtree {Object} - The subtree for the new cluster.
     * @param iteration {number} - The subtree for the new cluster.
     * @param subtree {maxNumIterations} - The subtree for the new cluster.
     */
    function computeDistances(subtree, iteration, maxNumIterations) {
        childInstance.computeDistances(subtree);
        clusteringInstance.remainingClusterNames.push(subtree.name);

        if (iteration === maxNumIterations-1)
            subtree.value = 0;

        outputData.remainingClusters.push(jQuery.extend(true, [], clusteringInstance.remainingClusterNames));  // for visualization
    }

    /**
     * Returns all algorithm output.
     * @return {Object} - Contains all output of the algorithm.
     */
    function getInput() {
        return inputData;
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
