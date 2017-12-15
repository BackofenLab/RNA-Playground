/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods ("getMatrixAsTable" is set static because creation of a full instance to get just the table would be too inefficient)
    namespace("bases.clustering", Clustering, getMatrixAsTable);

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
     * It is used by hierarchical clustering algorithms as superclass
     * to avoid code duplicates.
     * @param child {Object} - The child algorithm which inherits from this class.
     * @constructor
     */
    function Clustering(child) {
        clusteringInstance = this;
        newickEncoderInstance = new formats.newickEncoder.NewickEncoder();
        parserInstance = new formats.csvParser.CsvParser();

        // variables
        this.nameIndex = 0;  // only really needed, if getNextClusterName() function is used
        this.remainingClusterNames = [];
        this.currentClusterNames = [];
        this.lastCurrentClusterNames = [];
        this.removedKeys = [];
        this.treeParts = [];  // storing the current incurred tree parts

        // inheritance
        childInstance = child;

        // public class methods
        this.setInput = setInput;
        this.getInput = getInput;
        this.compute = compute;
        this.initializeStructs = initializeStructs;
        this.determineMatrixMinimum = determineMatrixMinimum;
        this.mergeClusters = mergeClusters;
        this.appendToTree = appendToTree;
        this.getMatrixKeys = getMatrixKeys;
        this.getMatrixValue = getMatrixValue;
        this.getMatrixAsTable = getMatrixAsTable;  // yes, it should be possible to execute this function from class
        this.getPositionByName = getPositionByName;
        this.getOutput = getOutput;
        this.setIO = setIO;
        this.getNewickEncoder = getNewickEncoder;
    }

    /**
     * Sets the algorithm input for an appropriate algorithm
     * which is using the inputViewmodel properties in its computations.
     * @param inputViewmodel {Object} - The InputViewmodel of an appropriate algorithm.
     */
    function setInput(inputViewmodel) {
        var noError = inputViewmodel.errorInput() === SYMBOLS.EMPTY;

        inputData.clusteringSubalgorithm = inputViewmodel.selectedApproach()[0];  // because it's array from which we choose

        if (inputData.clusteringSubalgorithm !== CLUSTERING_ALGORITHMS.NEIGHBOUR_JOINING)  // to work with correct child instance
            childInstance = new agglomerativeClustering.AgglomerativeClustering();

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
        if (inputData.clusteringSubalgorithm === CLUSTERING_ALGORITHMS.NEIGHBOUR_JOINING)
            return computeNeighbourJoining();

        var distanceMatrixCopy = jQuery.extend(true, {}, outputData.distanceMatrix);  // deep copy
        var numOfIterations = inputData.numOfStartClusters - 1;  // always lower by one in fundamental hierarchical clustering algorithms

        initializeStructs();
        if (inputData.clusteringSubalgorithm === CLUSTERING_ALGORITHMS.COMPLETE_LINKAGE
            || inputData.clusteringSubalgorithm === CLUSTERING_ALGORITHMS.SINGLE_LINKAGE
            || inputData.clusteringSubalgorithm === CLUSTERING_ALGORITHMS.UPGMA)
            initializeClusters();

        for (var i = 0; i < numOfIterations; i++) {
            var minimum = determineMatrixMinimum(outputData.distanceMatrix);
            var newClusterName = mergeClusters(minimum.cluster1Name, minimum.cluster2Name);
            var subtree = appendToTree(minimum.cluster1Name, minimum.cluster2Name, newClusterName, minimum.distance / 2, minimum.distance / 2);
            computeDistances(subtree, i, numOfIterations);
        }

        getMatrixKeys(outputData.distanceMatrix, true);  // only for visualization called again, to store also the last matrix
        outputData.distanceMatrix = distanceMatrixCopy;  // write-back
        outputData.newickString = newickEncoderInstance.getEncoding(outputData.treeBranches[outputData.treeBranches.length - 1]);
        return [inputData, outputData];
    }

    /**
     * Computes Neighbour-Joining output.
     * @return {[input,output]} - The input data and output data.
     */
    function computeNeighbourJoining() {
        var algorithm = new neighbourJoining.NeighbourJoining();
        algorithm.setIO(inputData, outputData);
        return algorithm.compute();
    }

    /**
     * Initializes structs used in the algorithm.
     */
    function initializeStructs() {
        clusteringInstance.remainingClusterNames = outputData.clusterNames.slice();  // shallow copy (because they won't be changed)
        clusteringInstance.currentClusterNames = outputData.clusterNames.slice();  // stores the current clusters of the distance matrix
        clusteringInstance.lastCurrentClusterNames = outputData.clusterNames.slice();
        clusteringInstance.removedKeys = [];
        clusteringInstance.treeParts = [];

        outputData.cardinalities = {};  // needed for distance computations (for example in UPGMA)
        outputData.clusterMembers = {};  // needed for distance computation (for example in Single Linkage)

        // stores the created tree branches in the order
        // in which they were created to avoid a repeated tree traversal during the progressive alignment
        // and for possible step by step visualizations of tree-growthment with libraries
        outputData.treeBranches = [];

        outputData.allClusterNames = outputData.clusterNames.slice();
        outputData.remainingClusters = [jQuery.extend(true, [], clusteringInstance.remainingClusterNames)];
        outputData.distanceMatrices = [];
        outputData.keys = [];
        outputData.minimums = [];

        // neighbour-joining algorithm
        outputData.neighbourJoiningMatrices = [];
        outputData.matrixTables = []; // to avoid a recomputation
        outputData.totalDistancesPerRound = [];
        outputData.lastComputedDistance = 0;  // stores the last computed distance (needed to get last distance)
        outputData.lastComputedClusterNames = [];  // stores the name of the last cluster
    }

    /**
     * Initializes the size-parameters of the clusters.
     */
    function initializeClusters() {
        clusteringInstance.nameIndex = inputData.initialNamingIndex;  // do not change that!

        for (var i = 0; i < clusteringInstance.nameIndex; i++) {
            outputData.cardinalities[outputData.clusterNames[i]] = 1;
            outputData.clusterMembers[outputData.clusterNames[i]] = outputData.clusterNames[i];
        }
    }

    /**
     * Determines the entry with the lowest distance
     * in the remaining entries.
     * @param distanceMatrix {Object} - The distance matrix in which it is searched for the minimum.
     * @return {Object} - The matrix-entry which contains the minimum value.
     * @see: It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function determineMatrixMinimum(distanceMatrix) {
        var minKey = SYMBOLS.EMPTY;
        var minValue = Number.POSITIVE_INFINITY;

        var keys = getMatrixKeys(distanceMatrix, true);

        // searching for minimum by going over all keys
        for (var i = 0; i < keys.length; i++) {
            var key = keys[i].split(SYMBOLS.COMMA);  // string to array: Javascript stores keys as strings
            var value = distanceMatrix[key];

            if (value < minValue) {  // hint: diagonals were not computed because they are zero
                minKey = key;  // name: [a,b] with a,b cluster names
                minValue = value;  // distance
            }
        }

        outputData.minimums.push(Math.round(minValue * 10) / 10);  // only for visualization

        // create structure for better understandable access
        var minimum = {};
        minimum.cluster1Name = minKey[0];
        minimum.cluster2Name = minKey[1];
        minimum.distance = minValue;
        return minimum;
    }

    /**
     * Returns the remaining acceptable keys from the matrix.
     * @param distanceMatrix {Object} - The associative array
     * @param storeMatrix {boolean} - Tells if the matrix should be stored for visualization purposes.
     * from which the cluster names should be determined.
     * @return {Array} - Array of cluster names.
     */
    function getMatrixKeys(distanceMatrix, storeMatrix) {
        var keys = Object.keys(distanceMatrix);
        var remainingKeys = [];

        for (var i = 0; i < keys.length; i++) {
            if (clusteringInstance.removedKeys.indexOf(keys[i]) === -1)  // if (not contained)
                remainingKeys.push(keys[i]);
        }

        if (storeMatrix) {
            outputData.distanceMatrices.push(jQuery.extend(true, {}, outputData.distanceMatrix));  // only for visualization the matrix of each round is stored
            outputData.keys.push(remainingKeys);

            if (inputData.clusteringSubalgorithm === CLUSTERING_ALGORITHMS.NEIGHBOUR_JOINING)
                outputData.neighbourJoiningMatrices.push(jQuery.extend(true, {}, distanceMatrix));
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
        var newClusterName = createNewCluster(cluster1Name, cluster2Name);
        clusteringInstance.lastCurrentClusterNames = clusteringInstance.currentClusterNames.slice();  // shallow copy
        removeEntriesWith(cluster1Name, cluster2Name);
        clusteringInstance.currentClusterNames.push(newClusterName);
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

        if (inputData.clusteringSubalgorithm === CLUSTERING_ALGORITHMS.UPGMA) {
            var firstClusterCardinality = outputData.cardinalities[cluster1Name];
            var SecondClusterCardinality = outputData.cardinalities[cluster2Name];
            outputData.cardinalities[newClusterName] = firstClusterCardinality + SecondClusterCardinality;
        } else if (inputData.clusteringSubalgorithm === CLUSTERING_ALGORITHMS.SINGLE_LINKAGE
            || inputData.clusteringSubalgorithm === CLUSTERING_ALGORITHMS.COMPLETE_LINKAGE) {  // store names of members
            outputData.clusterMembers[newClusterName]
                = outputData.clusterMembers[cluster1Name].concat(outputData.clusterMembers[cluster2Name]);
        }

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
     * @param clusterName {string} - The name which should be removed.
     */
    function removeFromRemaining(clusterName) {
        var index = clusteringInstance.remainingClusterNames.indexOf(clusterName);

        if (index >= 0)
            clusteringInstance.remainingClusterNames.splice(index, 1);

        index = clusteringInstance.currentClusterNames.indexOf(clusterName);

        if (index >= 0)
            clusteringInstance.currentClusterNames.splice(index, 1);
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
     * @param distance1 {number} - The distance value for edge above cluster 1.
     * @param distance2 {number} - The distance value for edge above cluster 2.
     * @return {Object} - The new tree part.
     */
    function appendToTree(cluster1Name, cluster2Name, newClusterName, distance1, distance2) {
        // create node
        var node = {};
        node.rightChild = getNode(cluster1Name, distance1);
        node.leftChild = getNode(cluster2Name, distance2);
        node.name = newClusterName;
        node.value = distance2;  // non-final evolutionary distance for edge above this node

        outputData.treeBranches.push(node);
        clusteringInstance.treeParts.push(node);
        return clusteringInstance.treeParts[clusteringInstance.treeParts.length - 1];
    }

    /**
     * Creates a new node with the given name or uses an existing if it exists.
     * @param name {string} - The name of the node which have to be created.
     * @param value {number} - The computed distance value for the node.
     * @return {Object} - The node with the given parameters.
     */
    function getNode(name, value) {
        var node = {};

        // search for an existing node
        for (var i = 0; i < clusteringInstance.treeParts.length; i++) {
            if (clusteringInstance.treeParts[i].name === name) {
                node = clusteringInstance.treeParts.splice(i, 1)[0];  // removes and returns the removed element

                if (inputData.clusteringSubalgorithm === CLUSTERING_ALGORITHMS.NEIGHBOUR_JOINING)
                    node.value = value;  // do not subtract in neighbour-joining, because already final value
                else
                    node.value = value - node.value;  // computing final evolutionary distance for edge above the node
                return node;
            }
        }

        // create node
        node.rightChild = undefined;
        node.leftChild = undefined;
        node.name = name;
        node.value = value;  // the value of a edge above a node

        return node;
    }

    /**
     * Computes the distance of the new cluster to the other clusters.
     * @example:
     * dist(c, k = i union j) = [|i|* dist(c, i) + |j|* dist(c, j)] / [|i|+|j|]
     * @param subtree {Object} - The subtree for the new cluster.
     * @param iteration {number} - The subtree for the new cluster.
     * @param maxNumIterations {number} - The maximum number of iterations.
     */
    function computeDistances(subtree, iteration, maxNumIterations) {
        childInstance.computeDistances(subtree);
        clusteringInstance.remainingClusterNames.push(subtree.name);

        if (iteration === maxNumIterations - 1)
            subtree.value = 0;

        outputData.remainingClusters.push(jQuery.extend(true, [], clusteringInstance.remainingClusterNames));  // for visualization
    }

    /**
     * Returns the distance matrix value from the given entry.
     * Hint: Only one half of the matrix is filled.
     * But the other half is just a mirrored version.
     * This is why this function is needed.
     * @param distanceMatrix {Array} - The array from which you want the values.
     * @param cluster1Name {string} - The name of the first cluster.
     * @param cluster2Name {string} - The name of the second cluster.
     * @return {number} - The value from an entry.
     */
    function getMatrixValue(distanceMatrix, cluster1Name, cluster2Name) {
        var value1 = distanceMatrix[[cluster1Name, cluster2Name]];
        var value2 = distanceMatrix[[cluster2Name, cluster1Name]];

        if (isNaN(value1))
            return value2;
        return value1;
    }

    /**
     * Converts the given distance matrix (associative array) into an two-dimensional array.
     * Hint: "Associative arrays" do not have a defined order (browser-dependant).
     * @param distanceMatrix {Object} - The distance-value which are converted into a two-dimensional array.
     * @param distanceMatrixLength {number} - The number of clusters in the distance matrix.
     * @param remainingClusterNames {Array} - The existent cluster names in the matrix.
     * @param matrixKeys {Array} - The keys (tuples) from the distance matrix (associative array).
     * @param fillBoth {boolean} - Tells if both halves of the matrix should be filled or not. If false only upper half is filled.
     * @return {Array} - The matrix as an array.
     */
    function getMatrixAsTable(distanceMatrix, distanceMatrixLength, remainingClusterNames, matrixKeys, fillBoth) {
        var matrix = createMatrix(distanceMatrixLength);

        if (matrixKeys === undefined)
            matrixKeys = Object.keys(distanceMatrix);  // argument possibilities {a,b}, {a,c}, ...

        // fill diagonals with zero
        for (var i = 0; i < matrix.length; i++) {
            for (var j = 0; j < matrix.length; j++) {
                if (i === j)
                    matrix[i][j] = 0;
            }
        }

        // fill right upper half and left lower half
        for (var j = 0; j < matrixKeys.length; j++) {
            var key = matrixKeys[j].split(SYMBOLS.COMMA);
            var cluster1Position = getPositionByName(key[0], remainingClusterNames);
            var cluster2Position = getPositionByName(key[1], remainingClusterNames);
            var value = distanceMatrix[key];

            matrix[cluster1Position][cluster2Position] = value;  // upper right half
            if (fillBoth)
                matrix[cluster2Position][cluster1Position] = value;  // lower left half
        }

        return matrix;
    }

    /**
     * Returns for a cluster-name, its position in the distance matrix.
     * @param clusterName {string} - The name of the cluster.
     * @param remainingClusterNames {Array} - The remaining cluster names after execution of UPGMA.
     * @return {number} - The position of a name in the remaining clusters.
     */
    function getPositionByName(clusterName, remainingClusterNames) {
        var position = -1;

        for (var i = 0; i < remainingClusterNames.length; i++) {
            if (clusterName === remainingClusterNames[i]) {
                position = i;
                break;
            }
        }

        return position;
    }

    /**
     * Creates a matrix with the given size.
     * @param size {number} - The width and height of the matrix.
     */
    function createMatrix(size) {
        var matrix = new Array(size);

        for (var i = 0; i < size; i++) {
            matrix[i] = [];
        }

        return matrix;
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

    /**
     * Returns a Newick format encoder.
     * @return {Object} - The Newick encoder.
     */
    function getNewickEncoder() {
        return newickEncoderInstance;
    }
}());
