/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("neighbourJoining", NeighbourJoining);

    // instances
    var neighbourJoiningInstance;
    var clusteringInstance;

    // shared variables
    var inputData = {};  // stores the input of the algorithm
    var outputData = {};  // stores the output of the algorithm

    /**
     * Computes an agglomerative clustering with a hierarchical approach.
     * It is executed by the class for agglomerative clustering.
     * Hint: This algorithm earns an own class, because it is very different
     * to the common agglomerative clustering algorithms like Single Linkage (Nearest Neighbour),
     * UPGMA (Average Linkage), WPGMA (Simple Average), Complete Linkage (Furthest Neighbour)
     * and many more.
     * @constructor
     * @augments Clustering
     *
     * Saitou, Naruya, and Masatoshi Nei.
     * "The neighbor-joining method: a new method for reconstructing phylogenetic trees."
     * Molecular biology and evolution 4.4 (1987): 406-425.
     *                  and
     * Studier, James A., and Karl J. Keppler.
     * "A note on the neighbor-joining algorithm of Saitou and Nei."
     * Molecular biology and evolution 5.6 (1988): 729-731.
     */
    function NeighbourJoining() {
        neighbourJoiningInstance = this;

        // variables
        this.type = ALGORITHMS.NEIGHBOUR_JOINING;

        // inheritance
        clusteringInstance = new bases.clustering.Clustering(this);

        this.setInput = setInput;
        this.compute = compute;
        this.getOutput = getOutput;

        this.setIO = setIO;

        // public class methods
        this.computeDistances = computeDistances;
        
        this.getSuperclass = getSuperclass;
    }

    /**
     * Sets the algorithm input for an appropriate algorithm
     * which is using the inputViewmodel properties in its computations.
     * @param inputViewmodel {Object} - The InputViewmodel of an appropriate algorithm.
     */
    function setInput(inputViewmodel) {
        clusteringInstance.setIO(inputData, outputData);
        clusteringInstance.setInput(inputViewmodel);
    }

    /**
     * Starts the computation.
     * Hint: Because the distance matrix is changing during the following procedure
     * and it has to be visualized later on, a copy is made which is written back.
     */
    function compute() {
        var distanceMatrixCopy = jQuery.extend(true, {}, outputData.distanceMatrix);  // deep copy
        var numOfIterations = inputData.numOfStartClusters - 2;  // because of the formula with N-2 (special case)

        clusteringInstance.initializeStructs();
        clusteringInstance.initializeCardinalities();

        for (var i = 0; i < numOfIterations; i++) {
            // Step 1 from Unit-Test.
            var neighbourJoiningMatrix
                = getNeighbourJoiningMatrix(outputData.distanceMatrix, getTotalDistances(outputData.distanceMatrix));

            // Step 2 from Unit-Test
            var minimum = clusteringInstance.determineMatrixMinimum(neighbourJoiningMatrix);
            var newClusterName = clusteringInstance.mergeClusters(minimum.cluster1Name, minimum.cluster2Name);

            // Step 3 from Unit-Test
            var subtree  // Step 3.1
                = clusteringInstance.computeBranchLengths(minimum.cluster1Name, minimum.cluster2Name, newClusterName);
            computeDistances(subtree, i, numOfIterations);  // Step 3.2
        }

        // step to create correct newick string

        clusteringInstance.getMatrixKeys(outputData.distanceMatrix);  // only for visualization called again, to store also the last matrix
        outputData.distanceMatrix = distanceMatrixCopy;  // write-back

        // outputData.newickString = newickEncoderInstance.getEncoding(outputData.treeBranches[outputData.treeBranches.length-1]);
        return [inputData, outputData];
    }

    /**
     * Computes the total distances (sum of values in every line)
     * in the distance matrix and returns them.
     * Hint: Step 1.1 from Unit-Test.
     * @param distanceMatrix {Object} - The distance matrix which is used to compute total distances.
     * @return {Array} - The total distance for each line.
     */
    function getTotalDistances(distanceMatrix) {
        // execute getMatrixAsTable

        // after returned, search minima with determineMatrixMinimum
        return undefined;
    }

    /**
     * Returns the given distance matrix into an two-dimensional array.
     * @param distanceMatrix {Object} - The distance-value which are converted into a two-dimensional array.
     * @return {Array} - The matrix as an array.
     */
    function getMatrixAsTable(distanceMatrix) {
        return undefined;
    }

    /**
     * Computes a neighbour-joining matrix out of a distance matrix and the total distances.
     * Hint: Step 1.2 from Unit-Test.
     * The
     */
    function getNeighbourJoiningMatrix(distanceMatrix, totalDistances) {
        // after returned, search minima with determineMatrixMinimum
        return undefined;
    }

    /**
     * Computes the branch length between to clusters and a new cluster.
     * Hint: Step 3.1 from Unit-Test.
     * @param cluster1Name {string} - The name of the first cluster.
     * @param cluster2Name {string} - The name of the second cluster.
     * @param newClusterName {string} - The name of the new cluster.
     */
    function computeBranchLengths(cluster1Name, cluster2Name, newClusterName) {
        // compute the two distances
        
        // execute appendToTree
    }

    /**
     * Computes the distances of the new cluster to the other clusters.
     * Hint: Step 3.2 from Unit-Test.
     * @param subtree {Object} - The subtree for the new cluster.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     */
    function computeDistances(subtree, inputData, outputData) {
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
        return clusteringInstance;
    }
}());