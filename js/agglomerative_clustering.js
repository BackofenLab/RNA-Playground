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
    if (loaded === ALGORITHMS.AGGLOMERATIVE_CLUSTERING) {  // to avoid self execution on a script import
        agglomerativeClustering.startAgglomerativeClustering();
        loaded = ALGORITHMS.NONE;
    }
});

(function () {  // namespace
    // public methods
    namespace("agglomerativeClustering",
        AgglomerativeClustering, startAgglomerativeClustering, computeDistances, getSuperclass);

    // instances
    var agglomerativeClusteringInstance;
    var clusteringInstance;

    /**
     * Function managing objects.
     */
    function startAgglomerativeClustering() {
        var clusteringInterface = new interfaces.clusteringInterface.ClusteringInterface();
        clusteringInterface.startClustering(AgglomerativeClustering, ALGORITHMS.AGGLOMERATIVE_CLUSTERING);
    }

    /*---- ALGORITHM ----*/
    /**
     * Computes an agglomerative clustering with different hierarchical approaches.
     * @constructor
     * @augments Clustering
     *
     * Neighbour Joining (single linkage) by
     * Saitou, Naruya, and Masatoshi Nei.
     * "The neighbor-joining method: a new method for reconstructing phylogenetic trees."
     * Molecular biology and evolution 4.4 (1987): 406-425.
     *
     * UPGMA (average linkage) by
     * Sokal, Robert R.
     * "A statistical method for evaluating systematic relationship."
     * University of Kansas science bulletin 28 (1958): 1409-1438.
     *
     * WPGMA (simple average) by
     * McQuitty, Louis L.
     * "Single and multiple hierarchical classification by reciprocal pairs and rank order types."
     * Educational and Psychological Measurement 26.2 (1966): 253-265.
     */
    function AgglomerativeClustering() {
        agglomerativeClusteringInstance = this;

        // variables
        this.type = ALGORITHMS.AGGLOMERATIVE_CLUSTERING;

        // inheritance
        clusteringInstance = new bases.clustering.Clustering(this);

        this.setInput = clusteringInstance.setInput;
        this.compute = clusteringInstance.compute;
        this.getOutput = clusteringInstance.getOutput;

        this.setIO = clusteringInstance.setIO;

        // public class methods
        this.getSuperclass = getSuperclass;

        this.computeDistances = computeDistances;
    }

    /**
     * Computes the distance of the new cluster to the other clusters.
     * @param subtree {Object} - The subtree for the new cluster.
     */
    function computeDistances(subtree) {
        var outputData = clusteringInstance.getOutput();
        var inputData = clusteringInstance.getInput();

        if (inputData.clusteringSubalgorithm === CLUSTERING_ALGORITHMS.UPGMA)
            computeUpgmaDistance(subtree, outputData);
        else if (inputData.clusteringSubalgorithm === CLUSTERING_ALGORITHMS.WPGMA)
            computeWpgmaDistance(subtree, outputData);
    }

    /**
     * Computes the distance of the new cluster to the other clusters.
     * Hint: It is really UPGMA and not WPGMA!
     * Calculated distances in UPGMA are unweighted
     * with respect to the cluster-sizes. From this the "unweighted"-term results.
     * @example:
     * dist(c, k = i union j) = [|i|* dist(c, i) + |j|* dist(c, j)] / [|i|+|j|]
     * @param subtree {Object} - The subtree for the new cluster.
     * @param outputData {Object} - Contains all output data.
     */
    function computeUpgmaDistance(subtree, outputData) {
        // retrieve values
        var cluster1Name = subtree.leftChild.name;
        var cluster2Name = subtree.rightChild.name;
        var newClusterName = subtree.name;

        var cluster1Cardinality = outputData.cardinalities[cluster1Name];
        var cluster2Cardinality = outputData.cardinalities[cluster2Name];

        var clusterNames = clusteringInstance.remainingClusterNames;

        for (var i = 0; i < clusterNames.length; i++) {
            var product1 = cluster1Cardinality * getMatrixValue(outputData.distanceMatrix, clusterNames[i], cluster1Name);
            var product2 = cluster2Cardinality * getMatrixValue(outputData.distanceMatrix, clusterNames[i], cluster2Name);

            var dividendSum = product1 + product2;
            var divisorSum = cluster1Cardinality + cluster2Cardinality;

            var quotient = dividendSum / divisorSum;

            outputData.distanceMatrix[[clusterNames[i], newClusterName]] = quotient;  // hint: do not change order of arguments
        }
    }

    /**
     * Returns the distance matrix value from the given entry.
     * Hint: Only one half of the matrix is filled.
     * But the other half is just a mirrored version.
     * This is why this function is needed.
     * @param distanceMatrix {Array} - The array from which you want the values.
     * @param cluster1Name {string} - The name of the first cluster.
     * @param cluster2Name {string} - The name of the second cluster.
     * @return {number} - The value
     */
    function getMatrixValue(distanceMatrix, cluster1Name, cluster2Name) {
        var value1 = distanceMatrix[[cluster1Name, cluster2Name]];
        var value2 = distanceMatrix[[cluster2Name, cluster1Name]];

        if(isNaN(value1))
            return value2;
        return value1;
    }

    /**
     * Computes the distance of the new cluster to the other clusters.
     * Hint: It is really WPGMA and not UPGMA!
     * Calculated distances in UPGMA are unweighted
     * with respect to the cluster-sizes. From this the "unweighted"-term results.
     * @example:
     * dist(c, k = i union j) = dist(c, i) + dist(c, j) / 2
     * @param subtree {Object} - The subtree for the new cluster.
     * @param outputData {Object} - Contains all output data.
     */
    function computeWpgmaDistance(subtree, outputData) {
        // retrieve values
        var cluster1Name = subtree.leftChild.name;
        var cluster2Name = subtree.rightChild.name;
        var newClusterName = subtree.name;

        var clusterNames = clusteringInstance.remainingClusterNames;

        for (var i = 0; i < clusterNames.length; i++) {
            var summand1 = getMatrixValue(outputData.distanceMatrix, clusterNames[i], cluster1Name);
            var summand2 = getMatrixValue(outputData.distanceMatrix, clusterNames[i], cluster2Name);

            var dividendSum = summand1 + summand2;
            var quotient = dividendSum / 2;

            outputData.distanceMatrix[[clusterNames[i], newClusterName]] = quotient;  // hint: do not change order of arguments
        }
    }

    /**
     * Returns the superclass instance.
     * @return {Object} - Superclass instance.
     */
    function getSuperclass() {
        return clusteringInstance;
    }
}());