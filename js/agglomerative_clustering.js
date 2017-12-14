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
        AgglomerativeClustering, startAgglomerativeClustering);

    // instances
    var agglomerativeClusteringInstance;
    var clusteringInstance;

    /**
     * Function managing objects.
     */
    function startAgglomerativeClustering() {
        var clusteringInterface = new interfaces.clusteringInterface.ClusteringInterface();
        clusteringInterface.startClusteringInterface(AgglomerativeClustering, ALGORITHMS.AGGLOMERATIVE_CLUSTERING);
    }

    /*---- ALGORITHM ----*/
    /**
     * Computes an agglomerative clustering with different hierarchical approaches.
     * @constructor
     * @augments Clustering
     * @see: https://archive.org/details/cbarchive_33927_astatisticalmethodforevaluatin1902
     * and http://journals.sagepub.com/doi/abs/10.1177/001316446602600201
     *
     * Complete Linkage (Furthest Neighbour) by
     * Sørensen, Thorvald.
     * "A method of establishing groups of equal amplitude in plant sociology based on similarity
     * of species and its application to analyses of the vegetation on Danish commons."
     * Biol. Skr. 5 (1948): 1-34.
     *
     * Single Linkage (Nearest neighbour) by
     * Kazimierz Florek, Jozef Lukaszewicz, Julian Perkal, Hugo Steinhaus and Stefan Zubrzycki
     * "Sur la liaison et la division des points d'un ensemble fini."
     * Colloquium Mathematicae. Vol. 2. No. 3-4. 1951.
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
        this.computeDistances = computeDistances;

        this.getSuperclass = getSuperclass;
    }

    /**
     * Computes the distances of the new cluster to the other clusters.
     * @param subtree {Object} - The subtree for the new cluster.
     */
    function computeDistances(subtree) {
        var inputData = clusteringInstance.getInput();
        var outputData = clusteringInstance.getOutput();

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
        var cluster1Name = subtree.rightChild.name;
        var cluster2Name = subtree.leftChild.name;
        var newClusterName = subtree.name;

        var cluster1Cardinality = outputData.cardinalities[cluster1Name];
        var cluster2Cardinality = outputData.cardinalities[cluster2Name];

        var clusterNames = clusteringInstance.remainingClusterNames;

        for (var i = 0; i < clusterNames.length; i++) {
            var product1 = cluster1Cardinality * clusteringInstance.getMatrixValue(outputData.distanceMatrix, clusterNames[i], cluster1Name);
            var product2 = cluster2Cardinality * clusteringInstance.getMatrixValue(outputData.distanceMatrix, clusterNames[i], cluster2Name);

            var dividendSum = product1 + product2;
            var divisorSum = cluster1Cardinality + cluster2Cardinality;

            var quotient = dividendSum / divisorSum;

            outputData.distanceMatrix[[clusterNames[i], newClusterName]] = quotient;  // hint: do not change order of arguments
        }
    }

    /**
     * Computes the distance of the new cluster to the other clusters.
     * Hint: It is really WPGMA and not UPGMA!
     * Calculated distances in UPGMA are unweighted
     * with respect to the cluster-sizes. From this the "unweighted"-term results.
     * @example:
     * dist(c, k = i union j) = [dist(c, i) + dist(c, j)] / 2
     * @param subtree {Object} - The subtree for the new cluster.
     * @param outputData {Object} - Contains all output data.
     */
    function computeWpgmaDistance(subtree, outputData) {
        // retrieve values
        var cluster1Name = subtree.rightChild.name;
        var cluster2Name = subtree.leftChild.name;
        var newClusterName = subtree.name;

        var clusterNames = clusteringInstance.remainingClusterNames;

        for (var i = 0; i < clusterNames.length; i++) {
            var summand1 = clusteringInstance.getMatrixValue(outputData.distanceMatrix, clusterNames[i], cluster1Name);
            var summand2 = clusteringInstance.getMatrixValue(outputData.distanceMatrix, clusterNames[i], cluster2Name);

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