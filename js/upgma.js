/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("upgma", Upgma, computeDistances, getSuperclass);

    // instances
    var clusteringInstance;
    var upgmaInstance;

    /*---- ALGORITHM ----*/
    /**
     * Computes a clustering with the hierarchical approach
     * unweighted pair group method with arithmetic mean.
     * Hint: It is really UPGMA and not WPGMA!
     * @see
     * Sokal, Robert R. "A statistical method for evaluating systematic relationship."
     * University of Kansas science bulletin 28 (1958): 1409-1438.
     * @constructor
     * @augments Clustering
     */
    function Upgma() {
        upgmaInstance = this;

        // variables
        this.type = ALGORITHMS.UPGMA;

        // inheritance
        clusteringInstance = new bases.clustering.Clustering(this);

        this.setInput = clusteringInstance.setInput;
        this.compute = clusteringInstance.compute;
        this.getOutput = clusteringInstance.getOutput;

        this.setIO = clusteringInstance.setIO;

        // public class methods
        this.computeDistances = computeDistances;
    }

    /**
     * Computes the distance of the new cluster to the other clusters.
     * Hint: It is really UPGMA and not WPGMA!
     * Calculated distances in UPGMA are unweighted
     * with respect to the cluster-sizes. From this the "unweighted"-term results.
     * @example:
     * dist(c, k = i union j) = [|i|* dist(c, i) + |j|* dist(c, j)] / [|i|+|j|]
     * @param subtree {Object} - The subtree for the new cluster.
     */
    function computeDistances(subtree) {
        var outputData = clusteringInstance.getOutput();

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
     * Returns the superclass instance.
     * @return {Object} - Superclass instance.
     */
    function getSuperclass() {
        return clusteringInstance;
    }
}());