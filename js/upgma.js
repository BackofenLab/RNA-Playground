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
     * @see https://en.wikipedia.org/wiki/Robert_R._Sokal
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
     * dist(c, k = i union j) = [|i|*dist(c, i) + |j|*dist(c, j)] / [|i|+|j|]
     * @param subtree {Object} - The subtree for the new cluster.
     */
    function computeDistances(subtree) {
        var outputData = clusteringInstance.getOutput();

        // retrieve values
        var cluster1Name = subtree.leftChild;
        var cluster2Name = subtree.rightChild;
        var newClusterName = subtree.root;

        var cluster1Cardinality = outputData.cardinalities[cluster1Name];
        var cluster2Cardinality = outputData.cardinalities[cluster2Name];

        var remainingClusterNames = Object.keys(outputData.distanceMatrix);

        for (var i = 0; i < remainingClusterNames.length; i++) {
            var product1 = cluster1Cardinality * outputData.distanceMatrix[[remainingClusterNames[i], cluster1Name]];
            var product2 = cluster2Cardinality * outputData.distanceMatrix[[remainingClusterNames[i], cluster2Name]];

            var dividendSum = product1 + product2;
            var divisorSum = cluster1Cardinality + cluster2Cardinality;

            var quotient = dividendSum / divisorSum;

            outputData.distanceMatrix[[remainingClusterNames[i], newClusterName]] = quotient;
        }
    }

    /**
     * Returns the superclass instance.
     * @return {Object} - Superclass instance.
     */
    function getSuperclass() {
        return upgmaInstance;
    }
}());