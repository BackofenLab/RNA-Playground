/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("upgma", Upgma, getSuperclass);

    // instances
    var clusteringInstance;
    var upgmaInstance;

    /*---- ALGORITHM ----*/
    /**
     * Computes a clustering with the hierarchical approach
     * unweighted pair group method with arithmetic mean.
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
    }

    /**
     * Returns the superclass instance.
     * @return {Object} - Superclass instance.
     */
    function getSuperclass() {
        return upgmaInstance;
    }
}());