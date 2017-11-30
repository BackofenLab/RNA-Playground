/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("interfaces.clusteringInterface", ClusteringInterface);

    // instances
    var clusteringInterfaceInstance;

    /**
     * Is used to work with the input and output (the interface) of an alignment algorithm.
     * It contains the basic methods and the viewmodel for the output.
     * This class is used by the various interface scripts as superclass.
     * @constructor
     */
    function ClusteringInterface() {
        clusteringInterfaceInstance = this;

        // public class methods
        /*
        this.imports = imports;
        this.sharedInterfaceOperations = sharedInterfaceOperations;
        this.startProcessing = startProcessing;
        */
    }
}());