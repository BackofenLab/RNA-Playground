/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

TestCase("test_single_linkage", {
    /**
     * @see Test values are taken from https://en.wikipedia.org/wiki/Neighbor_joining
     * Hint: Also testing complete linkage, because exactly the same function is used,
     * just with Math.max and not Math.min as argument.
     */
    "test_1": function () {
        var algorithm = new agglomerativeClustering.AgglomerativeClustering();

        var inputData = {};
        var outputData = {};

        inputData.clusteringSubalgorithm = CLUSTERING_ALGORITHMS.SINGLE_LINKAGE;
        inputData.numOfStartClusters = 5;
        inputData.initialNamingIndex = 5;  // needed because of Feng-Doolittle

        outputData.distanceMatrix = {};
        outputData.distanceMatrix[["a", "b"]] = 5;
        outputData.distanceMatrix[["a", "c"]] = 9;
        outputData.distanceMatrix[["a", "d"]] = 9;
        outputData.distanceMatrix[["a", "e"]] = 8;

        outputData.distanceMatrix[["b", "c"]] = 10;
        outputData.distanceMatrix[["b", "d"]] = 10;
        outputData.distanceMatrix[["b", "e"]] = 9;

        outputData.distanceMatrix[["c", "d"]] = 8;
        outputData.distanceMatrix[["c", "e"]] = 7;

        outputData.distanceMatrix[["d", "e"]] = 3;

        outputData.clusterNames = ["a", "b", "c", "d", "e"];

        algorithm.setIO(inputData, outputData);

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals("(((e:1.5,d:1.5):2,c:3.5):0.5,(b:2.5,a:2.5):1.5);", outputData.newickString);
    }
});