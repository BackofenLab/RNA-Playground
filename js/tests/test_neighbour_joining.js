/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

TestCase("test_neighbour_joining", {
    /**
     * @see Test values are taken from https://en.wikipedia.org/wiki/Neighbor_joining
     */
    "test_1": function () {
        debugger;
        var algorithm = new neighbourJoining.NeighbourJoining();

        var inputData = {};
        var outputData = {};

        inputData.clusteringSubalgorithm = CLUSTERING_ALGORITHMS.NEIGHBOUR_JOINING;
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

        assertEquals("(((b:3,a:2):3,c:4):0,(e:1,d:2):2);", outputData.newickString);
    }
});