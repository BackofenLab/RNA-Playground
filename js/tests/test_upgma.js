/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

TestCase("test_upgma", {
    /**
     * @see Test values are taken from lecture Bioinformatics I.
     */
    "test_1": function () {
        debugger;
        var algorithm = new agglomerativeClustering.AgglomerativeClustering();

        var inputData = {};
        var outputData = {};

        inputData.clusteringSubalgorithm = CLUSTERING_ALGORITHMS.UPGMA;
        inputData.numOfStartClusters = 5;
        inputData.initialNamingIndex = 5;  // needed because of Feng-Doolittle

        outputData.distanceMatrix = {};
        outputData.distanceMatrix[["a", "b"]] = 6;
        outputData.distanceMatrix[["a", "c"]] = 10;
        outputData.distanceMatrix[["a", "d"]] = 10;
        outputData.distanceMatrix[["a", "e"]] = 10;

        outputData.distanceMatrix[["b", "c"]] = 10;
        outputData.distanceMatrix[["b", "d"]] = 10;
        outputData.distanceMatrix[["b", "e"]] = 10;

        outputData.distanceMatrix[["c", "d"]] = 2;
        outputData.distanceMatrix[["c", "e"]] = 6;

        outputData.distanceMatrix[["d", "e"]] = 6;

        outputData.clusterNames = ["a", "b", "c", "d", "e"];

        algorithm.setIO(inputData, outputData);

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals("(((d:1,c:1):2,e:3):2,(b:3,a:3):2);", outputData.newickString);
    },

    /**
     * @see Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_2": function () {
        debugger;
        var algorithm = new agglomerativeClustering.AgglomerativeClustering();

        var inputData = {};
        var outputData = {};

        inputData.clusteringSubalgorithm = CLUSTERING_ALGORITHMS.UPGMA;
        inputData.numOfStartClusters = 3;
        inputData.initialNamingIndex = 3;  // needed because of Feng-Doolittle

        outputData.distanceMatrix = {};
        outputData.distanceMatrix[["a", "b"]] = 2;
        outputData.distanceMatrix[["a", "c"]] = 3;
        outputData.distanceMatrix[["b", "c"]] = 3;

        outputData.clusterNames = ["a", "b", "c"];

        algorithm.setIO(inputData, outputData);

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals("((b:1,a:1):0.5,c:1.5);", outputData.newickString);
    },

    /**
     * @see Test values are taken from project Algorithms
     * for Bioninformatics of Alexander Mattheis
     * or from lecture Bioinformatics I.
     */
    "test_3": function () {
        debugger;
        var algorithm = new agglomerativeClustering.AgglomerativeClustering();

        var inputData = {};
        var outputData = {};

        inputData.clusteringSubalgorithm = CLUSTERING_ALGORITHMS.UPGMA;
        inputData.numOfStartClusters = 3;
        inputData.initialNamingIndex = 3;  // needed because of Feng-Doolittle

        outputData.distanceMatrix = {};
        outputData.distanceMatrix[["a", "b"]] = 3;
        outputData.distanceMatrix[["a", "c"]] = 6;
        outputData.distanceMatrix[["b", "c"]] = 6;

        outputData.clusterNames = ["a", "b", "c"];

        algorithm.setIO(inputData, outputData);

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals("((b:1.5,a:1.5):1.5,c:3);", outputData.newickString);
    }
});