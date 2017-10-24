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
        var algorithm = new upgma.Upgma();

        var inputData = {};
        inputData.numOfStartClusters = 5;

        inputData.distanceMatrix[["a", "b"]] = 6;
        inputData.distanceMatrix[["a", "c"]] = 10;
        inputData.distanceMatrix[["a", "d"]] = 10;
        inputData.distanceMatrix[["a", "e"]] = 10;

        inputData.distanceMatrix[["b", "c"]] = 10;
        inputData.distanceMatrix[["b", "d"]] = 10;
        inputData.distanceMatrix[["b", "e"]] = 10;

        inputData.distanceMatrix[["c", "d"]] = 2;
        inputData.distanceMatrix[["c", "e"]] = 6;

        inputData.distanceMatrix[["d", "e"]] = 6;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];
    }
});