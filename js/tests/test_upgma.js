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
        var algorithm = new upgma.Upgma();

        var inputData = {};
        var outputData = {};
        inputData.numOfStartClusters = 5;

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

        algorithm.setIO(inputData, outputData);

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        debugger;
        assertEquals("((a:3.0,b:3.0):2.0,(e:3.0,(e:1.0,d:1.0):2.0):2.0);", outputData.newickString);
    }
});