/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

TestCase("test_hirschberg", {
    /**
     * @see Test values are taken from lecture Bioinformatics I.
     */
    "test_1": function () {
        debugger;
        var algorithm = new hirschberg.Hirschberg();

        var inputData = {};
        inputData.sequenceB = "AACG";
        inputData.sequenceA = "AATCG";

        inputData.calculationType = "distance";

        inputData.deletion = 2;
        inputData.insertion = 2;
        inputData.match = -1;
        inputData.mismatch = 1;

        inputData.matrixHeight = inputData.sequenceA.length + 1;
        inputData.matrixWidth = inputData.sequenceB.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals("AATCG", outputData.firstSequences[0]);
        assertEquals("AACG", outputData.secondSequences[0]);

        assertEquals("AA", outputData.firstSequences[1]);
        assertEquals("AA", outputData.secondSequences[1]);
    }
});