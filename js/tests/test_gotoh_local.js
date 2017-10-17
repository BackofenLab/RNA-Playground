/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

TestCase("test_gotoh_local", {
    /**
     * Short sequences test.
     * @see Test sequences are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    /**
     * Initial sequences test.
     */
    "test_1": function () {
        var algorithm = new gotohLocal.GotohLocal();

        var inputData = {};
        inputData.sequenceA = "CCGA";
        inputData.sequenceB = "CG";

        inputData.calculationType = "similarity";

        inputData.baseCosts = -3;
        inputData.enlargement = -1;
        inputData.match = 1;
        inputData.mismatch = -1;

        inputData.matrixHeight = inputData.sequenceB.length + 1;
        inputData.matrixWidth = inputData.sequenceA.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(2, outputData.score);
        assertEquals("CG", outputData.alignments[0][0]);
        assertEquals("CG", outputData.alignments[0][2]);
    }
});