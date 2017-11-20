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

        // sequences
        assertEquals("AATCG", outputData.firstSequences[0]);
        assertEquals("AACG", outputData.secondSequences[0]);

        assertEquals("AA", outputData.firstSequences[1]);
        assertEquals("AA", outputData.secondSequences[1]);

        assertEquals("CG", outputData.firstSequences[2]);
        assertEquals("ACG", outputData.secondSequences[2]);

        // positions
        assertEquals([1,2,3,4,5], outputData.firstSequencePositions[0]);
        assertEquals([1,2,3,4], outputData.secondSequencePositions[0]);

        assertEquals([1,2], outputData.firstSequencePositions[1]);
        assertEquals([1,2], outputData.secondSequencePositions[1]);

        assertEquals([4,5], outputData.firstSequencePositions[2]);
        assertEquals([2,3,4], outputData.secondSequencePositions[2]);

        // sums
        assertEquals([7,2,0,2,7], outputData.addedRows[0]);
        assertEquals([3,-2,3], outputData.addedRows[1]);
        assertEquals([5,2,0,5], outputData.addedRows[2]);

        // local minimum
        assertEquals([3,2], outputData.minimum[0]);
        assertEquals([1,1], outputData.minimum[1]);
        assertEquals([1,2], outputData.minimum[2]);
    }
});