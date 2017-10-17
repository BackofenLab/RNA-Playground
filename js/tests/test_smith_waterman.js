/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

TestCase("test_smith_waterman", {
    /**
     * Initial sequences test.
     */
    "test_1": function () {
        var algorithm = new smithWaterman.SmithWaterman();

        var inputData = {};
        inputData.sequenceA = "AACG";
        inputData.sequenceB = "AATCG";

        inputData.calculationType = "similarity";

        inputData.deletion = -2;
        inputData.insertion = -2;
        inputData.match = 1;
        inputData.mismatch = -1;

        inputData.matrixHeight = inputData.sequenceB.length + 1;
        inputData.matrixWidth = inputData.sequenceA.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(2, outputData.score);
        assertEquals("AA", outputData.alignments[0][0]);
        assertEquals("AA", outputData.alignments[0][2]);

        assertEquals("CG", outputData.alignments[1][0]);
        assertEquals("CG", outputData.alignments[1][2]);
    },

    /**
     * Multi sequences test.
     * @see: Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_2": function () {
        debugger;
        var algorithm = new smithWaterman.SmithWaterman();

        var inputData = {};
        inputData.sequenceA = "TCCGA";
        inputData.sequenceB = "TACGCAGA";

        inputData.calculationType = "similarity";

        inputData.deletion = -1;
        inputData.insertion = -1;
        inputData.match = 1;
        inputData.mismatch = 0;

        inputData.matrixHeight = inputData.sequenceB.length + 1;
        inputData.matrixWidth = inputData.sequenceA.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(3, outputData.score);
        assertEquals("TCCG", outputData.alignments[0][0]);
        assertEquals("TACG", outputData.alignments[0][2]);

        assertEquals("TCCGA", outputData.alignments[1][0]);
        assertEquals("TACGC", outputData.alignments[1][2]);

        assertEquals("TCCG_A", outputData.alignments[2][0]);
        assertEquals("TACGCA", outputData.alignments[2][2]);

        assertEquals("CCGA", outputData.alignments[3][0]);
        assertEquals("CAGA", outputData.alignments[3][2]);
    },

    /**
     * Long sequences test.
     * @see: Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_3": function () {
        debugger;
        var algorithm = new smithWaterman.SmithWaterman();

        var inputData = {};
        inputData.sequenceA = "CCCCGCGACTCGGGTTCAAGGG";
        inputData.sequenceB = "GGGTGAGACCCCAGTTCAACCC";

        inputData.calculationType = "similarity";

        inputData.deletion = -2;
        inputData.insertion = -2;
        inputData.match = 4;
        inputData.mismatch = -1;

        inputData.matrixHeight = inputData.sequenceB.length + 1;
        inputData.matrixWidth = inputData.sequenceA.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(40, outputData.score);
        assertEquals("GCGACTCGGGTTCAA", outputData.alignments[0][0]);
        assertEquals("GAGACCCCAGTTCAA", outputData.alignments[0][2]);
    }
});
