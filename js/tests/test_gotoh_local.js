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
        inputData.sequenceB = "CCGA";
        inputData.sequenceA = "CG";

        inputData.calculationType = "similarity";

        inputData.baseCosts = -3;
        inputData.enlargement = -1;
        inputData.match = 1;
        inputData.mismatch = -1;

        inputData.matrixHeight = inputData.sequenceA.length + 1;
        inputData.matrixWidth = inputData.sequenceB.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(2, outputData.score);
        assertEquals("CG", outputData.alignments[0][2]);
        assertEquals("CG", outputData.alignments[0][0]);
    },

    /**
     * Simulating Smith-Waterman to test correctness.
     */
    "test_2": function () {
        var algorithm = new gotohLocal.GotohLocal();

        var inputData = {};
        inputData.sequenceB = "AACG";
        inputData.sequenceA = "AATCG";

        inputData.calculationType = "similarity";

        inputData.baseCosts = 0;
        inputData.enlargement = -2;
        inputData.match = 1;
        inputData.mismatch = -1;

        inputData.matrixHeight = inputData.sequenceA.length + 1;
        inputData.matrixWidth = inputData.sequenceB.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(2, outputData.score);
        assertEquals("AA", outputData.alignments[0][2]);
        assertEquals("AA", outputData.alignments[0][0]);

        assertEquals("CG", outputData.alignments[1][2]);
        assertEquals("CG", outputData.alignments[1][0]);
    },

    /**
     * Multi sequences test.
     * Simulating Smith-Waterman to test correctness.
     */
    "test_3": function () {
        var algorithm = new gotohLocal.GotohLocal();

        var inputData = {};
        inputData.sequenceB = "TCCGA";
        inputData.sequenceA = "TACGCAGA";

        inputData.calculationType = "similarity";

        inputData.baseCosts = 0;
        inputData.enlargement = -1;
        inputData.match = 1;
        inputData.mismatch = 0;

        inputData.matrixHeight = inputData.sequenceA.length + 1;
        inputData.matrixWidth = inputData.sequenceB.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(3, outputData.score);
        assertEquals("TCCG", outputData.alignments[0][2]);
        assertEquals("TACG", outputData.alignments[0][0]);

        assertEquals("TCCGA", outputData.alignments[1][2]);
        assertEquals("TACGC", outputData.alignments[1][0]);

        assertEquals("TCCG_A", outputData.alignments[2][2]);
        assertEquals("TACGCA", outputData.alignments[2][0]);

        assertEquals("CCGA", outputData.alignments[3][2]);
        assertEquals("CAGA", outputData.alignments[3][0]);
    },

    /**
     * Long sequences test.
     * Simulating Smith-Waterman to test correctness.
     * @see: Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_4": function () {
        var algorithm = new gotohLocal.GotohLocal();

        var inputData = {};
        inputData.sequenceB = "CCCCGCGACTCGGGTTCAAGGG";
        inputData.sequenceA = "GGGTGAGACCCCAGTTCAACCC";

        inputData.calculationType = "similarity";

        inputData.baseCosts = 0;
        inputData.enlargement = -2;
        inputData.match = 4;
        inputData.mismatch = -1;

        inputData.matrixHeight = inputData.sequenceA.length + 1;
        inputData.matrixWidth = inputData.sequenceB.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(40, outputData.score);
        assertEquals("GCGACTCGGGTTCAA", outputData.alignments[0][2]);
        assertEquals("GAGACCCCAGTTCAA", outputData.alignments[0][0]);
    }
});