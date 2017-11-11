/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

TestCase("test_needleman_wunsch", {
    /**
     * Short sequences test.
     * @see Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_1": function () {
        var algorithm = new needlemanWunsch.NeedlemanWunsch();

        var inputData = {};
        inputData.sequenceB = "AGTC";
        inputData.sequenceA = "ATC";

        inputData.calculationType = "similarity";

        inputData.deletion = -2;
        inputData.insertion = -2;
        inputData.match = 1;
        inputData.mismatch = -1;

        inputData.matrixHeight = inputData.sequenceA.length + 1;
        inputData.matrixWidth = inputData.sequenceB.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(1, outputData.score);
        assertEquals("AGTC", outputData.alignments[0][2]);
        assertEquals("A_TC", outputData.alignments[0][0]);
    },

    /**
     * Long sequences test.
     * @see Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_2": function () {
        var algorithm = new needlemanWunsch.NeedlemanWunsch();

        var inputData = {};
        inputData.sequenceB = "CCCCGCGACTCGGGTTCAAGGG";
        inputData.sequenceA = "GGGTGAGACCCCAGTTCAACCC";

        inputData.calculationType = "similarity";

        inputData.deletion = -2;
        inputData.insertion = -2;
        inputData.match = 4;
        inputData.mismatch = -1;

        inputData.matrixHeight = inputData.sequenceA.length + 1;
        inputData.matrixWidth = inputData.sequenceB.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(33, outputData.score);
        assertEquals("CCCCGCGACTCGGGTTCAAGGG", outputData.alignments[0][2]);
        assertEquals("GGGTGAGACCCCAGTTCAACCC", outputData.alignments[0][0]);
    },

    /**
     * Short sequences test.
     * @see Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_3": function () {
        var algorithm = new needlemanWunsch.NeedlemanWunsch();

        var inputData = {};
        inputData.sequenceB = "TCCGA";
        inputData.sequenceA = "TACGCGC";

        inputData.calculationType = "similarity";

        inputData.deletion = -1;
        inputData.insertion = -1;
        inputData.match = 1;
        inputData.mismatch = 0;

        inputData.matrixHeight = inputData.sequenceA.length + 1;
        inputData.matrixWidth = inputData.sequenceB.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(2, outputData.score);
        assertEquals("T_C_CGA", outputData.alignments[0][2]);
        assertEquals("TACGCGC", outputData.alignments[0][0]);
    },

    /**
     * Initial sequences test.
     */
    "test_4": function () {
        var algorithm = new needlemanWunsch.NeedlemanWunsch();

        var inputData = {};
        inputData.sequenceB = "AACG";
        inputData.sequenceA = "AATCG";

        inputData.calculationType = "similarity";

        inputData.deletion = -2;
        inputData.insertion = -2;
        inputData.match = 1;
        inputData.mismatch = -1;

        inputData.matrixHeight = inputData.sequenceA.length + 1;
        inputData.matrixWidth = inputData.sequenceB.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(2, outputData.score);
        assertEquals("AA_CG", outputData.alignments[0][2]);
        assertEquals("AATCG", outputData.alignments[0][0]);
    }
});
