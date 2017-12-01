/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

TestCase("test_hirschberg", {
    /**
     * Short sequences test.
     * @see Test values are taken from lecture Bioinformatics I.
     */
    "test_1": function () {
        var algorithm = new hirschberg.Hirschberg();

        var inputData = {};
        inputData.sequenceA = "AGTC";
        inputData.sequenceB = "ATC";

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

        assertEquals("AGTC", outputData.alignments[0][0]);
        assertEquals("A_TC", outputData.alignments[0][2]);
    },

    /**
     * Long sequences test.
     * @see Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_2": function () {
        var algorithm = new hirschberg.Hirschberg();

        var inputData = {};
        inputData.sequenceA = "GGGTGAGACCCCAGTTCAACCC";
        inputData.sequenceB = "CCCCGCGACTCGGGTTCAAGGG";

        inputData.calculationType = "distance";

        inputData.deletion = 2;
        inputData.insertion = 2;
        inputData.match = -4;
        inputData.mismatch = 1;

        inputData.matrixHeight = inputData.sequenceA.length + 1;
        inputData.matrixWidth = inputData.sequenceB.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals("GGGTGAGACCCCAGTTCAACCC", outputData.alignments[0][0]);
        assertEquals("CCCCGCGACTCGGGTTCAAGGG", outputData.alignments[0][2]);
    },

    /**
     * Short sequences test.
     * @see Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_3": function () {
        var algorithm = new hirschberg.Hirschberg();

        var inputData = {};
        inputData.sequenceA = "TCCGA";
        inputData.sequenceB = "TACGCGC";

        inputData.calculationType = "distance";

        inputData.deletion = 1;
        inputData.insertion = 1;
        inputData.match = -1;
        inputData.mismatch = 0;

        inputData.matrixHeight = inputData.sequenceA.length + 1;
        inputData.matrixWidth = inputData.sequenceB.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals("T_C_CGA", outputData.alignments[0][0]);
        assertEquals("TACGCGC", outputData.alignments[0][2]);
    },

    /**
     * Short sequences test.
     * @see Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_4": function () {
        var algorithm = new hirschberg.Hirschberg();

        var inputData = {};
        inputData.sequenceA = "AATCG";
        inputData.sequenceB = "AACG";

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

        assertEquals("AATCG", outputData.alignments[0][0]);
        assertEquals("AA_CG", outputData.alignments[0][2]);
    },

    /**
     * Multi traceback test.
     * @see Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_5": function () {
        debugger;
        var algorithm = new hirschberg.Hirschberg();

        var inputData = {};
        inputData.sequenceA = "AATCG";
        inputData.sequenceB = "ACG";

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

        assertEquals("AATCG", outputData.alignments[0][0]);
        assertEquals("_A_CG", outputData.alignments[0][2]);
    }
});