/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

TestCase("test_waterman_smith_beyer", {
    /**
     * Short sequences test.
     * @see Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_1": function () {
        var algorithm = new watermanSmithBeyer.WatermanSmithBeyer();

        var inputData = {};
        inputData.sequenceB = "TGGA";
        inputData.sequenceA = "GG";

        inputData.calculationType = "similarity";

        inputData.baseCosts = -3;
        inputData.enlargement = -1;
        inputData.match = 0;
        inputData.mismatch = -1;
        inputData.subadditiveFunction = SUBADDITIVE_FUNCTIONS.AFFINE;

        inputData.matrixHeight = inputData.sequenceA.length + 1;
        inputData.matrixWidth = inputData.sequenceB.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(-6, outputData.score);
        assertEquals("TGGA", outputData.alignments[0][2]);
        assertEquals("__GG", outputData.alignments[0][0]);

        assertEquals("TGGA", outputData.alignments[1][2]);
        assertEquals("GG__", outputData.alignments[1][0]);
    },

    /**
     * Initial sequences test.
     */
    "test_2": function () {
        var algorithm = new watermanSmithBeyer.WatermanSmithBeyer();

        var inputData = {};
        inputData.sequenceB = "CCGA";
        inputData.sequenceA = "CG";

        inputData.calculationType = "similarity";

        inputData.baseCosts = -3;
        inputData.enlargement = -1;
        inputData.match = 1;
        inputData.mismatch = -1;
        inputData.subadditiveFunction = SUBADDITIVE_FUNCTIONS.AFFINE;

        inputData.matrixHeight = inputData.sequenceA.length + 1;
        inputData.matrixWidth = inputData.sequenceB.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(-5, outputData.score);
        assertEquals("CCGA", outputData.alignments[0][2]);
        assertEquals("C__G", outputData.alignments[0][0]);

        assertEquals("CCGA", outputData.alignments[1][2]);
        assertEquals("CG__", outputData.alignments[1][0]);
    },

    /**
     * Multi sequences test.
     * @see Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis
     */
    "test_3": function () {
        var algorithm = new watermanSmithBeyer.WatermanSmithBeyer();

        var inputData = {};
        inputData.sequenceB = "TACGCAGA";
        inputData.sequenceA = "TCCGA";

        inputData.calculationType = "similarity";

        inputData.baseCosts = -4;
        inputData.enlargement = -1;
        inputData.match = 1;
        inputData.mismatch = 0;
        inputData.subadditiveFunction = SUBADDITIVE_FUNCTIONS.AFFINE;

        inputData.matrixHeight = inputData.sequenceA.length + 1;
        inputData.matrixWidth = inputData.sequenceB.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(-3, outputData.score);
        assertEquals("TACGCAGA", outputData.alignments[0][2]);
        assertEquals("T___CCGA", outputData.alignments[0][0]);

        assertEquals("TACGCAGA", outputData.alignments[1][2]);
        assertEquals("TCC___GA", outputData.alignments[1][0]);

        assertEquals("TACGCAGA", outputData.alignments[2][2]);
        assertEquals("TCCG___A", outputData.alignments[2][0]);
    },

    /**
     * Short sequences test.
     * @see Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_4": function () {
        var algorithm = new watermanSmithBeyer.WatermanSmithBeyer();

        var inputData = {};
        inputData.sequenceB = "ACCT";
        inputData.sequenceA = "CC";

        inputData.calculationType = "similarity";

        inputData.baseCosts = -4;
        inputData.enlargement = -1;
        inputData.match = 0;
        inputData.mismatch = -1;
        inputData.subadditiveFunction = SUBADDITIVE_FUNCTIONS.AFFINE;

        inputData.matrixHeight = inputData.sequenceA.length + 1;
        inputData.matrixWidth = inputData.sequenceB.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(-7, outputData.score);
        assertEquals("ACCT", outputData.alignments[0][2]);
        assertEquals("__CC", outputData.alignments[0][0]);

        assertEquals("ACCT", outputData.alignments[1][2]);
        assertEquals("CC__", outputData.alignments[1][0])
    },

    /**
     * Simulating Needleman-Wunsch to test correctness.
     * @see Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_5": function () {
        var algorithm = new watermanSmithBeyer.WatermanSmithBeyer();

        var inputData = {};
        inputData.sequenceB = "AGTC";
        inputData.sequenceA = "ATC";

        inputData.calculationType = "similarity";

        inputData.baseCosts = -2;
        inputData.enlargement = 0;
        inputData.match = 1;
        inputData.mismatch = -1;
        inputData.subadditiveFunction = SUBADDITIVE_FUNCTIONS.AFFINE;

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
     * Simulating Needleman-Wunsch to test correctness.
     * @see Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_6": function () {
        var algorithm = new watermanSmithBeyer.WatermanSmithBeyer();

        var inputData = {};
        inputData.sequenceB = "CCCCGCGACTCGGGTTCAAGGG";
        inputData.sequenceA = "GGGTGAGACCCCAGTTCAACCC";

        inputData.calculationType = "similarity";

        inputData.baseCosts = 0;
        inputData.enlargement = -2;
        inputData.match = 4;
        inputData.mismatch = -1;
        inputData.subadditiveFunction = SUBADDITIVE_FUNCTIONS.AFFINE;

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
     * Simulating Needleman-Wunsch to test correctness.
     * @see Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_7": function () {
        var algorithm = new watermanSmithBeyer.WatermanSmithBeyer();

        var inputData = {};
        inputData.sequenceB = "TCCGA";
        inputData.sequenceA = "TACGCGC";

        inputData.calculationType = "similarity";

        inputData.baseCosts = 0;
        inputData.enlargement = -1;
        inputData.match = 1;
        inputData.mismatch = 0;
        inputData.subadditiveFunction = SUBADDITIVE_FUNCTIONS.AFFINE;

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
     * Simulating Needleman-Wunsch to test correctness.
     * @see Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_8": function () {
        var algorithm = new watermanSmithBeyer.WatermanSmithBeyer();

        var inputData = {};
        inputData.sequenceB = "AACG";
        inputData.sequenceA = "AATCG";

        inputData.calculationType = "similarity";

        inputData.baseCosts = 0;
        inputData.enlargement = -2;
        inputData.match = 1;
        inputData.mismatch = -1;
        inputData.subadditiveFunction = SUBADDITIVE_FUNCTIONS.AFFINE;

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