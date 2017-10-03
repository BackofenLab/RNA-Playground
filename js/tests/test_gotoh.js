/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

TestCase("test_gotoh", {
    "test_1": function () {  // Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
        var algorithm = new gotoh.Gotoh();

        var inputData = {};
        inputData.sequenceA = "TGGA";
        inputData.sequenceB = "GG";

        inputData.calculationType = "similarity";

        inputData.baseCosts = -3;
        inputData.enlargement = -1;
        inputData.match = 0;
        inputData.mismatch = -1;

        inputData.matrixHeight = inputData.sequenceB.length + 1;
        inputData.matrixWidth = inputData.sequenceA.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(-6, outputData.score);
        assertEquals("TGGA", outputData.alignments[0][0]);
        assertEquals("__GG", outputData.alignments[0][2]);
        assertEquals("TGGA", outputData.alignments[1][0]);
        assertEquals("GG__", outputData.alignments[1][2]);
    },

    "test_2": function () {
        var algorithm = new gotoh.Gotoh();

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

        assertEquals(-5, outputData.score);
        assertEquals("CCGA", outputData.alignments[0][0]);
        assertEquals("C__G", outputData.alignments[0][2]);

        assertEquals("CCGA", outputData.alignments[1][0]);
        assertEquals("CG__", outputData.alignments[1][2]);
    },

    "test_3": function () {  // Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
        var algorithm = new gotoh.Gotoh();

        var inputData = {};
        inputData.sequenceA = "TACGCAGA";
        inputData.sequenceB = "TCCGA";

        inputData.calculationType = "similarity";

        inputData.baseCosts = -4;
        inputData.enlargement = -1;
        inputData.match = 1;
        inputData.mismatch = 0;

        inputData.matrixHeight = inputData.sequenceB.length + 1;
        inputData.matrixWidth = inputData.sequenceA.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(-3, outputData.score);
        assertEquals("TACGCAGA", outputData.alignments[0][0]);
        assertEquals("T___CCGA", outputData.alignments[0][2]);

        assertEquals("TACGCAGA", outputData.alignments[1][0]);
        assertEquals("TCC___GA", outputData.alignments[1][2]);

        assertEquals("TACGCAGA", outputData.alignments[2][0]);
        assertEquals("TCCG___A", outputData.alignments[2][2]);
    },

    "test_4": function () {  // Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
        var algorithm = new gotoh.Gotoh();

        var inputData = {};
        inputData.sequenceA = "ACCT";
        inputData.sequenceB = "CC";

        inputData.calculationType = "similarity";

        inputData.baseCosts = -4;
        inputData.enlargement = -1;
        inputData.match = 0;
        inputData.mismatch = -1;

        inputData.matrixHeight = inputData.sequenceB.length + 1;
        inputData.matrixWidth = inputData.sequenceA.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        assertEquals(-7, outputData.score);
        assertEquals("ACCT", outputData.alignments[0][0]);
        assertEquals("__CC", outputData.alignments[0][2]);

        assertEquals("ACCT", outputData.alignments[1][0]);
        assertEquals("CC__", outputData.alignments[1][2])
    }
});