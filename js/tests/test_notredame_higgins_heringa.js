/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

TestCase("test_notredame_higgins_heringa", {
    /**
     * Short sequences test.
     * @see Test sequences are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_1": function () {
        var algorithm = new notredameHigginsHeringa.NotredameHigginsHeringa();

        var inputData = {};
        inputData.sequences = ["ACGT", "AT", "GCT"];
        inputData.initialNamingIndex = 3;

        inputData.calculationType = "similarity";
        inputData.arrayPositionsOfRemovedSequences = [];

        inputData.baseCosts = 0;
        inputData.enlargement = -2;
        inputData.match = 1;
        inputData.mismatch = -1;

        inputData.globalAlignmentsPerSequencePair = 1;
        inputData.useLocalLibrary = false;
        inputData.localAlignmentsPerSequencePair = 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        // output: pairwise alignment
        assertEquals(4, outputData.alignmentLengths[0]);
        assertEquals(4, outputData.alignmentLengths[1]);
        assertEquals(3, outputData.alignmentLengths[2]);

        assertEquals(2, outputData.gapNumbers[0]);
        assertEquals(1, outputData.gapNumbers[1]);
        assertEquals(1, outputData.gapNumbers[2]);

        assertEquals(1, outputData.gapStarts[0]);
        assertEquals(1, outputData.gapStarts[1]);
        assertEquals(1, outputData.gapStarts[2]);

        assertEquals(-2, outputData.similarities[0]);
        assertEquals(-1, outputData.similarities[1]);
        assertEquals(-2, outputData.similarities[2]);

        var pairAB = outputData.sequencePairs[0];
        var pairAC = outputData.sequencePairs[1];
        var pairBC = outputData.sequencePairs[2];

        var alignmentAB = outputData.alignmentsAndScores[[pairAB[0], pairAB[1]]][0];
        var alignmentAC = outputData.alignmentsAndScores[[pairAC[0], pairAC[1]]][0];
        var alignmentBC = outputData.alignmentsAndScores[[pairBC[0], pairBC[1]]][0];

        assertEquals("ACGT", alignmentAB[0]);
        assertEquals("A__T", alignmentAB[2]);

        assertEquals("ACGT", alignmentAC[0]);
        assertEquals("GC_T", alignmentAC[2]);

        assertEquals("_AT", alignmentBC[0]);
        assertEquals("GCT", alignmentBC[2]);

        // output: primary weight library
        assertEquals(100, outputData.primaryWeightLib[["ACGT", "AT"]][[1, 1]]);
        assertEquals(100, outputData.primaryWeightLib[["ACGT", "AT"]][[4, 2]]);

        assertEquals(200 / 3, outputData.primaryWeightLib[["ACGT", "GCT"]][[1, 1]]);
        assertEquals(200 / 3, outputData.primaryWeightLib[["ACGT", "GCT"]][[2, 2]]);
        assertEquals(200 / 3, outputData.primaryWeightLib[["ACGT", "GCT"]][[4, 3]]);

        assertEquals(50, outputData.primaryWeightLib[["AT", "GCT"]][[1, 2]]);
        assertEquals(50, outputData.primaryWeightLib[["AT", "GCT"]][[2, 3]]);

        // output: extended weight library
        assertEquals(100, outputData.extendedWeightLib[["ACGT", "AT"]][[1, 1]]);
        assertEquals(50, outputData.extendedWeightLib[["ACGT", "AT"]][[2, 1]]);  // earlier zero-edge
        assertEquals(150, outputData.extendedWeightLib[["ACGT", "AT"]][[4, 2]]);

        assertEquals(200 / 3, outputData.extendedWeightLib[["ACGT", "GCT"]][[1, 1]]);
        assertEquals(50, outputData.extendedWeightLib[["ACGT", "GCT"]][[1, 2]]);  // earlier zero-edge
        assertEquals(200 / 3, outputData.extendedWeightLib[["ACGT", "GCT"]][[2, 2]]);
        assertEquals(350 / 3, outputData.extendedWeightLib[["ACGT", "GCT"]][[4, 3]]);

        assertEquals(200 / 3, outputData.extendedWeightLib[["AT", "GCT"]][[1, 1]]);  // earlier zero-edge
        assertEquals(50, outputData.extendedWeightLib[["AT", "GCT"]][[1, 2]]);
        assertEquals(350 / 3, outputData.extendedWeightLib[["AT", "GCT"]][[2, 3]]);

        // output: distances
        assertEquals(1, Math.round(outputData.distanceMatrix[["a", "b"]] * 10) / 10);
        assertEquals(1, Math.round(outputData.distanceMatrix[["a", "c"]] * 10) / 10);
        assertEquals(1.5, Math.round(outputData.distanceMatrix[["b", "c"]] * 10) / 10);

        // output: phylogenetic tree
        assertEquals("((b:0.4904,a:0.4904):0.136,c:0.6264);", outputData.newickString);

        // output: joinment
        // matrix a~b (is mirrored)
        assertEquals(100, outputData.groupMatrices["ab"][1][1]);
        assertEquals(100, outputData.groupMatrices["ab"][1][2]);

        assertEquals(100, outputData.groupMatrices["ab"][2][1]);
        assertEquals(100, outputData.groupMatrices["ab"][2][2]);

        assertEquals(100, outputData.groupMatrices["ab"][3][1]);
        assertEquals(100, outputData.groupMatrices["ab"][3][2]);

        assertEquals(100, outputData.groupMatrices["ab"][4][1]);
        assertEquals(250, outputData.groupMatrices["ab"][4][2]);

        // matrix ab~c (mirrored, so really "c~ab")
        debugger;
        assertEquals(200 / 3, outputData.groupMatrices["cab"][1][1]);
        assertEquals(200 / 3, outputData.groupMatrices["cab"][1][2]);
        assertEquals(200 / 3, outputData.groupMatrices["cab"][1][3]);
        assertEquals(200 / 3, outputData.groupMatrices["cab"][1][4]);

        assertEquals(200 / 3, outputData.groupMatrices["cab"][2][1]);
        assertEquals(100, outputData.groupMatrices["cab"][2][2]);
        assertEquals(100, outputData.groupMatrices["cab"][2][3]);
        assertEquals(100, outputData.groupMatrices["cab"][2][4]);

        assertEquals(200 / 3, outputData.groupMatrices["cab"][3][1]);
        assertEquals(100, outputData.groupMatrices["cab"][3][2]);
        assertEquals(100, outputData.groupMatrices["cab"][3][3]);
        assertEquals(Math.round(650 / 3), Math.round(outputData.groupMatrices["cab"][3][4]));

        // output: final
        assertEquals(-5, outputData.score);
        assertEquals("GC_T", outputData.progressiveAlignment[0]);
        assertEquals("ACGT", outputData.progressiveAlignment[1]);
        assertEquals("A__T", outputData.progressiveAlignment[2]);
    },

    /**
     * Short sequences test.
     * @see Test sequences are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_2": function () {
        var algorithm = new notredameHigginsHeringa.NotredameHigginsHeringa();

        var inputData = {};
        inputData.sequences = ["ACGT", "AT", "GCT"];
        inputData.initialNamingIndex = 3;

        inputData.calculationType = "similarity";
        inputData.arrayPositionsOfRemovedSequences = [];

        inputData.baseCosts = 0;
        inputData.enlargement = -2;
        inputData.match = 1;
        inputData.mismatch = -1;

        inputData.globalAlignmentsPerSequencePair = 1;
        inputData.useLocalLibrary = true;
        inputData.totalNumberAlignments = 3;
        inputData.localAlignmentsPerSequencePair = 1;

        inputData.baseCostsLocal = 0;
        inputData.enlargementLocal = -2;
        inputData.matchLocal = 1;
        inputData.mismatchLocal = -1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        // output: pairwise alignment
        assertEquals(4, outputData.alignmentLengths[0]);
        assertEquals(4, outputData.alignmentLengths[1]);
        assertEquals(3, outputData.alignmentLengths[2]);

        assertEquals(2, outputData.gapNumbers[0]);
        assertEquals(1, outputData.gapNumbers[1]);
        assertEquals(1, outputData.gapNumbers[2]);

        assertEquals(1, outputData.gapStarts[0]);
        assertEquals(1, outputData.gapStarts[1]);
        assertEquals(1, outputData.gapStarts[2]);

        assertEquals(-2, outputData.similarities[0]);
        assertEquals(-1, outputData.similarities[1]);
        assertEquals(-2, outputData.similarities[2]);

        var pairAB = outputData.sequencePairs[0];
        var pairAC = outputData.sequencePairs[1];
        var pairBC = outputData.sequencePairs[2];

        var alignmentAB = outputData.alignmentsAndScores[[pairAB[0], pairAB[1]]][0];
        var alignmentAC = outputData.alignmentsAndScores[[pairAC[0], pairAC[1]]][0];
        var alignmentBC = outputData.alignmentsAndScores[[pairBC[0], pairBC[1]]][0];

        assertEquals("ACGT", alignmentAB[0]);
        assertEquals("A__T", alignmentAB[2]);

        assertEquals("ACGT", alignmentAC[0]);
        assertEquals("GC_T", alignmentAC[2]);

        assertEquals("_AT", alignmentBC[0]);
        assertEquals("GCT", alignmentBC[2]);

        // output: primary weight library
        assertEquals(200, outputData.primaryWeightLib[["ACGT", "AT"]][[1, 1]]);
        assertEquals(100, outputData.primaryWeightLib[["ACGT", "AT"]][[4, 2]]);

        assertEquals(200 / 3, outputData.primaryWeightLib[["ACGT", "GCT"]][[1, 1]]);
        assertEquals(Math.round(500 / 3), Math.round(outputData.primaryWeightLib[["ACGT", "GCT"]][[2, 2]]));
        assertEquals(200 / 3, outputData.primaryWeightLib[["ACGT", "GCT"]][[4, 3]]);

        assertEquals(50, outputData.primaryWeightLib[["AT", "GCT"]][[1, 2]]);
        assertEquals(150, outputData.primaryWeightLib[["AT", "GCT"]][[2, 3]]);

        // output: extended weight library
        assertEquals(200, outputData.extendedWeightLib[["ACGT", "AT"]][[1, 1]]);
        assertEquals(Math.round(500 / 3), Math.round(outputData.extendedWeightLib[["ACGT", "AT"]][[4, 2]]));

        assertEquals(200 / 3, outputData.extendedWeightLib[["ACGT", "GCT"]][[1, 1]]);
        assertEquals(Math.round(500 / 3), Math.round(outputData.extendedWeightLib[["ACGT", "GCT"]][[2, 2]]));
        assertEquals(Math.round(500 / 3), Math.round(outputData.extendedWeightLib[["ACGT", "GCT"]][[4, 3]]));

        assertEquals(50, outputData.extendedWeightLib[["AT", "GCT"]][[1, 2]]);
        assertEquals(Math.round(650 / 3), Math.round(outputData.extendedWeightLib[["AT", "GCT"]][[2, 3]]));

        // output: joinment
        // matrix a~b (is mirrored)
        assertEquals(200, outputData.groupMatrices["ab"][1][1]);
        assertEquals(200, outputData.groupMatrices["ab"][1][2]);

        assertEquals(200, outputData.groupMatrices["ab"][2][1]);
        assertEquals(200, outputData.groupMatrices["ab"][2][2]);

        assertEquals(200, outputData.groupMatrices["ab"][3][1]);
        assertEquals(200, outputData.groupMatrices["ab"][3][2]);

        assertEquals(200, outputData.groupMatrices["ab"][4][1]);
        assertEquals(1100 / 3, outputData.groupMatrices["ab"][4][2]);

        // matrix ab~c (mirrored, so really "c~ab")
        assertEquals(200 / 3, outputData.groupMatrices["cab"][1][1]);
        assertEquals(200 / 3, outputData.groupMatrices["cab"][1][2]);
        assertEquals(200 / 3, outputData.groupMatrices["cab"][1][3]);
        assertEquals(200 / 3, outputData.groupMatrices["cab"][1][4]);

        assertEquals(200 / 3, outputData.groupMatrices["cab"][2][1]);
        assertEquals(150, outputData.groupMatrices["cab"][2][2]);
        assertEquals(150, outputData.groupMatrices["cab"][2][3]);
        assertEquals(150, outputData.groupMatrices["cab"][2][4]);

        assertEquals(200 / 3, outputData.groupMatrices["cab"][3][1]);
        assertEquals(150, outputData.groupMatrices["cab"][3][2]);
        assertEquals(150, outputData.groupMatrices["cab"][3][3]);
        assertEquals(Math.round(1025 / 3), Math.round(outputData.groupMatrices["cab"][3][4]));

        // output: final
        assertEquals(-5, outputData.score);
        assertEquals("GC_T", outputData.progressiveAlignment[0]);
        assertEquals("ACGT", outputData.progressiveAlignment[1]);
        assertEquals("A__T", outputData.progressiveAlignment[2]);
    },

    /**
     * Two alignments test.
     */
    "test_3": function () {
        var algorithm = new notredameHigginsHeringa.NotredameHigginsHeringa();

        var inputData = {};
        inputData.sequences = ["ACGT", "AT", "GCT"];
        inputData.initialNamingIndex = 3;

        inputData.calculationType = "similarity";
        inputData.arrayPositionsOfRemovedSequences = [];

        inputData.baseCosts = 0;
        inputData.enlargement = -2;
        inputData.match = 1;
        inputData.mismatch = -1;

        inputData.globalAlignmentsPerSequencePair = 2;  // two alignments
        inputData.useLocalLibrary = false;
        inputData.totalNumberAlignments = 3;
        inputData.localAlignmentsPerSequencePair = 1;

        inputData.baseCostsLocal = 0;
        inputData.enlargementLocal = -2;
        inputData.matchLocal = 1;
        inputData.mismatchLocal = -1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        // output: primary weight library
        assertEquals(200, outputData.primaryWeightLib[["ACGT", "AT"]][[1, 1]]);
        assertEquals(200, outputData.primaryWeightLib[["ACGT", "AT"]][[4, 2]]);

        assertEquals(200 / 3, outputData.primaryWeightLib[["ACGT", "GCT"]][[1, 1]]);
        assertEquals(200 / 3, outputData.primaryWeightLib[["ACGT", "GCT"]][[2, 2]]);
        assertEquals(200 / 3, outputData.primaryWeightLib[["ACGT", "GCT"]][[4, 3]]);

        assertEquals(50, outputData.primaryWeightLib[["AT", "GCT"]][[1, 1]]);
        assertEquals(50, outputData.primaryWeightLib[["AT", "GCT"]][[1, 2]]);
        assertEquals(100, outputData.primaryWeightLib[["AT", "GCT"]][[2, 3]]);
    },

    /**
     * Two local alignments test.
     */
    "test_4": function () {
        var algorithm = new notredameHigginsHeringa.NotredameHigginsHeringa();

        var inputData = {};
        inputData.sequences = ["ACGT", "AT", "GCT"];
        inputData.initialNamingIndex = 3;

        inputData.calculationType = "similarity";
        inputData.arrayPositionsOfRemovedSequences = [];

        inputData.baseCosts = 0;
        inputData.enlargement = -2;
        inputData.match = 1;
        inputData.mismatch = -1;

        inputData.globalAlignmentsPerSequencePair = 1;  // two alignments
        inputData.useLocalLibrary = true;
        inputData.totalNumberAlignments = 3;
        inputData.localAlignmentsPerSequencePair = 2;

        inputData.baseCostsLocal = 0;
        inputData.enlargementLocal = -2;
        inputData.matchLocal = 1;
        inputData.mismatchLocal = -1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        // output: primary weight library
        assertEquals(200, outputData.primaryWeightLib[["ACGT", "AT"]][[1, 1]]);
        assertEquals(200, outputData.primaryWeightLib[["ACGT", "AT"]][[4, 2]]);

        assertEquals(200 / 3, outputData.primaryWeightLib[["ACGT", "GCT"]][[1, 1]]);
        assertEquals(Math.round(500 / 3), Math.round(outputData.primaryWeightLib[["ACGT", "GCT"]][[2, 2]]));
        assertEquals(100, outputData.primaryWeightLib[["ACGT", "GCT"]][[3, 1]]);
        assertEquals(200 / 3, outputData.primaryWeightLib[["ACGT", "GCT"]][[4, 3]]);

        assertEquals(50, outputData.primaryWeightLib[["AT", "GCT"]][[1, 2]]);
        assertEquals(150, outputData.primaryWeightLib[["AT", "GCT"]][[2, 3]]);
    }
});