/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

TestCase("test_notredame_higgins_heringa", {
    /**
     * Short sequences test.
     * @see Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_1": function () {
        debugger;

        var algorithm = new notredameHigginsHeringa.NotredameHigginsHeringa();

        var inputData = {};
        inputData.sequences = ["ACGT", "AT", "GCT"];

        inputData.calculationType = "similarity";
        inputData.arrayPositionsOfRemovedSequences = [];

        inputData.baseCosts = 0;
        inputData.enlargement = -2;
        inputData.match = 1;
        inputData.mismatch = -1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        debugger;
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
        assertEquals(100, outputData.primaryWeightLib[["ACGT","AT"]][[1,1]]);
        assertEquals(100, outputData.primaryWeightLib[["ACGT","AT"]][[4,2]]);

        assertEquals(200/3, outputData.primaryWeightLib[["ACGT","GCT"]][[1,1]]);
        assertEquals(200/3, outputData.primaryWeightLib[["ACGT","GCT"]][[2,2]]);
        assertEquals(200/3, outputData.primaryWeightLib[["ACGT","GCT"]][[4,3]]);

        assertEquals(50, outputData.primaryWeightLib[["AT","GCT"]][[1,2]]);
        assertEquals(50, outputData.primaryWeightLib[["AT","GCT"]][[2,3]]);

        // output: extended weight library
        assertEquals(100, outputData.extendedWeightLib[["ACGT","AT"]][[1,1]]);
        assertEquals(150, outputData.extendedWeightLib[["ACGT","AT"]][[4,2]]);

        assertEquals(200/3, outputData.extendedWeightLib[["ACGT","GCT"]][[1,1]]);
        assertEquals(200/3, outputData.extendedWeightLib[["ACGT","GCT"]][[2,2]]);
        assertEquals(350/3, outputData.extendedWeightLib[["ACGT","GCT"]][[4,3]]);

        assertEquals(50, outputData.extendedWeightLib[["AT","GCT"]][[1,2]]);
        assertEquals(350/3, outputData.extendedWeightLib[["AT","GCT"]][[2,3]]);

        // output: distances
        assertEquals(1, Math.round(outputData.distanceMatrix[["a", "b"]] * 10) / 10);
        assertEquals(1, Math.round(outputData.distanceMatrix[["a", "c"]] * 10) / 10);
        assertEquals(1.5, Math.round(outputData.distanceMatrix[["b", "c"]] * 10) / 10);

        // output: phylogenetic tree
        assertEquals("(c:0.6264,(a:0.4904,b:0.4904):0.136);", outputData.newickString);

        debugger;
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
        assertEquals(100/3, outputData.groupMatrices["cab"][1][1]);
        assertEquals(100/3, outputData.groupMatrices["cab"][1][2]);
        assertEquals(100/3, outputData.groupMatrices["cab"][1][3]);
        assertEquals(100/3, outputData.groupMatrices["cab"][1][4]);

        assertEquals(100/3, outputData.groupMatrices["cab"][2][1]);
        assertEquals(200/3, outputData.groupMatrices["cab"][2][2]);
        assertEquals(200/3, outputData.groupMatrices["cab"][2][3]);
        assertEquals(200/3, outputData.groupMatrices["cab"][2][4]);

        assertEquals(100/3, outputData.groupMatrices["cab"][3][1]);
        assertEquals(200/3, outputData.groupMatrices["cab"][3][2]);
        assertEquals(200/3, outputData.groupMatrices["cab"][3][3]);
        assertEquals(550/3, outputData.groupMatrices["cab"][3][4]);

        // output: final
        assertEquals(-5, outputData.score);
        assertEquals("GC_T", outputData.progressiveAlignment[0]);
        assertEquals("ACGT", outputData.progressiveAlignment[1]);
        assertEquals("A__T", outputData.progressiveAlignment[2]);
    }
});