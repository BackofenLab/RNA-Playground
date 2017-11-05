/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

TestCase("test_feng_doolittle", {
    /**
     * Short sequences test.
     * @see Test values are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_1": function () {
        var algorithm = new fengDoolittle.FengDoolittle();

        var inputData = {};
        inputData.sequences = ["ACGT", "AT", "GCC"];

        inputData.calculationType = "similarity";
        inputData.arrayPositionsOfRemovedSequences = [];

        inputData.baseCosts = 0;
        inputData.enlargement = -2;
        inputData.match = 1;
        inputData.mismatch = -1;

        algorithm.setIO(inputData, {});

        debugger;
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
        assertEquals(-3, outputData.similarities[1]);
        assertEquals(-4, outputData.similarities[2]);

        var pairAB = outputData.sequencePairs[0];
        var pairAC = outputData.sequencePairs[1];
        var pairBC = outputData.sequencePairs[2];

        var alignmentAB = outputData.alignmentsAndScores[[pairAB[0], pairAB[1]]][0];
        var alignmentAC = outputData.alignmentsAndScores[[pairAC[0], pairAC[1]]][0];
        var alignmentBC = outputData.alignmentsAndScores[[pairBC[0], pairBC[1]]][0];

        assertEquals("ACGT", alignmentAB[0]);
        assertEquals("A__T", alignmentAB[2]);

        assertEquals("ACGT", alignmentAC[0]);
        assertEquals("GC_C", alignmentAC[2]);

        assertEquals("_AT", alignmentBC[0]);
        assertEquals("GCC", alignmentBC[2]);

        // output: distances
        assertEquals(1, Math.round(outputData.distanceMatrix[["a", "b"]]));
        assertEquals(3, Math.round(outputData.distanceMatrix[["a", "c"]]));
        assertEquals(9, Math.round(outputData.distanceMatrix[["b", "c"]]));

        // output: phylogenetic tree
        assertEquals("(c:2.8547,(a:0.4904,b:0.4904):2.3642);", outputData.newickString);

        // output: final
        assertEquals(-9, outputData.score);
        assertEquals("GC_C", outputData.progressiveAlignment[0]);
        assertEquals("ACGT", outputData.progressiveAlignment[1]);
        assertEquals("A__T", outputData.progressiveAlignment[2]);
    }
});