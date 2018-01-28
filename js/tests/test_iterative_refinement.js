/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

TestCase("test_iterative_refinement", {
    /**
     * @see Test values are taken partially from lecture Bioinformatics I.
     */
    "test_1": function () {
        var algorithm = new iterativeRefinement.IterativeRefinement();

        var inputData = {};
        inputData.iterativeRefinementSubalgorithm = ITERATIVE_REFINEMENT_STRATEGIES.MIN_DISTANCE_PAIR;
        inputData.roundRobinOrder = ITERATIVE_REFINEMENT_ORDERS.RIGHT_FIRST;

        inputData.sequences = ["ACGT", "AT", "GCT", "GC"];
        inputData.arrayPositionsOfRemovedSequences = [];
        inputData.initialNamingIndex = 4;

        inputData.calculationType = "similarity";

        inputData.baseCosts = -1;
        inputData.enlargement = -3;
        inputData.match = 1;
        inputData.mismatch = 0;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        // output: final
        assertEquals(-26, outputData.refinedScore);
        assertEquals("ACGT", outputData.refinedProgressiveAlignment[0]);  // a
        assertEquals("GC_T", outputData.refinedProgressiveAlignment[1]);  // c
        assertEquals("GC__", outputData.refinedProgressiveAlignment[2]);  // d
        assertEquals("A__T", outputData.refinedProgressiveAlignment[3]);  // b
    }
});