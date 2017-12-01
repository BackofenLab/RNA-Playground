/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

TestCase("test_arslan_egecioglu_pevzner", {
    /**
     * Long sequences test.
     * @see Test sequences are taken from lecture Bioinformatics I.
     */
    "test_1": function () {
        var algorithm = new arslanEgeciougluPevzner.ArslanEgeciougluPevzner();

        var inputData = {};
        inputData.sequenceB = "GCATTUGCCUU";
        inputData.sequenceA = "CTTGACCATU";

        inputData.calculationType = "similarity";

        inputData.deletion = -2;
        inputData.insertion = -2;
        inputData.match = 3;
        inputData.mismatch = -1;
        inputData.length = 10;

        inputData.matrixHeight = inputData.sequenceA.length + 1;
        inputData.matrixWidth = inputData.sequenceB.length + 1;

        algorithm.setIO(inputData, {});

        var ioData = algorithm.compute();
        var outputData = ioData[1];

        // checking values in every iteration
        // iterationData[i][j] = [score, length, lambda, deletion, insertion, match, mismatch, alignments, matrix]

        // [first possibility]
        // [1]
        var firstIterationData = outputData.iterationData[0][0];

        var score = firstIterationData[0];
        var alignmentLength = firstIterationData[1];
        var lambda = firstIterationData[2];
        var deletion = firstIterationData[3];
        var insertion = firstIterationData[4];
        var match = firstIterationData[5];
        var mismatch = firstIterationData[6];
        var alignments = firstIterationData[7];

        assertEquals(12, score);
        assertEquals(15, alignmentLength);
        assertEquals(0.48, lambda);

        assertEquals(-2.48, deletion);
        assertEquals(-2.48, insertion);
        assertEquals(2.04, match);
        assertEquals(-1.96, mismatch);

        assertEquals("CATTUG_CC", alignments[0][2]);
        assertEquals("C_TT_GACC", alignments[0][0]);

        assertEquals("CATTUG_CC_UU", alignments[1][2]);
        assertEquals("C_TT_GACCATU", alignments[1][0]);

        assertEquals("CATTUG_CCU_U", alignments[2][2]);
        assertEquals("C_TT_GACCATU", alignments[2][0]);

        // [2]
        var secondIterationData = outputData.iterationData[0][1];

        score = secondIterationData[0];
        alignmentLength = secondIterationData[1];
        lambda = secondIterationData[2];
        deletion = secondIterationData[3];
        insertion = secondIterationData[4];
        match = secondIterationData[5];
        mismatch = secondIterationData[6];
        alignments = secondIterationData[7];

        assertEquals(9, score);
        assertEquals(6, alignmentLength);
        assertEquals(0.5625, lambda);

        assertEquals(-2.5625, deletion);
        assertEquals(-2.5625, insertion);
        assertEquals(1.875, match);
        assertEquals(-2.125, mismatch);

        assertEquals("CAT", alignments[0][0]);
        assertEquals("CAT", alignments[0][2]);

        // [3]
        var thirdIterationData = outputData.iterationData[0][2];

        score = thirdIterationData[0];
        alignmentLength = thirdIterationData[1];
        lambda = thirdIterationData[2];
        deletion = thirdIterationData[3];
        insertion = thirdIterationData[4];
        match = thirdIterationData[5];
        mismatch = thirdIterationData[6];
        alignments = thirdIterationData[7];

        assertEquals(9, score);
        assertEquals(6, alignmentLength);
        assertEquals(0.5625, lambda);

        assertEquals(-2.5625, deletion);
        assertEquals(-2.5625, insertion);
        assertEquals(1.875, match);
        assertEquals(-2.125, mismatch);

        assertEquals("CAT", alignments[0][2]);
        assertEquals("CAT", alignments[0][0]);

        /*
        // [All path test-cases below!]
        // [second/third possibility]
        // [4,7]
        var fourthIterationData = outputData.iterationData[1][0];
        var seventhIterationData = outputData.iterationData[2][0];

        score = fourthIterationData[0] === seventhIterationData[0];
        alignmentLength = fourthIterationData[1] === seventhIterationData[1];
        lambda = fourthIterationData[2] === seventhIterationData[2];
        deletion = fourthIterationData[3] === seventhIterationData[3];
        insertion = fourthIterationData[4] === seventhIterationData[4];
        match = fourthIterationData[5] === seventhIterationData[5];
        mismatch = fourthIterationData[6] === seventhIterationData[6];
        alignmentsFourth = fourthIterationData[7];
        alignmentsSeventh = seventhIterationData[7];

        assertEquals(true, score);
        assertEquals(true, alignmentLength);
        assertEquals(true, lambda);

        assertEquals(true, deletion);
        assertEquals(true, insertion);
        assertEquals(true, match);
        assertEquals(true, mismatch);

        assertEquals("CATTUG_CC", alignmentsFourth[0][2]);
        assertEquals("C_TT_GACC", alignmentsFourth[0][0]);

        assertEquals("CATTUG_CC_UU", alignmentsFourth[1][2]);
        assertEquals("C_TT_GACCATU", alignmentsFourth[1][0]);

        assertEquals("CATTUG_CCU_U", alignmentsFourth[2][2]);
        assertEquals("C_TT_GACCATU", alignmentsFourth[2][0]);

        assertEquals("CATTUG_CC", alignmentsSeventh[0][2]);
        assertEquals("C_TT_GACC", alignmentsSeventh[0][0]);

        assertEquals("CATTUG_CC_UU", alignmentsSeventh[1][2]);
        assertEquals("C_TT_GACCATU", alignmentsSeventh[1][0]);

        assertEquals("CATTUG_CCU_U", alignmentsSeventh[2][2]);
        assertEquals("C_TT_GACCATU", alignmentsSeventh[2][0]);

        // [5,8]
        var fifthIterationData = outputData.iterationData[1][1];
        var eighthIterationData = outputData.iterationData[2][1];

        score = fifthIterationData[0] === eighthIterationData[0];
        alignmentLength = fifthIterationData[1] === eighthIterationData[1];
        lambda = fifthIterationData[2] === eighthIterationData[2];
        deletion = fifthIterationData[3] === eighthIterationData[3];
        insertion = fifthIterationData[4] === eighthIterationData[4];
        match = fifthIterationData[5] === eighthIterationData[5];
        mismatch = fifthIterationData[6] === eighthIterationData[6];
        alignmentsFifth = fifthIterationData[7];
        alignmentsEighth = eighthIterationData[7];

        assertEquals(true, score);
        assertEquals(true, alignmentLength);
        assertEquals(true, lambda);

        assertEquals(true, deletion);
        assertEquals(true, insertion);
        assertEquals(true, match);
        assertEquals(true, mismatch);

        assertEquals("CAT", alignmentsFifth[0][2]);
        assertEquals("CAT", alignmentsFifth[0][0]);

        assertEquals("CAT", alignmentsEighth[0][2]);
        assertEquals("CAT", alignmentsEighth[0][0]);

        // [6,9]
        var sixthIterationData = outputData.iterationData[1][2];
        var ninthIterationData = outputData.iterationData[2][2];

        score = sixthIterationData[0] === ninthIterationData[0];
        alignmentLength = sixthIterationData[1] === ninthIterationData[1];
        lambda = sixthIterationData[2] === ninthIterationData[2];
        deletion = sixthIterationData[3] === ninthIterationData[3];
        insertion = sixthIterationData[4] === ninthIterationData[4];
        match = sixthIterationData[5] === ninthIterationData[5];
        mismatch = sixthIterationData[6] === ninthIterationData[6];
        alignmentsSixth = sixthIterationData[7];
        alignmentsNinth = ninthIterationData[7];

        assertEquals(true, score);
        assertEquals(true, alignmentLength);
        assertEquals(true, lambda);

        assertEquals(true, deletion);
        assertEquals(true, insertion);
        assertEquals(true, match);
        assertEquals(true, mismatch);

        assertEquals("CAT", alignmentsSixth[0][2]);
        assertEquals("CAT", alignmentsSixth[0][0]);

        assertEquals("CAT", alignmentsNinth[0][2]);
        assertEquals("CAT", alignmentsNinth[0][0]);
        */
    }
});