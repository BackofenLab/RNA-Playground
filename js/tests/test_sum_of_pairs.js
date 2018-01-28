/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

TestCase("test_sum_of_pairs", {
    /**
     * Gotoh alignment to test correctness.
     * Short sequences test.
     * @see Test sequences are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_1": function () {
        var inputData = {};
        var alignment1 = ["TGGA", "__GG"];
        var alignment2 = ["TGGA", "GG__"];

        inputData.baseCosts = -3;
        inputData.enlargement = -1;
        inputData.match = 0;
        inputData.mismatch = -1;

        var multiSequenceAlignment = new bases.multiSequenceAlignment.MultiSequenceAlignment(undefined);
        var score1 = multiSequenceAlignment.getAffineSumOfPairsScore(inputData, alignment1);
        var score2 = multiSequenceAlignment.getAffineSumOfPairsScore(inputData, alignment2);
        assertEquals(-6, score1);
        assertEquals(-6, score2);
    },

    /**
     * Gotoh alignment to test correctness.
     */
    "test_2": function () {
        var inputData = {};
        var alignment1 = ["CCGA", "C__G"];
        var alignment2 = ["CCGA", "CG__"];

        inputData.baseCosts = -3;
        inputData.enlargement = -1;
        inputData.match = 1;
        inputData.mismatch = -1;

        var multiSequenceAlignment = new bases.multiSequenceAlignment.MultiSequenceAlignment(undefined);
        var score1 = multiSequenceAlignment.getAffineSumOfPairsScore(inputData, alignment1);
        var score2 = multiSequenceAlignment.getAffineSumOfPairsScore(inputData, alignment2);
        assertEquals(-5, score1);
        assertEquals(-5, score2);
    },

    /**
     * Gotoh alignment to test correctness.
     * @see Test sequences are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_3": function () {
        var inputData = {};
        var alignment1 = ["TACGCAGA", "T___CCGA"];
        var alignment2 = ["TACGCAGA", "TCC___GA"];
        var alignment3 = ["TACGCAGA", "TCCG___A"];

        inputData.baseCosts = -4;
        inputData.enlargement = -1;
        inputData.match = 1;
        inputData.mismatch = 0;

        var multiSequenceAlignment = new bases.multiSequenceAlignment.MultiSequenceAlignment(undefined);
        var score1 = multiSequenceAlignment.getAffineSumOfPairsScore(inputData, alignment1);
        var score2 = multiSequenceAlignment.getAffineSumOfPairsScore(inputData, alignment2);
        var score3 = multiSequenceAlignment.getAffineSumOfPairsScore(inputData, alignment3);
        assertEquals(-3, score1);
        assertEquals(-3, score2);
        assertEquals(-3, score3);
    },

    /**
     * Needleman-Wunsch alignment to test correctness.
     * @see Test sequences are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_4": function () {
        var inputData = {};
        var alignment1 = ["CCCCGCGACTCGGGTTCAAGGG", "GGGTGAGACCCCAGTTCAACCC"];

        inputData.baseCosts = 0;
        inputData.enlargement = -2;
        inputData.match = 4;
        inputData.mismatch = -1;

        var multiSequenceAlignment = new bases.multiSequenceAlignment.MultiSequenceAlignment(undefined);
        var score1 = multiSequenceAlignment.getAffineSumOfPairsScore(inputData, alignment1);
        assertEquals(33, score1);
    },

    /**
     * Needleman-Wunsch alignment to test correctness..
     * @see Test sequences are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_5": function () {
        var inputData = {};
        var alignment1 = ["T_C_CGA", "TACGCGC"];

        inputData.baseCosts = 0;
        inputData.enlargement = -1;
        inputData.match = 1;
        inputData.mismatch = 0;

        var multiSequenceAlignment = new bases.multiSequenceAlignment.MultiSequenceAlignment(undefined);
        var score1 = multiSequenceAlignment.getAffineSumOfPairsScore(inputData, alignment1);
        assertEquals(2, score1);
    },

    /**
     * Needleman-Wunsch alignment to test correctness.
     * @see Test sequences are taken from project Algorithms for Bioninformatics of Alexander Mattheis.
     */
    "test_6": function () {
        var inputData = {};
        var alignment1 = ["AA_CG", "AATCG"];

        inputData.baseCosts = 0;
        inputData.enlargement = -2;
        inputData.match = 1;
        inputData.mismatch = -1;

        var multiSequenceAlignment = new bases.multiSequenceAlignment.MultiSequenceAlignment(undefined);
        var score1 = multiSequenceAlignment.getAffineSumOfPairsScore(inputData, alignment1);
        assertEquals(2, score1);
    }
});