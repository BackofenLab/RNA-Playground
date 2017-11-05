/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("formats.scoringFunctions", getAffineSumOfPairsScore);

    /**
     * Computes the affine sum-of-pairs score for a given alignment.
     * The Gotoh-score between every alignment is computed.
     * These pairwise scores are summed up.
     * Hint: The score for overlying gaps is zero.
     * @param inputData {Object} - Contains all input data.
     * @param alignment {Object} - The alignment for which the sum-of-pairs score is computed.
     * @return {number} - The sum-of-pairs score.
     */
    function getAffineSumOfPairsScore(inputData, alignment) {
        debugger;
        var score = 0;

        // every alignment with every other alignment
        for (var i = 1; i < alignment.length; i++) {
            for (var j = 0; j < i; j++) {

                var sequences = removeOverlyingGaps(alignment[j], alignment[i]);

                var sequenceA = sequences[0];
                var sequenceB = sequences[1];

                for (var k = 0; k < sequenceA.length; k++) {
                    if (sequenceA[k] === SYMBOLS.GAP) {
                        var gapSize = getGapSize(k, sequenceA);
                        score += inputData.baseCosts + gapSize * inputData.enlargement;
                        k += gapSize-1;
                    } else if (sequenceB[k] === SYMBOLS.GAP) {
                        var gapSize = getGapSize(k, sequenceB);
                        score += inputData.baseCosts + gapSize * inputData.enlargement;
                        k += gapSize-1;
                    } else
                        score += sequenceA[k] === sequenceB[k] ? inputData.match : inputData.mismatch;
                }
            }
        }

        return score;
    }

    /**
     * Removes overlying gaps in sequences, because they are not counted.
     * @param sequence1 {string} - The first sequence.
     * @param sequence2 {string} - The second sequence.
     * @return {[string, string]} - The alignment without overlying gaps.
     */
    function removeOverlyingGaps(sequence1, sequence2) {
        var newSequence1 = SYMBOLS.EMPTY;
        var newSequence2 = SYMBOLS.EMPTY;

        for (var i = 0; i < sequence1.length; i++) {
            if (sequence1[i] === sequence2[i] && sequence1[i] === SYMBOLS.GAP)  // if overlying gaps
                continue;

            newSequence1 += sequence1[i];
            newSequence2 += sequence2[i];
        }

        return [newSequence1, newSequence2];
    }

    /**
     * Returns the size of gap.
     * @param gapStartPosition {number} - The start position.
     * @param sequence {string} - The sequence with the gap.
     * @return {number} - The gap size.
     */
    function getGapSize(gapStartPosition, sequence) {
        var gapSize = 0;

        var i = gapStartPosition;

        while (i < sequence.length && sequence[i] === SYMBOLS.GAP) {
            gapSize++;
            i++;
        }

        return gapSize;
    }
}());
