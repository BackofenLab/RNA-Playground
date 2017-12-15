/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("formats.csvParser", CsvParser);

    /**
     * Creates a "comma-separated values" parser.
     * @constructor
     */
    function CsvParser() {
        // public class methods
        this.checkInput = checkInput;
        this.getCSVData = getCSVData;
        this.getMatrix = getMatrix;
        this.getNumOfTableLines = getNumOfTableLines;
    }

    /**
     * Checks if the number of columns equals the number of non-empty rows
     * and checks if the number of column separators is the same in every line.
     * If not, it returns the error output.
     * @param csvData {string} - The csv string which has to be converted into a cluster algorithm distance matrix.
     * @return {string} - The error output, if it exists.
     */
    function checkInput(csvData) {
        var lines = getNonEmptyLines(csvData.split(SYMBOLS.NEW_LINE));

        var lastNumberOfSeparators = -1;

        for (var i = 0; i < lines.length; i++) {
            var line = lines[i];

            if (line.replace(MULTI_SYMBOLS.SPACE, SYMBOLS.EMPTY) !== SYMBOLS.EMPTY) {  // ignore empty lines
                var numberOfSeparators = (line.match(MULTI_SYMBOLS.SEPARATORS) || []).length;

                if (lastNumberOfSeparators !== -1 && lastNumberOfSeparators !== numberOfSeparators)
                    return ERRORS.WRONG_NUMBER_OF_COLUMNS_IN_ROW + (i + 1);  // "+1" because humans start counting with 1

                lastNumberOfSeparators = numberOfSeparators;
            }
        }

        if (lines.length !== lastNumberOfSeparators + 1)  // #rows !== #columns
            return ERRORS.DIFFERENT_NUMBER_OF_COLUMNS_AND_ROWS;

        return SYMBOLS.EMPTY;
    }

    /**
     * Removes empty lines.
     * @param lines {Array} - The
     */
    function getNonEmptyLines(lines) {
        var nonEmptyLines = [];

        for (var i = 0; i < lines.length; i++) {
            if (lines[i].replace(MULTI_SYMBOLS.SPACE, SYMBOLS.EMPTY) !== SYMBOLS.EMPTY)
                nonEmptyLines.push(lines[i]);
        }

        return nonEmptyLines;
    }

    /**
     * Returns the CSV string of the given matrix.
     * @param matrix {Array} - The two-dimensional matrix.
     * @param leftString {String} - The symbols string on the left side of the matrix (for example characters of NW matrix).
     * @return {string} - The CSV-string.
     */
    function getCSVData(matrix, leftString) {
        var string = SYMBOLS.EMPTY;

        for (var i = 0; i < matrix.length; i++) {
            if (i === 0)
                string += SYMBOLS.COMMA;
            else
                string += leftString.charAt(i - 1) + SYMBOLS.COMMA;

            string += round(matrix[i]) + SYMBOLS.NEW_LINE;  // Hint: it is allowed to have a line break in the last line
        }

        return string;
    }

    /**
     * Rounds values to four decimal places if it is possible.
     * @param row {number} - The row of which values are rounded.
     * @return {Array} - Row with rounded values and original values.
     */
    function round(row) {
        var matrixRow = [];

        for (var i = 0; i < row.length; i++) {
            if (typeof row[i] === "number")
                matrixRow.push(Math.round(row[i] * 10000) / 10000);
            else
                matrixRow.push(row[i]);
        }

        return matrixRow;
    }

    /**
     * Returns a matrix object.
     * @param csvData {string} - The csv string which has to be converted into a matrix.
     * @param getNames {Function} - The function which retrieves the names for the columns/rows.
     * @return {Object} - The matrix.
     */
    function getMatrix(csvData, getNames) {
        var distances = [];
        var lines = getNonEmptyLines(csvData.split(SYMBOLS.NEW_LINE));

        // read in all distances
        for (var i = 0; i < lines.length; i++) {
            var line = lines[i];
            var lineEntries = line.split(SYMBOLS.SEMICOLON);

            for (var j = i + 1; j < lineEntries.length; j++) {  // "+1" because diagonal an every below it not saved
                var lineEntry = lineEntries[j];

                distances.push(Number(lineEntry));
            }
        }

        // create a distance matrix object
        var clusterNames = getNames(lines.length);
        var distanceMatrix = {};

        var k = 0;  // current distance value index
        for (var i = 0; i < clusterNames.length; i++) {
            var clusterName1 = clusterNames[i];

            for (var j = i + 1; j < clusterNames.length; j++) {  // "+1" because diagonal an every below it not saved
                var clusterName2 = clusterNames[j];

                distanceMatrix[[clusterName1, clusterName2]] = distances[k++];
            }
        }

        return distanceMatrix;
    }

    /**
     * Returns the number of columns.
     * @param csvData {string} - The csv string which has to be converted into a matrix.
     * @return {Object} - The matrix.
     */
    function getNumOfTableLines(csvData) {
        return getNonEmptyLines(csvData.split(SYMBOLS.NEW_LINE)).length;
    }
}());