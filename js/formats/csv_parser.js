/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("formats.csvParser", checkInput, getCSVData, getMatrix);

    /**
     * Checks with the help of a CSV-parser,
     * if the input is correct or not and returns the error output.
     * @param csvData {string} - The csv string which has to be converted into a cluster algorithm distance matrix.
     * @return {string} - The error output.
     */
    function checkInput(csvData) {

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
                string += leftString.charAt(i-1) + SYMBOLS.COMMA;

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
                matrixRow.push(Math.round(row[i]*10000)/10000);
            else
                matrixRow.push(row[i]);
        }

        return matrixRow;
    }

    /**
     * Returns the distance matrix object created by a CSV-parser.
     * @param csvData {string} - The csv string which has to be converted into a cluster algorithm distance matrix.
     * @return {Object} - The distance matrix.
     */
    function getMatrix(csvData) {
    }
}());