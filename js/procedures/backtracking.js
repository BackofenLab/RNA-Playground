/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("procedures.backtracking", Vector, Backtracking, getNeighboured, getMultiNeighboured);

    /**
     * Creates a vector with a label, which is by the default the default matrix label.
     * @param i {number} - First matrix axis.
     * @param j {number} - Second matrix axis.
     * @constructor
     */
    function Vector(i, j) {
        this.i = i;
        this.j = j;
        this.label = MATRICES.DEFAULT;
    }

    /**
     * Changes the label of a vector. It is used to create a vector with a specific label.
     * @exmaple
     * create(new Vector(i,j), MATRICES.HORIZONTAL)
     * @param Vector
     * @param label {string} - The string label from "defaults.js".
     * @return {Vector} - The vector with the changed label.
     */
    function create(Vector, label) {
        Vector.label = label;
        return Vector;
    }

    /**
     * Contains backtracking neighbouring functions shared by the alignment algorithms.
     * @constructor
     */
    function Backtracking() {
        this.getNeighboured = getNeighboured;
        this.getMultiNeighboured = getMultiNeighboured;
    }

    /**
     * Returns the neighbours to which you can go from the current cell position used as input.
     * @param position {Vector} - Current cell position in matrix.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     * @param algorithm {Object} - Contains an alignment algorithm.
     * @return {Array} - Contains neighboured positions as Vector-objects.
     */
    function getNeighboured(position, inputData, outputData, algorithm) {
        var neighboured = [];

        var left = position.j - 1;
        var up = position.i - 1;

        // retrieve values
        var aChar = left >= 0 ? inputData.sequenceA[left] : SYMBOLS.EMPTY;
        var bChar = up >= 0 ? inputData.sequenceB[up] : SYMBOLS.EMPTY;

        var currentValue = outputData.matrix[position.i][position.j];

        var matchOrMismatch = aChar === bChar ? inputData.match : inputData.mismatch;

        var diagonalValue = left >= 0 && up >= 0 ? outputData.matrix[up][left] : NaN;
        var upValue = up >= 0 ? outputData.matrix[up][position.j] : NaN;
        var leftValue = left >= 0 ? outputData.matrix[position.i][left] : NaN;

        // check
        var isMatchMismatch = currentValue === (diagonalValue + matchOrMismatch);
        var isDeletion = currentValue === (upValue + inputData.deletion);
        var isInsertion = currentValue === (leftValue + inputData.insertion);

        if (algorithm.type === ALGORITHMS.SMITH_WATERMAN) {
            isMatchMismatch = isMatchMismatch || currentValue === 0 && up >= 0 && left >= 0;
            isDeletion = isDeletion || currentValue === 0 && up >= 0;
            isInsertion = isInsertion || currentValue === 0 && left >= 0;
        }

        // add
        if (isMatchMismatch)
            neighboured.push(new Vector(up, left));

        if (isDeletion)
            neighboured.push(new Vector(up, position.j));

        if (isInsertion)
            neighboured.push(new Vector(position.i, left));

        return neighboured;
    }

    /**
     * Returns the neighbours to which you can go from the current cell position used as input.
     * @param position {Vector} - Current cell position in matrix.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     * @return {Array} - Contains neighboured positions as Vector-objects.
     */
    function getMultiNeighboured(position, inputData, outputData) {
        debugger;
        var neighboured = [];

        if (position.label === MATRICES.VERTICAL)
            return getVerticalNeighboured(position, inputData, outputData);
        else if (position.label === MATRICES.HORIZONTAL)
            return getHorizontalNeighboured(position, inputData, outputData);

        var left = position.j - 1;
        var up = position.i - 1;

        // retrieve values
        var aChar = left >= 0 ? inputData.sequenceA[left] : SYMBOLS.EMPTY;
        var bChar = up >= 0 ? inputData.sequenceB[up] : SYMBOLS.EMPTY;

        var currentValue = outputData.matrix[position.i][position.j];

        var matchOrMismatch = aChar === bChar ? inputData.match : inputData.mismatch;

        var diagonalValue = left >= 0 && up >= 0 ? outputData.matrix[up][left] : NaN;
        var verticalValue = up >= 0 ? outputData.verticalGaps[position.i][position.j] : NaN;
        var horizontalValue = left >= 0 ? outputData.horizontalGaps[position.i][position.j] : NaN;

        var upValue = up >= 0 && position.j === 0 ? outputData.matrix[up][position.j] : Number.NaN;
        var leftValue = left >= 0 && position.i === 0 ? outputData.matrix[position.i][left] : Number.NaN;

        // check
        var isMatchMismatch = currentValue === (diagonalValue + matchOrMismatch);
        var isChangeToP = currentValue === verticalValue;
        var isChangeToQ = currentValue === horizontalValue;

        var isDeletion = currentValue === upValue + inputData.enlargement;
        var isInsertion = currentValue === leftValue + inputData.enlargement;

        // add
        if (isMatchMismatch)
            neighboured.push(create(new Vector(up, left), MATRICES.DEFAULT));

        if (isChangeToP)
            neighboured.push(create(new Vector(position.i, position.j), MATRICES.VERTICAL));

        if (isChangeToQ)
            neighboured.push(create(new Vector(position.i, position.j), MATRICES.HORIZONTAL));

        if (isInsertion)
            neighboured.push(create(new Vector(position.i, left), MATRICES.DEFAULT));

        if (isDeletion)
            neighboured.push(create(new Vector(up, position.j), MATRICES.DEFAULT));

        if (!(isMatchMismatch || isChangeToP || isChangeToQ || isInsertion || isDeletion)
            && (position.i !== 0 || position.j !== 0))
            neighboured.push(create(new Vector(0, 0), MATRICES.DEFAULT));

        return neighboured;
    }

    /**
     * Returns the neighbours to which you can go from the current cell position in the matrix for vertical gap costs.
     * @param position {Vector} - Current cell position in matrix.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     * @return {Array} - Contains neighboured positions as Vector-objects.
     */
    function getVerticalNeighboured(position, inputData, outputData) {
        var neighboured = [];

        var up = position.i - 1;

        // retrieve values
        var currentValue = outputData.verticalGaps[position.i][position.j];

        var pUpValue = Number.NaN;
        var xUpValue = Number.NaN;

        if (position.i >= 0 && up >= 0) {
            pUpValue = outputData.verticalGaps[up][position.j];
            xUpValue = outputData.matrix[up][position.j];
        }

        // check
        var isUpInP = currentValue === pUpValue + inputData.enlargement;
        var isUpInX = currentValue === xUpValue + inputData.baseCosts + inputData.enlargement;

        // add
        if (isUpInP)
            neighboured.push(create(new Vector(up, position.j), MATRICES.VERTICAL));

        if (isUpInX)
            neighboured.push(create(new Vector(up, position.j), MATRICES.DEFAULT));

        return neighboured;
    }

    /**
     * Returns the neighbours to which you can go from the current cell position in the matrix for horizontal gap costs.
     * @param position {Vector} - Current cell position in matrix.
     * @param inputData {Object} - Contains all input data.
     * @param outputData {Object} - Contains all output data.
     * @return {Array} - Contains neighboured positions as Vector-objects.
     */
    function getHorizontalNeighboured(position, inputData, outputData) {
        var neighboured = [];

        var left = position.j - 1;

        // retrieve values
        var currentValue = outputData.horizontalGaps[position.i][position.j];

        var qLeftValue = Number.NaN;
        var xLeftValue = Number.NaN;

        if (position.i >= 0 && left >= 0) {
            qLeftValue = outputData.horizontalGaps[position.i][left];
            xLeftValue = outputData.matrix[position.i][left];
        }

        // check
        var isLeftInQ = currentValue === qLeftValue + inputData.enlargement;
        var isLeftInX = currentValue === xLeftValue + inputData.baseCosts + inputData.enlargement;

        // add
        if (isLeftInQ)
            neighboured.push(create(new Vector(position.i, left), MATRICES.HORIZONTAL));

        if (isLeftInX)
            neighboured.push(create(new Vector(position.i, left), MATRICES.DEFAULT));

        return neighboured;
    }
}());
