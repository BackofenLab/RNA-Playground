/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("procedures.backtracking", Vector, Backtracking, backtrace);

    function Vector(i, j) {
        this.i = i;
        this.j = j;
        this.label = MATRICES.DEFAULT;
    }

    function create(Vector, label) {
        Vector.label = label;
        return Vector;
    }

    function Backtracking() {
        this.backtrace = backtrace;
    }

    /**
     *
     * @param pathLength - Number of edges (arrows) in the path.
     * @return {Array}
     */
    function backtrace(algorithm, path, inputData, outputData, pathLength) {
        var paths = [];

        switch (algorithm.type) {
            case ALGORITHMS.GOTOH:
                /*
                 It is based on the code of Alexander Mattheis
                 in project Algorithms for Bioninformatics.
                 */
                this.gotohTraceback = function (path, inputData, outputData, pathLength) {
                    var currentPosition = path[path.length - 1];
                    var neighboured = getMultiNeighboured(currentPosition, inputData, outputData);

                    for (var i = 0; i < neighboured.length; i++) {
                        if ((neighboured[i].i === 0 && neighboured[i].j === 0)
                            || (pathLength !== -1 && path.length >= pathLength)) {

                            path.push(neighboured[i]);
                            paths.push(path.slice());  // creating a shallow copy
                            path.pop();
                        } else {
                            path.push(neighboured[i]);
                            this.gotohTraceback(path, inputData, outputData, pathLength);
                            path.pop();
                        }
                    }
                };

                this.gotohTraceback(path, inputData, outputData, pathLength);
                return paths;

            case ALGORITHMS.NEEDLEMAN_WUNSCH:
                /*
                 It is based on the code of Alexander Mattheis
                 in project Algorithms for Bioninformatics.
                 */
                this.needlemanWunschTraceback = function (path, inputData, outputData, pathLength) {
                    var currentPosition = path[path.length - 1];
                    var neighboured = getNeighboured(currentPosition, inputData, outputData, algorithm);

                    for (var i = 0; i < neighboured.length; i++) {
                        if ((neighboured[i].i === 0 && neighboured[i].j === 0)
                            || (pathLength !== -1 && path.length >= pathLength)) {

                            path.push(neighboured[i]);
                            paths.push(path.slice());  // creating a shallow copy
                            path.pop();
                        } else {
                            path.push(neighboured[i]);
                            this.needlemanWunschTraceback(path, inputData, outputData, pathLength);
                            path.pop();
                        }
                    }
                };

                this.needlemanWunschTraceback(path, inputData, outputData, pathLength);
                return paths;

            case ALGORITHMS.SMITH_WATERMAN:
                /*
                 It is based on the code of Alexander Mattheis
                 in project Algorithms for Bioninformatics.
                 */
                this.smithWatermanTraceback = function (path, inputData, outputData, pathLength) {
                    var currentPosition = path[path.length - 1];
                    var neighboured = getNeighboured(currentPosition, inputData, outputData, algorithm);

                    for (var i = 0; i < neighboured.length; i++) {
                        if (outputData.matrix[neighboured[i].i][neighboured[i].j] === 0
                            || (pathLength !== -1 && path.length >= pathLength)) {

                            path.push(neighboured[i]);
                            paths.push(path.slice());  // creating a shallow copy
                            path.pop();
                        } else {
                            path.push(neighboured[i]);
                            this.smithWatermanTraceback(path, inputData, outputData, pathLength);
                            path.pop();
                        }
                    }
                };

                this.smithWatermanTraceback(path, inputData, outputData, pathLength);
                return paths;
        }
    }

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
}());
