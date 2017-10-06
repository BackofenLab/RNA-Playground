/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("postProcessing.visualizer", Visualizer,
        shareInformation, showFlow,
        showTraceback, highlight, downloadTable, replaceInfinityStrings, redrawAllLines, removeAllContents);

    // instances
    var visualizerInstance;

    /**
     * Contains functions for visualization
     * and helper functions which easify the access on information.
     * @constructor
     */
    function Visualizer() {
        visualizerInstance = this;

        // variables
        this.algorithm = {};

        this.cellLines = [];
        this.container = document.getElementById("overlay");

        this.lastFlows = [];
        this.lastPath = [];
        this.lastRowNumber = -1;

        this.input = {};
        this.output = {};

        this.svg = createSVG();

        // bindings
        ko.bindingHandlers.drawChar = {
            update: function (element, valueAccessor) {
                var values = ko.unwrap(valueAccessor());
                var character = values[0];
                var index = values[1];

                if (character !== undefined)
                    element.innerHTML = character.toUpperCase() + SUB.START_TAG + index + SUB.END_TAG;
            }
        };

        // public methods (linking)
        this.shareInformation = shareInformation;
        this.showFlow = showFlow;
        this.showTraceback = showTraceback;
        this.highlight = highlight;
        this.downloadTable = downloadTable;
        this.replaceInfinityStrings = replaceInfinityStrings;
        this.redrawAllLines = redrawAllLines;
        this.removeAllContents = removeAllContents;
    }

    function createSVG() {
        var svg = document.createElementNS(SVG.NAME_SPACE, "svg");
        svg.appendChild(createEndMarker(SVG.TRACEBACK_LONG_ARROW_COLOR, SVG.MARKER.ID_TRACEBACK));
        svg.appendChild(createEndMarker(SVG.FLOW_LONG_ARROW_COLOR, SVG.MARKER.ID_FLOW));
        return svg;
    }

    /* created with the help of <marker> definition https://developer.mozilla.org/de/docs/Web/SVG/Element/marker */
    function createEndMarker(color, id) {
        // create triangle
        var trianglePath = document.createElementNS(SVG.NAME_SPACE, "path");
        trianglePath.setAttribute("d", SVG.TRIANGLE.D);  // M: move to, L: line to
        trianglePath.setAttribute("fill", color);

        // create marker object using path defined above
        var markerEnd = document.createElementNS(SVG.NAME_SPACE, "marker");
        markerEnd.setAttribute("id", id);
        markerEnd.setAttribute("orient", SVG.MARKER.ORIENT);
        markerEnd.setAttribute("refX", SVG.MARKER.BOUNDS.REF_X);
        markerEnd.setAttribute("refY", SVG.MARKER.BOUNDS.REF_Y);
        markerEnd.setAttribute("markerWidth", SVG.MARKER.BOUNDS.WIDTH);
        markerEnd.setAttribute("markerHeight", SVG.MARKER.BOUNDS.HEIGHT);
        markerEnd.setAttribute("viewBox", SVG.MARKER.VIEW_BOX);
        markerEnd.appendChild(trianglePath);

        return markerEnd;
    }

    function shareInformation(algorithm, input, output) {
        visualizerInstance.algorithm = algorithm;
        visualizerInstance.input = input;
        visualizerInstance.output = output;
    }

    function showFlow(cellCoordinates, calculationVerticalTable, table, calculationHorizontalTable) {
        var flows = visualizerInstance.algorithm.getTraces([cellCoordinates], visualizerInstance.input, visualizerInstance.output, 1);

        for (i = 0; i < visualizerInstance.lastFlows.length; i++)
            demarkCells(visualizerInstance.lastFlows[i], calculationVerticalTable, table, calculationHorizontalTable, i, true);

        for (var i = 0; i < flows.length; i++)
            markCells(flows[i].reverse(), calculationVerticalTable, table, calculationHorizontalTable, i, true, true);

        visualizerInstance.lastFlows = flows;
    }

    function demarkCells(path, calculationVerticalTable, table, calculationHorizontalTable, colorClass, flowMode) {
        flowMode = flowMode || false;

        var currentTable;

        if (path.length > 0) {
            for (var j = 0; j < path.length; j++) {
                currentTable = getRightTable(path, j, calculationVerticalTable, table, calculationHorizontalTable);

                var posI = path[j].i + 1;
                var posJ = path[j].j + 1;

                if (currentTable.rows[posI].cells[posJ] !== undefined) {  // if table has shrinked

                    switch (colorClass) {
                        case -1:
                            currentTable.rows[posI].cells[posJ].classList.remove("selected");
                            break;
                        default:
                            removeColors(currentTable, posI, posJ);
                    }
                }

                removeArrow(currentTable, posI, posJ);
            }
        }

        removeAllLines();  // below last flows/paths redrawn

        if (flowMode)  // redraw last traceback
            markCells(visualizerInstance.lastPath, calculationVerticalTable, table, calculationHorizontalTable, -1, true, false);
        else {  // redraw last flow
            var lastFlows = visualizerInstance.lastFlows;
            for (var i = 0; i < lastFlows.length; i++)
                markCells(lastFlows[i], calculationVerticalTable, table, calculationHorizontalTable, i, true, true);
        }
    }

    function getRightTable(path, j, calculationVerticalTable, table, calculationHorizontalTable) {
        if (path[j].label === MATRICES.VERTICAL)
            return calculationVerticalTable;
        else if (path[j].label === MATRICES.DEFAULT)
            return table;
        // else if (path[j].label === MATRICES.HORIZONTAL)
        return calculationHorizontalTable;
    }

    function removeColors(table, posI, posJ) {
        table.rows[posI].cells[posJ].classList.remove("selected_light_red");
        table.rows[posI].cells[posJ].classList.remove("selected_very_light_red");
        table.rows[posI].cells[posJ].classList.remove("selected_red");
        table.rows[posI].cells[posJ].classList.remove("selected_green");
    }

    function removeArrow(table, i, j) {
        var cell = table.rows[i].cells[j];
        var numChildren = cell.children.length;

        for (var k = 0; k < numChildren; k++)
            cell.children[0].outerHTML = SYMBOLS.EMPTY;
    }

    function removeAllLines() {
        var line;
        while ((line = visualizerInstance.cellLines.pop()) !== undefined) {
            if (visualizerInstance.svg.contains(line))
                visualizerInstance.svg.removeChild(line);
        }
    }

    function markCells(path, calculationVerticalTable, table, calculationHorizontalTable, colorClass, arrows, flowMode) {
        arrows = arrows || false;

        var lastPosI;
        var lastPosJ;
        var lastTable;

        var currentTable;

        for (var j = 0; j < path.length; j++) {
            currentTable = getRightTable(path, j, calculationVerticalTable, table, calculationHorizontalTable);

            var posI = path[j].i + 1;
            var posJ = path[j].j + 1;

            switch (colorClass) {
                case 0:
                    currentTable.rows[posI].cells[posJ].classList.add("selected_light_red");
                    break;
                case 1:
                    currentTable.rows[posI].cells[posJ].classList.add("selected_very_light_red");
                    break;
                case 2:
                    currentTable.rows[posI].cells[posJ].classList.add("selected_red");
                    break;
                default:
                    currentTable.rows[posI].cells[posJ].classList.add("selected");
            }

            if (j === path.length - 1 && colorClass !== -1) {
                removeColors(currentTable, posI, posJ);
                currentTable.rows[posI].cells[posJ].classList.add("selected_green");
            }

            if (arrows) {
                placeArrow(currentTable, posI, posJ, lastTable, lastPosI, lastPosJ, flowMode);
                lastPosI = posI;
                lastPosJ = posJ;
                lastTable = currentTable;
            }
        }
    }

    function placeArrow(table, i, j, lastTable, lastI, lastJ, flowMode) {
        if (lastI !== undefined && lastI !== undefined) {
            var cell = table.rows[i].cells[j];
            var lastCell = lastTable.rows[lastI].cells[lastJ];

            if (lastTable === table) {
                var isTop = i - lastI > 0;
                var isLeft = j - lastJ > 0;

                if (isTop && isLeft) {
                    if ($(cell).find(ARROWS.DIAGONAL_NAME).length !== 1)
                        $(cell).append(ARROWS.DIAGONAL);
                }
                else if (isLeft) {
                    if ($(cell).find(ARROWS.LEFT_NAME).length !== 1)
                        $(cell).append(ARROWS.LEFT);
                }
                else if (isTop) {
                    if ($(cell).find(ARROWS.TOP_NAME).length !== 1)
                        $(cell).append(ARROWS.TOP);
                }
            } else if (lastTable !== table) {
                var parentMatrix = getParentMatrix(cell);
                var lastParentMatrix = getParentMatrix(lastCell);

                var isPtoX = parentMatrix === MATRICES.VERTICAL && lastParentMatrix === MATRICES.DEFAULT;
                var isQtoX = parentMatrix === MATRICES.HORIZONTAL && lastParentMatrix === MATRICES.DEFAULT;
                var isXtoP = parentMatrix === MATRICES.DEFAULT && lastParentMatrix === MATRICES.VERTICAL;
                var isXtoQ = parentMatrix === MATRICES.DEFAULT && lastParentMatrix === MATRICES.HORIZONTAL;

                if (isPtoX)
                    drawLine(cell, lastCell, MOVE.P_TO_X, i, j, flowMode);
                else if (isQtoX)
                    drawLine(cell, lastCell, MOVE.Q_TO_X, i, j, flowMode);
                else if (isXtoP)
                    drawLine(cell, lastCell, MOVE.X_TO_P, i, j, flowMode);
                else if (isXtoQ)
                    drawLine(cell, lastCell, MOVE.X_TO_Q, i, j, flowMode);
            }
        }
    }

    function getParentMatrix(cell) {
        if (cell.parentNode.parentNode.parentNode.id === "calculation_horizontal")
            return MATRICES.HORIZONTAL;
        else if (cell.parentNode.parentNode.parentNode.id === "calculation_vertical")
            return MATRICES.VERTICAL;
        else if (cell.parentNode.parentNode.parentNode.id === "calculation")
            return MATRICES.DEFAULT;
    }
    
    function drawLine(cell, lastCell, move, i, j, flowMode) {
        var cellHeight = cell.offsetHeight;
        var cellWidth = cell.offsetWidth;

        var left;
        var top;
        var lastLeft;
        var lastTop;

        if (move === MOVE.X_TO_P) {
            left = (cell.offsetLeft + cellWidth * CELL_PERCENT).toString();
            top = (cell.offsetTop + cellHeight * CELL_PERCENT).toString();
            lastLeft = (lastCell.offsetLeft + cellWidth * CELL_PERCENT).toString();
            lastTop = (lastCell.offsetTop + cellHeight * (1 - CELL_PERCENT)).toString();
        } else if (move === MOVE.X_TO_Q) {
            left        =   (cell.offsetLeft        + cellWidth     * CELL_PERCENT).toString();
            top         =   (cell.offsetTop         + cellHeight    * (1-CELL_PERCENT)).toString();
            lastLeft    =   (lastCell.offsetLeft    + cellWidth     * CELL_PERCENT).toString();
            lastTop     =   (lastCell.offsetTop     + cellHeight    * CELL_PERCENT).toString();
        } else if (move === MOVE.Q_TO_X) {
            left        =   (cell.offsetLeft        + cellWidth     * CELL_PERCENT).toString();
            top         =   (cell.offsetTop         + cellHeight    * CELL_PERCENT).toString();
            lastLeft    =   (lastCell.offsetLeft    + cellWidth     * (1-CELL_PERCENT)).toString();
            lastTop     =   (lastCell.offsetTop     + cellHeight    * (1-CELL_PERCENT)).toString();
        } else if (move === MOVE.P_TO_X) {
            left        =   (cell.offsetLeft        + cellWidth     * (1-CELL_PERCENT)).toString();
            top         =   (cell.offsetTop         + cellHeight    * (1-CELL_PERCENT)).toString();
            lastLeft    =   (lastCell.offsetLeft    + cellWidth     * (1-CELL_PERCENT)).toString();
            lastTop     =   (lastCell.offsetTop     + cellHeight    * CELL_PERCENT).toString();
        }

        // define svg dimensions
        visualizerInstance.svg.setAttribute("width", document.body.offsetWidth);
        visualizerInstance.svg.setAttribute("height", document.body.offsetHeight);

        // create line with previously defined marker
        var line = document.createElementNS(SVG.NAME_SPACE, "line");
        if (flowMode) {
            line.setAttribute("marker-end", SVG.MARKER.URL_FLOW);
            line.setAttribute("stroke", SVG.FLOW_LONG_ARROW_COLOR);
        }
        else {
            line.setAttribute("marker-end", SVG.MARKER.URL_TRACEBACK);
            line.setAttribute("stroke", SVG.TRACEBACK_LONG_ARROW_COLOR);
        }
        line.setAttribute("x1", left);
        line.setAttribute("y1", top);
        line.setAttribute("x2", lastLeft);
        line.setAttribute("y2", lastTop);
        visualizerInstance.svg.appendChild(line);
        visualizerInstance.cellLines.push(line);

        if (visualizerInstance.container.childElementCount !== 1)
            visualizerInstance.container.appendChild(visualizerInstance.svg);
    }

    function showTraceback(traceNumber, calculationVerticalTable, calculationTable, calculationHorizontalTable) {
        debugger;
        var path = visualizerInstance.output.tracebackPaths[traceNumber];

        if (visualizerInstance.lastPath.length > 0) {
            var posI = visualizerInstance.lastPath[0].i + 1;
            var posJ = visualizerInstance.lastPath[0].j + 1;
            var tableCell = calculationTable.rows[posI].cells[posJ];

            if (path === visualizerInstance.lastPath
                && tableCell !== undefined
                && tableCell.classList.contains("selected")) {  // case: same path
                demarkCells(visualizerInstance.lastPath, calculationVerticalTable, calculationTable, calculationHorizontalTable, -1, false);
                visualizerInstance.lastPath = [];
            } else {  // case: different path
                demarkCells(visualizerInstance.lastPath, calculationVerticalTable, calculationTable, calculationHorizontalTable, -1, false);
                markCells(path, calculationVerticalTable, calculationTable, calculationHorizontalTable, -1, true, false);
                visualizerInstance.lastPath = path;
            }
        } else {  // case: first time selected
            markCells(path, calculationVerticalTable, calculationTable, calculationHorizontalTable, -1, true, false);
            visualizerInstance.lastPath = path;
        }
    }

    function highlight(rowNumber, table) {
        var start = 0;

        var cell = table.rows[rowNumber + start].cells[0];
        var lastCell = visualizerInstance.lastRowNumber >= 0
            ? table.rows[visualizerInstance.lastRowNumber + start].cells[0] : undefined;

        if (rowNumber === visualizerInstance.lastRowNumber
            && lastCell.classList.contains("selected")) {  // case: same cell
            lastCell.classList.remove("selected");
        } else if (lastCell !== undefined) {  // case: different cell
            lastCell.classList.remove("selected");
            cell.classList.add("selected");
        } else {  // case: first time clicked in table
            cell.classList.add("selected");
        }

        visualizerInstance.lastRowNumber = rowNumber;
    }

    /*
    Hint: " e.data.number" allows to distinguish between the different tables of an algorithm.
     */
    function downloadTable(e) {
        var number = e.data.number;

        var matrix = getMatrix(number);
        var upperString = visualizerInstance.input.sequenceA;
        var leftString = visualizerInstance.input.sequenceB;

        var tableCSV = tableToCSV(number, matrix, upperString, leftString);
        var tableFile = new File([tableCSV], {type: TABLE.TEXT_FILE_ENCODING});

        saveAs(tableFile, TABLE.DOWNLOAD_NAME);
    }

    function getMatrix(number) {
        switch (number) {
            case 0:
                return replaceInfinities(visualizerInstance.output.verticalGaps);
            case 2:
                return replaceInfinities(visualizerInstance.output.horizontalGaps);
        }

        return visualizerInstance.output.matrix;
    }

    function replaceInfinities(matrix) {
        for (var i = 0; i < matrix.length; i++) {
            for (var j = 0; j < matrix[0].length; j++) {
                if (matrix[i][j] === LATEX.NEGATIVE_INFINITY)
                    matrix[i][j] = SYMBOLS.INFINITY;
                else if (matrix[i][j] === LATEX.POSITIVE_INFINITY)
                    matrix[i][j] = SYMBOLS.NEGATIVE_INFINITY;
            }
        }

        return matrix;
    }

    /*
    CSV specification: https://www.ietf.org/rfc/rfc4180.txt
    Hint: It is allowed to have a line break in the last line.
     */
    function tableToCSV(number, matrix, upperString, leftString) {
        var string = SYMBOLS.EMPTY;

        switch (number) {
            case 0:
                string += MATRICES.VERTICAL + SYMBOLS.COMMA;
                break;
            case 1:
                string += MATRICES.DEFAULT + SYMBOLS.COMMA;
                break;
            case 2:
                string += MATRICES.HORIZONTAL + SYMBOLS.COMMA;
                break;
        };

        string += SYMBOLS.COMMA + upperString.split(SYMBOLS.EMPTY).toString() + SYMBOLS.NEW_LINE;

        for (var i = 0; i < matrix.length; i++) {
            if (i == 0)
                string += SYMBOLS.COMMA;
            else
                string += leftString.charAt(i-1) + SYMBOLS.COMMA;

            string += matrix[i] + SYMBOLS.NEW_LINE;
        }

        return string;
    }

    function replaceInfinityStrings(matrix) {
        for (var i = 0; i < matrix.length; i++) {
            for (var j = 0; j < matrix[0].length; j++) {
                if (matrix[i][j] === Number.POSITIVE_INFINITY)
                    matrix[i][j] = LATEX.POSITIVE_INFINITY;
                else if (matrix[i][j] === Number.NEGATIVE_INFINITY)
                    matrix[i][j] = LATEX.NEGATIVE_INFINITY;
            }
        }

        return matrix;
    }

    function redrawAllLines(e) {
        var calculationVerticalTable;
        var calculation = e.data.calculationTable[0];
        var calculationHorizontalTable;

        if (e.data.calculationVerticalTable !== undefined) {
            calculationVerticalTable = e.data.calculationVerticalTable[0];
            calculationHorizontalTable = e.data.calculationHorizontalTable[0];
        }

        removeAllLines();
        drawAllLines(calculationVerticalTable, calculation, calculationHorizontalTable);
    }

    function drawAllLines(calculationVerticalTable, table, calculationHorizontalTable) {
        drawArrowLine(visualizerInstance.lastPath, calculationVerticalTable, table, calculationHorizontalTable, false);

        var lastFlows = visualizerInstance.lastFlows;
        for (var i = 0; i < lastFlows.length; i++)
            drawArrowLine(lastFlows[i], calculationVerticalTable, table, calculationHorizontalTable, true);
    }

    function drawArrowLine(path, calculationVerticalTable, table, calculationHorizontalTable, flowMode) {
        var lastPosI;
        var lastPosJ;
        var lastTable;

        var currentTable;

        for (var j = 0; j < path.length; j++) {
            currentTable = getRightTable(path, j, calculationVerticalTable, table, calculationHorizontalTable);

            var posI = path[j].i + 1;
            var posJ = path[j].j + 1;

            placeArrow(currentTable, posI, posJ, lastTable, lastPosI, lastPosJ, flowMode);
            lastPosI = posI;
            lastPosJ = posJ;
            lastTable = currentTable;
        }
    }

    function removeAllContents() {
        debugger;
        removeAllLines();

        visualizerInstance.algorithm = {};

        visualizerInstance.cellLines = [];

        visualizerInstance.lastFlows = [];
        visualizerInstance.lastPath = [];
        visualizerInstance.lastRowNumber = -1;

        visualizerInstance.input = {};
        visualizerInstance.output = {};
    }
}());