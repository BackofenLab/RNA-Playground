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

        this.svg = createSVG(createEndMarker());

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

    function createSVG(endMarker) {
        var svg = document.createElementNS(SVG.NAME_SPACE, "svg");
        svg.appendChild(endMarker);
        return svg;
    }

    /* created with the help of <marker> definition https://developer.mozilla.org/de/docs/Web/SVG/Element/marker */
    function createEndMarker() {
        // create triangle
        var trianglePath = document.createElementNS(SVG.NAME_SPACE, "path");
        trianglePath.setAttribute("d", SVG.TRIANGLE.D);  // M: move to, L: line to
        trianglePath.setAttribute("fill", SVG.OBJECT_COLOR);

        // create marker object using path defined above
        var markerEnd = document.createElementNS(SVG.NAME_SPACE, "marker");
        markerEnd.setAttribute("id", SVG.MARKER.ID);
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
        var flows = new procedures.backtracking.Backtracking().backtrace(
            visualizerInstance.algorithm, [cellCoordinates], visualizerInstance.input, visualizerInstance.output, 1);

        for (i = 0; i < visualizerInstance.lastFlows.length; i++)
            demarkCells(visualizerInstance.lastFlows[i], calculationVerticalTable, table, calculationHorizontalTable, i, true);

        for (var i = 0; i < flows.length; i++)
            markCells(flows[i].reverse(), calculationVerticalTable, table, calculationHorizontalTable, i, true);

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
            markCells(visualizerInstance.lastPath, calculationVerticalTable, table, calculationHorizontalTable, -1, true);
        else {  // redraw last flow
            var lastFlows = visualizerInstance.lastFlows;
            for (var i = 0; i < lastFlows.length; i++)
                markCells(lastFlows[i], calculationVerticalTable, table, calculationHorizontalTable, i, true);
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

    function markCells(path, calculationVerticalTable, table, calculationHorizontalTable, colorClass, arrows) {
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
                placeArrow(currentTable, posI, posJ, lastTable, lastPosI, lastPosJ);
                lastPosI = posI;
                lastPosJ = posJ;
                lastTable = currentTable;
            }
        }
    }

    function placeArrow(table, i, j, lastTable, lastI, lastJ) {
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
                    drawLine(cell, lastCell, MOVE.P_TO_X, i, j);
                else if (isQtoX)
                    drawLine(cell, lastCell, MOVE.Q_TO_X, i, j);
                else if (isXtoP)
                    drawLine(cell, lastCell, MOVE.X_TO_P, i, j);
                else if (isXtoQ)
                    drawLine(cell, lastCell, MOVE.X_TO_Q, i, j);
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
    
    function drawLine(cell, lastCell, move, i, j) {
        var cellHeight = cell.offsetHeight;
        var cellWidth = cell.offsetWidth;

        var left;
        var top;
        var lastLeft;
        var lastTop;

        if (move === MOVE.X_TO_P) {
            left        =   (cell.offsetLeft        + cellWidth     * CELL_PERCENT).toString();
            top         =   (cell.offsetTop         + cellHeight    * CELL_PERCENT).toString();
            lastLeft    =   (lastCell.offsetLeft    + cellWidth     * CELL_PERCENT).toString();
            lastTop     =   (lastCell.offsetTop     + cellHeight    * (1-CELL_PERCENT)).toString();
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
        line.setAttribute("marker-end", SVG.MARKER.URL);
        line.setAttribute("stroke", SVG.OBJECT_COLOR);
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
                markCells(path, calculationVerticalTable, calculationTable, calculationHorizontalTable, -1, true);
                visualizerInstance.lastPath = path;
            }
        } else {  // case: first time selected
            markCells(path, calculationVerticalTable, calculationTable, calculationHorizontalTable, -1, true);
            visualizerInstance.lastPath = path;
        }
    }

    function highlight(rowNumber, table) {
        var start = 1;  // rows[0] contains the header and we want skip it

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
        var table = getTable(e);
        var tableCSV = htmlToCsv(table);
        var tableFile = new File([tableCSV], {type: TABLE.TEXT_FILE_ENCODING});

        saveAs(tableFile, TABLE.DOWNLOAD_NAME);
    }

    function getTable(e) {
        var number = e.data.number;

        var calculationVerticalTable;
        var calculationTable = e.data.calculationTable[0];
        var calculationHorizontalTable;

        if (e.data.calculationVerticalTable !== undefined) {
            calculationVerticalTable = e.data.calculationVerticalTable[0];
            calculationHorizontalTable = e.data.calculationHorizontalTable[0];
        }

        if (calculationVerticalTable !== undefined) {  // OR: calculationHor.. !== undefined (case: more than table)
            if (number === 0)
                return calculationVerticalTable;
            else if (number === 1)
                return calculationTable;

            return calculationHorizontalTable;
        }

        return calculationTable;
    }

    /*
    CSV specification: https://www.ietf.org/rfc/rfc4180.txt
    Hint: It is allowed to have a line break in the last line.
     */
    function htmlToCsv(table) {
        var a = SYMBOLS.EMPTY;
        for (var i = 0; i < table.rows.length - 2; i++) {
            var row = [];

            for (var j = 0; j < table.rows[i].cells.length; j++) {
                row[j] = table.rows[i].cells[j].innerHTML;
            }

            if (i < table.rows.length - 3)
                a = a + removeSubscripts(getString(row)) + NEW_LINE;
            else
                a = a + removeSubscripts(getString(row));
        }

        return a;
    }

    function getString(row) {
        var stringRow = row.toString();
        stringRow = stringRow.replace(SUB.START_TAGS, SYMBOLS.EMPTY);
        stringRow = stringRow.replace(SUB.END_TAGS, SYMBOLS.EMPTY);
        stringRow = stringRow.replace(MATH_JAX_TAGS, SYMBOLS.EMPTY);
        stringRow = stringRow.replace(DOUBLE_INFINITIES, SYMBOLS.INFINITY);
        return stringRow;
    }

    function removeSubscripts(string) {
        var subscriptPositions = [];
        var copyString = string.slice(0);

        var position = -1;
        while ((position = copyString.search(CHARACTER.BASE)) !== -1) {
            subscriptPositions.push(position);
            copyString = copyString.slice(0, position) + SYMBOLS.DUMMY + copyString.slice(position + 1);
        }

        for (var i = 0; i < subscriptPositions.length; i++)
            string = string.slice(0, subscriptPositions[i] + 1) + SYMBOLS.SPACE + string.slice(subscriptPositions[i] + 2);

        return string.replace(MULTI_SYMBOLS.SPACE, SYMBOLS.EMPTY);
    }

    function replaceInfinityStrings(matrix) {
        for (var i = 0; i < matrix.length; i++) {
            for (var j = 0; j < matrix[0].length; j++) {
                if (matrix[i][j] === Number.POSITIVE_INFINITY)
                    matrix[i][j] = SYMBOLS.POSITIVE_INFINITY;
                else if (matrix[i][j] === Number.NEGATIVE_INFINITY)
                    matrix[i][j] = SYMBOLS.NEGATIVE_INFINITY;
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
        drawArrowLine(visualizerInstance.lastPath, calculationVerticalTable, table, calculationHorizontalTable);

        var lastFlows = visualizerInstance.lastFlows;
        for (var i = 0; i < lastFlows.length; i++)
            drawArrowLine(lastFlows[i], calculationVerticalTable, table, calculationHorizontalTable);
    }

    function drawArrowLine(path, calculationVerticalTable, table, calculationHorizontalTable) {
        var lastPosI;
        var lastPosJ;
        var lastTable;

        var currentTable;

        for (var j = 0; j < path.length; j++) {
            currentTable = getRightTable(path, j, calculationVerticalTable, table, calculationHorizontalTable);

            var posI = path[j].i + 1;
            var posJ = path[j].j + 1;

            placeArrow(currentTable, posI, posJ, lastTable, lastPosI, lastPosJ);
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