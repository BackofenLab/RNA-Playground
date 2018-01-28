/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("postProcessing.visualizer", Visualizer);

    // instances
    var visualizerInstance;

    /**
     * Contains functions for visualization and creation of downloadable files.
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
        this.lastIterationNumber = -1;
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

                if (character !== undefined) {
                    element.innerHTML = character;

                    if (index !== undefined)
                        element.innerHTML += SUB.START_TAG + index + SUB.END_TAG;
                }
            }
        };

        // public class methods
        this.shareInformation = shareInformation;
        this.showFlow = showFlow;
        this.showTraceback = showTraceback;
        this.highlight = highlight;
        this.downloadTable = downloadTable;
        this.replaceInfinityStrings = replaceInfinityStrings;
        this.redrawOverlay = redrawOverlay;
        this.drawTree = drawTree;
        this.markMinima = markMinima;
        this.removeAllContents = removeAllContents;
    }

    /**
     * Creates fundamental SVG code, to draw SVGs on top of the webpage.
     * @return {Element} - SVG code which is used by drawLine()-function, to draw an arrow.
     */
    function createSVG() {
        var svg = document.createElementNS(SVG.NAME_SPACE, "svg");
        svg.appendChild(createEndMarker(SVG.TRACEBACK_LONG_ARROW_COLOR, SVG.MARKER.ID_TRACEBACK));
        svg.appendChild(createEndMarker(SVG.FLOW_LONG_ARROW_COLOR, SVG.MARKER.ID_FLOW));
        return svg;
    }

    /**
     * Creates the end piece of a line.
     * @param color {string} - The color of the marker.
     * @param id {string} - The id of the marker with which the marker can be accessed.
     * @return {Element} - The marker-Element (XML-code Object).
     * @see Created with the help of <marker> definition https://developer.mozilla.org/de/docs/Web/SVG/Element/marker
     */
    function createEndMarker(color, id) {
        // create triangle
        var trianglePath = document.createElementNS(SVG.NAME_SPACE, "path");
        trianglePath.setAttribute("d", SVG.TRIANGLE.D);  // triangle path
        trianglePath.setAttribute("fill", color);

        // create marker object using path defined above
        var markerEnd = document.createElementNS(SVG.NAME_SPACE, "marker");
        markerEnd.setAttribute("id", id);
        markerEnd.setAttribute("orient", SVG.MARKER.ORIENT);
        markerEnd.setAttribute("refX", SVG.MARKER.BOUNDS.REF_X);  // relative marker coordinate
        markerEnd.setAttribute("refY", SVG.MARKER.BOUNDS.REF_Y);  // relative marker coordinate
        markerEnd.setAttribute("markerWidth", SVG.MARKER.BOUNDS.WIDTH);
        markerEnd.setAttribute("markerHeight", SVG.MARKER.BOUNDS.HEIGHT);
        markerEnd.setAttribute("viewBox", SVG.MARKER.VIEW_BOX);
        markerEnd.appendChild(trianglePath);

        return markerEnd;
    }

    /**
     * Sharing the algorithm and its output and input with the visualizer such that the algorithm can work with the data.
     * @param algorithm {Object} - Contains an alignment algorithm.
     * @param input {Object} - Contains all input data.
     * @param output {Object} - Contains all output data.
     */
    function shareInformation(algorithm, input, output) {
        visualizerInstance.algorithm = algorithm;
        visualizerInstance.input = input;
        visualizerInstance.output = output;
    }

    /**
     * Shows to which cells a cell can backtrack.
     * @param cellCoordinates {Object} - The vector-position of the cell which have been clicked.
     * @param calculationVerticalTable {Element} - The table storing the vertical gap costs.
     * @param table {Element} - The default or main table.
     * @param calculationHorizontalTable {Element} - The table storing the horizontal gap costs.
     * @param iterationTablesArray {Array} - An array of tables.
     * @param mainOutput {Element} - The div containing only the calculation tables.
     * @param iteration {number} - The iteration from which table was selected.
     */
    function showFlow(cellCoordinates, calculationVerticalTable, table, calculationHorizontalTable,
                      iterationTablesArray, mainOutput, iteration) {

        var flows = getTraces(cellCoordinates, iteration);

        if (iterationTablesArray !== undefined) {  // if alignment algorithm with iterations
            removeAllFlows(table);

            for (var i = 0; i < flows.length; i++)
                markCells(flows[i].reverse(), calculationVerticalTable, table, calculationHorizontalTable, mainOutput, i, true, true);

        } else {  // in the case of a non-iterative alignment algorithm
            for (i = 0; i < visualizerInstance.lastFlows.length; i++)
                demarkCells(visualizerInstance.lastFlows[i], calculationVerticalTable, table, calculationHorizontalTable, mainOutput, i, true);

            for (var i = 0; i < flows.length; i++)
                markCells(flows[i].reverse(), calculationVerticalTable, table, calculationHorizontalTable, mainOutput, i, true, true);
        }

        visualizerInstance.lastFlows = flows;
    }

    /**
     * Returns the right traces for the right algorithm.
     * @param cellCoordinates {Object} - The vector-position of the cell which have been clicked.
     * @param iteration {number} - The iteration from which table was selected.
     * @return {Array} - The traces which have to be highlighted.
     */
    function getTraces(cellCoordinates, iteration) {
        var algorithm = visualizerInstance.algorithm;
        var superclass = algorithm.getSuperclass();
        var child = superclass.getLastChild();  // hint: has not to be algorithm

        var flows = [];
        algorithm.numberOfTracebacks = 0;  // to avoid counting and a cancellation after some reached limit
        child.numberOfTracebacks = 0;

        var input = getTraceComputationInput(algorithm, iteration);

        // calling the right neighborhood function
        if (algorithm.getNeighboured !== undefined) {
            if (GLOBAL_ALGORITHMS.indexOf(algorithm.type) >= 0)
                flows = superclass.getGlobalTraces([cellCoordinates], input, visualizerInstance.output, 1, algorithm.getNeighboured);
            else if (LOCAL_ALGORITHMS.indexOf(algorithm.type) >= 0)
                flows = superclass.getLocalTraces([cellCoordinates], input, visualizerInstance.output, 1, algorithm.getNeighboured);
        }
        else {
            if (GLOBAL_ALGORITHMS.indexOf(algorithm.type) >= 0)
                flows = superclass.getGlobalTraces([cellCoordinates], input, visualizerInstance.output, 1, superclass.getNeighboured);
            else if (LOCAL_ALGORITHMS.indexOf(algorithm.type) >= 0)
                flows = superclass.getLocalTraces([cellCoordinates], input, visualizerInstance.output, 1, superclass.getNeighboured);
        }

        return flows;
    }

    /**
     * Returns the input on which flow is computed.
     * @param algorithm {Object} - The algorithm for which flow is computed.
     * @param iteration {number} - The iteration from which table was selected.
     * @return {Object} - Input data for the algorithm which is used for backtracking.
     */
    function getTraceComputationInput(algorithm, iteration) {
        var input = jQuery.extend(true, {}, visualizerInstance.input);  // deep copy to avoid changes on original values

        if (algorithm.type === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER) {
            visualizerInstance.output.matrix = visualizerInstance.output.iterationData[0][iteration][8];
            if (iteration > 0) {  // for computed matrix Smith-Waterman is used with some lambda
                input.deletion = visualizerInstance.output.iterationData[0][iteration - 1][3];  // deletion value with lambda
                input.insertion = visualizerInstance.output.iterationData[0][iteration - 1][4];  // ....
                input.match = visualizerInstance.output.iterationData[0][iteration - 1][5];
                input.mismatch = visualizerInstance.output.iterationData[0][iteration - 1][6];
            }
            visualizerInstance.lastIterationNumber = iteration;
        } else {
            visualizerInstance.lastIterationNumber = -1;
        }

        return input;
    }

    /**
     * Removes all flows.
     * @param table {Element} - The table from which flows have to be removed.
     */
    function removeAllFlows(table) {
        for (var i = 1; i < table.rows.length; i++) {
            for (var j = 1; j < table.rows[i].cells.length; j++) {
                var cell = table.rows[i].cells[j];

                if (cell.classList.contains("selected")
                    && cell.innerText !== SMITH_WATERMAN_STOP)  // avoids removing tracebacks
                    removeArrows(table, i, j, true);
                else
                    removeArrows(table, i, j, false);

                removeFlowColors(table, i, j);
            }
        }
    }

    /**
     * Removes short sprite sheet arrows at a specific position in the table.
     * @param table {Element} - The table from which arrows have to be removed.
     * @param posI {number} - The first coordinate.
     * @param posJ {number} - The second coordinate.
     * @param keepLast {boolean} - Allows to keep last arrow or not.
     */
    function removeArrows(table, posI, posJ, keepLast) {
        var cell = table.rows[posI].cells[posJ];
        var numChildren = keepLast ? cell.children.length - 1 : cell.children.length;

        for (var k = 0; k < numChildren; k++)
            cell.children[0].outerHTML = SYMBOLS.EMPTY;
    }

    /**
     * Removes flow colors at a specific position in the table.
     * @param table {Element} - The table from which colors have to be removed.
     * @param posI {number} - The first coordinate.
     * @param posJ {number} - The second coordinate.
     */
    function removeFlowColors(table, posI, posJ) {
        table.rows[posI].cells[posJ].classList.remove("selected_light_red");
        table.rows[posI].cells[posJ].classList.remove("selected_very_light_red");
        table.rows[posI].cells[posJ].classList.remove("selected_red");
        table.rows[posI].cells[posJ].classList.remove("selected_green");
    }

    /**
     * Turns on cell highlights.
     * @param path {Array} - Array containing the first vector element from which on you want find a path.
     * @param calculationVerticalTable {Element} - The table storing the vertical gap costs.
     * @param table {Element} - The default or main table.
     * @param calculationHorizontalTable {Element} - The table storing the horizontal gap costs.
     * @param mainOutput {Element} - The div containing only the calculation tables.
     * @param colorClass {number} - The highlight which should be added to a cell.
     * @param arrows {boolean} - Tells if arrows should be drawn or not.
     * @param flowMode {boolean} - Tells if flows or traceback-paths are drawn.
     */
    function markCells(path, calculationVerticalTable, table, calculationHorizontalTable, mainOutput, colorClass, arrows, flowMode) {
        arrows = arrows || false;

        var lastPosI;
        var lastPosJ;
        var lastTable;

        var currentTable;

        if (path !== undefined) {
            // go over the whole path
            for (var j = 0; j < path.length; j++) {
                currentTable = getRightTable(path, j, calculationVerticalTable, table, calculationHorizontalTable);

                var posI = path[j].i + 1;
                var posJ = path[j].j + 1;

                switch (colorClass) {  // selecting by adding the right color class to the element
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

                if (j === path.length - 1 && colorClass !== -1) {  // start element should be green in a flow visualization
                    removeFlowColors(currentTable, posI, posJ);
                    currentTable.rows[posI].cells[posJ].classList.add("selected_green");
                }

                if (arrows) {  // draw arrows: YES or NO
                    placeArrow(currentTable, posI, posJ, mainOutput, lastTable, lastPosI, lastPosJ, flowMode);
                    lastPosI = posI;
                    lastPosJ = posJ;
                    lastTable = currentTable;
                }
            }
        }
    }

    /**
     * Allows to draw short and long arrows on top of tables.
     * @param table {Element} - The current table which is visited.
     * @param posI {number} - The first coordinate.
     * @param posJ {number} - The second coordinate.
     * @param mainOutput {Element} - The div containing only the calculation tables.
     * @param lastTable {Element} - The table that was visited before.
     * @param lastPosI {number} - The last first coordinate.
     * @param lastPosJ {number} - The last second coordinate.
     * @param flowMode {boolean} - Tells if flows or traceback-paths are drawn.
     */
    function placeArrow(table, posI, posJ, mainOutput, lastTable, lastPosI, lastPosJ, flowMode) {
        // function is executed only if one step in the table has already been done
        if (lastPosI !== undefined && lastPosI !== undefined) {
            var cell = table.rows[posI].cells[posJ];
            var lastCell = lastTable.rows[lastPosI].cells[lastPosJ];

            // drawings within the same table
            if (lastTable === table) {
                // case determination
                var isTop = posI - lastPosI === 1;
                var isLeft = posJ - lastPosJ === 1;
                var isVertical = posI - lastPosI > 1;
                var isHorizontal = posJ - lastPosJ > 1;

                // draw case
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
                } else if (isVertical) {
                    drawLine(cell, lastCell, MOVE.VERTICAL, mainOutput, flowMode);
                } else if (isHorizontal) {
                    drawLine(cell, lastCell, MOVE.HORIZONTAL, mainOutput, flowMode);
                }
            } else if (lastTable !== table) {  // drawings between two different tables
                var parentMatrix = getParentMatrix(cell);
                var lastParentMatrix = getParentMatrix(lastCell);

                // case determination
                var isPtoX = parentMatrix === MATRICES.VERTICAL && lastParentMatrix === MATRICES.DEFAULT;
                var isQtoX = parentMatrix === MATRICES.HORIZONTAL && lastParentMatrix === MATRICES.DEFAULT;
                var isXtoP = parentMatrix === MATRICES.DEFAULT && lastParentMatrix === MATRICES.VERTICAL;
                var isXtoQ = parentMatrix === MATRICES.DEFAULT && lastParentMatrix === MATRICES.HORIZONTAL;

                // draw case
                if (isPtoX)
                    drawLine(cell, lastCell, MOVE.P_TO_X, mainOutput, flowMode);
                else if (isQtoX)
                    drawLine(cell, lastCell, MOVE.Q_TO_X, mainOutput, flowMode);
                else if (isXtoP)
                    drawLine(cell, lastCell, MOVE.X_TO_P, mainOutput, flowMode);
                else if (isXtoQ)
                    drawLine(cell, lastCell, MOVE.X_TO_Q, mainOutput, flowMode);
            }
        }
    }

    /**
     * Returns the matrix identifier of a cell.
     * @param cell {Element} - The cell from which the matrix is determined.
     * @return {string} - The identifier of the matrix cell.
     */
    function getParentMatrix(cell) {
        if (cell.parentNode.parentNode.parentNode.className === "calculation_horizontal")
            return MATRICES.HORIZONTAL;
        else if (cell.parentNode.parentNode.parentNode.className === "calculation_vertical")
            return MATRICES.VERTICAL;
        else if (cell.parentNode.parentNode.parentNode.className === "calculation")
            return MATRICES.DEFAULT;
    }

    /**
     * Drawing a line from cell to the last cell. Hint: The visiting of cells is starting in the top left.
     * @param cell {Element} - The current cell which is visited.
     * @param lastCell {Element} - The cell that was visited before.
     * @param move {string} - The type of MOVE {P_TO_X, Q_TO_X, ...}.
     * @param mainOutput {Element} - The div containing only the calculation tables.
     * @param flowMode {boolean} - Tells if flows or traceback-paths are drawn.
     */
    function drawLine(cell, lastCell, move, mainOutput, flowMode) {
        var cellHeight = cell.offsetHeight;
        var cellWidth = cell.offsetWidth;

        var left;
        var top;
        var lastLeft;
        var lastTop;

        // with different moves, the long arrow has to be drawn different
        if (move === MOVE.HORIZONTAL) {
            left = (cell.offsetLeft + cellWidth * CELL_PERCENT.LINE).toString();
            top = (cell.offsetTop + cellHeight * CELL_PERCENT.LINE).toString();
            lastLeft = (lastCell.offsetLeft + cellWidth * (1 - CELL_PERCENT.LINE_HEAD_PENETRATION)).toString();
            lastTop = (lastCell.offsetTop + cellHeight * CELL_PERCENT.LINE).toString();
        } else if (move === MOVE.P_TO_X) {
            left = (cell.offsetLeft + cellWidth * (1 - CELL_PERCENT.LINE)).toString();
            top = (cell.offsetTop + cellHeight * (1 - CELL_PERCENT.LINE)).toString();
            lastLeft = (lastCell.offsetLeft + cellWidth * (1 - CELL_PERCENT.LINE)).toString();
            lastTop = (lastCell.offsetTop + cellHeight * CELL_PERCENT.LINE).toString();
        } else if (move === MOVE.Q_TO_X) {
            left = (cell.offsetLeft + cellWidth * CELL_PERCENT.LINE).toString();
            top = (cell.offsetTop + cellHeight * CELL_PERCENT.LINE).toString();
            lastLeft = (lastCell.offsetLeft + cellWidth * (1 - CELL_PERCENT.LINE)).toString();
            lastTop = (lastCell.offsetTop + cellHeight * (1 - CELL_PERCENT.LINE)).toString();
        } else if (move === MOVE.VERTICAL) {
            left = (cell.offsetLeft + cellWidth * CELL_PERCENT.LINE).toString();
            top = (cell.offsetTop + cellHeight * CELL_PERCENT.LINE).toString();
            lastLeft = (lastCell.offsetLeft + cellWidth * CELL_PERCENT.LINE).toString();
            lastTop = (lastCell.offsetTop + cellHeight * (1 - CELL_PERCENT.LINE_HEAD_PENETRATION)).toString();
        } else if (move === MOVE.X_TO_P) {
            left = (cell.offsetLeft + cellWidth * CELL_PERCENT.LINE_2).toString();
            top = (cell.offsetTop + cellHeight * CELL_PERCENT.LINE).toString();
            lastLeft = (lastCell.offsetLeft + cellWidth * CELL_PERCENT.LINE_2).toString();
            lastTop = (lastCell.offsetTop + cellHeight * (1 - CELL_PERCENT.LINE)).toString();
        } else if (move === MOVE.X_TO_Q) {
            left = (cell.offsetLeft + cellWidth * CELL_PERCENT.LINE).toString();
            top = (cell.offsetTop + cellHeight * (1 - CELL_PERCENT.LINE)).toString();
            lastLeft = (lastCell.offsetLeft + cellWidth * CELL_PERCENT.LINE).toString();
            lastTop = (lastCell.offsetTop + cellHeight * CELL_PERCENT.LINE).toString();
        }

        // define svg dimensions
        visualizerInstance.svg.setAttribute("width", mainOutput.offsetWidth);
        visualizerInstance.svg.setAttribute("height", mainOutput.offsetHeight);
        visualizerInstance.container.style.left = mainOutput.offsetLeft.toString() + SVG.PX;
        visualizerInstance.container.style.top = mainOutput.offsetTop.toString() + SVG.PX;

        // create line with previously defined marker
        var line = document.createElementNS(SVG.NAME_SPACE, "line");

        if (flowMode) {  // differentiate between drawing flows and tracebacks
            line.setAttribute("marker-end", SVG.MARKER.URL_FLOW);
            line.setAttribute("stroke", SVG.FLOW_LONG_ARROW_COLOR);
        }
        else {
            line.setAttribute("marker-end", SVG.MARKER.URL_TRACEBACK);
            line.setAttribute("stroke", SVG.TRACEBACK_LONG_ARROW_COLOR);
        }
        line.setAttribute("x1", left - mainOutput.offsetLeft - mainOutput.scrollLeft);
        line.setAttribute("y1", top - mainOutput.offsetTop - mainOutput.scrollTop);
        line.setAttribute("x2", lastLeft - mainOutput.offsetLeft - mainOutput.scrollLeft);
        line.setAttribute("y2", lastTop - mainOutput.offsetTop - mainOutput.scrollTop);
        line.setAttribute("stroke-dasharray", SVG.STROKE_DASHARRAY);
        visualizerInstance.svg.appendChild(line);
        visualizerInstance.cellLines.push(line);

        // creating an SVG overlay if there is not already one
        if (visualizerInstance.container.childElementCount !== 1)
            visualizerInstance.container.appendChild(visualizerInstance.svg);
    }

    /**
     * Turns off cell highlights.
     * @param path {Array} - Array containing the first vector element from which on you want find a path.
     * @param calculationVerticalTable {Element} - The table storing the vertical gap costs.
     * @param table {Element} - The default or main table.
     * @param calculationHorizontalTable {Element} - The table storing the horizontal gap costs.
     * @param mainOutput {Element} - The div containing only the calculation tables.
     * @param colorClass {number} - The highlight which should be deleted from the cell.
     * @param flowMode {boolean} - Tells if flows or traceback-paths were drawn.
     */
    function demarkCells(path, calculationVerticalTable, table, calculationHorizontalTable, mainOutput, colorClass, flowMode) {
        flowMode = flowMode || false;

        var currentTable;

        if (path.length > 0) {
            // go over the whole path
            for (var j = 0; j < path.length; j++) {
                currentTable = getRightTable(path, j, calculationVerticalTable, table, calculationHorizontalTable);

                var posI = path[j].i + 1;
                var posJ = path[j].j + 1;

                if (currentTable.rows[posI].cells[posJ] !== undefined) {  // if table has shrinked

                    switch (colorClass) {  // deselecting by removing the "selected" class from the element
                        case -1:
                            currentTable.rows[posI].cells[posJ].classList.remove("selected");
                            break;
                        default:
                            removeFlowColors(currentTable, posI, posJ);
                    }
                }

                removeArrows(currentTable, posI, posJ, false);
            }
        }

        removeAllLines();  // below last flows/paths redrawn

        if (flowMode)  // redraw last traceback
            markCells(visualizerInstance.lastPath, calculationVerticalTable, table, calculationHorizontalTable, mainOutput, -1, true, false);
        else {  // redraw last flow
            var lastFlows = visualizerInstance.lastFlows;
            for (var i = 0; i < lastFlows.length; i++)
                markCells(lastFlows[i], calculationVerticalTable, table, calculationHorizontalTable, mainOutput, i, true, true);
        }
    }

    /**
     * Returns the table of the path element on which is currently looked at.
     * @param path {Array} - Array containing the first vector element from which on you want find a path.
     * @param j {number} - Current path element from which the table must be determined.
     * @param calculationVerticalTable {Element} - The table storing the vertical gap costs.
     * @param table {Element} - The default or main table.
     * @param calculationHorizontalTable {Element} - The table storing the horizontal gap costs.
     * @return {Element} - Default table or table for vertical or horizontal gap costs are returned.
     */
    function getRightTable(path, j, calculationVerticalTable, table, calculationHorizontalTable) {
        if (path[j].label === MATRICES.VERTICAL)
            return calculationVerticalTable;
        else if (path[j].label === MATRICES.HORIZONTAL)
            return calculationHorizontalTable;
        else // if (path[j].label === MATRICES.DEFAULT)
            return table;
    }

    /**
     * Removes long SVG arrows in the table.
     */
    function removeAllLines() {
        var line;
        while ((line = visualizerInstance.cellLines.pop()) !== undefined) {
            if (visualizerInstance.svg.contains(line))
                visualizerInstance.svg.removeChild(line);
        }
    }

    /**
     * Highlights tracebacks in the matrix.
     * @param traceNumber {number} - The current path which should be drawn.
     * @param calculationVerticalTable {Element} - The table storing the vertical gap costs.
     * @param calculationTable {Element} - The default or main table.
     * @param calculationHorizontalTable {Element} - The table storing the horizontal gap costs.
     * @param iterationTablesArray {Array} - An array of tables..
     * @param mainOutput {Element} - The div containing only the calculation tables.
     */
    function showTraceback(traceNumber, calculationVerticalTable, calculationTable, calculationHorizontalTable, iterationTablesArray, mainOutput) {
        if (iterationTablesArray !== undefined) {  // if an iterative alignment algorithm
            // iterate over rounds
            for (var i = 0; i < visualizerInstance.output.iterationData[0].length; i++) {
                var firstPathOfRound = visualizerInstance.output.iterationData[0][i][9][0];
                markCells(firstPathOfRound, calculationVerticalTable, iterationTablesArray[i][0], calculationHorizontalTable,
                    mainOutput, -1, true, false);
            }
        } else {  // else if a non-iterative algorithm
            var path = visualizerInstance.output.tracebackPaths[traceNumber];

            // check if you want maybe disable (unhighlight) last drawn path
            if (visualizerInstance.lastPath.length > 0) {
                var posI = visualizerInstance.lastPath[0].i + 1;
                var posJ = visualizerInstance.lastPath[0].j + 1;
                var tableCell = calculationTable.rows[posI].cells[posJ];

                // check if you want disable "unhighlight" last drawn path (example: clicked second time on same path in results table)
                if (path === visualizerInstance.lastPath
                    && tableCell !== undefined
                    && tableCell.classList.contains("selected")) {  // case: same path
                    demarkCells(visualizerInstance.lastPath, calculationVerticalTable, calculationTable, calculationHorizontalTable, mainOutput, -1, false);
                    visualizerInstance.lastPath = [];
                } else {  // case: different path (click on a new path)
                    demarkCells(visualizerInstance.lastPath, calculationVerticalTable, calculationTable, calculationHorizontalTable, mainOutput, -1, false);
                    markCells(path, calculationVerticalTable, calculationTable, calculationHorizontalTable, mainOutput, -1, true, false);
                    visualizerInstance.lastPath = path;
                }
            } else {  // case: first time selected
                markCells(path, calculationVerticalTable, calculationTable, calculationHorizontalTable, mainOutput, -1, true, false);
                visualizerInstance.lastPath = path;
            }
        }
    }

    /**
     * Highlights a selected entry for example in the results table.
     * @param rowNumber {number} - The row number of the row which was clicked.
     * @param table {Element} - The default or main table.
     */
    function highlight(rowNumber, table) {
        var start = 0;  // the row number in which highlighting starts

        if (table.rows.length > 0) {
            // determining cells
            var cell = table.rows[rowNumber + start].cells[0];
            var lastCell = visualizerInstance.lastRowNumber >= 0
                ? table.rows[visualizerInstance.lastRowNumber + start].cells[0] : undefined;

            // check if the row has to be selected or deselected
            if (rowNumber === visualizerInstance.lastRowNumber
                && lastCell.classList.contains("selected")) {  // case: same cell (clicked second time) -> deselect
                lastCell.classList.remove("selected");
            } else if (lastCell !== undefined) {  // case: different cell -> deselect last cell and select current cell
                lastCell.classList.remove("selected");
                cell.classList.add("selected");
            } else {  // case: first time clicked in table -> select current cell
                cell.classList.add("selected");
            }

            visualizerInstance.lastRowNumber = rowNumber;
        }
    }

    /**
     * Allows downloading a table.
     * Hint: " e.data.number" allows to distinguish between the different tables of an algorithm.
     * @param e - Stores data relevant to the event called that function.
     */
    function downloadTable(e) {
        var number = e.data.number;

        var matrix = getMatrix(number);
        if (matrix !== undefined) {
            var upperString = visualizerInstance.input.sequenceB;
            var leftString = visualizerInstance.input.sequenceA;

            var tableCSV = tableToCSV(number, matrix, upperString, leftString);
            var tableFile = new File([tableCSV], {type: TABLE.TEXT_FILE_ENCODING});

            saveAs(tableFile, TABLE.DOWNLOAD_NAME);
        }
    }

    /**
     * Allows to select a preprocessed matrix by number identifier.
     * @param number - Table number, allows to select a table {DEFAULT, HORIZONTAL, VERTICAL, ITERATION_i}.
     * @return {matrix} - The appropriate matrix to the number which was passed.
     */
    function getMatrix(number) {
        switch (number) {
            case MATRICES.VERTICAL_NUMBER:
                return replaceInfinities(visualizerInstance.output.verticalGaps);
            case MATRICES.HORIZONTAL_NUMBER:
                return replaceInfinities(visualizerInstance.output.horizontalGaps);
            default:  // downloading matrices from AEP algorithm
                if (MATRICES.ITERATION_NUMBER_5 <= number && number <= MATRICES.ITERATION_NUMBER_1)
                    if (visualizerInstance.output.iterationData !== undefined
                        && visualizerInstance.output.iterationData[0] !== undefined
                        && visualizerInstance.output.iterationData[0][-(number + 1)] !== undefined
                        && visualizerInstance.output.iterationData[0][-(number + 1)][8] !== undefined)
                        return visualizerInstance.output.iterationData[0][-(number + 1)][8];  // iteration numbers are negative in "defaults.js"
        }

        return visualizerInstance.output.matrix;
    }

    /**
     * Replaces LaTeX-infinity-symbols with infinity-values by going through the whole table.
     * @param matrix {matrix} - The matrix in which you want replace LaTeX-infinity-symbols with infinity-values.
     * @return {matrix} - The appropriate matrix to the number which was passed.
     */
    function replaceInfinities(matrix) {
        if (matrix !== undefined) {
            for (var i = 0; i < matrix.length; i++) {
                for (var j = 0; j < matrix[0].length; j++) {
                    if (matrix[i][j] === Number.NEGATIVE_INFINITY)
                        matrix[i][j] = SYMBOLS.NEGATIVE_INFINITY;
                    else if (matrix[i][j] === Number.POSITIVE_INFINITY)
                        matrix[i][j] = SYMBOLS.INFINITY;
                }
            }

            return matrix;
        }
    }

    /**
     * Computes the string which is stored in a table download-file.
     * @param number - The type of matrix {0: MATRICES.VERTICAL, 1: MATRICES.DEFAULT, 2: MATRICES.HORIZONTAL}.
     * @param matrix - The matrix from which you get the values.
     * @param upperString - The string of the first input sequence.
     * @param leftString {string} - The string of the second input sequence.
     * @return {string} - CSV formatted data.
     * @see CSV specification: https://www.ietf.org/rfc/rfc4180.txt
     */
    function tableToCSV(number, matrix, upperString, leftString) {
        var string = SYMBOLS.EMPTY;

        // differentiate between the type of matrices you can download and add a symbol for the matrix type to the string
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
            default:
                string += MATRICES.DEFAULT + (-number) + SYMBOLS.COMMA;  // iteration numbers are negative
                break;
        }
        string += SYMBOLS.COMMA + upperString.split(SYMBOLS.EMPTY).toString() + SYMBOLS.NEW_LINE;

        // compute CSV
        var csvParser = new formats.csvParser.CsvParser();
        string += csvParser.getCSVData(matrix, leftString);
        return string;
    }

    /**
     * Post edits a matrix and replaces infinity-values with LaTeX-infinity-symbols.
     * @param matrix - The matrix in which you want replace infinity-values with LaTeX-symbols.
     * @return {Array} - The matrix in which symbols where replaced with LaTeX-symbols.
     */
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

    /**
     * Redraw overlay after a resize-, zoom-in- or scrolling-event it the browser window.
     * @param calculationVerticalTable {Element} - The table storing the vertical gap costs.
     * @param calculationTable {Element} - The default or main table.
     * @param calculationHorizontalTable {Element} - The table storing the horizontal gap costs.
     * @param mainOutput {Element} - The div containing only the calculation tables.
     */
    function redrawOverlay(calculationVerticalTable, calculationTable, calculationHorizontalTable, mainOutput) {
        removeAllLines();
        drawAllLines(calculationVerticalTable, calculationTable, calculationHorizontalTable, mainOutput);
    }

    /**
     * Draw all lines which were previously drawn (before some resize-, zoom-in- or scrolling-event).
     * @param calculationVerticalTable {Element} - The table storing the vertical gap costs.
     * @param table {Element} - The default or main table.
     * @param calculationHorizontalTable {Element} - The table storing the horizontal gap costs.
     * @param mainOutput {Element} - The div containing only the calculation tables.
     */
    function drawAllLines(calculationVerticalTable, table, calculationHorizontalTable, mainOutput) {
        drawArrowLines(visualizerInstance.lastPath, calculationVerticalTable, table, calculationHorizontalTable, mainOutput, false);

        var lastFlows = visualizerInstance.lastFlows;
        for (var i = 0; i < lastFlows.length; i++)
            drawArrowLines(lastFlows[i], calculationVerticalTable, table, calculationHorizontalTable, mainOutput, true);
    }

    /**
     * Redrawing arrow lines.
     * @param path {Array} - Array containing the first vector element from which on you want find a path.
     * @param calculationVerticalTable {Element} - The table storing the vertical gap costs.
     * @param table {Element} - The default or main table.
     * @param calculationHorizontalTable {Element} - The table storing the horizontal gap costs.
     * @param mainOutput {Element} - The div containing only the calculation tables.
     * @param flowMode {boolean} - Tells if flows or traceback-paths are drawn.
     */
    function drawArrowLines(path, calculationVerticalTable, table, calculationHorizontalTable, mainOutput, flowMode) {
        var lastPosI;
        var lastPosJ;
        var lastTable;

        var currentTable;

        if (path !== undefined) {
            // going over the whole path and set right-positioned arrows by a recalculation
            for (var j = 0; j < path.length; j++) {
                currentTable = getRightTable(path, j, calculationVerticalTable, table, calculationHorizontalTable);

                var posI = path[j].i + 1;
                var posJ = path[j].j + 1;

                placeArrow(currentTable, posI, posJ, mainOutput, lastTable, lastPosI, lastPosJ, flowMode);
                lastPosI = posI;
                lastPosJ = posJ;
                lastTable = currentTable;
            }
        }
    }

    /**
     * Draws a phylogenetic tree from jsPhyloSVG-library.
     * @see: Bugfix-code for https://jsphylosvg.uservoice.com/forums/55902-general/suggestions/7252947-clear-and-reload-tree
     * was taken from
     * https://stackoverflow.com/questions/30667884/why-is-the-bottom-of-the-figure-cut-off
     * https://pastebin.com/9w4PXtLQ
     * and has been optimized, extended and commented.
     */
    function drawTree() {
        $("#phylogenetic_tree").remove();  // remove from container
        $(".tree_container").append(PHYLOGENETIC_TREE.SVG_CANVAS);  // add again

        var newick = visualizerInstance.output.newickString;

        if (visualizerInstance.output.newickString.length !== 1
            && newick.indexOf(SYMBOLS.MINUS) === -1) {  // if there is not only a ";" and if there are no negative values

            var numberOfClusters;

            if (visualizerInstance.algorithm.type === ALGORITHMS.AGGLOMERATIVE_CLUSTERING)
                numberOfClusters = visualizerInstance.input.initialNamingIndex;
            else  // if clustering algorithm
                numberOfClusters = visualizerInstance.input.sequences.length - visualizerInstance.input.arrayPositionsOfRemovedSequences.length;

            var svgHeight = numberOfClusters * PHYLOGENETIC_TREE.SVG_DIMENSION_FACTOR;  // make it dependant on the number of clusters

            visualizerInstance.phylogeneticTree
                = new Smits.PhyloCanvas(newick,
                PHYLOGENETIC_TREE.SVG_CANVAS_NAME,
                PHYLOGENETIC_TREE.SVG_WIDTH,
                svgHeight);

            var svgPaths = $("svg path");
            var svgTexts = $("svg text");

            var currentSvgY = parseInt($("svg")[0].attributes.getNamedItem("height").value);

            // search for maximum y in the defined SVG
            var definedMaxY = 0;

            // go through each path (in SVG figures defined as paths)
            for (var i = 0; i < svgPaths.length; i++) {
                var path = svgPaths[i].attributes.getNamedItem("d").value.split(",");

                // go through each element of the path
                for (var j = 1; j < path.length; j++) {
                    var currentY = parseInt(path[j].split("L")[0]);

                    if (definedMaxY < currentY)
                        definedMaxY = currentY;
                }
            }

            if (definedMaxY > currentSvgY) {  // if (SVG has moved)
                var heightAdjustment = definedMaxY - (currentSvgY - 20);  // compute difference between true value and desired value

                // adjust figures-heights
                // go through each path
                for (var i = 0; i < svgPaths.length; i++) {
                    var path = svgPaths[i].attributes.getNamedItem("d").value.split(SYMBOLS.COMMA);

                    var correctedPath = SYMBOLS.EMPTY + path[0];

                    // go through each element of the path and adjust the y-value
                    for (var j = 1; j < path.length; j++) {
                        var currentY = path[j].split("L");

                        // adjustment
                        if (currentY.length !== 1)
                            correctedPath = correctedPath + SYMBOLS.COMMA + (parseInt(currentY[0]) - heightAdjustment).toString() + "L" + currentY[1];
                        else
                            correctedPath = correctedPath + SYMBOLS.COMMA + (parseInt(currentY[0]) - heightAdjustment).toString();
                    }

                    svgPaths[i].setAttribute("d", correctedPath);
                }

                // adjust text-heights
                for (var i = 0; i < svgTexts.length; i++) {
                    svgTexts[i].setAttribute("y", (parseInt(svgTexts[i].attributes.getNamedItem("y").value) - heightAdjustment).toString());
                }
            }
        }
    }

    /**
     * Highlights minima in the matrix.
     * @param tracecellTable {Element} - The table in which the minimas has to be marked.
     */
    function markMinima(tracecellTable) {
        var tracecellLines = visualizerInstance.output.tracecellLines;
        var tracecellLinesKeys = Object.keys(tracecellLines);

        var height = visualizerInstance.input.matrixHeight;
        var width = visualizerInstance.input.matrixWidth;

        for (var i = 0; i < tracecellLinesKeys.length; i++) {
            var key = tracecellLinesKeys[i];
            var tracecell = tracecellLines[key];

            var posI = tracecell.i;  // hint: there can be empty lines
            var posJ = tracecell.j;

            // if in-between cell
            tracecellTable.rows[posI + 1].cells[posJ + 1].classList.add("selected");  // "+1" because in the first column/row, there is the header
        }

        tracecellTable.rows[1].cells[1].classList.add("selected_red");  // mark top left cell
        tracecellTable.rows[height].cells[width].classList.add("selected_green");  // mark bottom right cell
    }

    /**
     * Removes all contents stored in the visualizer
     * for example after the algorithm is changed
     * or before a recomputation is done.
     */
    function removeAllContents() {
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