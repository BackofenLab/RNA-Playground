/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("postProcessing.inputProcessor", InputProcessor);

    // instances
    var inputProcessorInstance;

    /**
     * Does the forwarding of input
     * to the algorithms and the visualizer.
     * In this connection, it post processes the input
     * by removing/updating wrong inputs
     * to the correct value.
     */
    function InputProcessor() {
        inputProcessorInstance = this;

        // variables
        this.inputUpdatesStarted = false;

        // public class methods
        this.activateInputUpdates = activateInputUpdates;
        this.inputUpdatesActivated = inputUpdatesActivated;
        this.linkElements = linkElements;
        this.postEdit = postEdit;
    }

    /**
     * Activates algorithm input updates after changes on the input.
     */
    function activateInputUpdates() {
        inputProcessorInstance.inputUpdatesStarted = true;
    }

    /**
     * Checks if input updates are activated or not.
     * @return {boolean} - "true" if updates are activated.
     */
    function inputUpdatesActivated() {
        return inputProcessorInstance.inputUpdatesStarted;
    }

    /**
     * Linking static and dynamic elements with a function.
     * @param algorithm {Object} - Algorithm used to update the user interface.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions.
     * @param processInput {Function} - Function from the algorithm which should process the input.
     * @param changeOutput {Function} - Function from the algorithm which should change the output after processing the input.
     */
    function linkElements(algorithm, viewmodels, processInput, changeOutput) {
        fixBrowserBugs();

        // retrieve parameters
        var algorithmInput = $("#algorithm_input");
        var functionParameters = algorithmInput.find(".fx_parameter");

        var mainOutput = $(".main_output");  // output containing the calculation tables
        var calculationTable = mainOutput.find(".calculation");
        var calculationHorizontalTable = mainOutput.find(".calculation_horizontal");
        var calculationVerticalTable = mainOutput.find(".calculation_vertical");

        var selectableEntryClass = ".selectable_entry";

        // linking (alphabetically sorted)
        linkBasicInputsBehaviour(algorithmInput, functionParameters);
        linkDownloadLinks(viewmodels.visual, calculationVerticalTable, calculationTable, calculationHorizontalTable);
        linkDynamicallyCreatedButtons(algorithm, viewmodels, processInput, changeOutput);  // do not exist from beginning
        linkGuiUpdater(algorithm, viewmodels, processInput, changeOutput);
        linkIterationTables(viewmodels.visual);
        linkOverlay(viewmodels.visual, calculationVerticalTable, calculationTable, calculationHorizontalTable, mainOutput);
        linkSelectables(viewmodels.visual, calculationVerticalTable, calculationTable, calculationHorizontalTable,
            mainOutput, selectableEntryClass);

        // highlighting
        doMinimaAndTracebackHighlighting(viewmodels.visual);
        doTracebackHighlighting(viewmodels.visual);

        // initial redrawing
        redrawInitialOverlay(viewmodels.visual, calculationVerticalTable, calculationTable, calculationHorizontalTable, mainOutput);
    }

    /**
     * Bug fixes for specific browsers or versions of browsers.
     */
    function fixBrowserBugs() {
        /*
         BUG-FIX for Mozilla:
         Inputs of type "number" doesn't get the focus
         in the Mozilla browser if one
         of the up- or down-buttons of the number-inputs
         is clicked.

         Detected: Mozilla Firefox 55.0b3 (32-Bit)
         */
        $(function () {
            $("input[type='number']").on("click", function () {
                $(this).focus();
            });
        });
    }

    /**
     * Linking inputs to get some special behaviour (removing non allowed bases).
     * @param algorithmInput {Element} - The input div in which the behaviour is changed.
     * @param functionParameters {Element} - The parameters which should have the given behaviour.
     */
    function linkBasicInputsBehaviour(algorithmInput, functionParameters) {
        var functionArguments = {"functionParameters": functionParameters};
        algorithmInput.find(".optimization_type").on("change", functionArguments, negateOptimizationParameters);

        algorithmInput.on("keyup", ".csv_data", removeNonAllowedCSVSymbols);
        algorithmInput.on("keyup", ".sequence", removeNonAllowedBases);
        algorithmInput.on("keyup", ".sequence_multi", removeNonAllowedBases);
        functionParameters.on("change", removeCriticalNumbers);
    }

    /**
     * Removes non-allowed symbols from a CSV-input.
     * Hint: The input has to be of class "sequence".
     */
    function removeNonAllowedCSVSymbols() {
        if (!CHARACTER.CSV_SYMBOLS.test(this.value)) {
            this.value = this.value.replace(CHARACTER.NON_CSV_SYMBOLS, SYMBOLS.EMPTY);
        }
    }

    /**
     * Removes non-english characters and special characters from an input-field.
     * Hint: The input has to be of class "sequence".
     */
    function removeNonAllowedBases() {
        if (!CHARACTER.BASES.test(this.value))
            this.value = this.value.replace(CHARACTER.NON_BASES, SYMBOLS.EMPTY);
    }

    /**
     * Negates the function parameters of an algorithm.
     * Hint: The input has to be of class "optimization_type".
     * @param e {Object} - Stores data relevant to the event called that function.
     */
    function negateOptimizationParameters(e) {
        var functionParameters = e.data.functionParameters;

        for (var i = 0; i < functionParameters.length; i++)
            functionParameters[i].value = -functionParameters[i].value;
    }

    /**
     * Removes values from an input-field
     * which can lead to problems with traceback or visualization.
     * It is especially used for the function parameters of an algorithm.
     * Hint: The input has to be of class "fx_parameter".
     * @param e {Object} - Stores data relevant to the event called that function.
     */
    function removeCriticalNumbers(e) {
        if (CHARACTER.NUMBERS.test(this.value))
            this.value = Math.round(this.value);

        if (e.type === "change") {
            if (this.id === "global_alignments_per_sequence_pair" || this.id === "local_alignments_per_sequence_pair") {
                this.value = this.value >= INPUT.ALIGNMENTS_MIN ? this.value : INPUT.ALIGNMENTS_MIN;
                this.value = this.value <= INPUT.ALIGNMENTS_MAX ? this.value : INPUT.ALIGNMENTS_MAX;
            } else if (this.id === "length") {
                this.value = this.value >= INPUT.LENGTH_MIN ? this.value : INPUT.LENGTH_MIN;
                this.value = this.value <= INPUT.LENGTH_MAX ? this.value : INPUT.LENGTH_MAX;
            }
            else {
                this.value = this.value >= INPUT.MIN ? this.value : INPUT.MIN;
                this.value = this.value <= INPUT.MAX ? this.value : INPUT.MAX;
            }
        }
    }

    /**
     * Linking table download links of tables with a download function.
     * @param visualViewmodel {Object} - Model which is used for example to highlight cells.
     * @param calculationTable {Element} - The default or main table.
     * @param calculationVerticalTable {Element} - The table storing the vertical gap costs.
     * @param calculationHorizontalTable {Element} - The table storing the horizontal gap costs.
     */
    function linkDownloadLinks(visualViewmodel, calculationVerticalTable, calculationTable, calculationHorizontalTable) {
        var tableDownload = $(".table_download");
        var tableVerticalDownload = $(".table_vertical_download");
        var tableHorizontalDownload = $(".table_horizontal_download");

        var functionArguments = {
            "calculationTable": calculationTable,
            "calculationHorizontalTable": calculationHorizontalTable,
            "calculationVerticalTable": calculationVerticalTable,
            "number": MATRICES.DEFAULT_NUMBER
        };
        tableDownload.on("click", functionArguments, visualViewmodel.downloadTable);

        if (MULTI_TABLE_ALGORITHMS.indexOf(visualViewmodel.algorithm.type) >= 0) {
            functionArguments = {
                "calculationTable": calculationTable,
                "calculationHorizontalTable": calculationHorizontalTable,
                "calculationVerticalTable": calculationVerticalTable,
                "number": MATRICES.VERTICAL_NUMBER
            };
            tableVerticalDownload.on("click", functionArguments, visualViewmodel.downloadTable);

            functionArguments = {
                "calculationTable": calculationTable,
                "calculationHorizontalTable": calculationHorizontalTable,
                "calculationVerticalTable": calculationVerticalTable,
                "number": MATRICES.HORIZONTAL_NUMBER
            };
            tableHorizontalDownload.on("click", functionArguments, visualViewmodel.downloadTable);
        }
    }

    /**
     * Dynamically created objects have to be bind on the whole document or they have to be rebinded.
     * @param algorithm {Object} - Algorithm used to update the user interface.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions.
     * @param processInput {Function} - Function from the algorithm which should process the input.
     * @param changeOutput {Function} - Function from the algorithm which should change the output after processing the input.
     */
    function linkDynamicallyCreatedButtons(algorithm, viewmodels, processInput, changeOutput) {
        if (MULTI_SEQUENCE_ALGORITHMS.indexOf(algorithm.type) >= 0) {  // only multi-sequence algorithms have dynamic buttons
            $(document).off("click", "**"); // remove delegated event handlers like the one below before reassignment

            $(document).on("click", ".add_remove", function () {
                update(algorithm, viewmodels, processInput, changeOutput);
            });
        }
    }

    /**
     * Links the user interface update-function which allows
     * to update the output after some kind of calculation with the input.
     * @param algorithm {Object} - Algorithm used to update the user interface.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions.
     * @param processInput {Function} - Function from the algorithm which should process the input.
     * @param changeOutput {Function} - Function from the algorithm which should change the output after processing the input.
     */
    function linkGuiUpdater(algorithm, viewmodels, processInput, changeOutput) {
        var input = $("#algorithm_input");

        input.on({
            change: function (event) {
                /*
                BUG-FIX for Knockout 3.4.2:
                Knockout fires an event, if the <select>-Tag is filled with Knockout:
                https://stackoverflow.com/questions/16521552/knockout-fires-change-event-when-select-list-initializing
                */
                if (event.cancelable !== undefined) {  // to filter out Knockout-events and let pass all other events
                    update(algorithm, viewmodels, processInput, changeOutput);
                }
            },

            keypress: function (e) {
                if (e.which === KEY_CODES.ENTER)
                    update(algorithm, viewmodels, processInput, changeOutput);
            }
        });
    }

    /**
     * Processes entered inputs with an algorithm to update the outputs.
     * @param algorithm {Object} - Algorithm used to update the user interface.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions.
     * @param processInput {Function} - Function from the algorithm which should process the input.
     * @param changeOutput {Function} - Function from the algorithm which should change the output after processing the input.
     */
    function update(algorithm, viewmodels, processInput, changeOutput) {
        // avoids using not updated values (especially in displayed formulas)
        // for example removeCriticalNumbers(e) needs to have enough time to be executed first
        // on a value change (uses same event: [..].on(change))
        setTimeout(function () {
            processInput(algorithm, inputProcessorInstance, viewmodels.input, viewmodels.visual);
            changeOutput(algorithm.getOutput(), inputProcessorInstance, viewmodels);
            postProcess(viewmodels);
        }, REUPDATE_TIMEOUT_MS);
    }

    /**
     * Does post processing after some kind of input by keyboard or mouse.
     * For example LaTeX-math is updated.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions.
     * @see MathJax has to be executed as last one!
     */
    function postProcess(viewmodels) {
        linkIterationTables(viewmodels.visual);  // iterative tables are not existing from the beginning and so they have to be relinked
        doMinimaAndTracebackHighlighting(viewmodels.visual);
        doTracebackHighlighting(viewmodels.visual);
        MathJax.Hub.Queue(["Typeset", MathJax.Hub]);  // reinterpret new LaTeX code
    }

    /**
     * Linking everything with tables for iterative algorithms.
     * @param visualViewmodel {Object} - Model which is used for example to highlight cells.
     * @see Used for AEP.
     */
    function linkIterationTables(visualViewmodel) {
        if (visualViewmodel.algorithm.type === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER) {
            var mainOutput = $(".main_output");  // output containing the calculation tables

            var calculationTable1 = mainOutput.find(".calculation_1");  // AEP iterations
            var calculationTable2 = mainOutput.find(".calculation_2");
            var calculationTable3 = mainOutput.find(".calculation_3");
            var calculationTable4 = mainOutput.find(".calculation_4");
            var calculationTable5 = mainOutput.find(".calculation_5");

            var iterationTablesArray = [calculationTable1, calculationTable2, calculationTable3, calculationTable4, calculationTable5];
            var selectableEntryClass = ".selectable_entry";

            doMultiTableHighlighting(visualViewmodel, mainOutput, iterationTablesArray);
            linkIterationDownloadLinks(visualViewmodel);

            for (var i = 0; i < iterationTablesArray.length; i++) {
                var functionArguments = {
                    "calculationHorizontalTable": [],
                    "calculationVerticalTable": [],
                    "iterationTablesArray": iterationTablesArray,
                    "mainOutput": mainOutput,
                    "number": -(i + 1),  // iteration numbers are negative in "defaults.js"
                    "selectableEntryClass": selectableEntryClass,
                    "visualViewmodel": visualViewmodel
                };

                iterationTablesArray[i].on("click", selectableEntryClass, functionArguments, selectCell);
            }
        }
    }

    /**
     * Some iterative algorithms has some initial highlighting.
     * @param visualViewmodel {Object} - Model which is used for example to highlight cells.
     * @param iterationTablesArray {Array} - An array of tables.
     * @param mainOutput {Element} - The div containing only the calculation tables.
     */
    function doMultiTableHighlighting(visualViewmodel, mainOutput, iterationTablesArray) {
        // the tables exist after just after the display has updated,
        // but the event is triggered before and so a delay is needed
        setTimeout(function () {
            visualViewmodel.showTraceback(0, undefined, undefined, undefined, iterationTablesArray, mainOutput[0]);
        }, REACTION_TIME_HIGHLIGHT);
    }

    /**
     * Does minima and traceback highlighting on two tables.
     * @param visualViewmodel {Object} - Model which is used for example to highlight cells.
     */
    function doMinimaAndTracebackHighlighting(visualViewmodel) {
        if (visualViewmodel.algorithm.type === ALGORITHMS.HIRSCHBERG) {
            var mainOutput = $(".main_output");  // output containing the calculation tables
            var traceTable = mainOutput.find(".trace_table");
            var tracecellsTable = mainOutput.find(".tracecells_table");

            // the table exist after just after the display has updated,
            // but the event is triggered before and so a delay is needed
            setTimeout(function () {
                visualViewmodel.showTraceback(0, undefined, traceTable[0], undefined, undefined, mainOutput[0]);
                visualViewmodel.markMinima(tracecellsTable[0]);
            }, REACTION_TIME_HIGHLIGHT);

            // fallback if something goes wrong and MathJax disrupts your SVG overlay
            var calculationTable = mainOutput.find(".calculation");
            var calculationHorizontalTable = mainOutput.find(".calculation_horizontal");
            var calculationVerticalTable = mainOutput.find(".calculation_vertical");

            setTimeout(function () {
                visualViewmodel.redrawOverlay(calculationVerticalTable[0], calculationTable[0], calculationHorizontalTable[0], mainOutput[0]);
            }, REACTION_TIME_REDRAWING.LONG_ARROWS.FIRST);  // it can happen that the rendering is disturbed by MathJax and in this case we wait for a while and redraw

            setTimeout(function () {
                visualViewmodel.redrawOverlay(calculationVerticalTable[0], calculationTable[0], calculationHorizontalTable[0], mainOutput[0]);
            }, REACTION_TIME_REDRAWING.LONG_ARROWS.SECOND);  // it can happen that the rendering is disturbed by MathJax and in this case we wait for a while and redraw
        }
    }

    /**
     * Does only traceback highlighting on two tables.
     * @param visualViewmodel {Object} - Model which is used for example to highlight cells.
     */
    function doTracebackHighlighting(visualViewmodel) {
        if (TABLE_INITIAL_HIGHLIGHT_ALGORITHMS.indexOf(visualViewmodel.algorithm.type) >= 0) {
            var mainOutput = $(".main_output");  // output containing the calculation tables
            var calculationTable = mainOutput.find(".calculation");
            var calculationHorizontalTable = mainOutput.find(".calculation_horizontal");
            var calculationVerticalTable = mainOutput.find(".calculation_vertical");
            var results = $(".results");

            setTimeout(function () {
                visualViewmodel.showTraceback(0, calculationVerticalTable[0], calculationTable[0], calculationHorizontalTable[0], undefined, mainOutput[0]);
                visualViewmodel.highlight(0, results[0]);
            }, REACTION_TIME_HIGHLIGHT);

            // fallback if something goes wrong and MathJax disrupts your SVG overlay
            if (SVG_ARROW_ALGORITHMS.indexOf(visualViewmodel.algorithm.type) >= 0) {
                setTimeout(function () {
                    visualViewmodel.redrawOverlay(calculationVerticalTable[0], calculationTable[0], calculationHorizontalTable[0], mainOutput[0]);
                }, REACTION_TIME_REDRAWING.LONG_ARROWS.FIRST);  // it can happen that the rendering is disturbed by MathJax and in this case we wait for a while and redraw

                setTimeout(function () {
                    visualViewmodel.redrawOverlay(calculationVerticalTable[0], calculationTable[0], calculationHorizontalTable[0], mainOutput[0]);
                }, REACTION_TIME_REDRAWING.LONG_ARROWS.SECOND);  // it can happen that the rendering is disturbed by MathJax and in this case we wait for a while and redraw
            }
        }
    }

    /**
     * Redraws the Overlay during loading once.
     * The problem is that MathJax disturbs the positioning of the overlay during the page loading.
     * Hint: window.onload or document.ready cannot be used,
     * because the table is generated after these events have fired with Knockout.
     * @param visualViewmodel {Object} - Model which is used for example to highlight cells.
     * @param calculationVerticalTable {Element} - The table storing the vertical gap costs.
     * @param calculationTable {Element} - The default or main table.
     * @param calculationHorizontalTable {Element} - The table storing the horizontal gap costs.
     * @param mainOutput {Element} - The div containing only the calculation tables.
     */
    function redrawInitialOverlay(visualViewmodel, calculationVerticalTable, calculationTable, calculationHorizontalTable, mainOutput) {
        if (TABLE_INITIAL_HIGHLIGHT_ALGORITHMS.indexOf(visualViewmodel.algorithm.type) >= 0 &&
            SVG_ARROW_ALGORITHMS.indexOf(visualViewmodel.algorithm.type) >= 0) {

            setTimeout(function () {
                visualViewmodel.redrawOverlay(calculationVerticalTable[0], calculationTable[0], calculationHorizontalTable[0], mainOutput[0]);
            }, REACTION_TIME_REDRAWING.FIRST);  // it can happen that the rendering is disturbed by MathJax and in this case we wait for a while and redraw

            setTimeout(function () {
                visualViewmodel.redrawOverlay(calculationVerticalTable[0], calculationTable[0], calculationHorizontalTable[0], mainOutput[0]);
            }, REACTION_TIME_REDRAWING.SECOND);  // it can happen that the rendering is disturbed by MathJax and in this case we wait for a while and redraw
        }
    }

    /**
     * Linking table download links of tables from several iterations with a download function.
     * @param visualViewmodel {Object} - Model which is used for example to highlight cells.
     */
    function linkIterationDownloadLinks(visualViewmodel) {
        for (var i = 0; i < MAX_NUMBER_ITERATIONS; i++) {
            var currentTableNumber = (i + 1);
            linkIterationDownloadLink(
                $(".table_download_" + currentTableNumber),
                $(".calculation_" + currentTableNumber),
                -currentTableNumber,  // iteration numbers are negative in "defaults.js"
                visualViewmodel);
        }
    }

    /**
     * Linking table download link of a table generated during an iterative algorithm with a download function.
     * @param download {Element} - The link-element <a>...</a> which is to be linked.
     * @param table {Element} - The table you want download.
     * @param number {number} - The table id.
     * @param visualViewmodel {Object} - Model which is used for example to highlight cells.
     */
    function linkIterationDownloadLink(download, table, number, visualViewmodel) {
        download.on("click", {"calculationTable": table, "number": number}, visualViewmodel.downloadTable);
    }

    /**
     * Allows to redraw an overlay on a scroll or resize event.
     * @param visualViewmodel {Object} - Model which is used for example to highlight cells.
     * @param calculationVerticalTable {Element} - The table storing the vertical gap costs.
     * @param calculationTable {Element} - The default or main table.
     * @param calculationHorizontalTable {Element} - The table storing the horizontal gap costs.
     * @param mainOutput {Element} - The div containing only the calculation tables.
     */
    function linkOverlay(visualViewmodel, calculationVerticalTable, calculationTable, calculationHorizontalTable, mainOutput) {
        if (SVG_ARROW_ALGORITHMS.indexOf(visualViewmodel.algorithm.type) >= 0) {  // algorithm uses SVG-arrows
            var browserWindow = $(window);

            var functionArguments = {
                "calculationTable": calculationTable,
                "calculationHorizontalTable": calculationHorizontalTable,
                "calculationVerticalTable": calculationVerticalTable,
                "mainOutput": mainOutput,
                "visualViewmodel": visualViewmodel
            };

            browserWindow.on("resize", functionArguments, reinitialize);
            mainOutput.on("scroll", functionArguments, reinitialize);
        }
    }

    /**
     * Reinitializes overlays after an event wih the browser window or a scrolling event of a table.
     * @param e {Object} - Stores data relevant to the event called that function.
     */
    function reinitialize(e) {
        var visualViewmodel = e.data.visualViewmodel;

        var mainOutput = e.data.mainOutput[0];
        var calculationVerticalTable;
        var calculationTable = e.data.calculationTable[0];
        var calculationHorizontalTable;

        if (e.data.calculationVerticalTable !== undefined) {
            calculationVerticalTable = e.data.calculationVerticalTable[0];
            calculationHorizontalTable = e.data.calculationHorizontalTable[0];
        }

        visualViewmodel.redrawOverlay(calculationVerticalTable, calculationTable, calculationHorizontalTable, mainOutput);
    }

    /**
     * Linking all elements which can be selected (results, cells, ...).
     * @param visualViewmodel {Object} - Model which is used for example to highlight cells.
     * @param calculationVerticalTable {Element} - The table storing the vertical gap costs.
     * @param calculationTable {Element} - The default or main table.
     * @param calculationHorizontalTable {Element} - The table storing the horizontal gap costs.
     * @param mainOutput {Element} - The div containing only the calculation tables.
     * @param selectableEntryClass {Object} - The class name of a selectable entry.
     */
    function linkSelectables(visualViewmodel, calculationVerticalTable, calculationTable, calculationHorizontalTable,
                             mainOutput, selectableEntryClass) {

        var results = $(".results");

        var functionArguments = {
            "calculationTable": calculationTable,
            "calculationHorizontalTable": calculationHorizontalTable,
            "calculationVerticalTable": calculationVerticalTable,
            "resultsTable": results,
            "mainOutput": mainOutput,
            "selectableEntryClass": selectableEntryClass,
            "visualViewmodel": visualViewmodel
        };
        results.on("click", selectableEntryClass, functionArguments, selectTableEntry);

        functionArguments = {
            "calculationTable": calculationTable,
            "calculationHorizontalTable": calculationHorizontalTable,
            "calculationVerticalTable": calculationVerticalTable,
            "mainOutput": mainOutput,
            "number": MATRICES.VERTICAL_NUMBER,
            "selectableEntryClass": selectableEntryClass,
            "visualViewmodel": visualViewmodel
        };
        calculationVerticalTable.on("click", selectableEntryClass, functionArguments, selectCell);

        functionArguments = {
            "calculationTable": calculationTable,
            "calculationHorizontalTable": calculationHorizontalTable,
            "calculationVerticalTable": calculationVerticalTable,
            "mainOutput": mainOutput,
            "number": MATRICES.DEFAULT_NUMBER,
            "selectableEntryClass": selectableEntryClass,
            "visualViewmodel": visualViewmodel
        };
        calculationTable.on("click", selectableEntryClass, functionArguments, selectCell);

        functionArguments = {
            "calculationTable": calculationTable,
            "calculationHorizontalTable": calculationHorizontalTable,
            "calculationVerticalTable": calculationVerticalTable,
            "mainOutput": mainOutput,
            "number": MATRICES.HORIZONTAL_NUMBER,
            "selectableEntryClass": selectableEntryClass,
            "visualViewmodel": visualViewmodel
        };

        calculationHorizontalTable.on("click", selectableEntryClass, functionArguments, selectCell);
    }

    /**
     * Selects a table entry and triggers a function of the visualizer to highlight the entry.
     * @param e {Object} - Stores data relevant to the event called that function.
     */
    function selectTableEntry(e) {
        // retrieve data
        var mainOutput = e.data.mainOutput[0];
        var calculationVerticalTable;
        var calculationTable = e.data.calculationTable[0];
        var calculationHorizontalTable;

        if (e.data.calculationVerticalTable !== undefined) {
            calculationVerticalTable = e.data.calculationVerticalTable[0];
            calculationHorizontalTable = e.data.calculationHorizontalTable[0];
        }

        var resultsTable = e.data.resultsTable;
        var selectableEntries = resultsTable.find(e.data.selectableEntryClass);
        var visualViewmodel = e.data.visualViewmodel;

        // compute the selected position in the results table
        var selectedRow = -1;

        for (var i = 0; i < selectableEntries.length; i++) {
            if (this === selectableEntries[i])
                selectedRow = i;
        }

        visualViewmodel.highlight(selectedRow, resultsTable[0]);
        visualViewmodel.showTraceback(selectedRow, calculationVerticalTable, calculationTable, calculationHorizontalTable,
            undefined, mainOutput);
    }

    /**
     * Highlights a table cell and its neighbours by the help of a visualizer function
     * which shows the neighbours that were needed to compute its cell value.
     * @param e {Object} - Stores data relevant to the event called that function.
     */
    function selectCell(e) {
        // retrieve data
        var calculationHorizontalTable = e.data.calculationHorizontalTable;
        var calculationTable = e.data.calculationTable;
        var calculationVerticalTable = e.data.calculationVerticalTable;
        var visualViewmodel = e.data.visualViewmodel;

        var number = e.data.number;
        var mainOutput = e.data.mainOutput[0];
        var iterationTablesArray = e.data.iterationTablesArray;

        if (iterationTablesArray !== undefined)  // if table from an iteration
            calculationTable = iterationTablesArray[-(number + 1)];  // iteration number are negative in "defaults.js"

        var currentSelectedTable;
        var label;

        if (number === MATRICES.VERTICAL_NUMBER) {
            currentSelectedTable = calculationVerticalTable;
            label = MATRICES.VERTICAL;
        } else if (number === MATRICES.HORIZONTAL_NUMBER) {
            currentSelectedTable = calculationHorizontalTable;
            label = MATRICES.HORIZONTAL;
        } else {  // if (number === MATRICES.ITERATION_NUMBER_i || number === MATRICES.DEFAULT)
            currentSelectedTable = calculationTable;
            label = MATRICES.DEFAULT;
        }

        var selectableEntryClass = e.data.selectableEntryClass;
        var selectableEntries = currentSelectedTable.find(selectableEntryClass);

        // compute the selected position in the calculation table
        var selectedColumn = -1;
        var selectedRow = -1;

        var tableHeight = currentSelectedTable[0].rows.length - 1;
        var tableWidth = currentSelectedTable[0].rows[0].cells.length - 1;

        for (var i = 0; i < tableHeight; i++) {
            for (var j = 0; j < tableWidth; j++) {

                if (selectableEntries[i * tableWidth + j] === this) {
                    selectedColumn = j;
                    selectedRow = i;
                }
            }
        }

        // store selected position in a vector
        var cellCoordinates = new bases.alignment.Vector(selectedRow, selectedColumn);
        cellCoordinates.label = label;

        visualViewmodel.showFlow(cellCoordinates,
            calculationVerticalTable[0], calculationTable[0], calculationHorizontalTable[0], iterationTablesArray,
            mainOutput, -(number + 1));  // MATRICES.ITERATION_NUMBER_i are negative
    }

    /**
     * Post edits a matrix and replaces for example values with LaTeX-symbols.
     * @param matrix {matrix} - The matrix in which you want replace values with for example LaTeX-symbols.
     * @param visualViewmodel {Object} - The model for visualization which is replacing symbols.
     * @return {matrix} - The matrix in which symbols where replaced with LaTeX-symbols.
     */
    function postEdit(matrix, visualViewmodel) {
        return visualViewmodel.replaceInfinityStrings(matrix);
    }
}());