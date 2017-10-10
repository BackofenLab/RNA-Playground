/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods (declaration)
    namespace("postProcessing.inputProcessor",
        InputProcessor, activateInputUpdates, inputUpdatesActivated, linkElements, updateGUI, postEdit);

    // instances
    var inputProcessorInstance;

    /**
     * Does all the post-processing
     * like removing/updating wrong inputs
     * and it defines how different
     * input types behave.
     */
    function InputProcessor() {
        inputProcessorInstance = this;

        // variables
        this.inputUpdatesStarted = false;
        this.avoidFocusOutUpdate = false;

        // public methods (linking)
        this.activateInputUpdates = activateInputUpdates;
        this.inputUpdatesActivated = inputUpdatesActivated;
        this.linkElements = linkElements;
        this.updateGUI = updateGUI;
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
     * Linking static elements with a function.
     * @param visualViewmodel {Object} - Model which is used for example to highlight cells.
     */
    function linkElements(visualViewmodel) {
        fixBrowserBugs();
        changeDefaultKeyBehaviour();

        var algorithmInput = $("#algorithm_input");
        var functionParameters = algorithmInput.find(".fx_parameter");

        var mainOutput = $("#main_output");
        var calculation = mainOutput.find("#calculation");
        var calculationHorizontal = mainOutput.find("#calculation_horizontal");
        var calculationVertical = mainOutput.find("#calculation_vertical");
        var results = $("#results");

        var selectableEntryClass = ".selectable_entry";

        var tableDownload = $(".table_download");
        var tableVerticalDownload = $(".table_vertical_download");
        var tableHorizontalDownload = $(".table_horizontal_download");

        // 1st
        var functionArguments = {"functionParameters": functionParameters};
        algorithmInput.find(".optimization_type").on("change", functionArguments, negateOptimizationParameters);

        var browserWindow = $(window);

        // 2nd
        functionParameters.on("focusout keypress", removeCriticalNumbers);
        algorithmInput.find(".sequence").on("keyup", removeNonAllowedBases);

        // 3rd
        functionArguments = {
            "calculationTable": calculation,
            "calculationHorizontalTable": calculationHorizontal,
            "calculationVerticalTable": calculationVertical,
            "resultsTable": results,
            "selectableEntryClass": selectableEntryClass,
            "visualViewmodel": visualViewmodel
        };
        results.on("click", selectableEntryClass, functionArguments, selectTableEntry);

        // 4th
        functionArguments = {
            "calculationTable": calculation,
            "calculationHorizontalTable": calculationHorizontal,
            "calculationVerticalTable": calculationVertical,
            "selectableEntryClass": selectableEntryClass,
            "visualViewmodel": visualViewmodel,
            "number": MATRICES.VERTICAL_NUMBER
        };
        calculationVertical.on("click", selectableEntryClass, functionArguments, selectCell);

        functionArguments = {
            "calculationTable": calculation,
            "calculationHorizontalTable": calculationHorizontal,
            "calculationVerticalTable": calculationVertical,
            "selectableEntryClass": selectableEntryClass,
            "visualViewmodel": visualViewmodel,
            "number": MATRICES.DEFAULT_NUMBER
        };
        calculation.on("click", selectableEntryClass, functionArguments, selectCell);

        functionArguments = {
            "calculationTable": calculation,
            "calculationHorizontalTable": calculationHorizontal,
            "calculationVerticalTable": calculationVertical,
            "selectableEntryClass": selectableEntryClass,
            "visualViewmodel": visualViewmodel,
            "number": MATRICES.HORIZONTAL_NUMBER
        };
        calculationHorizontal.on("click", selectableEntryClass, functionArguments, selectCell);

        // 5th
        functionArguments = {
            "calculationTable": calculation,
            "calculationHorizontalTable": calculationHorizontal,
            "calculationVerticalTable": calculationVertical,
            "number": MATRICES.VERTICAL_NUMBER
        };
        tableVerticalDownload.on("click", functionArguments, visualViewmodel.downloadTable);

        functionArguments = {
            "calculationTable": calculation,
            "calculationHorizontalTable": calculationHorizontal,
            "calculationVerticalTable": calculationVertical,
            "number": MATRICES.DEFAULT_NUMBER
        };
        tableDownload.on("click", functionArguments, visualViewmodel.downloadTable);

        functionArguments = {
            "calculationTable": calculation,
            "calculationHorizontalTable": calculationHorizontal,
            "calculationVerticalTable": calculationVertical,
            "number": MATRICES.HORIZONTAL_NUMBER
        };
        tableHorizontalDownload.on("click", functionArguments, visualViewmodel.downloadTable);

        // 6th
        functionArguments = {
            "calculationTable": calculation,
            "calculationHorizontalTable": calculationHorizontal,
            "calculationVerticalTable": calculationVertical
        };
        browserWindow.on("resize", functionArguments, visualViewmodel.redrawOverlay);
    }

    /**
     * Bug fixes for specific browsers or versions of browsers.
     */
    function fixBrowserBugs() {
        /*
         BUG-FIX for Firefox:
         Inputs of type "number" doesn't get the focus
         in the Firefox browser if one
         of the up- or down-button of the inputs
         is clicked.

         Detected: Firefox 55.0b3 (32-Bit)
         */
        $(function () {
            $("input[type='number']").on("click", function () {
                $(this).focus();
            });
        });
    }

    /**
     * Redefines what for example should happen
     * when a specific key is pressed or on a specific action like pasting.
     */
    function changeDefaultKeyBehaviour() {
        var input = $("input");
        /*
         If "Enter" is pressed the UI should't update immediately.
         The time after which the UI is updated is specified
         with the global default "REUPDATE_TIMEOUT_MS".
         */
        input.keypress(function (e) {
            // Hint: Better would be a tabulator behaviour,
            // but the implementation of that is too complex.
            if (e.which === KEY_CODES.ENTER)
                e.preventDefault();
        });
    }

    /**
     * Negates the function parameters of an algorithm.
     * Hint: The input has to be of class "optimization_type".
     * @param e - Stores data relevant to the event called that function.
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
     * @param e - Stores data relevant to the event called that function.
     */
    function removeCriticalNumbers(e) {
        if (CHARACTER.NUMBERS.test(this.value))
            this.value = Math.round(this.value);

        if (e.which === KEY_CODES.ENTER || e.type === "focusout") {
            this.value = this.value >= INPUT.MIN ? this.value : INPUT.MIN;
            this.value = this.value <= INPUT.MAX ? this.value : INPUT.MAX;
        }
    }

    /**
     * Removes non-english characters from an input-field.
     * Hint: The input has to be of class "sequence".
     */
    function removeNonAllowedBases() {
        if (!CHARACTER.BASES.test(this.value))
            this.value = this.value.replace(CHARACTER.NON_BASES, SYMBOLS.EMPTY);
    }

    /**
     * Selects a table entry and triggers a function of the visualizer to highlight the entry.
     * @param e - Stores data relevant to the event called that function.
     */
    function selectTableEntry(e) {
        // retrieve data
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

        // some delay without it won't work properly
        setTimeout(function () {
            visualViewmodel.highlight(selectedRow, resultsTable[0]);
            visualViewmodel.showTraceback(selectedRow, calculationVerticalTable, calculationTable, calculationHorizontalTable);
        }, REACTION_TIME_HIGHLIGHT);
    }

    /**
     * Highlights a table cell and its neighbours by the help of a visualizer function
     * which shows the neighbours that were needed to compute its cell value.
     * @param e - Stores data relevant to the event called that function.
     */
    function selectCell(e) {
        // retrieve data
        var number = e.data.number;
        var calculationVerticalTable;
        var calculationTable = e.data.calculationTable;
        var calculationHorizontalTable;

        if (e.data.calculationVerticalTable !== undefined) {
            calculationVerticalTable = e.data.calculationVerticalTable;
            calculationHorizontalTable = e.data.calculationHorizontalTable;
        }

        var currentSelectedTable;
        var label;

        if (number === 0) {
            currentSelectedTable = calculationVerticalTable;
            label = MATRICES.VERTICAL;
        }
        else if (number === 1) {
            currentSelectedTable = calculationTable;
            label = MATRICES.DEFAULT;
        }
        else {  // if (number === 2)
            currentSelectedTable = calculationHorizontalTable;
            label = MATRICES.HORIZONTAL;
        }

        var selectableEntryClass = e.data.selectableEntryClass;
        var selectableEntries = currentSelectedTable.find(selectableEntryClass);
        var visualViewmodel = e.data.visualViewmodel;

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

        // some delay without it won't work properly
        setTimeout(function () {
            visualViewmodel.showFlow(cellCoordinates,
                calculationVerticalTable[0], calculationTable[0], calculationHorizontalTable[0]);
        }, REACTION_TIME_HIGHLIGHT);
    }

    /**
     * Updates the user interfaces of the loaded page.
     * @param algorithm {Object} - Algorithm used to update the user interface.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions.
     * @param processInput {Function} - Function from the algorithm which should process the input.
     * @param changeOutput {Function} - Function from the algorithm which should change the output after processing the input.
     */
    function updateGUI(algorithm, viewmodels, processInput, changeOutput) {
        var inputs = $("#algorithm_input").find("input");

        inputs.on({
            click: function () {
                if ($(this).is(":radio"))
                    updateAfterTimeout(algorithm, viewmodels, processInput, changeOutput);
            },

            focusout: function () {
                if (!inputProcessorInstance.avoidFocusOutUpdate) {
                    // Important: if "alert" is used in this function then it will fire this event forever in Chrome
                    // because the alert-box gets the focus and loses it when you click on "OK"
                    // Solution: use "console.log(..)"
                    if (!$(this).is(":radio"))  // a radio button should only fire on a click or you get bugs!
                        updateAfterTimeout(algorithm, viewmodels, processInput, changeOutput);
                }
                inputProcessorInstance.avoidFocusOutUpdate = false;
            },

            keypress: function () {
                if (e.which === KEY_CODES.ENTER)
                    updateAfterTimeout(algorithm, viewmodels, processInput, changeOutput);
            },

            mousedown: function () {
                if (this.className === "fx_parameter"
                    && this !== document.activeElement
                    && document.activeElement !== document.body)
                    inputProcessorInstance.avoidFocusOutUpdate = true;
                else
                    inputProcessorInstance.avoidFocusOutUpdate = false;
                postProcess();
            }
        });
    }

    /**
     * Does post processing after some kind of input by keyboard or mouse.
     * @example
     * LaTeX-math is updated or wrong characters are removed.
     */
    function postProcess() {
        removeWrongBases();
        MathJax.Hub.Queue(["Typeset", MathJax.Hub]);  // reinterpret new LaTeX code
    }

    /**
     * Removes non-english characters from an input-field.
     */
    function removeWrongBases() {
        var textElements = $(".sequence");

        for (var i = 0; i < textElements.length; i++)
            textElements[i].value = textElements[i].value.replace(CHARACTER.NON_BASES, SYMBOLS.EMPTY);
    }

    /**
     * Processes entered inputs with an algorithm to update the outputs.
     * @param algorithm {Object} - Algorithm used to update the user interface.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions.
     * @param processInput {Function} - Function from the algorithm which should process the input.
     * @param changeOutput {Function} - Function from the algorithm which should change the output after processing the input.
     */
    function updateAfterTimeout(algorithm, viewmodels, processInput, changeOutput) {
        setTimeout(function () {  // to avoid using not updated values
            processInput(algorithm, inputProcessorInstance, viewmodels.input, viewmodels.visual);
            changeOutput(algorithm.getOutput(), inputProcessorInstance, viewmodels);
            postProcess();
        }, REUPDATE_TIMEOUT_MS);
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