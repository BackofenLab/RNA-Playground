/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("interfaces.alignmentInterface", AlignmentInterface, imports, sharedInterfaceOperations, startProcessing);

    // instances
    var alignmentInterfaceInstance;

    function AlignmentInterface() {
        alignmentInterfaceInstance = this;
        this.imports = imports;
        this.sharedInterfaceOperations = sharedInterfaceOperations;
        this.startAlignmentAlgorithm = startAlignmentAlgorithm;
        this.startProcessing = startProcessing;
    }

    function startAlignmentAlgorithm(Algorithm) {
        imports();

        var inputViewmodel = new InputViewmodel();
        sharedInterfaceOperations(Algorithm, inputViewmodel, processInput, changeOutput);
    }

    /**
     * Manages all imported scripts.
     */
    function imports() {
        // third party libs
        MathJax.Hub.Queue(["Typeset", MathJax.Hub]);  // reinterpret new LaTeX code
        $.getScript(PATHS.LIBS.KNOCKOUT);  // to make knockout working whenever page is reloaded

        // design/controls logic
        $.getScript(PATHS.INPUT_PROCESSOR);
        $.getScript(PATHS.VISUALIZER);

        // algorithms logic
        $.getScript(PATHS.BACKTRACKING);
    }

    function sharedInterfaceOperations(Algorithm, inputViewmodel, processInput, changeOutput) {
        var visualViewmodel = new postProcessing.visualizer.Visualizer();

        var algorithm = new Algorithm();
        var inputProcessor = new postProcessing.inputProcessor.InputProcessor();
        processInput(algorithm, inputProcessor, inputViewmodel, visualViewmodel);
        var outputViewmodel = new OutputViewmodel(algorithm.type,
            edit(algorithm.type, inputProcessor, algorithm.getOutput(), visualViewmodel));

        var viewmodels = {
            input: inputViewmodel,
            visual: visualViewmodel,
            output: outputViewmodel
        };

        linkInputWithOutput(algorithm, viewmodels, inputProcessor, processInput, changeOutput);

        ko.applyBindings(viewmodels, document.getElementById("algorithm_view"));
    }

    /*---- INPUT ----*/
    /**
     * In the Model-View-Viewmodel, the view (HTML-page) is filled with data from
     * outside with the help of the viewmodel (here: InputViewmodel)
     * by getting data from a model (here: ALIGNMENT_DEFAULTS).
     * @see https://en.wikipedia.org/wiki/Model-view-viewmodel
     * @constructor
     */
    function InputViewmodel() {
        this.sequence1 = ko.observable(ALIGNMENT_DEFAULTS.SEQUENCE_1);
        this.sequence2 = ko.observable(ALIGNMENT_DEFAULTS.SEQUENCE_2);

        this.calculation = ko.observable(ALIGNMENT_DEFAULTS.CALCULATION);

        // function
        this.deletion = ko.observable(ALIGNMENT_DEFAULTS.FUNCTION.DELETION);
        this.insertion = ko.observable(ALIGNMENT_DEFAULTS.FUNCTION.INSERTION);
        this.match = ko.observable(ALIGNMENT_DEFAULTS.FUNCTION.MATCH);
        this.mismatch = ko.observable(ALIGNMENT_DEFAULTS.FUNCTION.MISMATCH);
    }

    function processInput(algorithm, inputProcessor, inputViewmodel, visualViewmodel) {
        visualViewmodel.removeAllContents();

        // when page was loaded the inputs have not to be updated or you get wrong inputs
        if (inputProcessor.inputUpdatesActivated()) {
            inputViewmodel.sequence1($("#sequence_1").val());
            inputViewmodel.sequence2($("#sequence_2").val());

            if ($("#distance").is(":checked"))
                inputViewmodel.calculation(ALIGNMENT_TYPES.DISTANCE);
            else
                inputViewmodel.calculation(ALIGNMENT_TYPES.SIMILARITY);

            inputViewmodel.deletion(Number($("#deletion").val()));
            inputViewmodel.insertion(Number($("#insertion").val()));
            inputViewmodel.match(Number($("#match").val()));
            inputViewmodel.mismatch(Number($("#mismatch").val()));
        } else
            inputProcessor.activateInputUpdates();

        startProcessing(algorithm, inputViewmodel, visualViewmodel);
    }

    function startProcessing(algorithm, inputViewmodel, visualViewmodel) {
        algorithm.setInput(inputViewmodel);
        var ioData = algorithm.compute();
        visualViewmodel.shareInformation(algorithm, ioData[0], ioData[1]);
    }

    function linkInputWithOutput(algorithm, viewmodels, inputProcessor, processInput, changeOutput) {
        inputProcessor.linkElements(viewmodels.visual);
        inputProcessor.updateGUI(algorithm, viewmodels, processInput, changeOutput);
    }

    function changeOutput(outputData, inputProcessor, viewmodels) {
        viewmodels.output.matrix(outputData.matrix);

        for (var i = 0; i < outputData.matrix.length; i++) {
            // new variables (rows) are not automatically functions
            // and so we have to convert new variables manually into functions
            // or we get the following error
            // 'Uncaught TypeError: inputOutputViewmodel.output.tableValues[i] is not a function'
            if (i > viewmodels.output.matrix.length)
                viewmodels.output.matrix[i] = new Function();

            viewmodels.output.matrix[i](outputData.matrix[i]);
        }

        viewmodels.output.alignments(outputData.alignments);
        viewmodels.output.score(outputData.score);
        viewmodels.output.moreTracebacks(outputData.moreTracebacks);
    }

    /*---- OUTPUT ----*/
    function OutputViewmodel(algorithmName, outputData) {
        var viewModel = this;

        this.matrix = ko.observableArray(outputData.matrix);

        for (var i = 0; i < outputData.matrix.length; i++)
            this.matrix[i] = ko.observableArray(outputData.matrix[i]);


        if (algorithmName === ALGORITHMS.GOTOH) {
            this.horizontalGaps = ko.observableArray(outputData.horizontalGaps);
            this.verticalGaps = ko.observableArray(outputData.verticalGaps);

            for (var i = 0; i < outputData.matrix.length; i++) {
                this.horizontalGaps[i] = ko.observableArray(outputData.horizontalGaps[i]);
                this.verticalGaps[i] = ko.observableArray(outputData.verticalGaps[i]);
            }
        }

        this.alignments = ko.observableArray(outputData.alignments);

        this.score = ko.observable(outputData.score);
        this.moreTracebacks = ko.observable(outputData.moreTracebacks);
    }

    function edit(algorithmName, inputProcessor, outputData, visualViewmodel) {
        if (algorithmName === ALGORITHMS.GOTOH) {
            outputData.horizontalGaps = inputProcessor.postEdit(outputData.horizontalGaps, visualViewmodel);
            outputData.verticalGaps = inputProcessor.postEdit(outputData.verticalGaps, visualViewmodel);
        }

        return outputData;
    }
}());
