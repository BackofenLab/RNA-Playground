/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("interfaces.affineAlignmentInterface", AffineAlignmentInterface, startAffineAlignmentAlgorithm);

    // instances
    var alignmentInterfaceInstance;
    var affineAlignmentInterfaceInstance;

    function AffineAlignmentInterface() {
        affineAlignmentInterfaceInstance = this;

        // inheritance
        alignmentInterfaceInstance = new interfaces.alignmentInterface.AlignmentInterface();

        // public methods (linking)
        this.startAffineAlignmentAlgorithm = startAffineAlignmentAlgorithm;
    }

    function startAffineAlignmentAlgorithm(Algorithm) {
        imports();

        var inputViewmodel = new InputViewmodel();
        alignmentInterfaceInstance.sharedInterfaceOperations(Algorithm, inputViewmodel, processInput, changeOutput);
    }

    /**
     * Manages all imported scripts.
     */
    function imports() {
        MathJax.Hub.Queue(["Typeset", MathJax.Hub]);  // reinterpret new LaTeX code
        alignmentInterfaceInstance.imports();

        // interfaces
        $.getScript(PATHS.ALIGNMENT_INTERFACE);
    }

    /*---- INPUT ----*/
    /**
     * In the Model-View-Viewmodel, the view (HTML-page) is filled with data from
     * outside with the help of the viewmodel (here: InputViewmodel)
     * by getting data from a model (here: AFFINE_ALIGNMENT_DEFAULTS).
     * @see https://en.wikipedia.org/wiki/Model-view-viewmodel
     * @constructor
     */
    function InputViewmodel() {
        this.sequence1 = ko.observable(AFFINE_ALIGNMENT_DEFAULTS.SEQUENCE_1);
        this.sequence2 = ko.observable(AFFINE_ALIGNMENT_DEFAULTS.SEQUENCE_2);

        this.calculation = ko.observable(AFFINE_ALIGNMENT_DEFAULTS.CALCULATION);

        // function
        this.baseCosts = ko.observable(AFFINE_ALIGNMENT_DEFAULTS.FUNCTION.BASE_COSTS);
        this.enlargement = ko.observable(AFFINE_ALIGNMENT_DEFAULTS.FUNCTION.ENLARGEMENT);
        this.match = ko.observable(AFFINE_ALIGNMENT_DEFAULTS.FUNCTION.MATCH);
        this.mismatch = ko.observable(AFFINE_ALIGNMENT_DEFAULTS.FUNCTION.MISMATCH);
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

            inputViewmodel.baseCosts(Number($("#base_costs").val()));
            inputViewmodel.enlargement(Number($("#enlargement").val()));
            inputViewmodel.match(Number($("#match").val()));
            inputViewmodel.mismatch(Number($("#mismatch").val()));
        } else
            inputProcessor.activateInputUpdates();

        alignmentInterfaceInstance.startProcessing(algorithm, inputViewmodel, visualViewmodel);
    }

    function changeOutput(outputData, inputProcessor, viewmodels) {
        viewmodels.output.matrix(outputData.matrix);
        viewmodels.output.horizontalGaps(edit(inputProcessor, outputData.horizontalGaps, viewmodels));
        viewmodels.output.verticalGaps(edit(inputProcessor, outputData.verticalGaps, viewmodels));

        for (var i = 0; i < outputData.matrix.length; i++) {
            // new variables (rows) are not automatically functions
            // and so we have to convert new variables manually into functions
            // or we get the following error
            // 'Uncaught TypeError: inputOutputViewmodel.output.tableValues[i] is not a function'
            if (i > viewmodels.output.matrix.length) {
                viewmodels.output.matrix[i] = new Function();
                viewmodels.output.verticalGaps[i] = new Function();
                viewmodels.output.horizontalGaps[i] = new Function();
            }

            viewmodels.output.matrix[i](outputData.matrix[i]);
            viewmodels.output.verticalGaps[i](outputData.verticalGaps[i]);
            viewmodels.output.horizontalGaps[i](outputData.horizontalGaps[i]);
        }

        viewmodels.output.alignments(outputData.alignments);
        viewmodels.output.score(outputData.score);
        viewmodels.output.moreTracebacks(outputData.moreTracebacks);
    }

    function edit(inputProcessor, matrix, viewmodels) {
        return inputProcessor.postEdit(matrix, viewmodels.visual);
    }
}());
