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

    /**
     * Is used to work with the input and output (the interface) of an affine alignment algorithm.
     * @constructor
     * @augments AlignmentInterface
     */
    function AffineAlignmentInterface() {
        affineAlignmentInterfaceInstance = this;

        // inheritance
        alignmentInterfaceInstance = new interfaces.alignmentInterface.AlignmentInterface();

        // public methods (linking)
        this.startAffineAlignmentAlgorithm = startAffineAlignmentAlgorithm;
    }

    /**
     * Function managing objects.
     */
    function startAffineAlignmentAlgorithm(Algorithm) {
        imports();

        var inputViewmodel = new InputViewmodel();
        alignmentInterfaceInstance.sharedInterfaceOperations(Algorithm, inputViewmodel, processInput, changeOutput);
    }

    /**
     * Handling imports.
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
        var viewmodel = this;

        this.sequence1 = ko.observable(AFFINE_ALIGNMENT_DEFAULTS.SEQUENCE_1);
        this.sequence2 = ko.observable(AFFINE_ALIGNMENT_DEFAULTS.SEQUENCE_2);

        this.calculation = ko.observable(AFFINE_ALIGNMENT_DEFAULTS.CALCULATION);

        // function
        this.baseCosts = ko.observable(AFFINE_ALIGNMENT_DEFAULTS.FUNCTION.BASE_COSTS);
        this.enlargement = ko.observable(AFFINE_ALIGNMENT_DEFAULTS.FUNCTION.ENLARGEMENT);
        this.match = ko.observable(AFFINE_ALIGNMENT_DEFAULTS.FUNCTION.MATCH);
        this.mismatch = ko.observable(AFFINE_ALIGNMENT_DEFAULTS.FUNCTION.MISMATCH);

        // displayed dynamic formulas
        this.gapStart = ko.computed(
            function () {
                return Number(viewmodel.baseCosts()) + Number(viewmodel.enlargement());
            }
        );

        this.formula = ko.computed(
            function getSelectedFormula() {
                // to fire LaTeX-Code reinterpretation after the selected formula was changed
                // HINT: only found solution which works on all browsers
                setTimeout(function () {
                    MathJax.Hub.Queue(["Typeset", MathJax.Hub])
                }, REUPDATE_TIMEOUT_MS);

                return getFormula(viewmodel, MATRICES.DEFAULT);
            }
        );

        this.formulaP = ko.computed(
            function getSelectedFormula() {
                // to fire LaTeX-Code reinterpretation after the selected formula was changed
                // HINT: only found solution which works on all browsers
                setTimeout(function () {
                    MathJax.Hub.Queue(["Typeset", MathJax.Hub])
                }, REUPDATE_TIMEOUT_MS);

                return getFormula(viewmodel, MATRICES.VERTICAL);
            }
        );

        this.formulaQ = ko.computed(
            function getSelectedFormula() {
                // to fire LaTeX-Code reinterpretation after the selected formula was changed
                // HINT: only found solution which works on all browsers
                setTimeout(function () {
                    MathJax.Hub.Queue(["Typeset", MathJax.Hub])
                }, REUPDATE_TIMEOUT_MS);

                return getFormula(viewmodel, MATRICES.HORIZONTAL);
            }
        );

        this.gapFunction = ko.computed(
            function getSelectedFormula() {
                // to fire LaTeX-Code reinterpretation after the selected formula was changed
                // HINT: only found solution which works on all browsers
                setTimeout(function () {
                    MathJax.Hub.Queue(["Typeset", MathJax.Hub])
                }, REUPDATE_TIMEOUT_MS);

                return getSubformula(viewmodel);
            }
        );
    }

    /**
     * Returns the LaTeX-code for formulas of affine algorithms.
     * @param viewmodel {InputViewmodel} - The viewmodel of the view displaying the formula.
     * @param matrix {string} - The matrix for which you want display the formula.
     * @return {string} - LaTeX code.
     */
    function getFormula(viewmodel, matrix) {
        var string = LATEX.MATH_REGION;  // starting LaTeX math region

            // differentiate between formulas and writing code for static one
            switch (matrix) {
                case MATRICES.VERTICAL:
                    string += LATEX.FORMULA.CURRENT_P;
                    break;
                case MATRICES.HORIZONTAL:
                    string += LATEX.FORMULA.CURRENT_Q;
                    break;
                default:
                    string += LATEX.FORMULA.CURRENT;
            }

            // look if we maximize or minimize
            if (viewmodel.calculation() === ALIGNMENT_TYPES.SIMILARITY)
                string += SYMBOLS.EQUAL + LATEX.MAX;
            else
                string += SYMBOLS.EQUAL + LATEX.MIN;

        // differentiate between formulas and writing code for dynamic one
            switch (matrix) {
                case MATRICES.VERTICAL:
                    string += LATEX.RECURSION.GOTOH_P;
                    string += getDynamicFormula(viewmodel, matrix);
                    break;
                case MATRICES.HORIZONTAL:
                    string += LATEX.RECURSION.GOTOH_Q;
                    string += getDynamicFormula(viewmodel, matrix);
                    break;
                default:
                    string += LATEX.RECURSION.GOTOH;
                    string += getDynamicFormula(viewmodel, matrix);
            }

        string += LATEX.MATH_REGION;  // stopping LaTeX math region
        return string;
    }

    /**
     * Returns the LaTeX-code for dynamic, input-dependant formulas of affine algorithms.
     * @param viewmodel {InputViewmodel} - The viewmodel of the view displaying the formula.
     * @param matrix {string} - The matrix for which you want display the formula.
     * @return {string} - LaTeX code.
     */
    function getDynamicFormula(viewmodel, matrix) {
        var string = SYMBOLS.EMPTY;

        var gapStart = viewmodel.gapStart();
        var gapExtension = viewmodel.enlargement();

        // look if we maximize or minimize
        if (viewmodel.calculation() === ALIGNMENT_TYPES.SIMILARITY)
            string += SYMBOLS.EQUAL + LATEX.MAX;
        else
            string += SYMBOLS.EQUAL + LATEX.MIN;

        string += LATEX.BEGIN_CASES;  // starting LaTeX case region

            // differentiate between formulas and writing code for dynamic one
            // Hint: -x and x with x is number should be right aligned -> LATEX.SPACE is added on positive numbers
            switch (matrix) {
                case MATRICES.VERTICAL:
                    string += LATEX.FORMULA.TOP + LATEX.ALIGNED_PLUS;
                    string += gapStart >= 0 ? LATEX.SPACE + gapStart : gapStart;
                    string += LATEX.NEW_LINE;

                    string += LATEX.FORMULA.TOP_P + LATEX.ALIGNED_PLUS;
                    string += gapExtension >= 0 ? LATEX.SPACE + gapExtension : gapExtension;
                    break;
                case MATRICES.HORIZONTAL:
                    string += LATEX.FORMULA.LEFT + LATEX.ALIGNED_PLUS;
                    string += gapStart >= 0 ? LATEX.SPACE + gapStart : gapStart;
                    string += LATEX.NEW_LINE;

                    string += LATEX.FORMULA.LEFT_Q + LATEX.ALIGNED_PLUS;
                    string += gapExtension >= 0 ? LATEX.SPACE + gapExtension : gapExtension;
                    break;
                default:
                    string += LATEX.FORMULA.DIAGONAL + LATEX.ALIGNED_PLUS;
                    string += viewmodel.match() >= 0 ? LATEX.SPACE + viewmodel.match() : viewmodel.match();
                    string += SYMBOLS.AND + LATEX.FORMULA.MATCH + LATEX.NEW_LINE;

                    string += LATEX.FORMULA.DIAGONAL + LATEX.ALIGNED_PLUS;
                    string += viewmodel.mismatch() >= 0 ? LATEX.SPACE + viewmodel.mismatch() : viewmodel.mismatch();
                    string += SYMBOLS.AND + LATEX.FORMULA.MISMATCH + LATEX.NEW_LINE;

                    string += LATEX.FORMULA.CURRENT_P + LATEX.NEW_LINE;

                    string += LATEX.FORMULA.CURRENT_Q;
            }

        string += LATEX.END_CASES;  // stopping LaTeX case region

        return string;
    }

    /**
     * Returns the LaTeX-code for sub-formulas like gap-functions of affine algorithms.
     * @param viewmodel {InputViewmodel} - The viewmodel of the view displaying the formula.
     * @return {string} - LaTeX code.
     */
    function getSubformula(viewmodel) {
        var string = LATEX.MATH_REGION;

            string += LATEX.RECURSION.GOTOH_GAP_FUNCTION;
            string += SYMBOLS.EQUAL + LATEX.MULTIPLIER;
            string += viewmodel.enlargement() >= 0
                ? viewmodel.enlargement()
                : SYMBOLS.BRACKET_LEFT + viewmodel.enlargement() + SYMBOLS.BRACKET_RIGHT;

            string += viewmodel.baseCosts() >= 0
                ? SYMBOLS.PLUS + viewmodel.baseCosts()
                : SYMBOLS.PLUS + SYMBOLS.BRACKET_LEFT + viewmodel.baseCosts() + SYMBOLS.BRACKET_RIGHT;

        string += LATEX.MATH_REGION;
        return string;
    }

    /**
     * Processing the input from the user.
     * This function is executed by the Input-Processor
     * and it is dependant on the algorithm.
     * @param algorithm {Object} - Algorithm used to update the user interface.
     * @param inputProcessor {Object} - The unit processing the input.
     * @param inputViewmodel {Object} - The InputViewmodel used to access inputs.
     * @param visualViewmodel {Object} - The VisualViewmodel used to access visualization functions.
     */
    function processInput(algorithm, inputProcessor, inputViewmodel, visualViewmodel) {
        visualViewmodel.removeAllContents();

        // when page was loaded the inputs have not to be updated or you get wrong inputs
        if (inputProcessor.inputUpdatesActivated()) {
            inputViewmodel.sequence1($("#sequence_1").val().toUpperCase());
            inputViewmodel.sequence2($("#sequence_2").val().toUpperCase());

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

    /**
     * Changes the output after processing the input.
     * @param outputData {Object} - Contains all output data.
     * @param inputProcessor {Object} - The unit processing the input.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions.
     */
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

    /**
     * Post edits a matrix and replaces for example values with LaTeX-symbols.
     * @param inputProcessor {Object} - The unit processing the input.
     * @param matrix {matrix} - The matrix in which you want replace values with for example LaTeX-symbols.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions.
     * @return {matrix} - The matrix in which symbols where replaced with LaTeX-symbols.
     * @augments AlignmentInterface.edit(inputProcessor, matrix, viewmodels)
     */
    function edit(inputProcessor, matrix, viewmodels) {
        return inputProcessor.postEdit(matrix, viewmodels.visual);
    }
}());
