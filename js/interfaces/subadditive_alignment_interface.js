/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("interfaces.subadditiveAlignmentInterface", SubadditiveAlignmentInterface);

    // instances
    var alignmentInterfaceInstance;
    var subadditiveAlignmentInterfaceInstance;

    /**
     * Is used to work with the input and output (the interface) of a subadditive alignment algorithm.
     * @constructor
     * @augments AlignmentInterface
     */
    function SubadditiveAlignmentInterface() {
        subadditiveAlignmentInterfaceInstance = this;

        // inheritance
        alignmentInterfaceInstance = new interfaces.alignmentInterface.AlignmentInterface();

        // public class methods
        this.startSubadditiveAlignmentAlgorithm = startSubadditiveAlignmentAlgorithm;
    }

    /**
     * Function managing objects.
     * @param Algorithm {Object} - The algorithm which is started.
     * @param algorithmName {string} - The name of the algorithm which is started.
     */
    function startSubadditiveAlignmentAlgorithm(Algorithm, algorithmName) {
        imports();

        var inputViewmodel = new InputViewmodel(algorithmName);
        alignmentInterfaceInstance.sharedInterfaceOperations(Algorithm, inputViewmodel, processInput, changeOutput);
    }

    /**
     * Handling imports.
     */
    function imports() {
        alignmentInterfaceInstance.imports();
        //$.getScript(PATHS.SUBADDITIVE_ALIGNMENT_INTERFACE);  // very important, because other interfaces are also using this class
    }

    /*---- INPUT ----*/
    /**
     * In the Model-View-Viewmodel, the view (HTML-page) is filled with data from
     * outside with the help of the viewmodel (here: InputViewmodel)
     * by getting data from a model (here: SUBADDITIVE_ALIGNMENT_DEFAULTS).
     * @param algorithmName {string} - The name of the algorithm.
     * @see https://en.wikipedia.org/wiki/Model-view-viewmodel
     * @constructor
     */
    function InputViewmodel(algorithmName) {
        var viewmodel = this;

        this.sequence1 = ko.observable(SUBADDITIVE_ALIGNMENT_DEFAULTS.SEQUENCE_1);
        this.sequence2 = ko.observable(SUBADDITIVE_ALIGNMENT_DEFAULTS.SEQUENCE_2);

        this.calculation = ko.observable(SUBADDITIVE_ALIGNMENT_DEFAULTS.CALCULATION);

        // function
        this.baseCosts = ko.observable(SUBADDITIVE_ALIGNMENT_DEFAULTS.FUNCTION.BASE_COSTS);
        this.enlargement = ko.observable(SUBADDITIVE_ALIGNMENT_DEFAULTS.FUNCTION.ENLARGEMENT);
        this.match = ko.observable(SUBADDITIVE_ALIGNMENT_DEFAULTS.FUNCTION.MATCH);
        this.mismatch = ko.observable(SUBADDITIVE_ALIGNMENT_DEFAULTS.FUNCTION.MISMATCH);

        // displayed dynamic formulas
        this.gapStart = ko.computed(
            function () {
                return Number(viewmodel.baseCosts()) + Number(viewmodel.enlargement());
            }
        );

        if (algorithmName === ALGORITHMS.WATERMAN_SMITH_BEYER) {
            this.subadditiveFunction = ko.observable(SUBADDITIVE_ALIGNMENT_DEFAULTS.GAP_FUNCTION);
        } else {  // if Gotoh or Gotoh (Local)
            this.formulaP = ko.computed(
                function getSelectedFormula() {
                    setTimeout(function () {
                        MathJax.Hub.Queue(["Typeset", MathJax.Hub])
                    }, REUPDATE_TIMEOUT_MS);

                    return getFormula(algorithmName, viewmodel, MATRICES.VERTICAL);
                }
            );

            this.formulaQ = ko.computed(
                function getSelectedFormula() {
                    setTimeout(function () {
                        MathJax.Hub.Queue(["Typeset", MathJax.Hub])
                    }, REUPDATE_TIMEOUT_MS);

                    return getFormula(algorithmName, viewmodel, MATRICES.HORIZONTAL);
                }
            );
        }

        this.formula = ko.computed(
            function getSelectedFormula() {
                // to fire LaTeX-Code reinterpretation after the selected formula was changed
                // HINT: only found solution which works on all browsers
                setTimeout(function () {
                    MathJax.Hub.Queue(["Typeset", MathJax.Hub])
                }, REUPDATE_TIMEOUT_MS);

                return getFormula(algorithmName, viewmodel, MATRICES.DEFAULT);
            }
        );

        this.gapFunction = ko.computed(
            function getSelectedFormula() {
                setTimeout(function () {
                    MathJax.Hub.Queue(["Typeset", MathJax.Hub])
                }, REUPDATE_TIMEOUT_MS);

                return getSubformula(algorithmName, viewmodel);
            }
        );
    }

    /**
     * Returns the LaTeX-code for formulas of subadditive algorithms.
     * @param algorithmName {string} - The name of the algorithm.
     * @param viewmodel {InputViewmodel} - The viewmodel of the view displaying the formula.
     * @param matrix {string} - The matrix for which you want display the formula.
     * @return {string} - LaTeX code.
     */
    function getFormula(algorithmName, viewmodel, matrix) {
        if (algorithmName === ALGORITHMS.WATERMAN_SMITH_BEYER)
            return getWatermanSmithBeyerFormula(viewmodel);

        return getGotohFormula(algorithmName, viewmodel, matrix);
    }

    /**
     * Returns the LaTeX-code for Waterman-Smith-Beyer formula.
     * @param viewmodel {InputViewmodel} - The viewmodel of the view displaying the formula.
     * @return {string} - LaTeX code.
     */
    function getWatermanSmithBeyerFormula(viewmodel) {
        var string = LATEX.MATH_REGION;  // starting LaTeX math region

        string += LATEX.FORMULA.CURRENT;

        // look if we maximize or minimize
        if (viewmodel.calculation() === ALIGNMENT_TYPES.SIMILARITY) {
            string += SYMBOLS.EQUAL + LATEX.MAX;
            string += LATEX.RECURSION.WATERMAN_SMITH_BEYER_MAX;
        }
        else {
            string += SYMBOLS.EQUAL + LATEX.MIN;
            string += LATEX.RECURSION.WATERMAN_SMITH_BEYER_MIN;
        }

        string += getWatermanSmithBeyerDynamicFormula(viewmodel);

        string += LATEX.MATH_REGION;  // stopping LaTeX math region
        return string;
    }

    /**
     * Returns the dynamic, input-dependant LaTeX-code for Waterman-Smith-Beyer formula.
     * @param viewmodel {InputViewmodel} - The viewmodel of the view displaying the formula.
     * @return {string} - LaTeX code.
     */
    function getWatermanSmithBeyerDynamicFormula(viewmodel) {
        var string = SYMBOLS.EMPTY;

        // look if we maximize or minimize
        if (viewmodel.calculation() === ALIGNMENT_TYPES.SIMILARITY)
            string += SYMBOLS.EQUAL + LATEX.MAX;
        else
            string += SYMBOLS.EQUAL + LATEX.MIN;

        string += LATEX.BEGIN_CASES;  // starting LaTeX case region

        if (viewmodel.calculation() === ALIGNMENT_TYPES.SIMILARITY)
            string += LATEX.FORMULA.MAXIMIZE_HORIZONTAL + LATEX.NEW_LINE;
        else
            string += LATEX.FORMULA.MINIMIZE_HORIZONTAL + LATEX.NEW_LINE;

        string += LATEX.FORMULA.DIAGONAL + LATEX.ALIGNED_PLUS;
        string += viewmodel.match() >= 0 ? LATEX.SPACE + viewmodel.match() : viewmodel.match();
        string += SYMBOLS.AND + LATEX.FORMULA.MATCH + LATEX.NEW_LINE;

        string += LATEX.FORMULA.DIAGONAL + LATEX.ALIGNED_PLUS;
        string += viewmodel.mismatch() >= 0 ? LATEX.SPACE + viewmodel.mismatch() : viewmodel.mismatch();
        string += SYMBOLS.AND + LATEX.FORMULA.MISMATCH + LATEX.NEW_LINE;

        if (viewmodel.calculation() === ALIGNMENT_TYPES.SIMILARITY)
            string += LATEX.FORMULA.MAXIMIZE_VERTICAL;
        else
            string += LATEX.FORMULA.MINIMIZE_VERTICAL;

        string += LATEX.END_CASES;  // stopping LaTeX case region

        return string;
    }

    /**
     * Returns the LaTeX-code for Gotoh formulas.
     * @param algorithmName {string} - The name of the algorithm.
     * @param viewmodel {InputViewmodel} - The viewmodel of the view displaying the formula.
     * @param matrix {string} - The matrix for which you want display the formula.
     * @return {string} - LaTeX code.
     */
    function getGotohFormula(algorithmName, viewmodel, matrix) {
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
                string += getGotohDynamicFormula(algorithmName, viewmodel, matrix);
                break;
            case MATRICES.HORIZONTAL:
                string += LATEX.RECURSION.GOTOH_Q;
                string += getGotohDynamicFormula(algorithmName, viewmodel, matrix);
                break;
            default:
                if (algorithmName === ALGORITHMS.GOTOH)
                    string += LATEX.RECURSION.GOTOH;
                else if (algorithmName === ALGORITHMS.GOTOH_LOCAL)
                    string += LATEX.RECURSION.GOTOH_LOCAL;

                string += getGotohDynamicFormula(algorithmName, viewmodel, matrix);
        }

        string += LATEX.MATH_REGION;  // stopping LaTeX math region

        if (algorithmName === ALGORITHMS.GOTOH_LOCAL)
            string = string.replace(MULTI_SYMBOLS.D_BIG, SYMBOLS.S_BIG);

        return string;
    }

    /**
     * Returns the dynamic, input-dependant LaTeX-code for Gotoh formulas.
     * @param algorithmName {string} - The name of the algorithm.
     * @param viewmodel {InputViewmodel} - The viewmodel of the view displaying the formula.
     * @param matrix {string} - The matrix for which you want display the formula.
     * @return {string} - LaTeX code.
     */
    function getGotohDynamicFormula(algorithmName, viewmodel, matrix) {
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

                if (algorithmName === ALGORITHMS.GOTOH_LOCAL)
                    string += LATEX.NEW_LINE + LATEX.FORMULA.ZERO;
        }

        string += LATEX.END_CASES;  // stopping LaTeX case region

        return string;
    }

    /**
     * Returns the LaTeX-code for sub-formulas like gap-functions of subadditive algorithms.
     * @param algorithmName {string} - The name of the algorithm.
     * @param viewmodel {InputViewmodel} - The viewmodel of the view displaying the formula.
     * @return {string} - LaTeX code.
     */
    function getSubformula(algorithmName, viewmodel) {
        var string = LATEX.MATH_REGION;

        string += getGapFunction(algorithmName, viewmodel);

        string += SYMBOLS.EQUAL;
        string += viewmodel.baseCosts() >= 0
            ? viewmodel.baseCosts() + SYMBOLS.PLUS
            : SYMBOLS.BRACKET_LEFT + viewmodel.baseCosts() + SYMBOLS.BRACKET_RIGHT + SYMBOLS.PLUS;

        string += viewmodel.enlargement() >= 0
            ? viewmodel.enlargement() + getMultiplier(algorithmName, viewmodel)
            : SYMBOLS.BRACKET_LEFT + viewmodel.enlargement() + SYMBOLS.BRACKET_RIGHT + getMultiplier(algorithmName, viewmodel);

        string += LATEX.MATH_REGION;
        return string;
    }

    /**
     * Returns the right gap function.
     * @param algorithmName {string} - The name of the algorithm.
     * @param viewmodel {InputViewmodel} - The viewmodel of the view displaying the formula.
     * @return {string} - LaTeX code.
     */
    function getGapFunction(algorithmName, viewmodel) {
        var finalFunction = LATEX.SUB_FORMULAS.GOTOH_GAP_FUNCTION;

        if (algorithmName === ALGORITHMS.WATERMAN_SMITH_BEYER) {
            switch (viewmodel.subadditiveFunction()) {
                case SUBADDITIVE_FUNCTIONS.LOGARITHMIC:
                    finalFunction = LATEX.FORMULA.GAP + SYMBOLS.EQUAL;
                    finalFunction += LATEX.ALPHA + SYMBOLS.PLUS + LATEX.BETA +
                        LATEX.DOT + LATEX.LN + SYMBOLS.BRACKET_LEFT + LATEX.FACTOR + SYMBOLS.BRACKET_RIGHT;
                    break;
                case SUBADDITIVE_FUNCTIONS.QUADRATIC:
                    finalFunction = LATEX.FORMULA.GAP + SYMBOLS.EQUAL;
                    finalFunction += LATEX.ALPHA + SYMBOLS.PLUS + LATEX.BETA + LATEX.DOT + LATEX.FACTOR + LATEX.POW2;
                    break;
            }
        }

        return finalFunction;
    }

    /**
     * Returns the right factor for the enlargement formula.
     * @param algorithmName {string} - The name of the algorithm.
     * @param viewmodel {InputViewmodel} - The viewmodel of the view displaying the formula.
     * @return {string} - LaTeX code.
     */
    function getMultiplier(algorithmName, viewmodel) {
        var finalFactor = LATEX.FACTOR;

        if (algorithmName === ALGORITHMS.WATERMAN_SMITH_BEYER) {
            switch (viewmodel.subadditiveFunction()) {
                case SUBADDITIVE_FUNCTIONS.LOGARITHMIC:
                    finalFactor = LATEX.LN + SYMBOLS.BRACKET_LEFT + LATEX.FACTOR + SYMBOLS.BRACKET_RIGHT;
                    break;
                case SUBADDITIVE_FUNCTIONS.QUADRATIC:
                    finalFactor += LATEX.POW2;
                    break;
            }
        }

        return LATEX.DOT + finalFactor;
    }

    /**
     * Processing the input from the user.
     * This function is executed by the Input-Processor
     * and it is dependant on the algorithm.
     * It is needed by the algorithm
     * to read in the current and not the last values of the inputs.
     * The problem is that processInput can be only executed before the observable changes
     * its values in the viewmodel and setting timeouts is not possible,
     * because it's very hardware dependant.
     * Also, some values have to be converted first for example
     * into a number. This is why all values are reseted with JQuery.
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

            if (algorithm.type === ALGORITHMS.WATERMAN_SMITH_BEYER) {
                if ($("#affine").is(":checked"))
                    inputViewmodel.subadditiveFunction(SUBADDITIVE_FUNCTIONS.AFFINE);
                else if ($("#logarithmic").is(":checked"))
                    inputViewmodel.subadditiveFunction(SUBADDITIVE_FUNCTIONS.LOGARITHMIC);
                else // if ($("#quadratic").is(":checked"))
                    inputViewmodel.subadditiveFunction(SUBADDITIVE_FUNCTIONS.QUADRATIC);
            }
        } else
            inputProcessor.activateInputUpdates();

        alignmentInterfaceInstance.startProcessing(algorithm, inputViewmodel, visualViewmodel);
    }

    /**
     * Changes the output after processing the input.
     * @param outputData {Object} - Contains all output data.
     * @param inputProcessor {Object} - The unit processing the input.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions and input.
     * @see Hint: The parameter inputProcessor is needed!
     */
    function changeOutput(outputData, inputProcessor, viewmodels) {
        if (viewmodels.input.subadditiveFunction !== undefined
            && viewmodels.input.subadditiveFunction() === SUBADDITIVE_FUNCTIONS.LOGARITHMIC)
            alignmentInterfaceInstance.roundValues(viewmodels.visual.algorithm.type, outputData);

        viewmodels.output.matrix(outputData.matrix);

        if (viewmodels.output.horizontalGaps !== undefined) {  // if Gotoh or Gotoh (Local)
            viewmodels.output.horizontalGaps(edit(inputProcessor, outputData.horizontalGaps, viewmodels));
            viewmodels.output.verticalGaps(edit(inputProcessor, outputData.verticalGaps, viewmodels));

            for (var i = 0; i < outputData.matrix.length; i++) {
                // new variables (rows) are not automatically functions
                // and so we have to convert new variables manually into functions
                // or we get the following error
                // 'Uncaught TypeError: viewmodels.output.matrix[i] is not a function'
                if (i > viewmodels.output.matrix.length) {
                    viewmodels.output.matrix[i] = new Function();
                    viewmodels.output.verticalGaps[i] = new Function();
                    viewmodels.output.horizontalGaps[i] = new Function();
                }

                viewmodels.output.matrix[i](outputData.matrix[i]);
                viewmodels.output.verticalGaps[i](outputData.verticalGaps[i]);
                viewmodels.output.horizontalGaps[i](outputData.horizontalGaps[i]);
            }
        } else {  // if Waterman-Smith-Beyer
            for (var i = 0; i < outputData.matrix.length; i++) {
                if (i > viewmodels.output.matrix.length)
                    viewmodels.output.matrix[i] = new Function();

                viewmodels.output.matrix[i](outputData.matrix[i]);
            }
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
