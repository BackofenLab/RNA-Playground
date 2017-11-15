/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("interfaces.linearAlignmentInterface", LinearAlignmentInterface, startLinearAlignmentAlgorithm);

    // instances
    var alignmentInterfaceInstance;
    var linearAlignmentInterfaceInstance;

    /**
     * Is used to work with the input and output (the interface) of a linear alignment algorithm.
     * @constructor
     * @augments AlignmentInterface
     */
    function LinearAlignmentInterface() {
        linearAlignmentInterfaceInstance = this;

        // inheritance
        alignmentInterfaceInstance = new interfaces.alignmentInterface.AlignmentInterface();

        // public class methods
        this.startLinearAlignmentAlgorithm = startLinearAlignmentAlgorithm;
    }

    /**
     * Function managing objects.
     */
    function startLinearAlignmentAlgorithm(Algorithm, algorithmName) {
        imports();

        var inputViewmodel = new InputViewmodel(algorithmName);
        alignmentInterfaceInstance.sharedInterfaceOperations(Algorithm, inputViewmodel, processInput, changeOutput);
    }

    /**
     * Handling imports.
     */
    function imports() {
        alignmentInterfaceInstance.imports();

        $.getScript(PATHS.ALIGNMENT_INTERFACE);  // very important, because other interfaces are also using this class
    }

    /*---- INPUT ----*/
    /**
     * In the Model-View-Viewmodel, the view (HTML-page) is filled with data from
     * outside with the help of the viewmodel (here: InputViewmodel)
     * by getting data from a model (here: ALIGNMENT_DEFAULTS).
     * @param algorithmName {string} - The name of the algorithm.
     * @see https://en.wikipedia.org/wiki/Model-view-viewmodel
     * @constructor
     */
    function InputViewmodel(algorithmName) {
        var viewmodel = this;

        if (algorithmName === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER) {
            this.sequence1 = ko.observable(NORMALIZED_ALIGNMENT_DEFAULTS.SEQUENCE_1);
            this.sequence2 = ko.observable(NORMALIZED_ALIGNMENT_DEFAULTS.SEQUENCE_2);

            this.calculation = ko.observable(NORMALIZED_ALIGNMENT_DEFAULTS.CALCULATION);  // needed to output correct formulas

            // function
            this.gap = ko.observable(NORMALIZED_ALIGNMENT_DEFAULTS.FUNCTION.GAP);
            this.match = ko.observable(NORMALIZED_ALIGNMENT_DEFAULTS.FUNCTION.MATCH);
            this.mismatch = ko.observable(NORMALIZED_ALIGNMENT_DEFAULTS.FUNCTION.MISMATCH);
            this.length = ko.observable(NORMALIZED_ALIGNMENT_DEFAULTS.LENGTH);
        } else {
            this.sequence1 = ko.observable(ALIGNMENT_DEFAULTS.SEQUENCE_1);
            this.sequence2 = ko.observable(ALIGNMENT_DEFAULTS.SEQUENCE_2);

            this.calculation = ko.observable(ALIGNMENT_DEFAULTS.CALCULATION);

            // function
            this.gap = ko.observable(ALIGNMENT_DEFAULTS.FUNCTION.GAP);
            this.match = ko.observable(ALIGNMENT_DEFAULTS.FUNCTION.MATCH);
            this.mismatch = ko.observable(ALIGNMENT_DEFAULTS.FUNCTION.MISMATCH);
        }

        this.formula = ko.computed(
            function getSelectedFormula() {
                // to fire LaTeX-Code reinterpretation after the selected formula was changed
                // HINT: only found solution which works on all browsers
                setTimeout(function () {
                    MathJax.Hub.Queue(["Typeset", MathJax.Hub])
                }, REUPDATE_TIMEOUT_MS);

                return getFormula(algorithmName, viewmodel);
            }
        );

        if (algorithmName === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER) {
            this.subFormula = ko.computed(
                function getSelectedSubFormula() {
                    // to fire LaTeX-Code reinterpretation after the selected formula was changed
                    // HINT: only found solution which works on all browsers
                    setTimeout(function () {
                        MathJax.Hub.Queue(["Typeset", MathJax.Hub])
                    }, REUPDATE_TIMEOUT_MS);

                    return getSubFormula();
                }
            );
        }
    }

    /**
     * Returns the LaTeX-code for formulas.
     * @param algorithmName {string} - The name of the algorithm.
     * @param viewmodel {InputViewmodel} - The viewmodel of the view displaying the formula.
     * @return {string} - LaTeX code.
     */
    function getFormula(algorithmName, viewmodel) {
        var string = LATEX.MATH_REGION;  // starting LaTeX math region

        if (viewmodel.calculation() === ALIGNMENT_TYPES.SIMILARITY)
            string += LATEX.FORMULA.CURRENT + SYMBOLS.EQUAL + LATEX.MAX;
        else
            string += LATEX.FORMULA.CURRENT + SYMBOLS.EQUAL + LATEX.MIN;

        if (algorithmName === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER)
            string += LATEX.RECURSION.SMITH_WATERMAN_MODIFIED;
        else if (algorithmName === ALGORITHMS.NEEDLEMAN_WUNSCH)
            string += LATEX.RECURSION.NEEDLEMAN_WUNSCH;
        else
            string += LATEX.RECURSION.SMITH_WATERMAN;

        if (viewmodel.calculation() === ALIGNMENT_TYPES.SIMILARITY)
            string += SYMBOLS.EQUAL + LATEX.MAX;
        else
            string += SYMBOLS.EQUAL + LATEX.MIN;

        string += LATEX.BEGIN_CASES;  // starting LaTeX case region
        // Hint: -x and x with x is number should be right aligned -> LATEX.SPACE is added on positive numbers
        string += LATEX.FORMULA.DIAGONAL + LATEX.ALIGNED_PLUS;
        string += viewmodel.match() >= 0 ? LATEX.SPACE + viewmodel.match() : viewmodel.match();
        string += SYMBOLS.AND + LATEX.FORMULA.MATCH + LATEX.NEW_LINE;

        string += LATEX.FORMULA.DIAGONAL + LATEX.ALIGNED_PLUS;
        string += viewmodel.mismatch() >= 0 ? LATEX.SPACE + viewmodel.mismatch() : viewmodel.mismatch();
        string += SYMBOLS.AND + LATEX.FORMULA.MISMATCH + LATEX.NEW_LINE;

        string += LATEX.FORMULA.TOP + LATEX.ALIGNED_PLUS;
        string += viewmodel.gap() >= 0 ? LATEX.SPACE + viewmodel.gap() : viewmodel.gap();
        string += SYMBOLS.AND + LATEX.FORMULA.DELETION + LATEX.NEW_LINE;

        string += LATEX.FORMULA.LEFT + LATEX.ALIGNED_PLUS;
        string += viewmodel.gap() >= 0 ? LATEX.SPACE + viewmodel.gap() : viewmodel.gap();
        string += SYMBOLS.AND + LATEX.FORMULA.INSERTION;

        if (algorithmName === ALGORITHMS.SMITH_WATERMAN || algorithmName === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER)
            string += LATEX.NEW_LINE + LATEX.FORMULA.ZERO;

        string += LATEX.END_CASES;  // stopping LaTeX case region

        if (algorithmName === ALGORITHMS.SMITH_WATERMAN || algorithmName === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER)
            string = string.replace(MULTI_SYMBOLS.D_BIG, SYMBOLS.S_BIG);

        string += LATEX.MATH_REGION;  // stopping LaTeX math region
        return string;
    }

    /**
     * Returns the LaTeX-code for sub-formulas like gap-functions of subadditive algorithms.
     * @return {string} - LaTeX code.
     */
    function getSubFormula() {
        var string = LATEX.MATH_REGION;

        string += LATEX.SUB_FORMULAS.ARSLAN_EGECIOGLE_PEVZNER_SCORING;

        string += LATEX.MATH_REGION;
        return string;
    }

    /**
     * Processing the input from the user.
     * This function is executed by the Input-Processor
     * and it is dependant on the algorithm.
     * It is needed by the algorithm
     * to read in the current and not the last values of the inputs.
     * The problem is that processInput can be only executed before the observable changes
     * its values in the viewmodel and setting timeouts is not possible,
     * because it's very hardware dependant when the view is fully reloaded.
     * This is why all values are reseted with JQuery.
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

            if (algorithm.type !== ALGORITHMS.SMITH_WATERMAN && algorithm.type !== ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER)
                if ($("#distance").is(":checked"))
                    inputViewmodel.calculation(ALIGNMENT_TYPES.DISTANCE);
                else
                    inputViewmodel.calculation(ALIGNMENT_TYPES.SIMILARITY);

            inputViewmodel.gap(Number($("#gap").val()));
            inputViewmodel.match(Number($("#match").val()));
            inputViewmodel.mismatch(Number($("#mismatch").val()));

            if (algorithm.type === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER)
                inputViewmodel.length(Number($("#length").val()));
        } else
            inputProcessor.activateInputUpdates();

        alignmentInterfaceInstance.startProcessing(algorithm, inputViewmodel, visualViewmodel);
    }

    /**
     * Changes the output after processing the input.
     * @param outputData {Object} - Contains all output data.
     * @param inputProcessor {Object} - The unit processing the input.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions.
     * @see Hint: The parameter inputProcessor is needed!
     */
    function changeOutput(outputData, inputProcessor, viewmodels) {
        alignmentInterfaceInstance.roundValues(viewmodels.visual.algorithm.type, outputData);

        if (outputData.iterationData !== undefined)  // if AEP
            changeAEPOutput(outputData, viewmodels);
        else if (outputData.matrix !== undefined) {
            viewmodels.output.matrix(outputData.matrix);

            for (var i = 0; i < outputData.matrix.length; i++) {
                // new variables (rows) are not automatically functions
                // and so we have to convert new variables manually into functions
                // or we get the following error
                // 'Uncaught TypeError: viewmodels.output.matrix[i] is not a function'
                if (i > viewmodels.output.matrix.length)
                    viewmodels.output.matrix[i] = new Function();

                viewmodels.output.matrix[i](outputData.matrix[i]);
            }

            viewmodels.output.alignments(outputData.alignments);
            viewmodels.output.score(outputData.score);
            viewmodels.output.moreTracebacks(outputData.moreTracebacks);
        }
    }

    /**
     * Changes the output of Arslan-Egecioglu-Pevzner algorithm after processing the input.
     * @param outputData {Object} - Contains all output data.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions.
     * @see Not nice without array, but the only found way it's working without any bugs!
     */
    function changeAEPOutput(outputData, viewmodels) {
        if (outputData.iterationData.length > 0 && outputData.iterationData[0].length > 0) {
            viewmodels.output.matrix1(outputData.iterationData[0][0][8]);

            for (var i = 0; i < outputData.iterationData[0][0][8].length; i++) {
                // new variables (rows) are not automatically functions
                // and so we have to convert new variables manually into functions
                // or we get the following error
                // 'Uncaught TypeError: viewmodels.output.matrix[i] is not a function'
                if (i > viewmodels.output.matrix1.length)
                    viewmodels.output.matrix1[i] = new Function();

                viewmodels.output.matrix1[i](outputData.iterationData[0][0][8][i]);
            }

            viewmodels.output.alignments1(outputData.iterationData[0][0][7]);

            viewmodels.output.score1(outputData.iterationData[0][0][0]);
            viewmodels.output.length1(outputData.iterationData[0][0][1]);
            viewmodels.output.lambda1(outputData.iterationData[0][0][2]);

            viewmodels.output.alignmentNumber1(outputData.iterationData[0][0][10]);
            viewmodels.output.moreTracebacks1(outputData.iterationData[0][0][11]);
        } else if (viewmodels.output.matrix1 !== undefined) {
            viewmodels.output.matrix1([]);
            viewmodels.output.matrix1[0]([]);
            viewmodels.output.alignments1([]);
            viewmodels.output.score1(0);
            viewmodels.output.length1(0);
            viewmodels.output.lambda1(0);
            viewmodels.output.alignmentNumber1(undefined);
            viewmodels.output.moreTracebacks1(false);
        }

        if (outputData.iterationData.length > 0 && outputData.iterationData[0].length > 1) {
            viewmodels.output.matrix2(outputData.iterationData[0][1][8]);

            for (var i = 0; i < outputData.iterationData[0][1][8].length; i++) {
                // new variables (rows) are not automatically functions
                // and so we have to convert new variables manually into functions
                // or we get the following error
                // 'Uncaught TypeError: viewmodels.output.matrix[i] is not a function'
                if (i > viewmodels.output.matrix2.length)
                    viewmodels.output.matrix2[i] = new Function();

                viewmodels.output.matrix2[i](outputData.iterationData[0][1][8][i]);
            }

            viewmodels.output.alignments2(outputData.iterationData[0][1][7]);

            viewmodels.output.score2(outputData.iterationData[0][1][0]);
            viewmodels.output.length2(outputData.iterationData[0][1][1]);
            viewmodels.output.lambda2(outputData.iterationData[0][1][2]);

            viewmodels.output.alignmentNumber2(outputData.iterationData[0][1][10]);
            viewmodels.output.moreTracebacks2(outputData.iterationData[0][1][11]);
        } else if (viewmodels.output.matrix2 !== undefined) {
            viewmodels.output.matrix2([]);
            viewmodels.output.matrix2[0]([]);
            viewmodels.output.alignments2([]);
            viewmodels.output.score2(0);
            viewmodels.output.length2(0);
            viewmodels.output.lambda2(0);
            viewmodels.output.alignmentNumber2(undefined);
            viewmodels.output.moreTracebacks2(false);
        }

        if (outputData.iterationData.length > 0 && outputData.iterationData[0].length > 2) {
            viewmodels.output.matrix3(outputData.iterationData[0][2][8]);

            for (var i = 0; i < outputData.iterationData[0][2][8].length; i++) {
                // new variables (rows) are not automatically functions
                // and so we have to convert new variables manually into functions
                // or we get the following error
                // 'Uncaught TypeError: viewmodels.output.matrix[i] is not a function'
                if (i > viewmodels.output.matrix3.length)
                    viewmodels.output.matrix3[i] = new Function();

                viewmodels.output.matrix3[i](outputData.iterationData[0][2][8][i]);
            }

            viewmodels.output.alignments3(outputData.iterationData[0][2][7]);

            viewmodels.output.score3(outputData.iterationData[0][2][0]);
            viewmodels.output.length3(outputData.iterationData[0][2][1]);
            viewmodels.output.lambda3(outputData.iterationData[0][2][2]);

            viewmodels.output.alignmentNumber3(outputData.iterationData[0][2][10]);
            viewmodels.output.moreTracebacks3(outputData.iterationData[0][2][11]);
        } else if (viewmodels.output.matrix3 !== undefined) {
            viewmodels.output.matrix3([]);
            viewmodels.output.matrix3[0]([]);
            viewmodels.output.alignments3([]);
            viewmodels.output.score3(0);
            viewmodels.output.length3(0);
            viewmodels.output.lambda3(0);
            viewmodels.output.alignmentNumber3(undefined);
            viewmodels.output.moreTracebacks3(false);
        }

        if (outputData.iterationData.length > 0 && outputData.iterationData[0].length > 3) {
            viewmodels.output.matrix4(outputData.iterationData[0][3][8]);

            for (var i = 0; i < outputData.iterationData[0][3][8].length; i++) {
                // new variables (rows) are not automatically functions
                // and so we have to convert new variables manually into functions
                // or we get the following error
                // 'Uncaught TypeError: viewmodels.output.matrix[i] is not a function'
                if (i > viewmodels.output.matrix4.length)
                    viewmodels.output.matrix4[i] = new Function();

                viewmodels.output.matrix4[i](outputData.iterationData[0][3][8][i]);
            }

            viewmodels.output.alignments4(outputData.iterationData[0][3][7]);

            viewmodels.output.score4(outputData.iterationData[0][3][0]);
            viewmodels.output.length4(outputData.iterationData[0][3][1]);
            viewmodels.output.lambda4(outputData.iterationData[0][3][2]);

            viewmodels.output.alignmentNumber4(outputData.iterationData[0][3][10]);
            viewmodels.output.moreTracebacks4(outputData.iterationData[0][3][11]);
        } else if (viewmodels.output.matrix4 !== undefined) {
            viewmodels.output.matrix4([]);
            viewmodels.output.matrix4[0]([]);
            viewmodels.output.alignments4([]);
            viewmodels.output.score4(0);
            viewmodels.output.length4(0);
            viewmodels.output.lambda4(0);
            viewmodels.output.alignmentNumber4(undefined);
            viewmodels.output.moreTracebacks4(false);
        }

        if (outputData.iterationData.length > 0 && outputData.iterationData[0].length > 4) {
            viewmodels.output.matrix5(outputData.iterationData[0][4][8]);

            for (var i = 0; i < outputData.iterationData[0][4][8].length; i++) {
                // new variables (rows) are not automatically functions
                // and so we have to convert new variables manually into functions
                // or we get the following error
                // 'Uncaught TypeError: viewmodels.output.matrix[i] is not a function'
                if (i > viewmodels.output.matrix5.length)
                    viewmodels.output.matrix5[i] = new Function();

                viewmodels.output.matrix5[i](outputData.iterationData[0][4][8][i]);
            }

            viewmodels.output.alignments5(outputData.iterationData[0][4][7]);

            viewmodels.output.score5(outputData.iterationData[0][4][0]);
            viewmodels.output.length5(outputData.iterationData[0][4][1]);
            viewmodels.output.lambda5(outputData.iterationData[0][4][2]);

            viewmodels.output.alignmentNumber5(outputData.iterationData[0][4][10]);
            viewmodels.output.moreTracebacks5(outputData.iterationData[0][4][11]);
        } else if (viewmodels.output.matrix5 !== undefined) {
            viewmodels.output.matrix5([]);
            viewmodels.output.matrix5[0]([]);
            viewmodels.output.alignments5([]);
            viewmodels.output.score5(0);
            viewmodels.output.length5(0);
            viewmodels.output.lambda5(0);
            viewmodels.output.alignmentNumber5(undefined);
            viewmodels.output.moreTracebacks5(false);
        }

        viewmodels.output.maxNumberIterations(outputData.maxNumberIterations);
    }
}());
