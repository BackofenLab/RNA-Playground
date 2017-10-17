/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("interfaces.alignmentInterface", AlignmentInterface, imports, sharedInterfaceOperations, startProcessing, roundValues);

    // instances
    var alignmentInterfaceInstance;

    /**
     * Is used to work with the input and output (the interface) of an alignment algorithm.
     * @constructor
     */
    function AlignmentInterface() {
        alignmentInterfaceInstance = this;

        // public class methods
        this.imports = imports;
        this.sharedInterfaceOperations = sharedInterfaceOperations;
        this.startAlignmentAlgorithm = startAlignmentAlgorithm;
        this.startProcessing = startProcessing;
        this.roundValues = roundValues;
    }

    /**
     * Function managing objects.
     */
    function startAlignmentAlgorithm(Algorithm, algorithmName) {
        imports();

        var inputViewmodel = new InputViewmodel(algorithmName);
        sharedInterfaceOperations(Algorithm, inputViewmodel, processInput, changeOutput);
    }

    /**
     * Handling imports.
     */
    function imports() {
        // third party libs
        $.getScript(PATHS.LIBS.KNOCKOUT);  // to make knockout working whenever page is reloaded

        // design/controls logic
        $.getScript(PATHS.INPUT_PROCESSOR);
        $.getScript(PATHS.VISUALIZER);
    }

    /**
     * Interface Operations that are shared between algorithms to initialize and start an algorithm.
     * @param Algorithm {Object} - The alignment algorithm which has to be initialized and started.
     * @param inputViewmodel {Object} - The InputViewmodel used to access inputs.
     * @param processInput {Function} - Function from the algorithm which should process the input.
     * @param changeOutput {Function} - Function from the algorithm which should change the output after processing the input.
     */
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
     * @param algorithmName {string} - The name of the algorithm.
     * @see https://en.wikipedia.org/wiki/Model-view-viewmodel
     * @constructor
     */
    function InputViewmodel(algorithmName) {
        var viewmodel = this;

        if (algorithmName === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER) {
            this.sequence1 = ko.observable(NORMALIZED_ALIGNMENT_DEFAULTS.SEQUENCE_1);
            this.sequence2 = ko.observable(NORMALIZED_ALIGNMENT_DEFAULTS.SEQUENCE_2);

            this.calculation = ko.observable(NORMALIZED_ALIGNMENT_DEFAULTS.CALCULATION);

            // function
            this.deletion = ko.observable(NORMALIZED_ALIGNMENT_DEFAULTS.FUNCTION.DELETION);
            this.insertion = ko.observable(NORMALIZED_ALIGNMENT_DEFAULTS.FUNCTION.INSERTION);
            this.match = ko.observable(NORMALIZED_ALIGNMENT_DEFAULTS.FUNCTION.MATCH);
            this.mismatch = ko.observable(NORMALIZED_ALIGNMENT_DEFAULTS.FUNCTION.MISMATCH);
            this.length = ko.observable(NORMALIZED_ALIGNMENT_DEFAULTS.LENGTH);
        } else {
            this.sequence1 = ko.observable(ALIGNMENT_DEFAULTS.SEQUENCE_1);
            this.sequence2 = ko.observable(ALIGNMENT_DEFAULTS.SEQUENCE_2);

            this.calculation = ko.observable(ALIGNMENT_DEFAULTS.CALCULATION);

            // function
            this.deletion = ko.observable(ALIGNMENT_DEFAULTS.FUNCTION.DELETION);
            this.insertion = ko.observable(ALIGNMENT_DEFAULTS.FUNCTION.INSERTION);
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

                    return getSubFormula(viewmodel);
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
                string += viewmodel.deletion() >= 0 ? LATEX.SPACE + viewmodel.deletion() : viewmodel.deletion();
                string += SYMBOLS.AND + LATEX.FORMULA.DELETION + LATEX.NEW_LINE;

                string += LATEX.FORMULA.LEFT + LATEX.ALIGNED_PLUS;
                string += viewmodel.insertion() >= 0 ? LATEX.SPACE + viewmodel.insertion() : viewmodel.insertion();
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
     * @param viewmodel {InputViewmodel} - The viewmodel of the view displaying the formula.
     * @return {string} - LaTeX code.
     */
    function getSubFormula(viewmodel) {
        var string = LATEX.MATH_REGION;

        string += LATEX.SUB_FORMULAS.ARSLAN_EGECIOGLE_PEVZNER_SCORING;

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

            if (algorithm.type !== ALGORITHMS.SMITH_WATERMAN && algorithm.type !== ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER)
                if ($("#distance").is(":checked"))
                    inputViewmodel.calculation(ALIGNMENT_TYPES.DISTANCE);
                else
                    inputViewmodel.calculation(ALIGNMENT_TYPES.SIMILARITY);

            inputViewmodel.deletion(Number($("#deletion").val()));
            inputViewmodel.insertion(Number($("#insertion").val()));
            inputViewmodel.match(Number($("#match").val()));
            inputViewmodel.mismatch(Number($("#mismatch").val()));

            if (algorithm.type === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER)
                inputViewmodel.length(Number($("#length").val()));
        } else
            inputProcessor.activateInputUpdates();

        startProcessing(algorithm, inputViewmodel, visualViewmodel);
    }

    /**
     * Start processing the input from the user by computing the algorithm output.
     * @param algorithm {Object} - Algorithm used to update the user interface.
     * @param inputViewmodel {Object} - The InputViewmodel used to access inputs.
     * @param visualViewmodel {Object} - The VisualViewmodel used to access visualization functions.
     */
    function startProcessing(algorithm, inputViewmodel, visualViewmodel) {
        algorithm.setInput(inputViewmodel);
        var ioData = algorithm.compute();
        visualViewmodel.shareInformation(algorithm, ioData[0], ioData[1]);
    }

    /**
     * Binding Viewmodel-functions to InputProcessor elements.
     * This allows for example to highlight a selected entry from a table.
     * @param algorithm {Object} - The algorithm used to update the user interface.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions.
     * @param inputProcessor {Object} - The unit processing the input.
     * @param processInput {Function} - Function from the algorithm which should process the input.
     * @param changeOutput {Function} - Function from the algorithm which should change the output after processing the input.
     */
    function linkInputWithOutput(algorithm, viewmodels, inputProcessor, processInput, changeOutput) {
        inputProcessor.linkElements(viewmodels.visual);
        inputProcessor.updateGUI(algorithm, viewmodels, processInput, changeOutput);
    }

    /**
     * Changes the output after processing the input.
     * @param outputData {Object} - Contains all output data.
     * @param inputProcessor {Object} - The unit processing the input.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions.
     * @see Hint: The parameter inputProcessor is needed!
     */
    function changeOutput(outputData, inputProcessor, viewmodels) {
        roundValues(outputData);

        if (outputData.iterationData !== undefined)
            changeAEPOutput(outputData, viewmodels);
        else if (outputData.matrix !== undefined) {
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
                // 'Uncaught TypeError: inputOutputViewmodel.output.tableValues[i] is not a function'
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
                // 'Uncaught TypeError: inputOutputViewmodel.output.tableValues[i] is not a function'
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
                // 'Uncaught TypeError: inputOutputViewmodel.output.tableValues[i] is not a function'
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
                // 'Uncaught TypeError: inputOutputViewmodel.output.tableValues[i] is not a function'
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
                // 'Uncaught TypeError: inputOutputViewmodel.output.tableValues[i] is not a function'
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

    /**
     * Rounds matrix values, scores and other parameters.
     * @param outputData {Object} - Output data which is modified.
     */
    function roundValues(outputData) {
        if (outputData.iterationData !== undefined) {  // AEP

            // every possibility
            for (var i = 0; i < outputData.iterationData.length; i++) {
                // every round
                for (var j = 0; j < outputData.iterationData[i].length; j++) {

                    // every matrix row
                    for (var l = 0; l < outputData.iterationData[i][j][8].length; l++) {

                        // every entry
                        for (var k = 0; k < outputData.iterationData[i][j][8][l].length; k++) {
                            outputData.iterationData[i][j][8][l][k]
                                = round(outputData.iterationData[i][j][8][l][k], 1);
                        }
                    }

                    outputData.iterationData[i][j][0] = round(outputData.iterationData[i][j][0], 4); // score
                    outputData.iterationData[i][j][2] = round(outputData.iterationData[i][j][2], 4); // lambda
                }
            }

        } else { // other algorithms
            for (var i = 0; i < outputData.matrix.length; i++)
                for (var j = 0; j < outputData.matrix[0].length; j++)
                    outputData.matrix[i][j] = round(outputData.matrix[i][j], 1);

            outputData.score = round(outputData.score, 1);
        }
    }

    /**
     * Rounds a value to one decimal place.
     * @param number {number} - The number which is rounded.
     * @param decimalPlaces {number} - The number of decimal places you want round to.
     * @return {number} - Rounded value.
     */
    function round(number, decimalPlaces) {
        var factor = Math.pow(10, decimalPlaces);
        return Math.round(number*factor)/factor;
    }

    /*---- OUTPUT ----*/
    /**
     * In the Model-View-Viewmodel, the view (HTML-page) is filled with data from
     * outside with the help of the viewmodel (here: OutputViewmodel)
     * by getting data from a model (here: outputData).
     * This OutputViewmodel is shared by the different alignment algorithms.
     * @param algorithmName {string} - The name of the algorithm which is executed.
     * @param outputData {Object} - Contains all output data.
     * @constructor
     * @see https://en.wikipedia.org/wiki/Model-view-viewmodel
     */
    function OutputViewmodel(algorithmName, outputData) {
        var viewmodel = this;
        roundValues(outputData);

        if (outputData.iterationData !== undefined && outputData.iterationData.length > 0) {  // AEP
            createAEPOutputViewmodel(viewmodel, outputData);
        } else if (outputData.matrix != undefined) {  // other algorithms
            this.matrix = ko.observableArray(outputData.matrix);

            for (var i = 0; i < outputData.matrix.length; i++) {
                this.matrix[i] = ko.observableArray(outputData.matrix[i]);
            }

            if (algorithmName === ALGORITHMS.GOTOH) {  // special cases regarding possible algorithms
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
    }

    /**
     * Creates the AEP OutputViewmodel.
     * @param viewmodel - The output viewmodel container which should be filled.
     * @param outputData - The data which is used to fill the viewmodel.
     * @see Not nice without array, but the only found way it's working without any bugs!
     */
    function createAEPOutputViewmodel(viewmodel, outputData) {
        if (outputData.iterationData[0].length > 0) {
            viewmodel.matrix1 = ko.observableArray(outputData.iterationData[0][0][8]);

            for (var i = 0; i < outputData.iterationData[0][0][8].length; i++) {
                viewmodel.matrix1[i] = ko.observableArray(outputData.iterationData[0][0][8][i]);
            }

            viewmodel.alignments1 = ko.observableArray(outputData.iterationData[0][0][7]);

            viewmodel.score1 = ko.observable(outputData.iterationData[0][0][0]);
            viewmodel.length1 = ko.observable(outputData.iterationData[0][0][1]);
            viewmodel.lambda1 = ko.observable(outputData.iterationData[0][0][2]);

            viewmodel.alignmentNumber1 = ko.observable(outputData.iterationData[0][0][10]);
            viewmodel.moreTracebacks1 = ko.observable(outputData.iterationData[0][0][11]);
        } else {
            viewmodel.matrix1 = ko.observableArray([]);
            viewmodel.matrix1[0] = ko.observableArray([]);
            viewmodel.alignments1 = ko.observableArray([]);
            viewmodel.score1 = ko.observable(undefined);
            viewmodel.length1 = ko.observable(undefined);
            viewmodel.lambda1 = ko.observable(undefined);
            viewmodel.alignmentNumber1 = ko.observable(undefined);
            viewmodel.moreTracebacks1 = ko.observable(false);
        }

        if (outputData.iterationData[0].length > 1) {
            viewmodel.matrix2 = ko.observableArray(outputData.iterationData[0][1][8]);

            for (var i = 0; i < outputData.iterationData[0][1][8].length; i++) {
                viewmodel.matrix2[i] = ko.observableArray(outputData.iterationData[0][1][8][i]);
            }

            viewmodel.alignments2 = ko.observableArray(outputData.iterationData[0][1][7]);

            viewmodel.score2 = ko.observable(outputData.iterationData[0][1][0]);
            viewmodel.length2 = ko.observable(outputData.iterationData[0][1][1]);
            viewmodel.lambda2 = ko.observable(outputData.iterationData[0][1][2]);

            viewmodel.alignmentNumber2 = ko.observable(outputData.iterationData[0][1][10]);
            viewmodel.moreTracebacks2 = ko.observable(outputData.iterationData[0][1][11]);
        } else {
            viewmodel.matrix2 = ko.observableArray([]);
            viewmodel.matrix2[0] = ko.observableArray([]);
            viewmodel.alignments2 = ko.observableArray([]);
            viewmodel.score2 = ko.observable(undefined);
            viewmodel.length2 = ko.observable(undefined);
            viewmodel.lambda2 = ko.observable(undefined);
            viewmodel.alignmentNumber2 = ko.observable(undefined);
            viewmodel.moreTracebacks2 = ko.observable(false);
        }

        if (outputData.iterationData[0].length > 2) {
            viewmodel.matrix3 = ko.observableArray(outputData.iterationData[0][2][8]);

            for (var i = 0; i < outputData.iterationData[0][2][8].length; i++) {
                viewmodel.matrix3[i] = ko.observableArray(outputData.iterationData[0][2][8][i]);
            }

            viewmodel.alignments3 = ko.observableArray(outputData.iterationData[0][2][7]);

            viewmodel.score3 = ko.observable(outputData.iterationData[0][2][0]);
            viewmodel.length3 = ko.observable(outputData.iterationData[0][2][1]);
            viewmodel.lambda3 = ko.observable(outputData.iterationData[0][2][2]);

            viewmodel.alignmentNumber3 = ko.observable(outputData.iterationData[0][2][10]);
            viewmodel.moreTracebacks3 = ko.observable(outputData.iterationData[0][2][11]);
        } else {
            viewmodel.matrix3 = ko.observableArray([]);
            viewmodel.matrix3[0] = ko.observableArray([]);
            viewmodel.alignments3 = ko.observableArray([]);
            viewmodel.score3 = ko.observable(undefined);
            viewmodel.length3 = ko.observable(undefined);
            viewmodel.lambda3 = ko.observable(undefined);
            viewmodel.alignmentNumber3 = ko.observable(undefined);
            viewmodel.moreTracebacks3 = ko.observable(false);
        }

        if (outputData.iterationData[0].length > 3) {
            viewmodel.matrix4 = ko.observableArray(outputData.iterationData[0][3][8]);

            for (var i = 0; i < outputData.iterationData[0][3][8].length; i++) {
                viewmodel.matrix4[i] = ko.observableArray(outputData.iterationData[0][3][8][i]);
            }

            viewmodel.alignments4 = ko.observableArray(outputData.iterationData[0][3][7]);

            viewmodel.score4 = ko.observable(outputData.iterationData[0][3][0]);
            viewmodel.length4 = ko.observable(outputData.iterationData[0][3][1]);
            viewmodel.lambda4 = ko.observable(outputData.iterationData[0][3][2]);

            viewmodel.alignmentNumber4 = ko.observable(outputData.iterationData[0][3][10]);
            viewmodel.moreTracebacks4 = ko.observable(outputData.iterationData[0][3][11]);
        } else {
            viewmodel.matrix4 = ko.observableArray([]);
            viewmodel.matrix4[0] = ko.observableArray([]);
            viewmodel.alignments4 = ko.observableArray([]);
            viewmodel.score4 = ko.observable(undefined);
            viewmodel.length4 = ko.observable(undefined);
            viewmodel.lambda4 = ko.observable(undefined);
            viewmodel.alignmentNumber4 = ko.observable(undefined);
            viewmodel.moreTracebacks4 = ko.observable(false);
        }

        if (outputData.iterationData[0].length > 4) {
            viewmodel.matrix5 = ko.observableArray(outputData.iterationData[0][4][8]);

            for (var i = 0; i < outputData.iterationData[0][4][8].length; i++) {
                viewmodel.matrix5[i] = ko.observableArray(outputData.iterationData[0][4][8][i]);
            }

            viewmodel.alignments5 = ko.observableArray(outputData.iterationData[0][4][7]);

            viewmodel.score5 = ko.observable(outputData.iterationData[0][4][0]);
            viewmodel.length5 = ko.observable(outputData.iterationData[0][4][1]);
            viewmodel.lambda5 = ko.observable(outputData.iterationData[0][4][2]);

            viewmodel.alignmentNumber5 = ko.observable(outputData.iterationData[0][4][10]);
            viewmodel.moreTracebacks5 = ko.observable(outputData.iterationData[0][4][11]);
        } else {
            viewmodel.matrix5 = ko.observableArray([]);
            viewmodel.matrix5[0] = ko.observableArray([]);
            viewmodel.alignments5 = ko.observableArray([]);
            viewmodel.score5 = ko.observable(undefined);
            viewmodel.length5 = ko.observable(undefined);
            viewmodel.lambda5 = ko.observable(undefined);
            viewmodel.alignmentNumber5 = ko.observable(undefined);
            viewmodel.moreTracebacks5 = ko.observable(false);
        }

        viewmodel.maxNumberIterations = ko.observable(outputData.maxNumberIterations);
    }

    /**
     * Post edits a matrix and replaces for example values with LaTeX-symbols.
     * @param algorithmName {string} - The name of the algorithm which is executed.
     * @param inputProcessor {Object} - The unit processing the input.
     * @param outputData {Object} - Contains all output data.
     * @param visualViewmodel {Object} - The VisualViewmodel used to access visualization functions.
     * @return outputData {Object} - Changed output data.
     */
    function edit(algorithmName, inputProcessor, outputData, visualViewmodel) {
        if (algorithmName === ALGORITHMS.GOTOH) {
            outputData.horizontalGaps = inputProcessor.postEdit(outputData.horizontalGaps, visualViewmodel);
            outputData.verticalGaps = inputProcessor.postEdit(outputData.verticalGaps, visualViewmodel);
        }

        return outputData;
    }
}());
