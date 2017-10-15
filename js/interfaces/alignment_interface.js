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
        $.getScript(PATHS.LIBS.MATH_JAX);  // in some cases (AEP-algorithm) the script have to be reloaded
        MathJax.Hub.Queue(["Typeset", MathJax.Hub]);  // reinterpret new LaTeX code
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

        this.sequence1 = ko.observable(ALIGNMENT_DEFAULTS.SEQUENCE_1);
        this.sequence2 = ko.observable(ALIGNMENT_DEFAULTS.SEQUENCE_2);

        this.calculation = ko.observable(ALIGNMENT_DEFAULTS.CALCULATION);

        // function
        this.deletion = ko.observable(ALIGNMENT_DEFAULTS.FUNCTION.DELETION);
        this.insertion = ko.observable(ALIGNMENT_DEFAULTS.FUNCTION.INSERTION);
        this.match = ko.observable(ALIGNMENT_DEFAULTS.FUNCTION.MATCH);
        this.mismatch = ko.observable(ALIGNMENT_DEFAULTS.FUNCTION.MISMATCH);

        if (algorithmName === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER)
            this.length = ko.observable(ALIGNMENT_DEFAULTS.LENGTH);

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
        if (viewmodels.output.iterationData !== undefined)
            changeAEPOutput(outputData, viewmodels);
        else {
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
     */
    function changeAEPOutput(outputData, viewmodels) {
        debugger;
        roundValues(outputData);
        viewmodels.output.iterationData(outputData.iterationData);

        for (var i = 0; i < outputData.iterationData.length; i++) {  // iteration over possibilities
            // new possibilities are not automatically functions
            // and so we have to convert new variables manually into functions
            // or we get the following error
            // 'Uncaught TypeError: ... is not a function'
            if (i > outputData.iterationData.length)
                viewmodels.output.iterationData[i] = new Function();

            viewmodels.output.iterationData[i](outputData.iterationData[i]);

            for (var j = 0; j < outputData.iterationData[i].length; j++) {  // iteration over rounds
                // new ... automatically functions
                if (j > outputData.iterationData[i].length)
                    viewmodels.output.iterationData[i][j] = new Function();

                viewmodels.output.iterationData[i][j](outputData.iterationData[i][j]);

                // iteration over [score, length, lambda, deletion, insertion, match, mismatch, alignments, matrix]
                for (var k = 0; k < outputData.iterationData[i][j].length; k++) {
                    viewmodels.output.iterationData[i][j][k](outputData.iterationData[i][j][k]);
                }

                // iteration over matrix rows
                for (var l = 0; l < outputData.iterationData[i][j][8].length; l++) {
                    // new ... not automatically functions
                    if (l > outputData.iterationData[i][j][8].length)
                        viewmodels.output.iterationData[i][j][8][l] = new Function();

                    viewmodels.output.iterationData[i][j][8][l](outputData.iterationData[i][j][8][l]);
                }
            }
        }
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
        roundValues(outputData);

        if (algorithmName === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER){
            this.iterationData = ko.observableArray(outputData.iterationData);

            for (var i = 0; i < outputData.iterationData.length; i++) {  // iteration over possibilities
                this.iterationData[i] = ko.observableArray(outputData.iterationData[i]);

                for (var j = 0; j < outputData.iterationData[i].length; j++) {  // iteration over rounds
                    this.iterationData[i][j] = ko.observableArray(outputData.iterationData[i][j]);

                    // iteration over [score, length, lambda, deletion, insertion, match, mismatch, alignments, matrix, tracebacks, alignmentNumber]
                    for (var k = 0; k < outputData.iterationData[i][j].length; k++) {
                        if (k === 7 || k === 8 || k === 9)
                            this.iterationData[i][j][k] = ko.observableArray(outputData.iterationData[i][j][k]);
                        else
                            this.iterationData[i][j][k] = ko.observable(outputData.iterationData[i][j][k]);
                    }

                    // iteration over matrix rows
                    for (var l = 0; l < outputData.iterationData[i][j][8].length; l++) {
                        this.iterationData[i][j][8][l] = ko.observableArray(outputData.iterationData[i][j][8][l]);
                    }
                }
            }
        } else {  // other algorithms
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
