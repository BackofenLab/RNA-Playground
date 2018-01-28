/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("interfaces.linearAlignmentInterface", LinearAlignmentInterface);

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
     * @param Algorithm {Object} - The algorithm which is started.
     * @param algorithmName {string} - The name of the algorithm which is started.
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
        //$.getScript(PATHS.LINEAR_ALIGNMENT_INTERFACE);  // very important, because other interfaces are also using this class
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
        var isHirschBerg = algorithmName === ALGORITHMS.HIRSCHBERG;

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

            this.calculation = ko.observable(isHirschBerg ? ALIGNMENT_DEFAULTS.CALCULATION_HIRSCHBERG : ALIGNMENT_DEFAULTS.CALCULATION);

            // function
            this.gap = ko.observable(isHirschBerg ? -ALIGNMENT_DEFAULTS.FUNCTION.GAP : ALIGNMENT_DEFAULTS.FUNCTION.GAP);
            this.match = ko.observable(isHirschBerg ? -ALIGNMENT_DEFAULTS.FUNCTION.MATCH : ALIGNMENT_DEFAULTS.FUNCTION.MATCH);
            this.mismatch = ko.observable(isHirschBerg ? -ALIGNMENT_DEFAULTS.FUNCTION.MISMATCH : ALIGNMENT_DEFAULTS.FUNCTION.MISMATCH);
        }

        this.formula = ko.computed(
            function getSelectedFormula() {
                // to fire LaTeX-Code reinterpretation after the selected formula was changed
                // HINT: only found solution which works on all browsers
                setTimeout(function () {
                    MathJax.Hub.Queue(["Typeset", MathJax.Hub])
                }, REUPDATE_TIMEOUT_MS);

                return getFormula(algorithmName, viewmodel, false);
            }
        );

        if (algorithmName === ALGORITHMS.HIRSCHBERG) {
            this.formulaBackward = ko.computed(
                function getSelectedFormula() {
                    // to fire LaTeX-Code reinterpretation after the selected formula was changed
                    // HINT: only found solution which works on all browsers
                    setTimeout(function () {
                        MathJax.Hub.Queue(["Typeset", MathJax.Hub])
                    }, REUPDATE_TIMEOUT_MS);

                    return getFormula(algorithmName, viewmodel, true);
                }
            );
        }

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
     * @param secondRecursion {boolean} - Tells if the function is called a second time. If it is so, then second formula selected.
     * @return {string} - LaTeX code.
     */
    function getFormula(algorithmName, viewmodel, secondRecursion) {
        var string = LATEX.MATH_REGION;  // starting LaTeX math region

        if (viewmodel.calculation() === ALIGNMENT_TYPES.SIMILARITY)
            string += LATEX.FORMULA.CURRENT + SYMBOLS.EQUAL + LATEX.MAX;
        else if (secondRecursion)
            string += LATEX.FORMULA.CURRENT_BACKWARD + SYMBOLS.EQUAL + LATEX.MIN;
        else
            string += LATEX.FORMULA.CURRENT + SYMBOLS.EQUAL + LATEX.MIN;

        if (algorithmName === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER)
            string += LATEX.RECURSION.SMITH_WATERMAN_MODIFIED;
        else if (algorithmName === ALGORITHMS.HIRSCHBERG)
            if (secondRecursion)
                string += LATEX.RECURSION.HIRSCHBERG_BACKWARD;
            else
                string += LATEX.RECURSION.HIRSCHBERG_FORWARD;
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
            string = string.replace(LATEX.FORMULA.D_BIG, LATEX.FORMULA.S_BIG);

        if (algorithmName === ALGORITHMS.HIRSCHBERG && secondRecursion) {
            // replace D with D'
            string = string.replace(LATEX.FORMULA.D_BIG_UNDERSCORE, LATEX.FORMULA.D_PRIME_UNDERSCORE);

            // replace i-1 with i+1
            string = string
                .replace(LATEX.FORMULA.I_MINUS_ONE, LATEX.FORMULA.I_PLUS_ONE)
                .replace(LATEX.FORMULA.J_MINUS_ONE, LATEX.FORMULA.J_PLUS_ONE);

            // replace a_i with a_{i+1} and b_j ...
            string = string
                .replace(LATEX.FORMULA.SEQ_A_I, LATEX.FORMULA.SEQ_A_I_PLUS_1)
                .replace(LATEX.FORMULA.SEQ_B_J, LATEX.FORMULA.SEQ_B_J_PLUS_1);

            // add initialization information
            string += LATEX.NEW_LINE + LATEX.RECURSION.HIRSCHBERG_INITIALIZATION;
        }

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
        var algorithmType = viewmodels.visual.algorithm.type;

        alignmentInterfaceInstance.roundValues(algorithmType, outputData);

        if (algorithmType === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER)  // if AEP
            changeAEPOutput(outputData, viewmodels);
        else if (algorithmType === ALGORITHMS.HIRSCHBERG) {
            changeHirschbergOutput(outputData, viewmodels);
        } else if (outputData.matrix !== undefined) {  // all other algorithms
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

    /**
     * Changes the output of Arslan-Egecioglu-Pevzner algorithm after processing the input.
     * @param outputData {Object} - Contains all output data.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions.
     */
    function changeHirschbergOutput(outputData, viewmodels) {
        var traceFunctionsData = alignmentInterfaceInstance.getLaTeXTraceFunctions(outputData);
        var rowData = alignmentInterfaceInstance.getRowData(outputData);
        var columnData = alignmentInterfaceInstance.getColumnData(outputData);
        var minimaData = alignmentInterfaceInstance.getMinimaData(rowData, columnData);
        var twoRowsData = alignmentInterfaceInstance.getTwoRowsSubmatricesData(outputData);

        // get matrices data
        var twoRowsMatrices = twoRowsData[0];
        var twoRowsCharacters = twoRowsData[1];
        var twoRowsCharactersPositions = twoRowsData[2];

        // divide data in forward and backward
        var forwardTwoRowsMatrices = twoRowsMatrices[0];
        var backwardTwoRowsMatrices = twoRowsMatrices[1];

        var forwardTwoRowsCharacters = twoRowsCharacters[0];
        var backwardTwoRowsCharacters = twoRowsCharacters[1];

        var forwardTwoRowsCharactersPositions = twoRowsCharactersPositions[0];
        var backwardTwoRowsCharactersPositions = twoRowsCharactersPositions[1];

        // main output
        viewmodels.output.forwardMatrices(outputData.forwardMatrices);
        viewmodels.output.backwardMatrices(outputData.backwardMatrices);

        // iteration over each matrix in forward matrices
        for (var i = 0; i < outputData.forwardMatrices.length; i++) {
            // new variables (rows) are not automatically functions...
            if (i >= viewmodels.output.forwardMatrices.length)
                viewmodels.output.forwardMatrices[i] = new Function();

            viewmodels.output.forwardMatrices[i](outputData.forwardMatrices[i]);

            // iteration over each row of the matrix
            for (var j = 0; j < outputData.forwardMatrices[i].length; j++) {
                // new variables (rows) are not automatically functions...
                if (j >= viewmodels.output.forwardMatrices[i].length)
                    viewmodels.output.forwardMatrices[i][j] = new Function();

                viewmodels.output.forwardMatrices[i][j](outputData.forwardMatrices[i][j]);
            }
        }

        // iteration over each matrix in backward matrices
        for (var i = 0; i < outputData.backwardMatrices.length; i++) {
            // new variables (rows) are not automatically functions...
            if (i >= viewmodels.output.backwardMatrices.length)
                viewmodels.output.backwardMatrices[i] = new Function();

            viewmodels.output.backwardMatrices[i](outputData.backwardMatrices[i]);

            // iteration over each row of the matrix
            for (var j = 0; j < outputData.backwardMatrices[i].length; j++) {
                // new variables (rows) are not automatically functions...
                if (j >= viewmodels.output.backwardMatrices[i].length)
                    viewmodels.output.backwardMatrices[i][j] = new Function();

                viewmodels.output.backwardMatrices[i][j](outputData.backwardMatrices[i][j]);
            }
        }

        viewmodels.output.alignments(outputData.alignments);

        // matrix of all minima
        viewmodels.output.tracecellLines(outputData.tracecellLines);
        viewmodels.output.globalMinima(minimaData);

        // header
        viewmodels.output.recursionNumbersContainer(outputData.recursionNumbersContainer);
        viewmodels.output.traceFunctions(traceFunctionsData);
        viewmodels.output.currentGlobalRow(rowData);

        // table header (to avoid a problem between Knockout and MathJax the LaTeX code is generated in viewmodel and not in the view)
        //viewmodels.output.matrixDLatex(alignmentInterfaceInstance.getLaTeXFormula(LATEX.FORMULA.D));
        //viewmodels.output.matrixDPrimeLatex(alignmentInterfaceInstance.getLaTeXFormula(LATEX.FORMULA.D_PRIME));
        //viewmodels.output.sumLatex(alignmentInterfaceInstance.getLaTeXFormula(LATEX.SUM));

        viewmodels.output.secondSequences(outputData.secondSequences);
        viewmodels.output.secondSequencePositions(outputData.secondSequencePositions);

        // addition table
        viewmodels.output.forwardRows(outputData.forwardRows);
        viewmodels.output.mirroredBackwardRows(outputData.mirroredBackwardRows);
        viewmodels.output.addedRows(outputData.addedRows);
        viewmodels.output.highlightPositions(outputData.relativeSplittingPoint);

        // generated two rows submatrices (intermediate steps)
        viewmodels.output.prefixTwoRowsCharacters(forwardTwoRowsCharacters);
        viewmodels.output.suffixTwoRowsCharacters(backwardTwoRowsCharacters);

        viewmodels.output.prefixTwoRowsCharactersPositions(forwardTwoRowsCharactersPositions);
        viewmodels.output.suffixTwoRowsCharactersPositions(backwardTwoRowsCharactersPositions);

        viewmodels.output.prefixTwoRowsMatrices(forwardTwoRowsMatrices);

        // iteration over each matrix (forward matrices)
        for (var i = 0; i < forwardTwoRowsMatrices.length; i++) {
            // new variables (rows) are not automatically functions...
            if (i >= viewmodels.output.prefixTwoRowsMatrices.length)
                viewmodels.output.prefixTwoRowsMatrices[i] = new Function();

            viewmodels.output.prefixTwoRowsMatrices[i](forwardTwoRowsMatrices[i]);

            // iteration over each row of the matrix
            for (var j = 0; j < forwardTwoRowsMatrices[i].length; j++) {
                // new variables (rows) are not automatically functions...
                if (j >= viewmodels.output.prefixTwoRowsMatrices[i].length)
                    viewmodels.output.prefixTwoRowsMatrices[i][j] = new Function();

                viewmodels.output.prefixTwoRowsMatrices[i][j](forwardTwoRowsMatrices[i][j]);
            }
        }

        viewmodels.output.suffixTwoRowsMatrices(backwardTwoRowsMatrices);

        // iteration over each matrix (backward matrices)
        for (var i = 0; i < backwardTwoRowsMatrices.length; i++) {
            // new variables (rows) are not automatically functions...
            if (i >= viewmodels.output.suffixTwoRowsMatrices.length)
                viewmodels.output.suffixTwoRowsMatrices[i] = new Function();

            viewmodels.output.suffixTwoRowsMatrices[i](backwardTwoRowsMatrices[i]);

            // iteration over each row of the matrix
            for (var j = 0; j < backwardTwoRowsMatrices[i].length; j++) {
                // new variables (rows) are not automatically functions...
                if (j >= viewmodels.output.suffixTwoRowsMatrices[i].length)
                    viewmodels.output.suffixTwoRowsMatrices[i][j] = new Function();

                viewmodels.output.suffixTwoRowsMatrices[i][j](backwardTwoRowsMatrices[i][j]);
            }
        }

        MathJax.Hub.Queue(["Typeset", MathJax.Hub]);  // reinterpret new LaTeX code of the trace functions
    }
}());
