/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("interfaces.alignmentInterface", AlignmentInterface);

    // instances
    var alignmentInterfaceInstance;
    var interfaceInstance;

    /**
     * Is used to work with the input and output (the interface) of an alignment algorithm.
     * It contains the basic methods and the viewmodel for the output.
     * This class is used by the various interface scripts as superclass.
     * @constructor
     */
    function AlignmentInterface() {
        alignmentInterfaceInstance = this;

        // inheritance (inherits the following methods to its childs)
        interfaceInstance = new interfaces.interface.Interface();

        this.startProcessing = interfaceInstance.startProcessing;
        this.roundValues = interfaceInstance.roundValues;
        this.getDistanceTables = interfaceInstance.getDistanceTables;
        this.getLaTeXFormula = interfaceInstance.getLaTeXFormula;

        // public class methods
        this.imports = imports;
        this.sharedInterfaceOperations = sharedInterfaceOperations;
        this.getLaTeXTraceFunctions = getLaTeXTraceFunctions;
        this.getRowData = getRowData;
        this.getColumnData = getColumnData;
        this.getMinimaData = getMinimaData;
        this.getTwoRowsSubmatricesData = getTwoRowsSubmatricesData;
        this.reorderGroupSequences = reorderGroupSequences;
        this.getLibrariesData = getLibrariesData;
        this.removeNeutralSymbols = removeNeutralSymbols;
        this.sortWithClusterTuples = sortWithClusterTuples;
        this.reorderFinalAlignments = reorderFinalAlignments;
    }

    /**
     * Handling imports.
     */
    function imports() {
        interfaceInstance.imports();
        $.getScript(PATHS.ALIGNMENT_INTERFACE);  // very important, because other interfaces are also using this class
    }

    /**
     * Interface Operations that are shared between algorithms to initialize and start an algorithm.
     * @param Algorithm {Object} - The alignment algorithm which has to be initialized and started.
     * @param inputViewmodel {Object} - The InputViewmodel used to access inputs.
     * @param processInput {Function} - Function from the algorithm which should process the input.
     * @param changeOutput {Function} - Function from the algorithm which should change the output after processing the input.
     * @augments Interface.sharedInterfaceOperations(..)
     */
    function sharedInterfaceOperations(Algorithm, inputViewmodel, processInput, changeOutput) {
        interfaceInstance.sharedInterfaceOperations(Algorithm, inputViewmodel, OutputViewmodel, processInput, changeOutput);
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

        if (GLOBAL_ALGORITHMS.indexOf(algorithmName) >= 0
            || LOCAL_ALGORITHMS.indexOf(algorithmName) >= 0) {  // if basic local or global algorithm
            interfaceInstance.roundValues(algorithmName, outputData);

            if (algorithmName === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER) {  // if AEP
                createAEPOutputViewmodel(viewmodel, outputData);
            } else if (algorithmName === ALGORITHMS.HIRSCHBERG) {
                createHirschbergOutputViewmodel(viewmodel, outputData);
            } else if (outputData.matrix !== undefined) {  // other algorithms
                createMainOutputViewmodel(algorithmName, viewmodel, outputData);
            }
        } else if (MULTI_SEQUENCE_ALGORITHMS.indexOf(algorithmName) >= 0) {  // if multi-sequence alignment algorithm
            if (algorithmName === ALGORITHMS.FENG_DOOLITTLE)
                createFengDoolittleOutputViewmodel(algorithmName, viewmodel, outputData);
            else if (algorithmName === ALGORITHMS.ITERATIVE_REFINMENT)
                createIterativeRefinementOutputViewmodel(algorithmName, viewmodel, outputData);
            else if (algorithmName === ALGORITHMS.NOTREDAME_HIGGINS_HERINGA)
                createTcoffeeOutputViewmodel(algorithmName, viewmodel, outputData);
        }
    }

    /**
     * Creates the AEP OutputViewmodel.
     * @param viewmodel {Object} - The output viewmodel container which should be filled.
     * @param outputData {Object} - The data which is used to fill the viewmodel.
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
     * Creates the OutputViewmodel for some local and global alignment algorithms.
     * @param viewmodel {Object} - The output viewmodel container which should be filled.
     * @param outputData {Object} - The data which is used to fill the viewmodel.
     */
    function createHirschbergOutputViewmodel(viewmodel, outputData) {
        var traceFunctionsData = getLaTeXTraceFunctions(outputData);
        var rowData = getRowData(outputData);
        var columnData = getColumnData(outputData);
        var minimaData = getMinimaData(rowData, columnData);
        var twoRowsData = getTwoRowsSubmatricesData(outputData);

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
        viewmodel.forwardMatrices = ko.observable(outputData.forwardMatrices).extend({deferred: true});
        viewmodel.backwardMatrices = ko.observable(outputData.backwardMatrices).extend({deferred: true});

        // iteration over each matrix
        for (var i = 0; i < outputData.forwardMatrices.length; i++) {
            viewmodel.forwardMatrices[i] = ko.observableArray(outputData.forwardMatrices[i]).extend({deferred: true});

            // iteration over each row of the matrix
            for (var j = 0; j < outputData.forwardMatrices[i].length; j++) {
                viewmodel.forwardMatrices[i][j] = ko.observableArray(outputData.forwardMatrices[i][j]).extend({deferred: true});
            }
        }

        for (var i = 0; i < outputData.backwardMatrices.length; i++) {
            viewmodel.backwardMatrices[i] = ko.observableArray(outputData.backwardMatrices[i]).extend({deferred: true});

            // iteration over each row of the matrix
            for (var j = 0; j < outputData.backwardMatrices[i].length; j++) {
                viewmodel.backwardMatrices[i][j] = ko.observableArray(outputData.backwardMatrices[i][j]).extend({deferred: true});
            }
        }

        viewmodel.alignments = ko.observableArray(outputData.alignments);

        // matrix of all minima
        viewmodel.tracecellLines = ko.observable(outputData.tracecellLines).extend({deferred: true});
        viewmodel.globalMinima = ko.observable(minimaData);

        // header
        viewmodel.recursionNumbersContainer = ko.observable(outputData.recursionNumbersContainer).extend({deferred: true});
        viewmodel.traceFunctions = ko.observable(traceFunctionsData);
        viewmodel.currentGlobalRow = ko.observable(rowData);

        // table header (to avoid a problem between Knockout and MathJax the LaTeX code is generated in viewmodel and not in the view)
        viewmodel.matrixDLatex = ko.observable(interfaceInstance.getLaTeXFormula(LATEX.FORMULA.D));
        viewmodel.matrixDPrimeLatex = ko.observable(interfaceInstance.getLaTeXFormula(LATEX.FORMULA.D_PRIME));
        viewmodel.sumLatex = ko.observable(interfaceInstance.getLaTeXFormula(LATEX.SUM));

        viewmodel.secondSequences = ko.observable(outputData.secondSequences);
        viewmodel.secondSequencePositions = ko.observable(outputData.secondSequencePositions).extend({deferred: true});

        // addition table
        viewmodel.forwardRows = ko.observable(outputData.forwardRows);
        viewmodel.mirroredBackwardRows = ko.observable(outputData.mirroredBackwardRows);
        viewmodel.addedRows = ko.observable(outputData.addedRows);
        viewmodel.highlightPositions = ko.observable(outputData.relativeSplittingPoint);

        // gimmicks/optimizations
        viewmodel.showMatrices = ko.observable(false);

        viewmodel.toggleVisibility = function () {
            viewmodel.showMatrices(!viewmodel.showMatrices());
        };

        viewmodel.toggleLinkText = ko.computed(
            function () {
                return viewmodel.showMatrices() ? TOGGLE_LINK_TEXT.HIDE : TOGGLE_LINK_TEXT.SHOW;
            }
        );

        // generated two rows submatrices (intermediate steps)
        viewmodel.prefixTwoRowsCharacters = ko.observable(forwardTwoRowsCharacters).extend({deferred: true});
        viewmodel.suffixTwoRowsCharacters = ko.observable(backwardTwoRowsCharacters).extend({deferred: true});

        viewmodel.prefixTwoRowsCharactersPositions = ko.observable(forwardTwoRowsCharactersPositions).extend({deferred: true});
        viewmodel.suffixTwoRowsCharactersPositions = ko.observable(backwardTwoRowsCharactersPositions).extend({deferred: true});

        viewmodel.prefixTwoRowsMatrices = ko.observableArray(forwardTwoRowsMatrices).extend({deferred: true});

        // iteration over each matrix (forward matrices)
        for (var i = 0; i < forwardTwoRowsMatrices.length; i++) {
            viewmodel.prefixTwoRowsMatrices[i] = ko.observableArray(forwardTwoRowsMatrices[i]).extend({deferred: true});

            // iteration over each row of the matrix
            for (var j = 0; j < forwardTwoRowsMatrices[i].length; j++) {
                viewmodel.prefixTwoRowsMatrices[i][j] = ko.observableArray(forwardTwoRowsMatrices[i][j]).extend({deferred: true});
            }
        }

        viewmodel.suffixTwoRowsMatrices = ko.observableArray(backwardTwoRowsMatrices).extend({deferred: true});

        // iteration over each matrix (backward matrices)
        for (var i = 0; i < backwardTwoRowsMatrices.length; i++) {
            viewmodel.suffixTwoRowsMatrices[i] = ko.observableArray(backwardTwoRowsMatrices[i]).extend({deferred: true});

            // iteration over each row of the matrix
            for (var j = 0; j < backwardTwoRowsMatrices[i].length; j++) {
                viewmodel.suffixTwoRowsMatrices[i][j] = ko.observableArray(backwardTwoRowsMatrices[i][j]).extend({deferred: true});
            }
        }

        MathJax.Hub.Queue(["Typeset", MathJax.Hub]);  // reinterpret new LaTeX code of the trace functions
    }

    /**
     * Returns the trace function for each round.
     * @param outputData {Object} - The data which is used to fill the viewmodel.
     * @return {Array} - The data which stores the trace function for each recursion round.
     */
    function getLaTeXTraceFunctions(outputData) {
        var traceFunctionData = [];

        for (var k = 0; k < outputData.recursionNumbersContainer.length; k++) {
            var string = LATEX.MATH_REGION;  // starting LaTeX math region
            string += LATEX.FORMULA.DPM + SYMBOLS.BRACKET_LEFT;  // starting function

            var firstSequence = outputData.firstSequences[k];
            var firstSequencePositions = outputData.firstSequencePositions[k];

            string += LATEX.SPACE_SMALL;

            for (var i = 0; i < firstSequence.length; i++) {
                if (i >= MAX_TRACE_FUNCTION_ARG_LEN) {  // cut after a certain number of arguments
                    string += SYMBOLS.SPACE + END_SO_ON;
                    break;
                }

                string += LATEX.TEXT_START;
                string += firstSequence[i];
                string += LATEX.CLOSE;
                if (i === 0 || i + 1 === firstSequence.length) {
                    string += LATEX.SUBORDINATE + LATEX.CURLY_BRACKET_LEFT + firstSequencePositions[i] + LATEX.CURLY_BRACKET_RIGHT;
                }
            }

            string += LATEX.SPACE_SMALL + SYMBOLS.VERTICAL_BAR + LATEX.SPACE_SMALL;

            var secondSequence = outputData.secondSequences[k];
            var secondSequencePositions = outputData.secondSequencePositions[k];

            for (var i = 0; i < secondSequence.length; i++) {
                if (i >= MAX_TRACE_FUNCTION_ARG_LEN) {  // cut after a certain number of arguments
                    string += SYMBOLS.SPACE + END_SO_ON;
                    break;
                }

                string += LATEX.TEXT_START;
                string += secondSequence[i];
                string += LATEX.CLOSE;
                if (i === 0 || i + 1 === secondSequence.length) {
                    string += LATEX.SUBORDINATE + LATEX.CURLY_BRACKET_LEFT + secondSequencePositions[i] + LATEX.CURLY_BRACKET_RIGHT;
                }
            }

            string += LATEX.SPACE_SMALL;
            string += SYMBOLS.BRACKET_RIGHT;  // stopping function
            string += LATEX.MATH_REGION;  // stopping LaTeX math region

            traceFunctionData.push(string);
        }

        return traceFunctionData;
    }

    /**
     * Returns the horizontal minimum position for each round.
     * @param outputData {Object} - The data which is used to fill the viewmodel.
     * @return {Array} - The data which stores the global minimum for each recursion round.
     */
    function getRowData(outputData) {
        var rows = [];

        // iterate over all rounds
        for (var k = 0; k < outputData.recursionNumbersContainer.length; k++) {
            rows.push(outputData.firstSequencePositions[k][outputData.relativeSplittingPoint[k][0] - 1]);
        }

        if (outputData.recursionNumbersContainer.length > 0 && !outputData.lastTracecellIsSource) {  // add last row minimum position
            var lastRound = outputData.relativeSplittingPoint.length - 1;
            rows.push(outputData.firstSequencePositions[0][outputData.relativeSplittingPoint[lastRound][0] - 1]);
        }

        return rows;
    }

    /**
     * Returns the vertical minimum position for each round.
     * @param outputData {Object} - The data which is used to fill the viewmodel.
     * @return {Array} - The data which stores the global minimum for each recursion round.
     */
    function getColumnData(outputData) {
        var columns = [];

        var column;
        // iterate over all rounds
        for (var k = 0; k < outputData.recursionNumbersContainer.length; k++) {
            column = outputData.secondSequencePositions[k][outputData.relativeSplittingPoint[k][1] - 1];

            if (column === undefined)
                column = outputData.secondSequencePositions[k][0] - 1;  // select first defined position and then "-1"

            columns.push(column);
        }

        if (outputData.recursionNumbersContainer.length > 0 && !outputData.lastTracecellIsSource) {  // add last row minimum position
            var lastRound = outputData.relativeSplittingPoint.length - 1;

            column = outputData.secondSequencePositions[0][outputData.relativeSplittingPoint[lastRound][1] - 1];

            if (column === undefined)
                column = outputData.secondSequencePositions[0][0] - 1;  // select first defined position and then "-1"

            columns.push(column);
        }

        return columns;
    }

    /**
     * Creates a string of minimum position tuples with the given data.
     * @param rows {Array} - The rows which contains minimum positions.
     * @param columns {Array} - The columns which contains minimum positions.
     * @return {Array} - String of minimum tuples.
     */
    function getMinimaData(rows, columns) {
        var minimaString = SYMBOLS.EMPTY;

        for (var k = 0; k < rows.length; k++) {  // k is the current round
            if (k !== 0 && k % 3 === 0)  // three entries per line
                minimaString += SYMBOLS.NEW_LINE;

            minimaString += SYMBOLS.BRACKET_LEFT + rows[k] + SYMBOLS.COMMA + columns[k] + SYMBOLS.BRACKET_RIGHT;

            if (k < rows.length - 1)
                minimaString += SYMBOLS.COMMA + SYMBOLS.SPACE;
        }

        return minimaString;
    }

    /**
     * Returns the data to create submatrices with two rows.
     * @return {Object} - The submatrices data.
     */
    function getTwoRowsSubmatricesData(outputData) {
        var forwardSubmatricesOfEachRound = [];
        var forwardSubstringsOfEachRound = [];
        var forwardSubpositionsOfEachRound = [];

        var backwardSubmatricesOfEachRound = [];
        var backwardSubstringsOfEachRound = [];
        var backwardSubpositionsOfEachRound = [];

        var forwardMatrices = outputData.forwardMatrices;
        var backwardMatrices = outputData.backwardMatrices;
        var splittingPoints = outputData.relativeSplittingPoint;
        var leftStrings = outputData.firstSequences;
        var leftStringsPositions = outputData.firstSequencePositions;

        if (forwardMatrices.length > 1) {
            for (var k = 0; k < forwardMatrices.length; k++) {  // or: backwardMatrices.length
                var splittingPosI = splittingPoints[k][0];

                var forwardData = getTwoRowsForwardData(forwardMatrices[k], leftStrings[k], leftStringsPositions[k], splittingPosI);
                var backwardData = getTwoRowsBackwardData(backwardMatrices[k], leftStrings[k], leftStringsPositions[k], splittingPosI);

                var twoRowsForwardMatrices = forwardData[0];
                var twoRowsBackwardMatrices = backwardData[0];

                var twoRowsForwardStrings = forwardData[1];
                var twoRowsBackwardStrings = backwardData[1];

                var twoRowsForwardPositions = forwardData[2];
                var twoRowsBackwardPositions = backwardData[2];

                forwardSubmatricesOfEachRound.push(twoRowsForwardMatrices);
                backwardSubmatricesOfEachRound.push(twoRowsBackwardMatrices);

                forwardSubstringsOfEachRound.push(twoRowsForwardStrings);
                backwardSubstringsOfEachRound.push(twoRowsBackwardStrings);

                forwardSubpositionsOfEachRound.push(twoRowsForwardPositions);
                backwardSubpositionsOfEachRound.push(twoRowsBackwardPositions);
            }
        }

        return [[forwardSubmatricesOfEachRound, backwardSubmatricesOfEachRound],
            [forwardSubstringsOfEachRound, backwardSubstringsOfEachRound],
            [forwardSubpositionsOfEachRound, backwardSubpositionsOfEachRound]];
    }

    /**
     * Returns the forward data in an array of two rows.
     * @param forwardMatrix {Array} - The forward matrix from which the two-rows submatrices are generated.
     * @param leftString {Array} - The string from the column on the left side with the characters
     * @param leftPositions {Array} - The character positions of the column string on the left in an imaginary matrix.
     * @param posI {number} - The row until which the two-rows submatrices are generated.
     * @return {[twoRowsMatrices, twoRowsCharacters, twoRowsLeftPositions]} - Data to visualize two-matrices.
     */
    function getTwoRowsForwardData(forwardMatrix, leftString, leftPositions, posI) {
        var twoRowsMatrices = [];
        var twoRowsCharacters = [];
        var twoRowsLeftPositions = [];

        var upperRow = [];
        var lowerRow = [];
        var upperChar = SYMBOLS.EMPTY;
        var lowerChar = SYMBOLS.EMPTY;
        var upperPos = -1;
        var lowerPos = -1;

        for (var i = 1; i <= posI; i++) {
            upperRow = forwardMatrix[i - 1];
            lowerRow = forwardMatrix[i];

            if (i - 1 === 0) {
                upperChar = SYMBOLS.EMPTY;
                lowerChar = leftString[i - 1];

                upperPos = SYMBOLS.EMPTY;
                lowerPos = leftPositions[i - 1];
            }
            else {
                upperChar = leftString[i - 2];
                lowerChar = leftString[i - 1];

                upperPos = leftPositions[i - 2];
                lowerPos = leftPositions[i - 1];
            }

            twoRowsMatrices.push([upperRow, lowerRow]);
            twoRowsCharacters.push([upperChar, lowerChar]);
            twoRowsLeftPositions.push([upperPos, lowerPos]);
        }

        return [twoRowsMatrices, twoRowsCharacters, twoRowsLeftPositions];
    }

    /**
     * Returns the backward data in an array of two rows.
     * @param backwardMatrix {Array} - The backward matrix from which the two-rows submatrices are generated.
     * @param leftString {Array} - The string from the column on the left side with the characters
     * @param leftPositions {Array} - The character positions of the column string on the left in an imaginary matrix.
     * @param posI {number} - The row until which the two-rows submatrices are generated.
     * @return {[twoRowsMatrices, twoRowsCharacters, twoRowsLeftPositions]} - Data to visualize two-matrices.
     */
    function getTwoRowsBackwardData(backwardMatrix, leftString, leftPositions, posI) {
        var twoRowsMatrices = [];
        var twoRowsCharacters = [];
        var twoRowsLeftPositions = [];

        var upperRow = [];
        var lowerRow = [];
        var upperChar = SYMBOLS.EMPTY;
        var lowerChar = SYMBOLS.EMPTY;
        var upperPos = -1;
        var lowerPos = -1;

        var rotatedBackwardMatrix = getRotatedMatrix(backwardMatrix);

        for (var i = rotatedBackwardMatrix.length - 1; i > posI; i--) {
            upperRow = rotatedBackwardMatrix[i - 1];
            lowerRow = rotatedBackwardMatrix[i];

            upperChar = leftString[i - 2];
            lowerChar = leftString[i - 1];

            upperPos = leftPositions[i - 2];
            lowerPos = leftPositions[i - 1];

            twoRowsMatrices.push([upperRow, lowerRow]);
            twoRowsCharacters.push([upperChar, lowerChar]);
            twoRowsLeftPositions.push([upperPos, lowerPos]);
        }

        return [twoRowsMatrices, twoRowsCharacters, twoRowsLeftPositions];
    }

    /**
     * Rotates the given matrix by 180 degrees.
     * @param backwardMatrix {Array} - An array of the several rows of the matrix.
     * @return {Array} - The rotated matrix.
     */
    function getRotatedMatrix(backwardMatrix) {
        var rotatedMatrix = [];

        for (var i = backwardMatrix.length - 1; i >= 0; i--)
            rotatedMatrix.push(backwardMatrix[i].reverse());

        return rotatedMatrix;
    }

    /**
     * Creates the OutputViewmodel for some local and global alignment algorithms.
     * @param algorithmName {string} - The name of the algorithm which is executed.
     * @param viewmodel {Object} - The output viewmodel container which should be filled.
     * @param outputData {Object} - The data which is used to fill the viewmodel.
     */
    function createMainOutputViewmodel(algorithmName, viewmodel, outputData) {
        viewmodel.matrix = ko.observableArray(outputData.matrix);

        for (var i = 0; i < outputData.matrix.length; i++) {
            viewmodel.matrix[i] = ko.observableArray(outputData.matrix[i]);
        }

        if (algorithmName === ALGORITHMS.GOTOH || algorithmName === ALGORITHMS.GOTOH_LOCAL) {  // special cases regarding possible algorithms
            viewmodel.horizontalGaps = ko.observableArray(outputData.horizontalGaps);
            viewmodel.verticalGaps = ko.observableArray(outputData.verticalGaps);

            for (var i = 0; i < outputData.matrix.length; i++) {
                viewmodel.horizontalGaps[i] = ko.observableArray(outputData.horizontalGaps[i]);
                viewmodel.verticalGaps[i] = ko.observableArray(outputData.verticalGaps[i]);
            }
        }

        viewmodel.alignments = ko.observableArray(outputData.alignments);

        viewmodel.score = ko.observable(outputData.score);
        viewmodel.moreTracebacks = ko.observable(outputData.moreTracebacks);
    }

    /**
     * Creates the OutputViewmodel for Feng-Doolittle.
     * @param algorithmName {string} - The name of the algorithm which is executed.
     * @param viewmodel {Object} - The output viewmodel container which should be filled.
     * @param outputData {Object} - The data which is used to fill the viewmodel.
     */
    function createFengDoolittleOutputViewmodel(algorithmName, viewmodel, outputData) {
        // distance matrices
        outputData.distanceMatrices = interfaceInstance.getDistanceTables(outputData, false, true);

        interfaceInstance.roundValues(algorithmName, outputData);

        viewmodel.distanceMatrices = ko.observableArray(outputData.distanceMatrices).extend({deferred: true});

        // iteration over each matrix
        for (var i = 0; i < outputData.distanceMatrices.length; i++) {
            viewmodel.distanceMatrices[i] = ko.observableArray(outputData.distanceMatrices[i]).extend({deferred: true});

            // iteration over each row of the matrix
            for (var j = 0; j < outputData.distanceMatrices[i].length; j++) {
                viewmodel.distanceMatrices[i][j] = ko.observableArray(outputData.distanceMatrices[i][j]).extend({deferred: true});
            }
        }

        viewmodel.remainingClusters = ko.observable(outputData.remainingClusters).extend({deferred: true});
        viewmodel.minimums = ko.observable(outputData.minimums).extend({deferred: true});

        // merge steps
        reorderGroupSequences(outputData);
        viewmodel.guideAlignments = ko.observable(outputData.guideAlignments);
        viewmodel.guideAlignmentsNames = ko.observable(outputData.guideAlignmentsNames);
        viewmodel.firstGroups = ko.observable(outputData.firstGroups);
        viewmodel.secondGroups = ko.observable(outputData.secondGroups);
        viewmodel.firstGroupsNames = ko.observable(outputData.firstGroupsNames);
        viewmodel.secondGroupsNames = ko.observable(outputData.secondGroupsNames);
        viewmodel.joinedGroups = ko.observable(outputData.joinedGroups);
        viewmodel.joinedGroupNames = ko.observable(outputData.joinedGroupNames);

        // tree and final output
        viewmodel.newickString = ko.observable(outputData.newickString);
        viewmodel.progressiveAlignment = ko.observable(outputData.progressiveAlignment);
        viewmodel.score = ko.observable(outputData.score);

        // pairwise data
        sortWithClusterTuples(outputData.sequencePairNames,
            [outputData.alignmentLengths, outputData.similarities, outputData.gapNumbers, outputData.gapStarts]);
        viewmodel.sequencePairNames = ko.observable(outputData.sequencePairNames);
        viewmodel.alignmentLengths = ko.observable(outputData.alignmentLengths);
        viewmodel.similarities = ko.observable(outputData.similarities);
        viewmodel.gapNumbers = ko.observable(outputData.gapNumbers);
        viewmodel.gapStarts = ko.observable(outputData.gapStarts);

        // gimmicks/optimizations
        viewmodel.showMatrices = ko.observable(false);

        viewmodel.toggleVisibility = function () {
            viewmodel.showMatrices(!viewmodel.showMatrices());
        };

        viewmodel.toggleLinkText = ko.computed(
            function () {
                return viewmodel.showMatrices() ? TOGGLE_LINK_TEXT.HIDE : TOGGLE_LINK_TEXT.SHOW;
            }
        );
    }

    /**
     * Reordering groups in alphabetical order for increasing readability.
     * @param outputData {Object} - The output on which reordering is applied.
     */
    function reorderGroupSequences(outputData) {
        if (outputData.joinedGroupNames.length > 0) {
            var finalGroupName = outputData.joinedGroupNames[outputData.joinedGroupNames.length - 1];
            var finalGroup = outputData.joinedGroups[outputData.joinedGroups.length - 1];

            var groupMemberNames = bases.multiSequenceAlignment.getIndividualSequenceNames(finalGroupName, true);
            var groupMemberRankings = getRankings(groupMemberNames, finalGroup);

            var reorderedGroups = [];
            var reorderedGroupNames = [];

            var reorderedFirstGroups = [];
            var reorderedFirstGroupsNames = [];

            var reorderedSecondGroups = [];
            var reorderedSecondGroupsNames = [];

            // iterate over all groups (result, group 1 and group 2)
            for (var i = 0; i < outputData.joinedGroups.length; i++) {
                var group = outputData.joinedGroups[i];
                var group1 = outputData.firstGroups[i];
                var group2 = outputData.secondGroups[i];

                var memberNames = bases.multiSequenceAlignment.getIndividualSequenceNames(outputData.joinedGroupNames[i], true);
                var member1Names = bases.multiSequenceAlignment.getIndividualSequenceNames(outputData.firstGroupsNames[i], true);
                var member2Names = bases.multiSequenceAlignment.getIndividualSequenceNames(outputData.secondGroupsNames[i], true);

                var sortedGroupAndNames = getSortedGroup(group, memberNames, groupMemberRankings);
                var sorted1GroupAndNames = getSortedGroup(group1, member1Names, groupMemberRankings);
                var sorted2GroupAndNames = getSortedGroup(group2, member2Names, groupMemberRankings);

                reorderedGroups.push(sortedGroupAndNames[0]);
                reorderedGroupNames.push(sortedGroupAndNames[1]);

                reorderedFirstGroups.push(sorted1GroupAndNames[0]);
                reorderedFirstGroupsNames.push(sorted1GroupAndNames[1]);

                reorderedSecondGroups.push(sorted2GroupAndNames[0]);
                reorderedSecondGroupsNames.push(sorted2GroupAndNames[1]);
            }

            outputData.joinedGroups = reorderedGroups;
            outputData.joinedGroupNames = reorderedGroupNames;

            outputData.firstGroups = reorderedFirstGroups;
            outputData.firstGroupsNames = reorderedFirstGroupsNames;

            outputData.secondGroups = reorderedSecondGroups;
            outputData.secondGroupsNames = reorderedSecondGroupsNames;

            outputData.progressiveAlignment = reorderedGroups[reorderedGroups.length - 1];
        }
    }

    /**
     * Returns the rankings of the individual members.
     * The ranking is the position within the cluster names.
     * Hint: memberNames.length <= outputData.clusterNames.length (because duplicate sequences are removed)
     * @param memberNames {Array} - The names of the used sequences (duplicate sequences are removed).
     * @param group {Array} - The group of the members.
     * @return {[rankings, highestRanking]} - The structure containing ranking and highest ranking.
     */
    function getRankings(memberNames, group) {
        var rankings = {};
        var highestRanking = Number.NEGATIVE_INFINITY;

        for (var i = 0; i < memberNames.length; i++) {
            var name = memberNames[i].split(SYMBOLS.COMMA);
            var character = name[0];
            var number = name.length > 1 ? (Number(name[1]) - 1) : 0;

            var characterPosition = CLUSTER_NAMES.indexOf(character);
            var sequence = group[i].replace(MULTI_SYMBOLS.GAP, SYMBOLS.EMPTY).replace(MULTI_SYMBOLS.NONE, SYMBOLS.EMPTY);
            rankings[sequence] = characterPosition + CLUSTER_NAMES.length * number;

            if (highestRanking < rankings[sequence])
                highestRanking = rankings[sequence];
        }

        return [rankings, highestRanking];
    }

    /**
     * Resorts the group and the group names alphabetically in linear time by using two arrays.
     * @param group {Array} - The group which is resorted.
     * @param groupMemberNames {Array} - The group names which are resorted.
     * @param groupMemberRankings {Array} - The rankings which are sued to sort elements.
     * @return {[finalSortedGroup, finalSortedGroupNames]} - The sorted group and names.
     */
    function getSortedGroup(group, groupMemberNames, groupMemberRankings) {
        var highestRanking = groupMemberRankings[1];

        var sortedGroup = new Array(highestRanking);  // with empty positions
        var sortedGroupNames = new Array(highestRanking);

        var finalSortedGroup = [];  // without empty positions
        var finalSortedGroupNames = SYMBOLS.EMPTY;  // without empty positions

        // going over non sorted group
        for (var i = 0; i < group.length; i++) {
            var sequence = group[i].replace(MULTI_SYMBOLS.GAP, SYMBOLS.EMPTY).replace(MULTI_SYMBOLS.NONE, SYMBOLS.EMPTY);
            var sequenceRanking = groupMemberRankings[0][sequence];

            sortedGroup[sequenceRanking] = group[i];
            sortedGroupNames[sequenceRanking] = groupMemberNames !== undefined ? groupMemberNames[i] : SYMBOLS.EMPTY;
        }

        // going over sorted array with empties (to remove the empty positions)
        for (var j = 0; j < sortedGroup.length; j++) {
            if (sortedGroup[j] !== undefined) {
                finalSortedGroup.push(sortedGroup[j]);
                finalSortedGroupNames += sortedGroupNames[j].replace(SYMBOLS.COMMA, SYMBOLS.EMPTY);
            }
        }

        return [finalSortedGroup, finalSortedGroupNames];
    }

    /**
     * Creates the OutputViewmodel for T-coffee.
     * @param algorithmName {string} - The name of the algorithm which is executed.
     * @param viewmodel {Object} - The output viewmodel container which should be filled.
     * @param outputData {Object} - The data which is used to fill the viewmodel.
     */
    function createTcoffeeOutputViewmodel(algorithmName, viewmodel, outputData) {
        outputData.librariesData = getLibrariesData(outputData);

        removeNeutralSymbols(outputData);
        interfaceInstance.roundValues(algorithmName, outputData);
        reorderGroupSequences(outputData);

        // final output
        viewmodel.progressiveAlignment = ko.observable(outputData.progressiveAlignment);
        viewmodel.score = ko.observable(outputData.score);

        // merge steps
        viewmodel.firstGroups = ko.observable(outputData.firstGroups);
        viewmodel.secondGroups = ko.observable(outputData.secondGroups);
        viewmodel.firstGroupsNames = ko.observable(outputData.firstGroupsNames);
        viewmodel.secondGroupsNames = ko.observable(outputData.secondGroupsNames);
        viewmodel.joinedGroups = ko.observable(outputData.joinedGroups);
        viewmodel.joinedGroupNames = ko.observable(outputData.joinedGroupNames);

        // tree
        viewmodel.newickString = ko.observable(outputData.newickString);

        // libraries
        viewmodel.sequencePairsNames = ko.observable(outputData.librariesData[0]);
        viewmodel.libPositionPairs = ko.observable(outputData.librariesData[1]);
        viewmodel.primLibValues = ko.observable(outputData.librariesData[2]);
        viewmodel.extendedLibValues = ko.observable(outputData.librariesData[3]);

        // alignments
        viewmodel.alignmentsGlobal = ko.observable(outputData.librariesData[4]).extend({deferred: true});
        viewmodel.alignmentsLocal = ko.observable(outputData.librariesData[5]).extend({deferred: true});
    }

    /**
     * Returns the data needed to display from primary and extended library (also alignment data to avoid second loop).
     * @param outputData {Object} - The data which is used to fill the viewmodel.
     * @return {[sequencePairsNames, positionPairs, primLibValues, extendedLibValues, alignmentsGlobal, alignmentsLocal]}
     * - The data from primary library, extended library and alignment.
     */
    function getLibrariesData(outputData) {
        var sequencePairsNames = [];
        var alignmentsGlobal = [];
        var alignmentsLocal = [];
        var positionPairs = [];
        var primLibValues = [];
        var extendedLibValues = [];

        var alignmentKeys = Object.keys(outputData.primaryWeightLib);

        // iterate overall alignments
        for (var i = 0; i < alignmentKeys.length; i++) {
            var alignmentKey = alignmentKeys[i];
            var alignmentAndScore = outputData.alignmentsAndScores[alignmentKey];
            var alignmentAndScoreLocal;

            if (outputData.alignmentsAndScoresLocal !== undefined)
                alignmentAndScoreLocal = outputData.alignmentsAndScoresLocal[alignmentKey];

            var positionKeys = Object.keys(outputData.primaryWeightLib[alignmentKey]);

            // split alignmentKey to get an array
            var splittedAlignmentKey = alignmentKey.split(SYMBOLS.COMMA);
            var sequence1Name = outputData.nameOfSequence[splittedAlignmentKey[0]];
            var sequence2Name = outputData.nameOfSequence[splittedAlignmentKey[1]];

            var tempPositionPairs = [];
            var tempPrimLibValues = [];
            var tempExtendedLibValues = [];

            // iterate overall positions in this alignments
            for (var j = 0; j < positionKeys.length; j++) {
                var positionKey = positionKeys[j];

                var valueL = outputData.primaryWeightLib[alignmentKey][positionKey];  // primary library value
                var valueEL = outputData.extendedWeightLib[alignmentKey][positionKey];  // extended library value

                // split positionKey to get an array
                var splittedPositionKey = positionKey.split(SYMBOLS.COMMA);

                if (valueEL !== 0) {
                    tempPositionPairs.push([splittedPositionKey[0], splittedPositionKey[1]]);
                    tempPrimLibValues.push(valueL);
                    tempExtendedLibValues.push(valueEL);
                }
            }

            if (tempPositionPairs.length !== 0) {  // don't show names of sequence pairs for which no L or EL exists
                sortWithNumberTuples(tempPositionPairs, [tempPrimLibValues, tempExtendedLibValues]);

                sequencePairsNames.push([sequence1Name, sequence2Name]);
                positionPairs.push(tempPositionPairs);
                primLibValues.push(tempPrimLibValues);
                extendedLibValues.push(tempExtendedLibValues);
            }

            alignmentsGlobal.push(alignmentAndScore[2]);
            alignmentsLocal.push(alignmentAndScoreLocal !== undefined ? alignmentAndScoreLocal[2] : []);
        }

        sortWithClusterTuples(sequencePairsNames, [positionPairs, primLibValues, extendedLibValues, alignmentsGlobal, alignmentsLocal]);
        return [sequencePairsNames, positionPairs, primLibValues, extendedLibValues, alignmentsGlobal, alignmentsLocal];
    }

    /**
     * Returns numerically sorted input arrays.
     * @param positionPairs {Array} - Array of number tupels which is sorted and used to sort the input arrays.
     * @param inputArrays {Array} - Input arrays.
     */
    function sortWithNumberTuples(positionPairs, inputArrays) {
        var switches = [];

        // documentation {sort} - https://developer.mozilla.org/de/docs/Web/JavaScript/Reference/Global_Objects/Array/sort
        positionPairs.sort(function (a, b) {
            var leftNumberA = Number(a[0]);
            var leftNumberB = Number(b[0]);

            var rightNumberA = Number(a[1]);
            var rightNumberB = Number(b[1]);

            return switchOrNotSwitch(leftNumberA, leftNumberB, rightNumberA, rightNumberB, switches);
        });

        for (var j = 0; j < inputArrays.length; j++) {
            var i = 0;

            inputArrays[j].sort(function (a, b) {  // sort with the sorting determined above
                return switches[i++];
            });
        }
    }

    /**
     * Returns a number which tells you if you have to switch two tuples or not.
     * @param leftNumberA {number} - The left number in a tuple A.
     * @param leftNumberB {number} - The left number in a tuple B.
     * @param rightNumberA {number} - The right number in a tuple A.
     * @param rightNumberB {number}  - The right number in a tuple B.
     * @param switches {Array} - Stores the switches and not-switches.
     */
    function switchOrNotSwitch(leftNumberA, leftNumberB, rightNumberA, rightNumberB, switches) {
        var value = 0;

        if (leftNumberA === leftNumberB) {
            value = rightNumberB > rightNumberA ? -1 : (rightNumberB < rightNumberA ? 1 : 0);
            switches.push(value);
            return value;
        }

        value = leftNumberA - leftNumberB;
        switches.push(value);
        return value;
    }

    /**
     * Returns cluster name sorted input arrays.
     * @param sequencePairsNames {Array} - Array of cluster-tupels which is sorted and used to sort the input arrays.
     * @param inputArrays {Array} - Input arrays.
     * @return {Array} - Sorted Input arrays.
     */
    function sortWithClusterTuples(sequencePairsNames, inputArrays) {
        var switches = [];

        // documentation {sort} - https://developer.mozilla.org/de/docs/Web/JavaScript/Reference/Global_Objects/Array/sort
        sequencePairsNames.sort(function (a, b) {
            var leftNumberA = getNumber(a[0]);
            var leftNumberB = getNumber(b[0]);

            var rightNumberA = getNumber(a[1]);
            var rightNumberB = getNumber(b[1]);

            return switchOrNotSwitch(leftNumberA, leftNumberB, rightNumberA, rightNumberB, switches);
        });

        for (var j = 0; j < inputArrays.length; j++) {
            var i = 0;

            inputArrays[j].sort(function (a, b) {  // sort with the sorting determined above
                return switches[i++];
            });
        }
    }

    /**
     * Translates a cluster-name into a number.
     * @example: (with 26 characters)
     * CLUSTER NAMES:
     * a, b, c, ..., z,         FIRST EPISODE   (0 <= index < 26)
     * a2, b2, c2, ..., z2,     SECOND EPISODE  (26 <= index < 52)
     * a3, b3, ...              THIRD ...       (52 <= index < 78)
     *
     * CALCULATION:
     * c = pos(c) + 0 * CLUSTER_NAMES.length = 2 + 0 = 2
     * a2 = pos(a) + (2-1) * CLUSTER_NAMES.length = 0 + 26 = 26
     * c2 = pos(a) + (2-1) * CLUSTER_NAMES.length = 2 + 26 = 28
     * @param cluster {string} - The name of the cluster which should be translated into a number.
     * @return {number} - The position of the cluster-name within all clusters.
     */
    function getNumber(cluster) {
        // hint: Number("[empty]") is replaced with 0 in JS
        var clusterNumber = Number(cluster.replace(MULTI_SYMBOLS.STRINGS, SYMBOLS.EMPTY));  // the episode
        var clusterPosition = CLUSTER_NAMES.indexOf(cluster.replace(MULTI_SYMBOLS.NUMBERS, SYMBOLS.EMPTY));  // position in alphabet

        clusterNumber = clusterNumber > 0 ? (clusterNumber - 1) : 0;  // see example above

        return clusterPosition + clusterNumber * CLUSTER_NAMES.length;
    }

    /**
     * T-Coffee does not really needs a neutral symbol.
     * It's just an implementation gimmick
     * (to use same functions as in Feng-Doolittle and to avoid code duplicates).
     * Because of this and to avoid confusion
     * with the algorithm it is removed from all created groups (for simplicity).
     * @param outputData {Object} - The data which is used to fill the viewmodel.
     */
    function removeNeutralSymbols(outputData) {
        var numJoinings = outputData.firstGroups.length;

        for (var i = 0; i < numJoinings; i++) {
            for (var j = 0; j < outputData.firstGroups[i].length; j++)
                outputData.firstGroups[i][j] = outputData.firstGroups[i][j].replace(MULTI_SYMBOLS.NONE, SYMBOLS.GAP);

            for (var j = 0; j < outputData.secondGroups[i].length; j++)
                outputData.secondGroups[i][j] = outputData.secondGroups[i][j].replace(MULTI_SYMBOLS.NONE, SYMBOLS.GAP);

            for (var j = 0; j < outputData.joinedGroups[i].length; j++)
                outputData.joinedGroups[i][j] = outputData.joinedGroups[i][j].replace(MULTI_SYMBOLS.NONE, SYMBOLS.GAP);
        }
    }

    /**
     * Creates the OutputViewmodel for Feng-Doolittle.
     * @param algorithmName {string} - The name of the algorithm which is executed.
     * @param viewmodel {Object} - The output viewmodel container which should be filled.
     * @param outputData {Object} - The data which is used to fill the viewmodel.
     */
    function createIterativeRefinementOutputViewmodel(algorithmName, viewmodel, outputData) {
        // final output
        reorderFinalAlignments(outputData);  // do not move down this function
        viewmodel.progressiveAlignment = ko.observable(outputData.progressiveAlignment);
        viewmodel.score = ko.observable(outputData.score);
        viewmodel.refinedProgressiveAlignment = ko.observable(outputData.refinedProgressiveAlignment);
        viewmodel.refinedScore = ko.observable(outputData.refinedScore);

        // realignment steps
        reorderGroupSequences(outputData);  // do not move up this function
        viewmodel.guideAlignments = ko.observable(outputData.guideAlignments);
        viewmodel.guideAlignmentsNames = ko.observable(outputData.guideAlignmentsNames);

        viewmodel.firstGroups = ko.observable(outputData.firstGroups);
        viewmodel.firstGroupsNames = ko.observable(outputData.firstGroupsNames);

        viewmodel.secondGroups = ko.observable(outputData.secondGroups);
        viewmodel.secondGroupsNames = ko.observable(outputData.secondGroupsNames);

        viewmodel.joinedGroups = ko.observable(outputData.joinedGroups);
        viewmodel.realignmentsScores = ko.observable(outputData.realignmentsScores);
        viewmodel.joinedGroupNames = ko.observable(outputData.joinedGroupNames);

        viewmodel.accepted = ko.observable(outputData.accepted);

        // tree
        viewmodel.newickString = ko.observable(outputData.newickString);

        // distance matrix
        viewmodel.remainingClusters = ko.observable(outputData.remainingClusters).extend({deferred: true});

        outputData.distanceMatrix
            = bases.clustering.getMatrixAsTable(outputData.distanceMatrix, outputData.distanceMatrixLength, outputData.remainingClusters[0], undefined, true);

        interfaceInstance.roundValues(algorithmName, outputData);

        viewmodel.distanceMatrix = ko.observableArray(outputData.distanceMatrix);

        for (var i = 0; i < outputData.distanceMatrix.length; i++) {
            viewmodel.distanceMatrix[i] = ko.observableArray(outputData.distanceMatrix[i]);
        }
    }

    /**
     * Reorders the sequences of final alignment in order they have been given as input.
     * @param outputData {Object} - The output on which reordering is applied.
     */
    function reorderFinalAlignments(outputData) {
        if (outputData.joinedGroupNames.length > 0) {
            // sort refined progressive alignment
            var finalGroupName = outputData.refinedProgressiveAlignmentName;
            var finalGroup = outputData.refinedProgressiveAlignment;

            var groupMemberNames = bases.multiSequenceAlignment.getIndividualSequenceNames(finalGroupName, true);
            var groupMemberRankings = getRankings(groupMemberNames, finalGroup);

            outputData.refinedProgressiveAlignment = getSortedGroup(finalGroup, groupMemberNames, groupMemberRankings)[0];

            // sort progressive alignment
            finalGroupName = outputData.progressiveAlignmentName;
            finalGroup = outputData.progressiveAlignment;
            groupMemberNames = bases.multiSequenceAlignment.getIndividualSequenceNames(finalGroupName, true);
            groupMemberRankings = getRankings(groupMemberNames, finalGroup);

            outputData.progressiveAlignment = getSortedGroup(finalGroup, groupMemberNames, groupMemberRankings)[0];
        }
    }
}());
