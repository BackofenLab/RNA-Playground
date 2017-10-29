/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("interfaces.alignmentInterface", AlignmentInterface,
        imports, sharedInterfaceOperations, roundValues, startProcessing);

    // instances
    var alignmentInterfaceInstance;

    /**
     * Is used to work with the input and output (the interface) of an alignment algorithm.
     * It contains the basic methods and the viewmodel for the output.
     * This class is used by the various interface scripts as superclass.
     * @constructor
     */
    function AlignmentInterface() {
        alignmentInterfaceInstance = this;

        // public class methods
        this.imports = imports;
        this.sharedInterfaceOperations = sharedInterfaceOperations;
        this.startProcessing = startProcessing;
        this.roundValues = roundValues;
    }

    /**
     * Handling imports.
     */
    function imports() {
        // third party libs
        $.getScript(PATHS.LIBS.KNOCKOUT);  // to make knockout working whenever page is reloaded

        // design/controls logic
        /*
        This two imports are very important!
        Without an import the classes are not reinitialized correctly for the next algorithm!
         */
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

    /**
     * Post edits a matrix and replaces for example values with LaTeX-symbols.
     * @param algorithmName {string} - The name of the algorithm which is executed.
     * @param inputProcessor {Object} - The unit processing the input.
     * @param outputData {Object} - Contains all output data.
     * @param visualViewmodel {Object} - The VisualViewmodel used to access visualization functions.
     * @return outputData {Object} - Changed output data.
     */
    function edit(algorithmName, inputProcessor, outputData, visualViewmodel) {
        if (algorithmName === ALGORITHMS.GOTOH || algorithmName === ALGORITHMS.GOTOH_LOCAL) {
            outputData.horizontalGaps = inputProcessor.postEdit(outputData.horizontalGaps, visualViewmodel);
            outputData.verticalGaps = inputProcessor.postEdit(outputData.verticalGaps, visualViewmodel);
        }

        return outputData;
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
        if (algorithmName !== ALGORITHMS.FENG_DOOLITTLE) {
            roundValues(outputData);

            if (algorithmName === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER) {
                createAEPOutputViewmodel(viewmodel, outputData);
            } else if (outputData.matrix !== undefined) {  // other algorithms
                this.matrix = ko.observableArray(outputData.matrix);

                for (var i = 0; i < outputData.matrix.length; i++) {
                    this.matrix[i] = ko.observableArray(outputData.matrix[i]);
                }

                if (algorithmName === ALGORITHMS.GOTOH || algorithmName === ALGORITHMS.GOTOH_LOCAL) {  // special cases regarding possible algorithms
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
     * Rounds a value with a given precision.
     * @param number {number} - The number which is rounded.
     * @param decimalPlaces {number} - The number of decimal places you want round to.
     * @return {number} - Rounded value.
     */
    function round(number, decimalPlaces) {
        var factor = Math.pow(10, decimalPlaces);
        return Math.round(number*factor)/factor;
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
     * Start processing the input from the user by computing the algorithm output.
     * @param algorithm {Object} - Algorithm used to update the user interface.
     * @param inputViewmodel {Object} - The InputViewmodel used to access inputs.
     * @param visualViewmodel {Object} - The VisualViewmodel used to access visualization functions.
     */
    function startProcessing(algorithm, inputViewmodel, visualViewmodel) {
        algorithm.setInput(inputViewmodel);
        var ioData = algorithm.compute();

        // deep copy of the output before rounding to avoid information loss
        visualViewmodel.shareInformation(algorithm, ioData[0], jQuery.extend(true, {}, ioData[1]));
    }
}());
