/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("interfaces.interface", Interface, sharedInterfaceOperations, startProcessing);

    // instances
    var interfaceInstance;

    /**
     * Is used to work with the input and output (the interface) of an alignment algorithm.
     * It contains the basic methods and the viewmodel for the output.
     * This class is used by the various interface scripts as superclass.
     * @constructor
     */
    function Interface() {
        interfaceInstance = this;

        // public class methods
        this.sharedInterfaceOperations = sharedInterfaceOperations;
        this.startProcessing = startProcessing;
    }

    /**
     * Interface Operations that are shared between algorithms to initialize and start an algorithm.
     * @param Algorithm {Object} - The alignment algorithm which has to be initialized and started.
     * @param inputViewmodel {Object} - The InputViewmodel used to access inputs.
     * @param OutputViewmodel {Object} - The OutputViewmodel used to output algorithm computation data.
     * @param processInput {Function} - Function from the algorithm which should process the input.
     * @param changeOutput {Function} - Function from the algorithm which should change the output after processing the input.
     */
    function sharedInterfaceOperations(Algorithm, inputViewmodel, OutputViewmodel, processInput, changeOutput) {
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

        inputProcessor.linkElements(algorithm, viewmodels, processInput, changeOutput);
        executeAlgorithmInterfaceCode(algorithm, viewmodels);

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
     * Executes code for specific algorithm interfaces.
     * @param algorithm {Object} - The algorithm for which interface specific code is executed.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions.
     */
    function executeAlgorithmInterfaceCode(algorithm, viewmodels) {
        debugger;
        if (TREE_ALGORITHMS.indexOf(algorithm.type) !== -1)
            viewmodels.visual.drawTree();
    }

    /**
     * Start processing the input from the user by computing the algorithm output.
     * @param algorithm {Object} - Algorithm used to update the user interface.
     * @param inputViewmodel {Object} - The InputViewmodel used to access inputs.
     * @param visualViewmodel {Object} - The VisualViewmodel used to access visualization functions.
     */
    function startProcessing(algorithm, inputViewmodel, visualViewmodel) {
        debugger;
        algorithm.setInput(inputViewmodel);
        var ioData = algorithm.compute();

        // deep copy of the output before rounding to avoid information loss
        visualViewmodel.shareInformation(algorithm, ioData[0], jQuery.extend(true, {}, ioData[1]));
    }
}());