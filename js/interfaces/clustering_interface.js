/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("interfaces.clusteringInterface", ClusteringInterface, startClusteringInterface);

    // instances
    var clusteringInterfaceInstance;
    var interfaceInstance;

    /**
     * Is used to work with the input and output (the interface) of a clustering algorithm.
     * It contains the basic methods and the viewmodel for the output and input.
     * This class should be used as a superclass for other clustering interfaces.
     * @constructor
     */
    function ClusteringInterface() {
        clusteringInterfaceInstance = this;

        // inheritance
        interfaceInstance = new interfaces.interface.Interface();

        // public methods
        this.startClusteringInterface = startClusteringInterface;
    }

    /**
     * Function managing objects.
     */
    function startClusteringInterface(Algorithm, algorithmName) {
        var inputViewmodel = new InputViewmodel(algorithmName);
        sharedInterfaceOperations(Algorithm, inputViewmodel, processInput, changeOutput);
    }

    /*---- INPUT ----*/
    /**
     * In the Model-View-Viewmodel, the view (HTML-page) is filled with data from
     * outside with the help of the viewmodel (here: InputViewmodel)
     * by getting data from a model (here: CLUSTERING_DEFAULTS).
     * @param algorithmName {string} - The name of the algorithm.
     * @see https://en.wikipedia.org/wiki/Model-view-viewmodel
     * @constructor
     */
    function InputViewmodel(algorithmName) {
        var viewmodel = this;

        this.availableApproaches = ko.observableArray(HIERARCHICAL_CLUSTERING_DEFAULTS.APPROACHES);
        this.selectedApproach = ko.observableArray(HIERARCHICAL_CLUSTERING_DEFAULTS.STANDARD_APPROACH);

        this.csvTable = ko.observable(HIERARCHICAL_CLUSTERING_DEFAULTS.CSV_TABLE);

        this.errorInput = ko.computed(function () {
            var parser = new formats.csvParser.CsvParser();
            return parser.checkInput(viewmodel.csvTable());
        });

        setTimeout(function () {
            MathJax.Hub.Queue(["Typeset", MathJax.Hub])
        }, REUPDATE_TIMEOUT_MS);
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
        debugger;
        visualViewmodel.removeAllContents();

        // when page was loaded the inputs have not to be updated or you get wrong inputs
        if (inputProcessor.inputUpdatesActivated()) {
            inputViewmodel.selectedApproach([$("#approach_selector option:selected").val()]);
            inputViewmodel.csvTable($("#csv_table").val());
        } else
            inputProcessor.activateInputUpdates();

        clusteringInterfaceInstance.startProcessing(algorithm, inputViewmodel, visualViewmodel);
    }

    /**
     * Changes the output after processing the input.
     * @param outputData {Object} - Contains all output data.
     * @param inputProcessor {Object} - The unit processing the input.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions.
     * @see Hint: The parameter inputProcessor is needed!
     */
    function changeOutput(outputData, inputProcessor, viewmodels) {
        viewmodels.output.newickString(outputData.newickString);
        viewmodels.visual.drawTree();
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
        // tree
        this.newickString = ko.observable(outputData.newickString);

        // distance matrix
        outputData.distanceMatrix
            = getDistanceTable(outputData.distanceMatrix, outputData.distanceMatrixLength, outputData.remainingClusters[0], undefined);

        viewmodel.distanceMatrix =  ko.observableArray(outputData.distanceMatrix);

        for (var i = 0; i < outputData.distanceMatrix.length; i++) {
            viewmodel.distanceMatrix[i] = ko.observableArray(outputData.distanceMatrix[i]);
        }
    }
}());