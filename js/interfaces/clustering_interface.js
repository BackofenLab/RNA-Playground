/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("interfaces.clusteringInterface", ClusteringInterface, startClustering);

    // instances
    var clusteringInterfaceInstance;
    var interfaceInstance;
    var parserInstance;

    /**
     * Is used to work with the input and output (the interface) of an alignment algorithm.
     * It contains the basic methods and the viewmodel for the output.
     * This class is used by the various interface scripts as superclass.
     * @constructor
     */
    function ClusteringInterface() {
        clusteringInterfaceInstance = this;

        // inheritance
        interfaceInstance = new interfaces.interface.Interface();
        parserInstance = new formats.csvParser.CsvParser();

        this.startProcessing = interfaceInstance.startProcessing;

        // public methods
        this.startClustering = startClustering;
    }

    /**
     * Function managing objects.
     */
    function startClustering(Algorithm, algorithmName) {
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
            return checkInput(viewmodel.csvTable());
        });

        this.distanceMatrix = ko.computed(function () {
           return getDistanceMatrix(viewmodel.errorInput() === SYMBOLS.EMPTY, viewmodel.csvTable());
        });


        this.numOfStartClusters = ko.computed(function () {
            return getNumOfStartClusters(viewmodel.errorInput() === SYMBOLS.EMPTY, viewmodel.csvTable());
        });

        this.clusterNames = ko.computed(function () {
            return getColumnNames(viewmodel.distanceMatrix());
        });

        setTimeout(function () {
            MathJax.Hub.Queue(["Typeset", MathJax.Hub])
        }, REUPDATE_TIMEOUT_MS);
    }

    /**
     * Checks with the help of a CSV-parser, 
     * if the input is correct or not and returns the error output.
     * @param csvData {string} - The csv string which has to be converted into a cluster algorithm distance matrix.
     * @return {string} - The error output, if it exists.
     */
    function checkInput(csvData) {
        return parserInstance.checkInput(csvData);
    }

    /**
     * Returns the distance matrix object created by a CSV-parser.
     * @param isCorrectData {boolean} - If false, it returns an empty object and else the distance matrix object.
     * @param csvData {string} - The csv string which has to be converted into a cluster algorithm distance matrix.
     * @return {Object} - The distance matrix.
     */
    function getDistanceMatrix(isCorrectData, csvData) {
        var distanceMatrix = {};

        if (isCorrectData)
            distanceMatrix = parserInstance.getMatrix(csvData, getClusterNames);

        return distanceMatrix;
    }

    /**
     * Returns names for clusters associated with the data.
     * Hint: After all characters are depleted,
     * a number is concatenated to the character
     * to make this function generic.
     * @param number {number} - The number of names you want create.
     * @example:
     * CLUSTER NAMES:
     * a, b, c, ..., z,         FIRST EPISODE
     * a2, b2, c2, ..., z2,     SECOND EPISODE
     * a3, b3, ...              THIRD ...
     * @return {Array} - The cluster names.
     */
    function getClusterNames(number) {
        var clusterNames = [];
        var currentEpisode = 1;

        // for every pairwise distance we need a symbol
        for (var i = 0; i < number; i++) {
            if (i < CLUSTER_NAMES.length)
                clusterNames.push(CLUSTER_NAMES[i]);  // add a, b, c, ..., z

            if (i >= CLUSTER_NAMES.length && i % CLUSTER_NAMES.length === 0)  // out of characters
                currentEpisode++;  // new episode

            // out of characters -> a2, b2, c2, ..., z2, a3, b3, ...
            if (i >= CLUSTER_NAMES.length)
                clusterNames.push(CLUSTER_NAMES[i % CLUSTER_NAMES.length] + SYMBOLS.EMPTY + currentEpisode);
        }

        return clusterNames;
    }

    /**
     * Returns the number of start clusters given the CSV data.
     * @param isCorrectData {boolean} - If false, it returns an empty object and else the distance matrix object.
     * @param csvData {string} - The csv string which has to be converted into a cluster algorithm distance matrix.
     * @return {Object} - The distance matrix.
     */
    function getNumOfStartClusters(isCorrectData, csvData) {
        var numOfStartClusters = 0;

        if (isCorrectData)
            numOfStartClusters = parserInstance.getNumOfTableColumns(csvData);

        return numOfStartClusters;
    }

    /**
     * Returns the names from a distance matrix.
     * @param distanceMatrix {Object} - The distance matrix object from which the column/row names are returned.
     */
    function getColumnNames(distanceMatrix) {
        var namePairs = Object.keys(distanceMatrix);

        var names = [];

        for (var i = 0; i < namePairs.length; i++) {
            var name1 = namePairs[i][0];
            var name2 = namePairs[i][1];

            if (names.indexOf(name1) === -1)
                names.push(name1);

            if (!names.indexOf(name2) === -1)
                names.push(name2);
        }

        return names;
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
    }
}());