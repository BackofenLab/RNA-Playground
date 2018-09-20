/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("interfaces.interface", Interface);

    // instances
    var interfaceInstance;

    /**
     * Is used to work with the input and output (the interface) of an algorithm.
     * This class is used by the various interface scripts as superclass.
     * @constructor
     */
    function Interface() {
        interfaceInstance = this;

        // public class methods
        this.imports = imports;
        this.sharedInterfaceOperations = sharedInterfaceOperations;
        this.startProcessing = startProcessing;
        this.roundValues = roundValues;
        this.getDistanceTables = getDistanceTables;
        this.getLaTeXFormula = getLaTeXFormula;
    }

    /**
     * Handling imports.
     */
    function imports() {
        $.getScript(PATHS.INTERFACE);  // very important, because other interfaces are also using this class
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
        algorithm.setInput(inputViewmodel);
        var ioData = algorithm.compute();

        // deep copy of the output before rounding to avoid information loss
        visualViewmodel.shareInformation(algorithm, ioData[0], jQuery.extend(true, {}, ioData[1]));
    }

    /**
     * Rounds matrix values, scores and other parameters.
     * @param algorithmName {string} - The name of the algorithm which is executed.
     * @param outputData {Object} - Output data which is modified.
     */
    function roundValues(algorithmName, outputData) {
        if (algorithmName === ALGORITHMS.AGGLOMERATIVE_CLUSTERING
            || algorithmName === ALGORITHMS.FENG_DOOLITTLE) {  // if Feng-Doolittle or Agglomerative
            // iterate over each distance matrix
            for (var k = 0; k < outputData.distanceMatrices.length; k++) {

                // iterate over each row
                for (var i = 0; i < outputData.distanceMatrices[k].length; i++) {

                    // iterate over each entry
                    for (var j = 0; j < outputData.distanceMatrices[k][i].length; j++) {
                        // if (j > i)  // only the values upper the diagonal
                        outputData.distanceMatrices[k][i][j]
                            = round(outputData.distanceMatrices[k][i][j], 1);
                    }
                }
            }

            if (algorithmName === ALGORITHMS.AGGLOMERATIVE_CLUSTERING) {
                for (var k = 0; k < outputData.neighbourJoiningMatrices.length; k++) {

                    // iterate over each row
                    for (var i = 0; i < outputData.neighbourJoiningMatrices[k].length; i++) {

                        // iterate over each entry
                        for (var j = 0; j < outputData.neighbourJoiningMatrices[k][i].length; j++) {
                            if (j > i)  // only the values upper the diagonal
                                outputData.neighbourJoiningMatrices[k][i][j]
                                    = round(outputData.neighbourJoiningMatrices[k][i][j], 1);
                        }
                    }
                }
            }
        } else if (algorithmName === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER) {
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
        } else if (algorithmName === ALGORITHMS.HIRSCHBERG) {
            // do nothing, because there is nothing to round
        } else if (algorithmName === ALGORITHMS.ITERATIVE_REFINMENT) {
            // iterate over each row
            for (var i = 0; i < outputData.distanceMatrix.length; i++) {
                // iterate over each entry
                for (var j = 0; j < outputData.distanceMatrix[i].length; j++)
                    outputData.distanceMatrix[i][j] = round(outputData.distanceMatrix[i][j], 1);
            }
        } else if (algorithmName === ALGORITHMS.NOTREDAME_HIGGINS_HERINGA) {
            var alignmentPairsCount = outputData.librariesData[0].length;
            var primLibValues = outputData.librariesData[2];
            var extendedLibValues = outputData.librariesData[3];

            for (var i = 0; i < alignmentPairsCount; i++) {
                var positionPairsCount = outputData.librariesData[1][i].length;

                for (var j = 0; j < positionPairsCount; j++) {
                    outputData.librariesData[2][i][j] = round(primLibValues[i][j], 1);
                    outputData.librariesData[3][i][j] = round(extendedLibValues[i][j], 1);
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
        return Math.round(number * factor) / factor;
    }

    /**
     * Converts the distances stored in associative array into a real distance matrix.
     * @param outputData {Object} - The output on which conversion is applied.
     * @param computeJoiningMatrices {boolean} - Tells if instead of distance matrices neighbour-joining matrices should be computed.
     * @param fillBothHalves {boolean} - Tells if both halves of the matrix should be filled with values or not.
     * @return {Object} - The outputData with converted distance matrices.
     */
    function getDistanceTables(outputData, computeJoiningMatrices, fillBothHalves) {
        var matrixLength = outputData.distanceMatrixLength;  // start length

        if (!computeJoiningMatrices) {
            // in each round the matrix gets smaller by one, because two matrices are merged
            for (var i = 0; i < outputData.distanceMatrices.length; i++) {
                outputData.distanceMatrices[i] = bases.clustering.getMatrixAsTable(outputData.distanceMatrices[i],
                    matrixLength - i, outputData.remainingClusters[i], outputData.keys[i], fillBothHalves);
            }

            return outputData.distanceMatrices;
        }

        // else compute Neighbour-Joining matrices
        for (var i = 0; i < outputData.neighbourJoiningMatrices.length; i++) {
            outputData.neighbourJoiningMatrices[i] = bases.clustering.getMatrixAsTable(outputData.neighbourJoiningMatrices[i],
                matrixLength - i, outputData.remainingClusters[i], outputData.keys[i], fillBothHalves);
        }

        return outputData.neighbourJoiningMatrices;
    }

    /**
     * Returns a LaTeX enclosed formula.
     * @param formula {string} - The string which has to be enclosed in LaTeX math mode.
     * @return {string} - The LaTeX math mode enclosed formula.
     */
    function getLaTeXFormula(formula) {
        return LATEX.MATH_REGION + formula + LATEX.MATH_REGION;
    }
}());