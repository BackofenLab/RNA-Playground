/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("interfaces.multiSequenceInterface", MultiSequenceInterface, startMultiSequenceInterface);

    // instances
    var alignmentInterfaceInstance;
    var multiSequenceInterfaceInstance;

    /**
     * Is used to work with the input and output (the interface) of a multi-sequence algorithm.
     * @constructor
     * @augments AlignmentInterface
     */
    function MultiSequenceInterface() {
        multiSequenceInterfaceInstance = this;

        // inheritance
        alignmentInterfaceInstance = new interfaces.alignmentInterface.AlignmentInterface();

        // public class methods
        this.startMultiSequenceInterface = startMultiSequenceInterface;
    }

    /**
     * Function managing objects.
     * @param Algorithm - The algorithm which is started.
     * @param algorithmName - The name of the algorithm which is started.
     */
    function startMultiSequenceInterface(Algorithm, algorithmName) {
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
     * by getting data from a model (here: SUBADDITIVE_ALIGNMENT_DEFAULTS).
     * @param algorithmName {string} - The name of the algorithm.
     * @see https://en.wikipedia.org/wiki/Model-view-viewmodel
     * @constructor
     */
    function InputViewmodel(algorithmName) {
        var viewmodel = this;

        this.sequences = ko.observableArray(MULTI_SEQUENCE_DEFAULTS.SEQUENCES);

        this.calculation = ko.observable(MULTI_SEQUENCE_DEFAULTS.CALCULATION);

        // function
        this.baseCosts = ko.observable(MULTI_SEQUENCE_DEFAULTS.FUNCTION.BASE_COSTS);
        this.enlargement = ko.observable(MULTI_SEQUENCE_DEFAULTS.FUNCTION.ENLARGEMENT);
        this.match = ko.observable(MULTI_SEQUENCE_DEFAULTS.FUNCTION.MATCH);
        this.mismatch = ko.observable(MULTI_SEQUENCE_DEFAULTS.FUNCTION.MISMATCH);

        this.clusterNames = ko.computed(
            function () {
                return getClusterNames(viewmodel.sequences().length);
            }
        );

        this.addRow = function() {
            setTimeout(function () {  // to reinterpret in next statement dynamically created LaTeX-code
                MathJax.Hub.Queue(["Typeset", MathJax.Hub])
            }, REUPDATE_TIMEOUT_MS);

            viewmodel.sequences.push(SYMBOLS.EMPTY);
        };

        this.removeRow = function() {
            setTimeout(function () {  // to reinterpret in next statement dynamically created LaTeX-code
                MathJax.Hub.Queue(["Typeset", MathJax.Hub])
            }, REUPDATE_TIMEOUT_MS);

            viewmodel.sequences.pop();
        };

        // displayed dynamic formulas
        this.gapStart = ko.computed(
            function () {
                return Number(viewmodel.baseCosts()) + Number(viewmodel.enlargement());
            }
        );

        this.gapFunction = ko.computed(
            function getSelectedFormula() {
                setTimeout(function () {  // to reinterpret in next statement dynamically created LaTeX-code
                    MathJax.Hub.Queue(["Typeset", MathJax.Hub])
                }, REUPDATE_TIMEOUT_MS);

                return getSubformula(viewmodel);
            }
        );
    }

    /**
     * Returns LaTeX-formatted names for clusters.
     * Hint: After all characters are depleted,
     * a number is concatenated to the character
     * to make this function generic.
     * @param number {number} - The number of cluster names which should be generated.
     * @example:
     * CLUSTER NAMES:
     * a, b, c, ..., z,         FIRST EPISODE
     * a2, b2, c2, ..., z2,     SECOND EPISODE
     * a3, b3, ...              THIRD ...
     * @return {Array} - The LaTeX-formatted cluster names.
     */
    function getClusterNames(number) {
        var clusterNames = [];
        var currentEpisode = 1;

        // for every pairwise distance we need a symbol
        for (var i = 0; i < number; i++) {
            if (i < CLUSTER_NAMES.length)
                clusterNames.push(LATEX.MATH_REGION + CLUSTER_NAMES[i] + LATEX.MATH_REGION);  // add a, b, c, ..., z

            if (i >= CLUSTER_NAMES.length && i % CLUSTER_NAMES.length === 0)  // out of characters
                currentEpisode++;  // new episode

            if (i >= CLUSTER_NAMES.length)  // out of characters -> a2, b2, c2, ..., z2, a3, b3, ...
                clusterNames.push(
                    LATEX.MATH_REGION
                    + CLUSTER_NAMES[i % CLUSTER_NAMES.length] + LATEX.SUBORDINATE + currentEpisode +
                    LATEX.MATH_REGION);
        }

        return clusterNames;
    }

    /**
     * Returns the LaTeX-code for sub-formulas like gap-functions of subadditive algorithms.
     * @param viewmodel {InputViewmodel} - The viewmodel of the view displaying the formula.
     * @return {string} - LaTeX code.
     */
    function getSubformula(viewmodel) {
        var string = LATEX.MATH_REGION;

        string += LATEX.SUB_FORMULAS.GOTOH_GAP_FUNCTION;

        string += SYMBOLS.EQUAL;
        string += viewmodel.baseCosts() >= 0
            ? viewmodel.baseCosts() + SYMBOLS.PLUS
            : SYMBOLS.BRACKET_LEFT + viewmodel.baseCosts() + SYMBOLS.BRACKET_RIGHT + SYMBOLS.PLUS;

        string += viewmodel.enlargement() >= 0
            ? viewmodel.enlargement() + LATEX.DOT + LATEX.FACTOR
            : SYMBOLS.BRACKET_LEFT + viewmodel.enlargement() + SYMBOLS.BRACKET_RIGHT + LATEX.DOT + LATEX.FACTOR;

        string += LATEX.MATH_REGION;
        return string;
    }

    /**
     * Processing the input from the user.
     * This function is executed by the Input-Processor
     * and it is dependant on the algorithm.
     * It is needed by the algorithm
     * to read in the current and not the last values.
     * @param algorithm {Object} - Algorithm used to update the user interface.
     * @param inputProcessor {Object} - The unit processing the input.
     * @param inputViewmodel {Object} - The InputViewmodel used to access inputs.
     * @param visualViewmodel {Object} - The VisualViewmodel used to access visualization functions.
     */
    function processInput(algorithm, inputProcessor, inputViewmodel, visualViewmodel) {
        visualViewmodel.removeAllContents();

        // when page was loaded the inputs have not to be updated or you get wrong inputs
        if (inputProcessor.inputUpdatesActivated()) {
            inputViewmodel.sequences(getSequencesArray(inputViewmodel));

            inputViewmodel.baseCosts(Number($("#base_costs").val()));
            inputViewmodel.enlargement(Number($("#enlargement").val()));
            inputViewmodel.match(Number($("#match").val()));
            inputViewmodel.mismatch(Number($("#mismatch").val()));
        } else
            inputProcessor.activateInputUpdates();

        alignmentInterfaceInstance.startProcessing(algorithm, inputViewmodel, visualViewmodel);
    }

    /**
     * Returns the sequences of dynamically created inputs.
     * @param inputViewmodel {Object} - The InputViewmodel used to access inputs.
     * @return {Array} - The array of dynamically created inputs.
     */
    function getSequencesArray(inputViewmodel) {
        var sequenceArray = [];

        for (var i = 0; i < inputViewmodel.sequences().length; i++) {
            sequenceArray.push($(".sequence")[i].value.toUpperCase());
        }

        /* bug-fix for a Knockout-problem -> dynamically generated inputs get wrong values after typing in something */
        MULTI_SEQUENCE_DEFAULTS.SEQUENCES = MULTI_SEQUENCE_DEFAULTS.SEQUENCES_COPY.slice();
        inputViewmodel.sequences.removeAll();  // avoids changing on the as constant defined value

        return sequenceArray;
    }

    /**
     * Changes the output after processing the input.
     * @param outputData {Object} - Contains all output data.
     * @param inputProcessor {Object} - The unit processing the input.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions and input.
     */
    function changeOutput(outputData, inputProcessor, viewmodels) {
        // creates a visually representable distance matrix
        outputData.distanceMatrix
            = alignmentInterfaceInstance.getDistanceTable(outputData.distanceMatrix, outputData.distanceMatrixLength,
            outputData.remainingClusters[0], undefined);

        // distance matrix
        viewmodels.output.distanceMatrix(outputData.distanceMatrix);

        for (var i = 0; i < outputData.distanceMatrix.length; i++) {
            // new variables (rows) are not automatically functions
            // and so we have to convert new variables manually into functions
            // or we get the following error
            // 'Uncaught TypeError: viewmodels.output.distanceMatrix[i] is not a function'
            if (i > viewmodels.output.distanceMatrix.length)
                viewmodels.output.distanceMatrix[i] = new Function();

            viewmodels.output.distanceMatrix[i](outputData.distanceMatrix[i]);
        }

        // distance matrices
        outputData.distanceMatrices = alignmentInterfaceInstance.getDistanceTables(outputData);

        viewmodels.output.distanceMatrices(outputData.distanceMatrices);

        // iteration over each matrix
        for (var i = 0; i < outputData.distanceMatrices.length; i++) {
            // new variables (rows) are not automatically functions...
            if (i >= viewmodels.output.distanceMatrices.length)
                viewmodels.output.distanceMatrices[i] = new Function();

            viewmodels.output.distanceMatrices[i](outputData.distanceMatrices[i]);

            // iteration over each row of the matrix
            for (var j = 0; j < outputData.distanceMatrices[i].length; j++) {
                // new variables (rows) are not automatically functions...
                if (j >= viewmodels.output.distanceMatrices[i].length)
                    viewmodels.output.distanceMatrices[i][j] = new Function();

                viewmodels.output.distanceMatrices[i][j](outputData.distanceMatrices[i][j]);
            }
        }

        // merge steps
        viewmodels.output.guideAlignments(outputData.guideAlignments);
        viewmodels.output.firstGroups(outputData.firstGroups);
        viewmodels.output.secondGroups(outputData.secondGroups);
        viewmodels.output.joinedGroups(outputData.joinedGroups);
        viewmodels.output.joinedGroupNames(outputData.joinedGroupNames);

        // tree and final output
        viewmodels.output.newickString(outputData.newickString);
        viewmodels.output.progressiveAlignment(outputData.progressiveAlignment);

        viewmodels.visual.drawTree();

        // pairwise data
        viewmodels.output.sequencePairNames(outputData.sequencePairNames);
        viewmodels.output.alignmentLengths(outputData.alignmentLengths);
        viewmodels.output.similarities(outputData.similarities);
        viewmodels.output.gapNumbers(outputData.gapNumbers);
        viewmodels.output.gapStarts(outputData.gapStarts);
    }
}());
