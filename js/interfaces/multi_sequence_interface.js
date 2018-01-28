/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("interfaces.multiSequenceInterface", MultiSequenceInterface);

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

        // flags
        this.sequencesNumberChanged = false;
        this.lastNumberOfSequences = 0;

        // public class methods
        this.startMultiSequenceInterface = startMultiSequenceInterface;
        this.getMaxNumberOfAlignments = getMaxNumberOfAlignments;
    }

    /**
     * Function managing objects.
     * @param Algorithm {Object} - The algorithm which is started.
     * @param algorithmName {string} - The name of the algorithm which is started.
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
        //$.getScript(PATHS.MULTI_SEQUENCE_INTERFACE);  // very important, because other interfaces are also using this class
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
        var isTcoffee = algorithmName === ALGORITHMS.NOTREDAME_HIGGINS_HERINGA;
        var isIterativeRefinement = algorithmName === ALGORITHMS.ITERATIVE_REFINMENT;

        this.sequences = ko.observableArray(isIterativeRefinement ? ITERATIVE_SEQUENCE_DEFAULTS.SEQUENCES : MULTI_SEQUENCE_DEFAULTS.SEQUENCES);

        this.calculation = ko.observable(MULTI_SEQUENCE_DEFAULTS.CALCULATION);  // equal for all used MSA approaches

        // function
        this.baseCosts = ko.observable(isIterativeRefinement
            ? ITERATIVE_SEQUENCE_DEFAULTS.FUNCTION.BASE_COSTS : MULTI_SEQUENCE_DEFAULTS.FUNCTION.BASE_COSTS);

        this.enlargement = ko.observable(isIterativeRefinement
            ? ITERATIVE_SEQUENCE_DEFAULTS.FUNCTION.ENLARGEMENT : MULTI_SEQUENCE_DEFAULTS.FUNCTION.ENLARGEMENT);

        this.match = ko.observable(isIterativeRefinement
            ? ITERATIVE_SEQUENCE_DEFAULTS.FUNCTION.MATCH : MULTI_SEQUENCE_DEFAULTS.FUNCTION.MATCH);

        this.mismatch = ko.observable(isIterativeRefinement
            ? ITERATIVE_SEQUENCE_DEFAULTS.FUNCTION.MISMATCH : MULTI_SEQUENCE_DEFAULTS.FUNCTION.MISMATCH);

        if (isTcoffee) {
            multiSequenceInterfaceInstance.lastNumberOfSequences = viewmodel.sequences().length;

            this.baseCostsLocal = ko.observable(MULTI_SEQUENCE_DEFAULTS.FUNCTION.BASE_COSTS_LOCAL);
            this.enlargementLocal = ko.observable(MULTI_SEQUENCE_DEFAULTS.FUNCTION.ENLARGEMENT_LOCAL);
            this.matchLocal = ko.observable(MULTI_SEQUENCE_DEFAULTS.FUNCTION.MATCH_LOCAL);
            this.mismatchLocal = ko.observable(MULTI_SEQUENCE_DEFAULTS.FUNCTION.MISMATCH_LOCAL);

            this.globalAlignmentsPerSequencePair = ko.observable(MULTI_SEQUENCE_DEFAULTS.GLOBAL_ALIGNMENTS_PER_SEQUENCE);

            this.useLocalLibrary = ko.observable(MULTI_SEQUENCE_DEFAULTS.USE_LOCAL_LIBRARY);
            this.totalNumberAlignments = ko.observable(getMaxNumberOfAlignments(viewmodel));
            this.localAlignmentsPerSequencePair = ko.observable(MULTI_SEQUENCE_DEFAULTS.LOCAL_ALIGNMENTS_PER_SEQUENCE);

            // displayed dynamic formulas
            this.gapStartLocal = ko.computed(
                function () {
                    return Number(viewmodel.baseCostsLocal()) + Number(viewmodel.enlargementLocal());
                }
            );

            this.gapFunctionLocal = ko.computed(
                function getSelectedFormula() {
                    setTimeout(function () {  // to reinterpret in next statement dynamically created LaTeX-code
                        MathJax.Hub.Queue(["Typeset", MathJax.Hub])
                    }, REUPDATE_TIMEOUT_MS);

                    return getSubformula(viewmodel, true);
                }
            );
        } else if (isIterativeRefinement) {
            this.availableApproaches = ko.observableArray(ITERATIVE_SEQUENCE_DEFAULTS.APPROACHES);
            this.selectedApproach = ko.observableArray(ITERATIVE_SEQUENCE_DEFAULTS.STANDARD_APPROACH);

            this.availableOrders = ko.observableArray(ITERATIVE_SEQUENCE_DEFAULTS.ORDERS);
            this.selectedOrder = ko.observableArray(ITERATIVE_SEQUENCE_DEFAULTS.STANDARD_ORDER);
        }

        this.clusterNames = ko.computed(
            function () {
                return getClusterNames(viewmodel.sequences().length);
            }
        );

        this.addRow = function () {
            setTimeout(function () {  // to reinterpret in next statement dynamically created LaTeX-code
                MathJax.Hub.Queue(["Typeset", MathJax.Hub])
            }, REUPDATE_TIMEOUT_MS);

            viewmodel.sequences.push(SYMBOLS.EMPTY);
        };

        this.removeRow = function () {
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

                return getSubformula(viewmodel, false);
            }
        );
    }

    /**
     * Computes the maximum number of alignments.
     * @param inputViewmodel {Object} - The viewmodel which is used to find out the number of alignments.
     * @return {number} - The number of existing alignments.
     */
    function getMaxNumberOfAlignments(inputViewmodel) {
        var sequences = inputViewmodel.sequences();
        var hashTable = {};

        var uniqueSequences = sequences.filter(function (sequence) {
            if (hashTable.hasOwnProperty(sequence))
                return false;  // remove/filter

            hashTable[sequence] = true;

            return true;  // do not remove/ do not filter
        });

        var n = uniqueSequences.length;
        return (n * (n - 1)) / 2;  // Gauss-formula
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
                    + CLUSTER_NAMES[i % CLUSTER_NAMES.length] + LATEX.SUBORDINATE + LATEX.CURLY_BRACKET_LEFT + currentEpisode + LATEX.CURLY_BRACKET_RIGHT +
                    LATEX.MATH_REGION);
        }

        return clusterNames;
    }

    /**
     * Returns the LaTeX-code for sub-formulas like gap-functions of subadditive algorithms.
     * @param viewmodel {InputViewmodel} - The viewmodel of the view displaying the formula.
     * @param local {boolean} - Tells if local or global parameters should be used.
     * @return {string} - LaTeX code.
     */
    function getSubformula(viewmodel, local) {
        var string = LATEX.MATH_REGION;

        string += LATEX.SUB_FORMULAS.GOTOH_GAP_FUNCTION;

        string += SYMBOLS.EQUAL;

        if (local) {
            string += viewmodel.baseCostsLocal() >= 0
                ? viewmodel.baseCostsLocal() + SYMBOLS.PLUS
                : SYMBOLS.BRACKET_LEFT + viewmodel.baseCostsLocal() + SYMBOLS.BRACKET_RIGHT + SYMBOLS.PLUS;

            string += viewmodel.enlargementLocal() >= 0
                ? viewmodel.enlargementLocal() + LATEX.DOT + LATEX.FACTOR
                : SYMBOLS.BRACKET_LEFT + viewmodel.enlargementLocal() + SYMBOLS.BRACKET_RIGHT + LATEX.DOT + LATEX.FACTOR;
        } else {
            string += viewmodel.baseCosts() >= 0
                ? viewmodel.baseCosts() + SYMBOLS.PLUS
                : SYMBOLS.BRACKET_LEFT + viewmodel.baseCosts() + SYMBOLS.BRACKET_RIGHT + SYMBOLS.PLUS;

            string += viewmodel.enlargement() >= 0
                ? viewmodel.enlargement() + LATEX.DOT + LATEX.FACTOR
                : SYMBOLS.BRACKET_LEFT + viewmodel.enlargement() + SYMBOLS.BRACKET_RIGHT + LATEX.DOT + LATEX.FACTOR;
        }

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
            inputViewmodel.sequences(getSequencesArray(inputViewmodel));

            inputViewmodel.baseCosts(Number($("#base_costs").val()));
            inputViewmodel.enlargement(Number($("#enlargement").val()));
            inputViewmodel.match(Number($("#match").val()));
            inputViewmodel.mismatch(Number($("#mismatch").val()));

            if (algorithm.type === ALGORITHMS.ITERATIVE_REFINMENT) {
                inputViewmodel.selectedApproach([$("#approach_selector option:selected").val()]);
                inputViewmodel.selectedOrder([$("#order_selector option:selected").val()]);
            } else if (algorithm.type === ALGORITHMS.NOTREDAME_HIGGINS_HERINGA) {
                inputViewmodel.baseCostsLocal(Number($("#base_costs_local").val()));
                inputViewmodel.enlargementLocal(Number($("#enlargement_local").val()));
                inputViewmodel.matchLocal(Number($("#match_local").val()));
                inputViewmodel.mismatchLocal(Number($("#mismatch_local").val()));

                inputViewmodel.globalAlignmentsPerSequencePair(Number($("#global_alignments_per_sequence_pair").val()));
                inputViewmodel.localAlignmentsPerSequencePair(Number($("#local_alignments_per_sequence_pair").val()));

                var maxNumberOfAlignments = getMaxNumberOfAlignments(inputViewmodel);

                if (multiSequenceInterfaceInstance.sequencesNumberChanged) { // set to maximum
                    inputViewmodel.totalNumberAlignments(maxNumberOfAlignments);
                } else { // use user specified value, if possible
                    var alignmentNumber = Number($("#total_number_alignments").val());

                    if (alignmentNumber < 0)
                        inputViewmodel.totalNumberAlignments(0);
                    else if (alignmentNumber > maxNumberOfAlignments)
                        inputViewmodel.totalNumberAlignments(maxNumberOfAlignments);
                    else
                        inputViewmodel.totalNumberAlignments(alignmentNumber);
                }
            }
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

        // setting flag
        var currentNumberOfSequences = inputViewmodel.sequences().length;
        var lastNumberOfSequences = multiSequenceInterfaceInstance.lastNumberOfSequences;

        multiSequenceInterfaceInstance.sequencesNumberChanged = currentNumberOfSequences !== lastNumberOfSequences;
        multiSequenceInterfaceInstance.lastNumberOfSequences = currentNumberOfSequences;  // change to current value

        // change value
        for (var i = 0; i < currentNumberOfSequences; i++) {
            sequenceArray.push($(".sequence_multi")[i].value.toUpperCase());
        }

        /* bug-fix for a Knockout-problem -> dynamically generated inputs get wrong values after typing in something */
        MULTI_SEQUENCE_DEFAULTS.SEQUENCES = MULTI_SEQUENCE_DEFAULTS.SEQUENCES_COPY.slice();
        ITERATIVE_SEQUENCE_DEFAULTS.SEQUENCES = ITERATIVE_SEQUENCE_DEFAULTS.SEQUENCES_COPY.slice();
        inputViewmodel.sequences.removeAll();  // avoids changing on the as constant defined value

        return sequenceArray;
    }

    /**
     * Changes the output after processing the input.
     * @param outputData {Object} - Contains all output data.
     * @param inputProcessor {Object} - The unit processing the input.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions and input.
     * @see Hint: The parameter inputProcessor is needed!
     */
    function changeOutput(outputData, inputProcessor, viewmodels) {
        if (viewmodels.visual.algorithm.type === ALGORITHMS.FENG_DOOLITTLE)
            changeFengDoolittleOutput(outputData, viewmodels);
        else if (viewmodels.visual.algorithm.type === ALGORITHMS.ITERATIVE_REFINMENT)
            changeIterativeRefinementOutput(outputData, viewmodels);
        else if (viewmodels.visual.algorithm.type === ALGORITHMS.NOTREDAME_HIGGINS_HERINGA)
            changeTcoffeeOutput(outputData, viewmodels);
    }

    /**
     * Changes the output of Feng-Doolittle algorithm.
     * @param outputData {Object} - Contains all output data.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions and input.
     */
    function changeFengDoolittleOutput(outputData, viewmodels) {
        // distance matrices
        outputData.distanceMatrices = alignmentInterfaceInstance.getDistanceTables(outputData, false, true);

        alignmentInterfaceInstance.roundValues(viewmodels.visual.algorithm.type, outputData);

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

        viewmodels.output.remainingClusters(outputData.remainingClusters);
        viewmodels.output.minimums(outputData.minimums);

        // merge steps
        alignmentInterfaceInstance.reorderGroupSequences(outputData);
        viewmodels.output.guideAlignments(outputData.guideAlignments);
        viewmodels.output.guideAlignmentsNames(outputData.guideAlignmentsNames);
        viewmodels.output.firstGroups(outputData.firstGroups);
        viewmodels.output.secondGroups(outputData.secondGroups);
        viewmodels.output.firstGroupsNames(outputData.firstGroupsNames);
        viewmodels.output.secondGroupsNames(outputData.secondGroupsNames);
        viewmodels.output.joinedGroups(outputData.joinedGroups);
        viewmodels.output.joinedGroupNames(outputData.joinedGroupNames);

        // tree and final output
        viewmodels.output.newickString(outputData.newickString);
        viewmodels.output.progressiveAlignment(outputData.progressiveAlignment);
        viewmodels.output.score(outputData.score);

        viewmodels.visual.drawTree();

        // pairwise data
        alignmentInterfaceInstance.sortWithClusterTuples(outputData.sequencePairNames,
            [outputData.alignmentLengths, outputData.similarities, outputData.gapNumbers, outputData.gapStarts]);
        viewmodels.output.sequencePairNames(outputData.sequencePairNames);
        viewmodels.output.alignmentLengths(outputData.alignmentLengths);
        viewmodels.output.similarities(outputData.similarities);
        viewmodels.output.gapNumbers(outputData.gapNumbers);
        viewmodels.output.gapStarts(outputData.gapStarts);
    }

    /**
     * Changes the output of Notredame-Higgins-Heringa algorithm.
     * @param outputData {Object} - Contains all output data.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions and input.
     */
    function changeTcoffeeOutput(outputData, viewmodels) {
        outputData.librariesData = alignmentInterfaceInstance.getLibrariesData(outputData);

        alignmentInterfaceInstance.removeNeutralSymbols(outputData);
        alignmentInterfaceInstance.roundValues(viewmodels.visual.algorithm.type, outputData);
        alignmentInterfaceInstance.reorderGroupSequences(outputData);

        // final output
        viewmodels.output.progressiveAlignment(outputData.progressiveAlignment);
        viewmodels.output.score(outputData.score);

        // merge steps
        viewmodels.output.firstGroups(outputData.firstGroups);
        viewmodels.output.secondGroups(outputData.secondGroups);
        viewmodels.output.firstGroupsNames(outputData.firstGroupsNames);
        viewmodels.output.secondGroupsNames(outputData.secondGroupsNames);
        viewmodels.output.joinedGroups(outputData.joinedGroups);
        viewmodels.output.joinedGroupNames(outputData.joinedGroupNames);

        // tree
        viewmodels.output.newickString(outputData.newickString);
        viewmodels.visual.drawTree();

        // alignments (using index of sequencePairsNames, so they have to stay above)
        viewmodels.output.alignmentsGlobal(outputData.librariesData[4]);
        viewmodels.output.alignmentsLocal(outputData.librariesData[5]);

        // libraries
        viewmodels.output.sequencePairsNames(outputData.librariesData[0]);
        viewmodels.output.libPositionPairs(outputData.librariesData[1]);
        viewmodels.output.primLibValues(outputData.librariesData[2]);
        viewmodels.output.extendedLibValues(outputData.librariesData[3]);
    }

    /**
     * Changes the output of an iterative refinement algorithm.
     * @param outputData {Object} - Contains all output data.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions and input.
     */
    function changeIterativeRefinementOutput(outputData, viewmodels) {
        // final output
        alignmentInterfaceInstance.reorderFinalAlignments(outputData);  // do not move down this function
        viewmodels.output.progressiveAlignment(outputData.progressiveAlignment);
        viewmodels.output.score(outputData.score);
        viewmodels.output.refinedProgressiveAlignment(outputData.refinedProgressiveAlignment);
        viewmodels.output.refinedScore(outputData.refinedScore);

        // realignment steps
        alignmentInterfaceInstance.reorderGroupSequences(outputData);   // do not move up this function
        viewmodels.output.guideAlignments(outputData.guideAlignments);
        viewmodels.output.guideAlignmentsNames(outputData.guideAlignmentsNames);

        viewmodels.output.firstGroups(outputData.firstGroups);
        viewmodels.output.firstGroupsNames(outputData.firstGroupsNames);

        viewmodels.output.secondGroups(outputData.secondGroups);
        viewmodels.output.secondGroupsNames(outputData.secondGroupsNames);

        viewmodels.output.joinedGroups(outputData.joinedGroups);
        viewmodels.output.realignmentsScores(outputData.realignmentsScores);
        viewmodels.output.joinedGroupNames(outputData.joinedGroupNames);

        viewmodels.output.accepted(outputData.accepted);

        // tree
        viewmodels.output.newickString(outputData.newickString);
        viewmodels.visual.drawTree();

        // distance matrix
        viewmodels.output.remainingClusters(outputData.remainingClusters);

        outputData.distanceMatrix
            = bases.clustering.getMatrixAsTable(outputData.distanceMatrix, outputData.distanceMatrixLength, outputData.remainingClusters[0], undefined, true);

        alignmentInterfaceInstance.roundValues(viewmodels.visual.algorithm.type, outputData);

        viewmodels.output.distanceMatrix(outputData.distanceMatrix);

        // iteration over each row of the matrix
        for (var i = 0; i < outputData.distanceMatrix.length; i++) {
            // new variables (rows) are not automatically functions...
            if (i >= viewmodels.output.distanceMatrix.length)
                viewmodels.output.distanceMatrix[i] = new Function();

            viewmodels.output.distanceMatrix[i](outputData.distanceMatrix[i]);
        }
    }
}());
