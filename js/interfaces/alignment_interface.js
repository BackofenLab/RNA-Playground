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
        imports, sharedInterfaceOperations, roundValues, getDistanceTable, getDistanceTables,
        reorderGroupSequences, getLibrariesData, startProcessing);

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
        this.roundValues = roundValues;
        this.getDistanceTable = getDistanceTable;
        this.getDistanceTables = getDistanceTables;
        this.reorderGroupSequences = reorderGroupSequences;
        this.getLibrariesData = getLibrariesData;
        this.startProcessing = startProcessing;
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
        if (algorithm.type === ALGORITHMS.FENG_DOOLITTLE || algorithm.type === ALGORITHMS.NOTREDAME_HIGGINS_HERINGA)
            viewmodels.visual.drawTree();
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
            roundValues(algorithmName, outputData);

            if (algorithmName === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER) {  // if AEP
                createAEPOutputViewmodel(viewmodel, outputData);
            } else if (outputData.matrix !== undefined) {  // other algorithms
                createOutputViewmodel(algorithmName, viewmodel, outputData);
            }
        } else if (MULTI_SEQUENCE_ALGORITHMS.indexOf(algorithmName) >= 0) {  // if multi-sequence alignment algorithm
            if (algorithmName === ALGORITHMS.FENG_DOOLITTLE)
                createFengDoolittleOutputViewmodel(algorithmName, viewmodel, outputData);
            else
                createTcoffeeOutputViewmodel(algorithmName, viewmodel, outputData);
        }
    }

    /**
     * Rounds matrix values, scores and other parameters.
     * @param algorithmName {string} - The name of the algorithm which is executed.
     * @param outputData {Object} - Output data which is modified.
     */
    function roundValues(algorithmName, outputData) {
        if (algorithmName === ALGORITHMS.ARSLAN_EGECIOGLU_PEVZNER) {
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
        } else if (algorithmName === ALGORITHMS.FENG_DOOLITTLE) {  // if Feng-Doolittle or UPGMA
            // iterate over each distance matrix
            for (var k = 0; k < outputData.distanceMatrices.length; k++) {

                // iterate over each row
                for (var i = 0; i < outputData.distanceMatrices[k].length; i++) {

                    // iterate over each entry
                    for (var j = 0; j < outputData.distanceMatrices[k][i].length; j++) {
                        if (j > i)  // only the values upper the diagonal
                            outputData.distanceMatrices[k][i][j]
                                = round(outputData.distanceMatrices[k][i][j], 1);
                    }
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
     * @param algorithmName {string} - The name of the algorithm which is executed.
     * @param viewmodel {Object} - The output viewmodel container which should be filled.
     * @param outputData {Object} - The data which is used to fill the viewmodel.
     */
    function createOutputViewmodel(algorithmName, viewmodel, outputData) {
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
        // distance matrix
        outputData.distanceMatrix
            = getDistanceTable(outputData.distanceMatrix, outputData.distanceMatrixLength, outputData.remainingClusters[0], undefined);

        viewmodel.distanceMatrix =  ko.observableArray(outputData.distanceMatrix);

        for (var i = 0; i < outputData.distanceMatrix.length; i++) {
            viewmodel.distanceMatrix[i] = ko.observableArray(outputData.distanceMatrix[i]);
        }

        // distance matrices
        outputData.distanceMatrices = getDistanceTables(outputData);

        roundValues(algorithmName, outputData);

        viewmodel.distanceMatrices = ko.observableArray(outputData.distanceMatrices).extend({ deferred: true });

        // iteration over each matrix
        for (var i = 0; i < outputData.distanceMatrices.length; i++) {
            viewmodel.distanceMatrices[i] = ko.observableArray(outputData.distanceMatrices[i]).extend({ deferred: true });

            // iteration over each row of the matrix
            for (var j = 0; j < outputData.distanceMatrices[i].length; j++) {
                viewmodel.distanceMatrices[i][j] = ko.observableArray(outputData.distanceMatrices[i][j]).extend({ deferred: true });
            }
        }

        viewmodel.remainingClusters = ko.observable(outputData.remainingClusters).extend({ deferred: true });
        viewmodel.minimums = ko.observable(outputData.minimums).extend({ deferred: true });

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
        viewmodel.sequencePairNames = ko.observable(outputData.sequencePairNames);
        viewmodel.alignmentLengths = ko.observable(outputData.alignmentLengths);
        viewmodel.similarities = ko.observable(outputData.similarities);
        viewmodel.gapNumbers = ko.observable(outputData.gapNumbers);
        viewmodel.gapStarts = ko.observable(outputData.gapStarts);

        viewmodel.showMatrices = ko.observable(false);

        // gimmicks/optimizations
        viewmodel.toggleVisibility = function() {
            viewmodel.showMatrices(!viewmodel.showMatrices());
        };

        viewmodel.toggleLinkText = ko.computed(
            function () {
                return viewmodel.showMatrices() ? TOGGLE_LINK_TEXT.HIDE : TOGGLE_LINK_TEXT.SHOW;
            }
        );
    }

    /**
     * Converts the distances stored in associative array into a real distance matrix.
     * @param outputData - The output on which conversion is applied.
     * @return {Object} - The outputData with converted distance matrices.
     */
    function getDistanceTables(outputData) {
        var matrixLength = outputData.distanceMatrixLength;  // start length

        // in each round the matrix gets smaller by one, because two matrices are merged
        for (var i = 0; i < outputData.distanceMatrices.length; i++) {
            outputData.distanceMatrices[i]
                = getDistanceTable(outputData.distanceMatrices[i], matrixLength-i, outputData.remainingClusters[i], outputData.keys[i]);
        }

        return outputData.distanceMatrices;
    }

    /**
     * The distance matrix is an "associative array" and this has to be converted
     * into a 2D-array which is displayable.
     * Hint: "Associative arrays" do not have a defined order (browser-dependant).
     */
    function getDistanceTable(distanceMatrix, distanceMatrixLength, remainingClusters, matrixKeys) {
        var matrix = createMatrix(distanceMatrixLength);
        if (matrixKeys === undefined)
            matrixKeys = Object.keys(distanceMatrix);  // argument possibilities {a,b}, {a,c}, ...

        // fill diagonals with zero
        for (var i = 0; i < matrix.length; i++) {
            for (var j = 0; j < matrix.length; j++) {
                if (i === j)
                    matrix[i][j] = 0;
            }
        }

        // fill right upper half
        for (var j = 0; j < matrixKeys.length; j++) {
            var key = matrixKeys[j].split(SYMBOLS.COMMA);
            var cluster1Position = getPositionByName(key[0], remainingClusters);
            var cluster2Position = getPositionByName(key[1], remainingClusters);
            var value = distanceMatrix[key];

            matrix[cluster1Position][cluster2Position] = value;
        }

        return matrix;
    }

    /**
     * Creates a matrix with the given size.
     * @param size - The width and height of the matrix.
     */
    function createMatrix (size) {
        var matrix = new Array(size);

        for (var i = 0; i < size; i++) {
            matrix[i] = [];
        }

        return matrix;
    }

    /**
     * Returns for a cluster-name, its position in the distance matrix.
     * @param clusterName {string} - The name of the cluster.
     */
    function getPositionByName(clusterName, remainingClusterNames) {
        var position = -1;

        for (var i = 0; i < remainingClusterNames.length; i++) {
            if (clusterName === remainingClusterNames[i]) {
                position = i;
                break;
            }
        }

        return position;
    }

    /**
     * Reordering groups in alphabetical order for increasing readability.
     * @param outputData - The output on which conversion is applied.
     */
    function reorderGroupSequences(outputData) {
        if (outputData.joinedGroupNames.length > 0) {
            var finalGroupName = outputData.joinedGroupNames[outputData.joinedGroupNames.length-1];
            var finalGroup = outputData.joinedGroups[outputData.joinedGroups.length-1];

            var groupMemberNames = getIndividualElementNames(finalGroupName);
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

                var memberNames = getIndividualElementNames(outputData.joinedGroupNames[i]);
                var member1Names = getIndividualElementNames(outputData.firstGroupsNames[i]);
                var member2Names = getIndividualElementNames(outputData.secondGroupsNames[i]);

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
     * Returns the individual names of the group members,
     * where a name character separated by a comma from the name number.
     * @param groupName - The group name from which the names extracted.
     * @return {Array} - The array with the individual names.
     */
    function getIndividualElementNames(groupName) {
        var names = [];

        for (var i = 0; i < groupName.length; i++) {
            var character = groupName[i];
            var number = SYMBOLS.EMPTY;

            while (i + 1 < groupName.length && groupName[i + 1].match(CHARACTER.NUMBER)) {
                number += groupName[i + 1];
                i++;
            }
            names.push(number.length > 0 ? character + SYMBOLS.COMMA + number : character);
        }

        return names;
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
            var number = name.length > 1 ? Number(name[1] - 1) : 0;

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
     * @return {[sortedGroupMembers, sortedGroupMemberNames]} - The sorted group and names.
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

        roundValues(algorithmName, outputData);
        alignmentInterfaceInstance.reorderGroupSequences(outputData);

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
    }

    /**
     * Returns the data needed to display from primary and extended library.
     * @param outputData {Object} - The data which is used to fill the viewmodel.
     * @return {[sequencePairsNames, positionPairs, primLibValues, extendedLibValues]} - The data from primary and extended library.
     */
    function getLibrariesData(outputData) {
        var sequencePairsNames = [];
        var positionPairs = [];
        var primLibValues = [];
        var extendedLibValues = [];

        var alignmentKeys = Object.keys(outputData.primaryWeightLib);

        // iterate overall alignments
        for (var i = 0; i < alignmentKeys.length; i++) {
            var alignmentKey = alignmentKeys[i];
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

                tempPositionPairs.push([splittedPositionKey[0], splittedPositionKey[1]]);
                tempPrimLibValues.push(valueL);
                tempExtendedLibValues.push(valueEL);
            }

            if (tempPositionPairs.length !== 0) {  // don't show names of sequence pairs for which no L or EL exists
                sort(tempPositionPairs, tempPrimLibValues, tempExtendedLibValues);
                
                sequencePairsNames.push([sequence1Name, sequence2Name]);
                positionPairs.push(tempPositionPairs);
                primLibValues.push(tempPrimLibValues);
                extendedLibValues.push(tempExtendedLibValues);
            }
        }

        return [sequencePairsNames, positionPairs, primLibValues, extendedLibValues];
    }

    /**
     * Returns numerically sorted input arrays.
     * @param positionPairs {Array} - Array of tuples which is sorted and used to sort the other two arrays.
     * @param primLibValues {Array} - Array of primary lib values.
     * @param extendedLibValues {Array} - Array of extended lib values.
     * @return {[sortedPositionPairs, sortedPrimLibValues, sortedExtendedLibValues]} - The sorted arrays as a triple.
     */
    function sort(positionPairs, primLibValues, extendedLibValues) {
        var switches = [];

        // documentation {sort} - https://developer.mozilla.org/de/docs/Web/JavaScript/Reference/Global_Objects/Array/sort
        var sortedPositionPairs = positionPairs.sort(function (a,b) {
            var leftNumberA = Number(a[0]);
            var leftNumberB = Number(b[0]);

            var rightNumberA = Number(a[1]);
            var rightNumberB = Number(b[1]);

            var value = 0;

            if (leftNumberB === leftNumberA) {
                value = rightNumberB > rightNumberA ? -1 : (rightNumberB > rightNumberA ? 1 : 0);
                switches.push(value);
                return value;
            }

            var value = leftNumberA - leftNumberB;
            switches.push(value);
            return value;
        });

        var i = 0;

        var sortedPrimLibValues = primLibValues.sort(function (a,b) {
            return switches[i++];
        });

        var i = 0;

        var sortedExtendedLibValuess = extendedLibValues.sort(function (a,b) {
            return switches[i++];
        });
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
