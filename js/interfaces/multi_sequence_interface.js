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
    var subadditiveAlignmentInterfaceInstance;

    /**
     * Is used to work with the input and output (the interface) of an affine alignment algorithm.
     * @constructor
     * @augments AlignmentInterface
     */
    function MultiSequenceInterface() {
        subadditiveAlignmentInterfaceInstance = this;

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

        // displayed dynamic formulas
        this.gapStart = ko.computed(
            function () {
                return Number(viewmodel.baseCosts()) + Number(viewmodel.enlargement());
            }
        );

        this.gapFunction = ko.computed(
            function getSelectedFormula() {
                setTimeout(function () {
                    MathJax.Hub.Queue(["Typeset", MathJax.Hub])
                }, REUPDATE_TIMEOUT_MS);

                return getSubformula(viewmodel);
            }
        );
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
     * @param algorithm {Object} - Algorithm used to update the user interface.
     * @param inputProcessor {Object} - The unit processing the input.
     * @param inputViewmodel {Object} - The InputViewmodel used to access inputs.
     * @param visualViewmodel {Object} - The VisualViewmodel used to access visualization functions.
     */
    function processInput(algorithm, inputProcessor, inputViewmodel, visualViewmodel) {
        visualViewmodel.removeAllContents();

        // when page was loaded the inputs have not to be updated or you get wrong inputs
        if (inputProcessor.inputUpdatesActivated()) {
            // read out dynamically created inputs

            inputViewmodel.baseCosts(Number($("#base_costs").val()));
            inputViewmodel.enlargement(Number($("#enlargement").val()));
            inputViewmodel.match(Number($("#match").val()));
            inputViewmodel.mismatch(Number($("#mismatch").val()));
        } else
            inputProcessor.activateInputUpdates();

        alignmentInterfaceInstance.startProcessing(algorithm, inputViewmodel, visualViewmodel);
    }

    /**
     * Changes the output after processing the input.
     * @param outputData {Object} - Contains all output data.
     * @param inputProcessor {Object} - The unit processing the input.
     * @param viewmodels {Object} - The viewmodels used to access visualization functions and input.
     */
    function changeOutput(outputData, inputProcessor, viewmodels) {
    }
}());
