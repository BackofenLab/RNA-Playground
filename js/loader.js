/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

var loaded = ALGORITHMS.NONE;  // tells globally which algorithm was loaded

(function () {  // namespace
    // public methods
    namespace("loader", updateDocumentView);

    /**
     * Loads the HTML-page of the algorithm into current document
     * by removing non-characters and replacing symbols in the algorithm name.
     * @param algorithm {string} - The algorithm name without extension you want to load.
     * @param view {Element} - The view in which you want load the page.
     */
    function updateDocumentView(algorithm, view) {
        reinitialize();
        removeOverlay($("#overlay"));

        var htmlName = algorithm
            .replace(MULTI_SYMBOLS.BRACKET_LEFT, SYMBOLS.EMPTY)
            .replace(MULTI_SYMBOLS.BRACKET_RIGHT, SYMBOLS.EMPTY)
            .replace(MULTI_SYMBOLS.G_LITTLE_SPECIAL, SYMBOLS.G_LITTLE)
            .replace(MULTI_SYMBOLS.SPACE, SYMBOLS.HYPHEN);

        var javascriptName = algorithm.toLowerCase()
            .replace(MULTI_SYMBOLS.BRACKET_LEFT, SYMBOLS.EMPTY)
            .replace(MULTI_SYMBOLS.BRACKET_RIGHT, SYMBOLS.EMPTY)
            .replace(MULTI_SYMBOLS.DELIMITER, SYMBOLS.SEPARATOR)
            .replace(MULTI_SYMBOLS.G_LITTLE_SPECIAL, SYMBOLS.G_LITTLE)
            .replace(MULTI_SYMBOLS.SPACE, SYMBOLS.SEPARATOR);

        loaded = javascriptName;

        view.load(PATHS.MAIN.PAGES + htmlName + FILE_EXTENSIONS.HYPERTEXT_MARKUP_LANGUAGE, function () {
            $.getScript(PATHS.MAIN.SCRIPTS + javascriptName + FILE_EXTENSIONS.JAVASCRIPT);
        });
    }

    /**
     * Handling reimports to initialize libraries.
     */
    function reinitialize() {
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
     * Removes an overlay.
     * @param overlay - The overlay you want to empty.
     */
    function removeOverlay(overlay) {
        overlay.empty();
    }
}());