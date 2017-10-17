/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

var loaded = ALGORITHMS.NONE;  // tells globally which algorithm was loaded

/**
 * Defines tasks after page-loading.
 */
$(document).ready(function () {
    // to avoid the execution of the drop down menu
    if (typeof ALIGNMENT_WEBTITLE !== "undefined" && document.title === ALIGNMENT_WEBTITLE) {
    	loader.startLoader();
    } 
});

(function () {  // namespace
    // public methods
    namespace("loader", startLoader, updateDocumentView);

    /**
     * Function managing objects.
     */
    function startLoader() {
        linkElements();
    }

    /**
     * Linking elements with a function.
     */
    function linkElements() {
        var algorithmMenu = $("#algorithm_menu");
        var dropDownButton = $(".dropdown_button");

        var showHideMenuArguments = {
            dropDownName: "dropdown_button",
            menu: algorithmMenu,
            newDropDownName: "dropdown_reversed"
        };

        dropDownButton.on("click", showHideMenuArguments, showHideMenu);

        var loadAlgorithmArguments = {
            displayToUpdate: $("#algorithm_display"),
            dropDown: dropDownButton,
            dropDownName: "dropdown_button",
            menu: algorithmMenu,
            newDropDownName: "dropdown_reversed",
            overlay: $("#overlay"),
            viewToUpdate: $("#algorithm_view")
        };

        algorithmMenu.find("ul li").on("click", loadAlgorithmArguments, loadAlgorithm);
    }

    /**
     * Shows a hidden menu or closes an opened menu.
     * It updates the menu button during this process.
     * @param e - Stores data relevant to the event called that function.
     */
    function showHideMenu(e) {
        var menu = e.data.menu;
        var dropDownName = e.data.dropDownName;
        var newDropDownName = e.data.newDropDownName;

        if (menu.css("display") === "none") {
            menu.show();
            $(this).removeClass(dropDownName).addClass(newDropDownName);
        } else {
            menu.hide();
            $(this).removeClass(newDropDownName).addClass(dropDownName);
        }
    }

    /**
     * Loads an algorithm with its HTML-page
     * and updates the display showing current selected algorithm.
     * Then it closes the algorithm menu again.
     * @param e - Stores data relevant to the event called that function.
     */
    function loadAlgorithm(e) {
        var algorithm = $(this).text();
        var display = e.data.displayToUpdate;
        var dropDown = e.data.dropDown;
        var dropDownName = e.data.dropDownName;
        var menu = e.data.menu;
        var newDropDownName = e.data.newDropDownName;
        var overlay = e.data.overlay;
        var view = e.data.viewToUpdate;

        updateDisplay(display, algorithm);
        updateDocumentView(algorithm, view);

        menu.hide();
        dropDown.addClass(dropDownName).removeClass(newDropDownName);
    }

    /**
     * Updates a display placed on the HTML document with the text of the caller object.
     * @param {string} display - The display to update.
     * @param {string} displayText - The new text for the display.
     */
    function updateDisplay(display, displayText) {
        display.text(displayText);
    }

    /**
     * Loads the HTML-page of the algorithm into current document
     * by removing non-characters and replacing symbols in the algorithm name.
     * @param algorithm - The algorithm name without extension you want to load.
     * @param view - The view in which you want load the page.
     */
    function updateDocumentView(algorithm, view) {
        removeOverlay($("#overlay"));

        var htmlName = algorithm
            .replace(MULTI_SYMBOLS.BRACKET_LEFT, SYMBOLS.EMPTY)
            .replace(MULTI_SYMBOLS.BRACKET_RIGHT, SYMBOLS.EMPTY)
            .replace(MULTI_SYMBOLS.G_LITTLE_SPECIAL, SYMBOLS.G_LITTLE)
            .replace(MULTI_SYMBOLS.SPACE, SYMBOLS.EMPTY);

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
     * Removes an overlay.
     * @param overlay - The overlay you want to remove.
     */
    function removeOverlay(overlay) {
        overlay.empty();
    }

}());