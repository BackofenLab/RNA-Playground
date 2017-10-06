/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

$(document).ready(startIndex);

/**
 * Function managing objects.
 */
function startIndex() {
    index.linkElements();
}

(function () {  // namespace
    // public methods
    namespace("index", linkElements);

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
     * And it updates the menu button during this process.
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

        removeOverlay(overlay);
        updateDisplay(display, algorithm);
        updateDocumentView(algorithm, view);

        menu.hide();
        dropDown.addClass(dropDownName).removeClass(newDropDownName);
    }

    /**
     * Removes an overlay.
     * @param overlay - The overlay you want to remove.
     */
    function removeOverlay(overlay) {
        overlay.empty();
    }

    /**
     * Updates a display placed on the HTML document with the text of the caller object.
     * @param {String} display - The display to update.
     * @param {String} displayText - The new text for the display.
     */
    function updateDisplay(display, displayText) {
        display.text(displayText);
    }

    /**
     * Loads the HTML-page of the algorithm into current document.
     * @param algorithm - The HTML-filename without extension you want to load.
     * Hint: "-" are replaced with "_".
     * @param view - The view in which you want load the page.
     */
    function updateDocumentView(algorithm, view) {
        debugger;
        var algorithmName = algorithm.toLowerCase().replace(MULTI_SYMBOLS.DELIMITER, SYMBOLS.SEPARATOR)
            .replace(MULTI_SYMBOLS.SPACE, SYMBOLS.SEPARATOR);  // "-" and " " are replaced with "_"

        view.load(PATHS.MAIN.PAGES + algorithm + FILE_EXTENSIONS.HYPERTEXT_MARKUP_LANGUAGE, function () {
            $.getScript(PATHS.MAIN.SCRIPTS + algorithmName + FILE_EXTENSIONS.JAVASCRIPT);
        });
    }
}());