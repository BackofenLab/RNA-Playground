/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

/**
 * Defines tasks after page-loading.
 */
$(document).ready(function () {
    dropDownMenu.startDropDownMenu();
});

(function () {  // namespace
    // public methods
    namespace("dropDownMenu", startDropDownMenu);

    /**
     * Function managing objects.
     */
    function startDropDownMenu() {
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
            viewToUpdate: $("#algorithm_view")
        };

        algorithmMenu.find("ul li").on("click", loadAlgorithmArguments, loadAlgorithm);
    }

    /**
     * Shows a hidden menu or closes an opened menu.
     * It updates the menu button during this process.
     * @param e {Object} - Stores data relevant to the event called that function.
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
     * @param e {Object} - Stores data relevant to the event called that function.
     */
    function loadAlgorithm(e) {
        var algorithm = $(this).text();
        var display = e.data.displayToUpdate;
        var dropDown = e.data.dropDown;
        var dropDownName = e.data.dropDownName;
        var menu = e.data.menu;
        var newDropDownName = e.data.newDropDownName;
        var view = e.data.viewToUpdate;

        updateDisplay(display, algorithm);
        loader.updateDocumentView(algorithm, view);

        menu.hide();
        dropDown.addClass(dropDownName).removeClass(newDropDownName);
    }

    /**
     * Updates a display placed on the HTML document with the text of the caller object.
     * @param {Element} display - The display to update.
     * @param {string} displayText - The new text for the display.
     */
    function updateDisplay(display, displayText) {
        display.text(displayText);
    }
}());