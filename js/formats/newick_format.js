/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("formats.newickFormat", getEncoding);

    /**
     * Returns a Newick-Format representation of a phylogenetic tree
     * by a (depth-first-search) post order traversal of the tree.
     * @param tree - The phylogenetic tree from which you want get the representation.
     * @return {string} - The Newick string of the given tree.
     */
    function getEncoding(tree) {
        var newickString = SYMBOLS.BRACKET_LEFT;
        postOrder(tree, newickString, false);
        newickString += SYMBOLS.SEMICOLON;
        return newickString;
    }

    /**
     * Does a post-order traversal to create the Newick-format string.
     * @param node {Object} - The node from which on traversal takes place. At the beginning it is the root node.
     * @param newickString {string} - The string which has to be filled up. At the beginning it is empty.
     * @param isLeftChild {boolean} - Tells if the node is a left child or not. It is needed to create the Newick-formatted string.
     * @see: It is based on the code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function postOrder(node, newickString, isLeftChild) {
        if (node === undefined)
            return;

        postOrder(node.leftChild, newickString, true);
        postOrder(node.rightChild, newickString, false);

        if (isLeftChild)
            newickString += SYMBOLS.BRACKET_LEFT;

        newickString += node.name + SYMBOLS.COLON + Math.round(node.value * 10) / 10;  // rounding to one digit after decimal point

        if (isLeftChild)
            newickString += SYMBOLS.COMMA;
        else
            newickString += SYMBOLS.BRACKET_RIGHT;
    }
}());