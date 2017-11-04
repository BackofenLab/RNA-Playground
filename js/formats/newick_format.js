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

    var newickString;

    /**
     * Returns a Newick-Format representation of a phylogenetic tree
     * by a (depth-first-search) post order traversal of the tree.
     * @param tree - The phylogenetic tree from which you want get the representation.
     * @return {string} - The Newick string of the given tree.
     */
    function getEncoding(tree) {
        newickString = SYMBOLS.EMPTY;
        postOrder(tree, false);

        // hint: this is working, because value from last cluster is always 0
        // example: ..[CLUSTER_NAME]:0) -to-> ..[CLUSTER_NAME]
        newickString = newickString.slice(0, newickString.length-3) + SYMBOLS.SEMICOLON;
        return newickString;
    }

    /**
     * Does a post-order traversal to create the Newick-format string.
     * @param node {Object} - The node from which on traversal takes place. At the beginning it is the root node.
     * @param isLeftChild {boolean} - Tells if the node is a left child or not. It is needed to create the Newick-formatted string.
     * @see: Idea to use post-order-traversal came from code of Alexander Mattheis in project Algorithms for Bioninformatics.
     */
    function postOrder(node, isLeftChild) {
        if (node === undefined)
            return;

        if (isLeftChild)  // whenever you go left, you have to set a left bracket
            newickString += SYMBOLS.BRACKET_LEFT;

        postOrder(node.leftChild, true);
        postOrder(node.rightChild, false);

        var isLeaf = node.leftChild === undefined && node.rightChild === undefined;
        newickString += isLeaf ? node.name: SYMBOLS.EMPTY;
        newickString +=  SYMBOLS.COLON + Math.round(node.value * 10) / 10;  // rounding to one digit after decimal point

        if (isLeftChild) // whenever you wrote down a node name, you have to set a comma if it is a left child and else a right bracket
            newickString += SYMBOLS.COMMA;
        else
            newickString += SYMBOLS.BRACKET_RIGHT;
    }
}());