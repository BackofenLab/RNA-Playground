/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

"use strict";

(function () {  // namespace
    // public methods
    namespace("formats.newickEncoder", NewickEncoder);

    var newickEncoderInstance;

    /**
     * Creates an encoder which can encode binary trees
     * of phylogenetic data into the newick format.
     * @constructor
     */
    function NewickEncoder() {
        newickEncoderInstance = this;

        // variables
        this.newickString = SYMBOLS.EMPTY;

        // public class methods
        this.getEncoding = getEncoding;
    }

    /**
     * Returns a Newick-Format representation of a phylogenetic tree
     * by a (depth-first-search) post order traversal of the tree.
     * @param tree {Object} - The phylogenetic tree from which you want get the representation.
     * @return {string} - The Newick string of the given tree.
     */
    function getEncoding(tree) {
        newickEncoderInstance.newickString = SYMBOLS.EMPTY;
        postOrder(tree, false);

        // hint: this is working, because value from last cluster is always 0
        // example: ..[CLUSTER_NAME]:0) -to-> ..[CLUSTER_NAME]
        newickEncoderInstance.newickString
            = newickEncoderInstance.newickString.slice(0, newickEncoderInstance.newickString.length - 3) + SYMBOLS.SEMICOLON;

        return newickEncoderInstance.newickString;
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
            newickEncoderInstance.newickString += SYMBOLS.BRACKET_LEFT;

        postOrder(node.leftChild, true);
        postOrder(node.rightChild, false);

        var isLeaf = node.rightChild === undefined && node.leftChild === undefined;
        newickEncoderInstance.newickString += isLeaf ? node.name : SYMBOLS.EMPTY;
        // Hint: Do not change the rounding or the displayed tree could get wrong, because of floating-point errors!
        newickEncoderInstance.newickString += SYMBOLS.COLON + Math.round(node.value * 10000) / 10000;  // rounded to four digits after decimal point

        if (isLeftChild) // whenever you wrote down a node name, you have to set a comma if it is a left child and else a right bracket
            newickEncoderInstance.newickString += SYMBOLS.COMMA;
        else
            newickEncoderInstance.newickString += SYMBOLS.BRACKET_RIGHT;
    }
}());