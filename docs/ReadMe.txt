/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

/*-------------------------*/
/*---------READ ME---------*/
/*-------------------------*/

/*-----------------------------------------------------DEVELOPMENT----------------------------------------------------*/
The DNA-bionformatics-part of
the project was developed by Alexander Mattheis (supervised by Martin Raden)
with the cross-platform IDE "PhpStorm 2017.2.4".

/*----------------------------------------------EXTENDING VISUALIZATION-----------------------------------------------*/
All functions which are used for visualization like highlighting, drawing arrows
and functions which are used to access information are stored in "js/post_processing/visualizer.js".
This class is working closely together with the "js/post_processing/input_processor.js" class
which is handling the input.

/*-----------------------------------------------------INHERITANCE----------------------------------------------------*/
In this project simulated inheritance is used because ECMAScript 5
does not support real, Java-like inheritance like in ECMAScript 6 (from 2015).
But ECMAScript 5 is also fully supported by older browsers and tablet/smartphone browsers.

Examples for simulated inheritance you find under "js/procedures/bases/":
from the classes there, other classes like Needleman-Wunsch inherit functions
by exchanging instances.
So, the parent class is executing the child-class by an instance
and the child class can exchange information with the parent class
through an instance of the parent class (with the related functions: getInput, getOutput, setIO).

/*---------------------------------------------------INPUT HANDLING---------------------------------------------------*/
The "InputProcessor"-class under "js/post_processing/input_processor.js"
is used for forwarding input to the algorithms and the visualizer.
It removes wrong inputs before data is passed to a algorithm
and it defines the behaviour of different input types.

Example:    Result selection
-   after a click on a result,
    the input processor is selecting a table entry of the results table (id: "results")
-   the table and the entry is then sent from this class to the visualizer
    which highlights the entry

/*-----------------------------------------------------INTERFACES-----------------------------------------------------*/
Interfaces under "interfaces/" can be used to offload the algorithm input and output (the interface).
Therefore the "updateDocumentView"-function in the "alignment.js" have to be extended.

This functions tells you which interface (input + output) have to be loaded with which algorithm.
The interface contents will be loaded into the div-element "<div id="algorithm_interface"></div>"
of the algorithm page.

For example: "gotoh.html" is loading "affine_alignment_interface.html" into this div-element

An interface is always controlled
by a javascript file under "js/interfaces/", which have the same name like the related page under "interfaces/",
but the file extension is ".js".

/*------------------------------------IMPLEMENTING NEW ALGORITHM FOR VISUALIZATION------------------------------------*/
-   create a new page "example_algorithm.html" in "/"
    (use underscores for the page name to make "alignment.js" recognizing the new page)
-   create a new script "example_algorithm.js" in "js/"
-   go into "index.html" and enter your algorithm name (without underscores) in the algorithm menu area:

    <div id="algorithm_menu">
        <ul>
            <li><a href="#">Example Algorithm</a></li>
        </ul>
    </div>
	
/*------------------------------------------------------LOADING-------------------------------------------------------*/
The loading of an algorithm-page into the "alignment.html" is done in the "alignment.js"-file.
In this file is a function "updateDocumentView(algorithm, view)" which loads the page into the "div" with
id = "algorithm_view".

/*-------------------------------------------------------PATHS--------------------------------------------------------*/
All Javascript-paths are stored under "js/defaults.js". 

Hint:	If you want change file-paths, then you have also to change paths in the "alignment.html".

/*-----------------------------------------------------PROCEDURES-----------------------------------------------------*/
To save loading time, the number of files which are sequentially loaded is limited
by storing all algorithms of same type in one file "js/procedures/".
For example all algorithms for backtracking are stored in "js/procedures/backtracking.js"
and all algorithms for clustering are stored in "js/procedures/clustering.js".

/*-------------------------------------------------------TESTS--------------------------------------------------------*/
Unit-Tests are stored under "js/tests/". The file "UnitTests.jstd" is specifying Unit-Test file paths.
The tool/plugin which was used to create the Unit-Tests is "JsTestDriver".
