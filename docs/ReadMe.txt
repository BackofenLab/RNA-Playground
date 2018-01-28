/*
University of Freiburg WS 2017/2018
Chair for Bioinformatics
Supervisor: Martin Raden
Author: Alexander Mattheis
*/

/*-------------------------*/
/*---------READ ME---------*/
/*-------------------------*/

/*-------------------------------------------------------COMMENTS-----------------------------------------------------*/
In comments use this " (symbol) and not this ' (symbol) to emphasize something.

/*------------------------------------------------------CONSTANTS-----------------------------------------------------*/
Nearly all constants are stored in the "defaults.js".
Namespaces, HTML class-names and events like "mouse-over" etc. are not stored in the "defaults.js"
because else it would slow down the development process.
If you have multiple constants of same "logic" (for example with similar naming)
then you store them in a structure.

/*-----------------------------------------------------DEVELOPMENT----------------------------------------------------*/
The "Bioinformatics Algorithms" part of
the project (see "alignment.html") was developed by Alexander Mattheis (supervised by Martin Raden)
with the cross-platform IDE "PhpStorm 2017.2.4".

/*----------------------------------------------EXTENDING VISUALIZATION-----------------------------------------------*/
All functions which are used for visualization like highlighting, drawing arrows
and functions which are used to create downloadable files are stored in "js/post_processing/visualizer.js".
This class is working closely together with the "js/post_processing/input_processor.js" class
which is handling the input.

Hint:
Algorithm names have to be eventually added
into the lists of algorithms (SVG_ARROW_ALGORITHMS, MULTI_TABLE_ALGORITHMS, ...) in the "default.js"
to activate automatically an algorithm visualization
(use same HTML-element naming (like for example in the "Needleman-Wunsch.html") to get a visualization.

/*------------------------------------IMPLEMENTING NEW ALGORITHM FOR VISUALIZATION------------------------------------*/
-   create a new page "Example-Algorithm.html" in "/"
    (use hyphens for the page name i.e. "Needleman-Wunsch.html" or "Agglomerative-Clustering.html")
-   create a new script "example_algorithm.js" in "js/"
    (use underscores for the script name and lowercase letters)
-   go into "index.html" and enter your algorithm name (without underscores) in the algorithm menu area:

    <div id="algorithm_menu">
        <ul>
            <li><a href="#">Example-Algorithm</a></li>
        </ul>
    </div>

/*-----------------------------------------------------INHERITANCE----------------------------------------------------*/
In this project simulated inheritance is used because ECMAScript 5
does not support real, Java-like inheritance like in ECMAScript 6 (from 2015).
But, ECMAScript 5 is also fully supported by older browsers and tablet/smartphone browsers.

Examples for simulated inheritance you find under "js/procedures/bases/":
from the classes there, other classes like Needleman-Wunsch inherit functions
by exchanging instances.
So, the parent class is executing the child-class by an instance
and the child class can exchange information with the parent class
through an instance of the parent class (with the related functions: getInput, getOutput, setIO).

/*--------------------------------------------------INPUT PROCESSING--------------------------------------------------*/
The "InputProcessor"-class under "js/post_processing/input_processor.js"
is used for forwarding input to the algorithms and output to the visualizer.
It removes wrong inputs before data is passed to an algorithm
and it defines the behaviour of different input types.

Example:    Result selection
-   after a click on a result,
    the input processor is selecting a table entry of the results table (id: "results")
-   the table and the entry row number is then sent from this class to the visualizer
    which highlights the entry in the table

/*-----------------------------------------------------INTERFACES-----------------------------------------------------*/
Interfaces under "js/interfaces/" are used to work with the input and output (the interface) of the algorithm.
Such interfaces are bound to specific algorithms.
The "InputProcessor" is executing interfaces to get an algorithm output.

/*--------------------------------------------------------HINTS-------------------------------------------------------*/
# HINT 1
Some tables maybe consist out of three separate tables
(instead of one table with header and footer),
because of a Firefox-Bug (bottom right edge disappears).

see: https://github.com/PerfectionMaschine/RNA-Playground/issues/1
<<  Das war wohl ein Firefox-Bug.
    In der Developer-Version (also, vermutlich auch in der nächsten Firefox-Version)
    tritt er nicht auf und auch in keinem anderen Browser.
    Habe jetzt die Tabellen umgeschrieben und den Footer
    von der eigentlichen Tabelle separiert. Jetzt tritt das Problem nimmer auf. >>

translated:
<<  That was probably a Firefox bug.
    In the Developer version (so, probably in the next Firefox version),
    it does not appear and also in no other browser.
    I have now rewritten the tables and separated the footer
    from the actual table. Now, the problem never occurs. >>

/*------------------------------------------------------KEYWORDS------------------------------------------------------*/
 /*---- "this"-keyword ----*/
The "this"-keyword is not used in this project to reference the method-object,
instead an instance on the object is created in the constructor.
The idea behind this? In Javascript "this" points on the
invoker-object which has called the method. So, with the strict-mode it is undefined
in methods which are called by other methods.
This can lead to bugs and so the "this"-keyword
is only used in constructors of this project and for invokers which are not the class-object.

/*--- "use strict"; ----*/
To write good Javascript-code the strict-mode has to be used.
The strict-command has to be written
into the second line of every Javascript-file after the comment-header.

/*------------------------------------------------------LOADING-------------------------------------------------------*/
The loading of an algorithm-page into a webpage is done in the "loader.js"-file.
In this file is a function "updateDocumentView(algorithm, view)" which loads the page into a "div".
If you want add a new algorithm, you have probably to extend this file
and add the Javascript algorithm name into the structure ALGORITHMS in the file "defaults.js".

Hint: Do also have a look into EXTENDING VISUALIZATION!

/*----------------------------------------------------OFFLINE USAGE---------------------------------------------------*/
MathJax
1. 	Download MathJax v2.7.1 source code from https://github.com/mathjax/MathJax/releases
2. 	Extract it somewhere. 
3. 	Rename the extracted folder to "MathJax".
4. 	Copy extracted folder into the downloaded RNA-Playground folder
5. 	Replace your MathJax-Import 
	<script type="text/javascript" async 
		src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-MML-AM_CHTML">
	</script>
	
	with local
	
	<script type="text/javascript" src="MathJax/MathJax.js?config=TeX-MML-AM_CHTML"></script>

Corporate-Design
1. 	Download http://rna.informatik.uni-freiburg.de/resources/corporate-design.css
2. 	Copy the file to the other css-files of the downloaded RNA-Playground folder 
3. 	Replace your CSS-Import
	<link rel="stylesheet" type="text/css" href="http://rna.informatik.uni-freiburg.de/resources/corporate-design.css"/> 
	
	with local
	
	<link rel="stylesheet" type="text/css" href="css/corporate-design.css"/>

/*-----------------------------------------------ORDER IN CONSTRUCTORS------------------------------------------------*/
1.  variables:
    this.inputUpdatesStarted = false;

2.  bindings:
    ko.bindingHandlers.drawChar = {
        ...
    };

3.  public methods (linking):
    this.activateInputUpdates = activateInputUpdates;
    this.inputUpdatesActivated = inputUpdatesActivated;
    this.linkElements = linkElements;
    this.updateGUI = updateGUI;

/*---------------------------------------------ORDER OF PROGRAM ELEMENTS----------------------------------------------*/
1.  namespace with its static methods:
    namespace("needlemanWunsch", startNeedlemanWunsch);

2.  class instances:
    var needlemanWunschInstance;

3.  shared variables
    var inputData = {};
    var outputData = {};

4.  start-function:
    startNeedlemanWunsch() {...}

5.  imports:
    imports()

6.1 Input constructor (with its used methods below):
    function InputViewmodel() {...}

    HINT:   Methods are defined directly below in the order
            in which they are called in the parent method.
            With this design you can read a class like a book or an article
            (speeds up workflow).

            Example:
            function parent() {
                child1();
                child2();
            }

            function child1() {
                someMethod();
                ...
            }

            function someMethod() {
                ...
            }

            function child2() {
                ...
            }

6.2 Algorithm constructor (with its methods):
    NeedlemanWunsch() {...}

6.3 Output constructor (with its methods):
    function OutputViewmodel(outputData) {...}

/*-------------------------------------------------------PATHS--------------------------------------------------------*/
Many Javascript-paths are stored under "js/defaults.js".

Hint:	If you want change file-paths, then you have also to change paths in the "alignment.html" of your webpage.

/*--------------------------------------------------PUBLIC METHODS----------------------------------------------------*/
Namespaces are storing methods, which are visible out of the namespace:
    namespace("bases.multiSequenceAlignment", MultiSequenceAlignment, getAffineSumOfPairsScore);

and constructors storing methods which are visible out of the namespace over an instance of the class:
    function Interface() {
        ...

        // public class methods
        this.sharedInterfaceOperations = sharedInterfaceOperations;
        this.startProcessing = startProcessing;
        this.roundValues = roundValues;
        this.getDistanceTables = getDistanceTables;
        this.getDistanceTable = getDistanceTable;
    }

Hint: You have not to put the public class methods in the namespace!

/*-------------------------------------------------SHARED VARIABLES---------------------------------------------------*/
Variables which are stored outside a constructor.

/*-------------------------------------------------------TESTS--------------------------------------------------------*/
Unit-Tests are stored under "js/tests/". The file "UnitTests.jstd" is specifying Unit-Test file paths.
The tool/plugin which was used to create the Unit-Tests is called "JsTestDriver".

/*-------------------------------------------------UPLOAD ON SERVER---------------------------------------------------*/
Look into alignment.html for information.

/*----------------------------------------------------VARIABLES-------------------------------------------------------*/
Variables have to be always stored inside the constructor
to avoid variable conflicts with other classes in the same namespace.
An exception are "shared variables" which are used between the "classes".

/*----------------------------------------------------VIDEO LINK------------------------------------------------------*/
Presentation (without licensed pictures) explaining the architecture, 
objected oriented programming within the project and the Unit-Tests:
https://youtu.be/mRXmV0J1nAA (available after presentation as not listed video)
