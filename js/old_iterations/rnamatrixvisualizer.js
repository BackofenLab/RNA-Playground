//////////////////////////////////////////////////////////////////////////////////////
// RNA MATRIX VISUALIZER                                                            //
// Author: Muhammad Numair Mansur                                                   //
// Institute: Department of BioInformatics, University of Freiburg, Germany         //
// Last Updated: August 18th, 2015 (11:00 AM)                                       //
// Author's email: numair.mansur@gmail.com                                          //
//////////////////////////////////////////////////////////////////////////////////////


/**
 * @Main file for Visualiztions of RNA-algorithms-JS project.
 This file contain the following main items:
 -Boxes
 -other important stuff
 * @author Muhammad Numair Mansur (numair.mansur@gmail.com).
 */


/**
 * Generate Data Values array.
 * @param {object} Matrix Object
 * @param {string} Entered RNA sequence
 * @returns {array} An Array of data values to be put in the visual matrix
 */

function datagenerator(matrix_object, sequence) {
    var data = [];
    var length = sequence.length + 1;
    for (var i = 0; i < length; i++) {
        for (var j = 0; j < length; j++) {
            //get the matrxi value for row i and column j
            data.push(availableAlgorithms[formula_clicked].computeMatrix(sequence, looplength).getValue(i, j));
        }
    }
    //set the first row values of the matrix to the sequence values
    data[0] = ' ';
    for (var i = 0; i < sequence.length; i++) {
        data[i + 1] = sequence[i];
    }
    return data;
}

/**
 * Generate and return one box object.
 * @param {int} X-coordinate of the box
 * @param {int} Y-coordinate of the box
 * @param {string} String value to be displayed inside the box
 * @param {string} Color string. Either color name or the hexadecimal color code.
 * @param {int} Width of the box.
 * @param {int} height of the box.
 * @returns {object} A box array with the properties set by the user.
 */

function box(xcord, ycord, value, color, width, height) {

    var boxx = new joint.shapes.basic.Rect({
        position: {x: xcord, y: ycord},
        size: {width: width, height: height},
        attrs: {rect: {fill: color}, text: {text: value, fill: 'black'}}
    });
    return boxx;
}

/**
 * Generate an array of matrix box objects
 * @param {arrat} array of integer matrix values
 * @param {int} Length of the input RNA sequence
 * @returns {array} An Array of matrix box objects
 */

function matrix_box_combine(data, n, x_location, y_location) {
    var counter = 0;
    var boxgen = [];
    var xcord = x_location + 40;
    var ycord = y_location;
    //generate the boxes and put values in them
    for (var i = 0; i <= n; i++) {
        for (var j = 0; j <= n; j++) {
            boxgen[counter] = box(xcord, ycord, data[counter], 'white', 40, 30);
            counter = counter + 1;
            xcord = xcord + 40;
        }
        ycord = ycord + 30;
        xcord = x_location + 40;
    }
    for (var i = 0; i < n + 1; i++) {
        boxgen[i].attr('rect/fill', '#f2f2f0'); //setting the color of the first line of boxes to grey
    }
    return boxgen;
}

/**
 * Generates the side boxes in the matrix visual, which contains the sequence.
 * @param {arrat} array of integer matrix values
 * @param {int} Length of the input RNA sequence
 * @returns {array} An Array of matrix box objects
 */

function sideboxgenerator(sequence, x_location, y_location) {
    var side_boxes = [];
    var ycord = y_location;
    var xcord = x_location;
    side_boxes[0] = box(xcord, ycord, ' ', '#f2f2f0', 40, 30);
    ycord = ycord + 30;
    for (var i = 1; i <= sequence.length; i++) {
        var text = sequence[i - 1];
        side_boxes[i] = box(0, ycord, sequence[i - 1], '#f2f2f0', 40, 30);
        ycord = ycord + 30;
    }
    return side_boxes;
}


/**
 * Generates the invisible base pair boxes above the matrix.
 * @param {string} RNA Sequence
 * @param {int} X-coordinate of the top left corner of the matrix
 * @param {int} Y-coordinate of the top left corner of the matrix
 * @returns {array} An Array of base pair box objects.
 */

function basepairboxesgenerator(sequence, x_location, y_location) {
    var basepairboxes = [];
    var ycord_base = y_location;
    var xcord_base = x_location + 80; //80 because it has to be 2 boxes ahead on the x-axis
    for (var i = 0; i < sequence.length; i++) {
        basepairboxes[i] = new joint.shapes.basic.Rect({
            position: {x: xcord_base, y: ycord_base},
            size: {width: 40, height: 30},
            attrs: {rect: {stroke: 'none', 'fill-opacity': 0}, text: {text: ' ', fill: 'red', 'font-size': 25}}
        });
        xcord_base = xcord_base + 40;
    }
    return basepairboxes;
}


/**
 * Visalize the matrix object on the screen
 * @param {string} Entered RNA Sequence
 * @param {object} Matrix object
 * @param {object} graph object
 * @param {object} paper object
 * @param {int} X-Coordinate of the top left corner of the matrix
 * @param {int} Y-Coordinate of the top left corner of the matrix
 * @returns {visual_matrix} Generate a matrix on the screen
 */

function visualize_matrix(matrix_object, sequence, graph, paper, x_location, y_location) {
    //Generate an array of data values to be put into the boxes
    var data = datagenerator(matrix_object, sequence);
    //Generate an array of matrix box objects
    var boxes = matrix_box_combine(data, sequence.length, x_location, y_location + 30); //+30 because of base pair boxes
    var sideboxes = sideboxgenerator(sequence, x_location, y_location + 30); //+30 because of base pair boxes
    var headerbasepairboxes = basepairboxesgenerator(sequence, x_location, y_location);
    graph.addCells(boxes);
    graph.addCells(sideboxes);
    graph.addCells(headerbasepairboxes);
    rerendermath();
    return [boxes, headerbasepairboxes];
}


/**
 * Visualize the information box on the top.
 * @param {object} JOINTJS graph object
 * @param {object} JOINTJS paper object
 * @param {int} x-coordinate of the information box
 * @param {int} y-coordinate of the information box
 * @param {int} width
 * @param {int} height
 * @returns {object} box object of the information box
 */

function visualize_information_box(graph, paper, x_cord, y_cord, width, height, color) {
    var informationbox = [];
    informationbox[0] = box(x_cord, y_cord, 'Multi-click on a cell for its traces', color, width, height)

    graph.addCells(informationbox);
    return informationbox;
}


/**
 * Determines where you clicked on the work area.
 * @param {id} box_clicked_id
 * @param {array} boxes
 * @param {id} single_optimal_traceback_button_id
 * @param {id} all_optimal_traceback_button_id
 * @returns {int} An interger value depending on where you clicked in the work area
 */

function type_of_box_clicked_function(box_clicked_id, boxes, single_optimal_traceback_button_id, all_optimal_traceback_button_id) {
    if (box_clicked_id == single_optimal_traceback_button_id) {
        return 2;
    }
    else if (box_clicked_id == all_optimal_traceback_button_id) {
        return 3;
    }
    else {
        for (var i = 0; i < boxes.length; i++) {
            if (box_clicked_id == boxes[i].id) {
                return 1;
            }
        }
    }
}

/**
 * Every box is an object with its own id. This function takes a box id and return the index number of the box in the boxes array.
 * @param {id} id of a box
 * @param {array} boxes
 * @returns {int} box number in the boxes array
 */
function boxnumberfromid(id, boxes) {
    for (var i = 0; i <= boxes.length; i++) {
        if (id == boxes[i].id) {
            return i;
        }
    }
}

/**
 * Takes as input the box number and returns the box id.
 * @param {int} box number
 * @param {array} boxes
 * @returns {id} ID of the box number
 */
var idfromboxnumber = function (n, boxes) {
    for (i = 0; i <= boxes.length; i++) {
        if (n == i) {
            return boxes[i].id;
        }
    }
}

/**
 * Takes the box number clicked and returns the coordinates of the box in the form of a string in brackets.
 * @param {int} boxnumber
 * @param {int} sequence length
 * @returns {string} Coordinates in the form of a string.
 */
function coordinates_from_box_number(boxnumber, sequence_length) {
    var counter = 0;
    for (var i = 0; i <= sequence_length; i++) {
        for (var j = 0; j <= sequence_length; j++) {
            if (counter == boxnumber) {
                var string = "(" + i + "," + j + ")";
                return (string);
            }
            counter = counter + 1;
        }
    }
}

/**
 * Takes the box number clicked and returns the coordinates of the box in the form of a string in brackets.
 * @param {int} boxnumber
 * @param {int} sequence length
 * @returns {string} Coordinates in the form of a string.
 */
function coordinate_values_from_box_numbers(boxnumber, sequence_length) {
    var counter = 0;
    var coordinate_array = [];
    for (var i = 0; i <= sequence_length; i++) {
        for (var j = 0; j <= sequence_length; j++) {
            if (counter == boxnumber) {
                coordinate_array.push(i);
                coordinate_array.push(j);
                return coordinate_array;
            }
            counter = counter + 1;
        }
    }
}

/**
 * Determines which link should be generated next based on the number of parents, and the number of times the user have clicked on the box,
 * @param {array} givemelinksforarray
 * @param {id} boxclicked
 * @param {id} last_clicked_box
 * @param {int} numberofways
 * @returns {array} givemelinksforarray
 */
function givemelinksfor(givemelinksforarray, boxclicked, last_clicked_box, numberofways) {
    if (boxclicked == last_clicked_box) {
        if (givemelinksforarray[1] < numberofways - 1) {
            givemelinksforarray[1]++;
        }
        else if (numberofways > -1) {
            givemelinksforarray[1] = 0;
        }
        else {
            givemelinksforarray[1] = -1;
        }
    }
    else {
        givemelinksforarray[0] = boxclicked;
        if (numberofways >= 0) {
            givemelinksforarray[1] = 0;
        }
        else {
            givemelinksforarray[1] = -1;
        }
    }
    return givemelinksforarray;
}

/**
 * Takes coordinates and returns the box number
 * @param {int} x-Coordinate
 * @param {int} y-Coordinate
 * @param {int} Length
 * @returns {int} box number
 */
function boxnumberfromcoordinates(x, y, length) // length is the length of the total data values
{
    var counter = 0;
    for (var i = 0; i < Math.sqrt(length); i++) {
        for (var j = 0; j < Math.sqrt(length); j++) {

            if (i == x && j == y) {
                return counter
            }
            counter = counter + 1;
        }
    }

}

/**
 * Plus one box generator
 * @param {array} boxes
 * @param {int} boxclicked
 * @param {int} box_clicked_value
 * @param {array} parent_boxes_values
 * @param {int} graph
 * @returns {array} An array containg information plus one box object and a flag weather it is being displayed or not.
 */
var plusone = function (boxes, boxclicked, box_clicked_value, parent_boxes_values, graph) {
    var parents_sum = 0;
    var xpos = parseInt(boxes[boxclicked].position().x, 10);
    var ypos = parseInt(boxes[boxclicked].position().y, 10);
    for (var i = 0; i < parent_boxes_values.length; i++) {
        parents_sum = parents_sum + parent_boxes_values[i];
    }

    if (parents_sum != box_clicked_value) {
        var plusonebox = [];
        plusonebox[0] = box(xpos + 25, ypos - 20, "+1", 'yellow', 40, 30);
        graph.addCells(plusonebox);
        return [plusonebox[0], 1];
    }
    else {
        return [0, 0];
    }

}


/**
 * Refresh the work area.
 */
function refresher(graph, boxes, sequence, live_links, plusonebox, headerbasepairboxes, information_box,
                   single_optimal_traceback_button, all_optimal_traceback_button, all_optimal_traceback_select_box, live_links_sel) {	//Refesh all the boxes to
    //console.log("inside refresher function", live_links_sel);
    for (var i = 0; i < boxes.length; i++) {
        boxes[i].attr('rect/fill', 'white');
    }
    //First row should be gray
    for (var i = 0; i < sequence.length + 1; i++) {
        boxes[i].attr('rect/fill', '#f2f2f0');
    }
    for (var i = 0; i < live_links.length; i++) // LINK REFRESHER
    {
        graph.getCell(live_links[i]).remove();

    }

    for (var i = 0; i < live_links_sel.length; i++) // LINK REFRESHER
    {
        //console.log("Before splice", live_links_sel);
        graph.getCell(live_links_sel[i]).remove();
        live_links_sel.splice(i, 1);
        //console.log("after splice", live_links_sel);
        i = i - 1;

    }

    if (plusonebox[1] == 1)      // PLUS ONE BOX REFRESHER
    {
        if (graph.getCell(plusonebox[0]) != undefined) {
            graph.getCell(plusonebox[0]).remove();
        }
    }
    for (i = 0; i < headerbasepairboxes.length; i++) {
        headerbasepairboxes[i].attr('text/text', ' ');
    }
    information_box[0].attr('text/text', 'Multi-click on a cell for its traces');
    information_box[0].attr('rect/fill', '#D8D8D8');
    information_box[0].attr('text/fill', 'black');
    single_optimal_traceback_button.attr('rect/fill', '#D8D8D8');
    single_optimal_traceback_button.attr('text/fill', 'black');
    all_optimal_traceback_button.attr('rect/fill', '#D8D8D8');
    all_optimal_traceback_button.attr('text/fill', 'black');
    // REMOVE THE ALL TRACE BACK SELECT BOX
    if (graph.getCell(all_optimal_traceback_select_box[0]) != undefined) {
        graph.getCell(all_optimal_traceback_select_box[0]).remove();
        $("#selecterbox").remove();
    }

}

/**
 * Plus one box generator
 * @param {object} graph
 * @param {object} paper
 * @param {int} x_cord
 * @param {int} y_cord
 * @param {string} text
 * @param {int} width
 * @param {int} height
 * @param {hex} color
 * @returns {object} Trace back information box
 */
function visualize_traceback_info_box(graph, paper, x_cord, y_cord, text, width, height, color) {
    var tracebackinfobox = [];
    tracebackinfobox[0] = box(x_cord, y_cord, text, color, width, height)
    graph.addCells(tracebackinfobox);
    return tracebackinfobox;
}

/**
 * Function that gets called when you select an option from the drop down menu below the optimal tracebacks.
 * @param {object} Select object
 */

function selectorfunction(a) {
    //console.log("selector input:", a)

    $("traces").on('click', 'li', function(e) {
        $(this).parent().find('li.active').removeClass('active');
        $(this).addClass('active');
    });

    for (var i = 0; i < live_links.length; i++) // LINK REFRESHER
    {
        graph.getCell(live_links[i]).remove();

    }

    for (var i = 0; i < boxes.length; i++) {
        boxes[i].attr('rect/fill', 'white');
    }


    for (var i = 0; i < live_links_sel.length; i++) // LINK REFRESHER
    {
        //console.log("Before splice", live_links_sel);
        graph.getCell(live_links_sel[i]).remove();
        live_links_sel.splice(i, 1);
        //console.log("after splice", live_links_sel);
        i = i - 1;

    }

    //console.log("You clicked in the selector thingie");
    //console.log("YOU CLICKED VALUE NUMBER " + a);
    var base_pair_string = allOpt[a].structure;
    //console.log("YOU CLICKED VALUE NUMBER " + a);
    var optionTexts = [];
    $("#traces li").each(function() { optionTexts.push($(this).text()) });

    for(var t=0; t<optionTexts.length; t++){
        lid = "#list" + t;
        if(t==a) $(lid).css("background-color", '#567DEA');
        else $(lid).css("background-color", '#D8D8D8');
    }

    var selected_value_index = a;

    //console.log("Base pair string value = " + base_pair_string[1]);
    for (var i = 0; i < base_pair_string.length; i++) {
        matrix_boxes[1][i].attr('text/text', base_pair_string[i]);
    }


    var boxes_sel = matrix_boxes[0];
    var singe_traceback_trace_information_sel = allOpt[selected_value_index].traces;
    //console.log("SINGLE TRACE ", JSON.stringify(singe_traceback_trace_information_sel));
    var colorarray = ['#a292bc', '#a868c0', '#40b868', '#f0d048', '#3088f0', '#f05868', '#a27dfa', '#259073', '#3cbeb7', '#f5ac78', '#9db7f5', '#f05868',
        '#a27dfa', '#259073', '#3cbeb7', '#f5ac78', '#9db7f5', '#a292bc', '#a868c0', '#40b868', '#f0d048', '#3088f0', '#f05868', '#a27dfa', '#259073', '#3cbeb7',
        '#f5ac78', '#9db7f5', '#f05868', '#a27dfa', '#259073', '#3cbeb7', '#f5ac78', '#9db7f5', '#a292bc', '#a868c0', '#40b868', '#f0d048', '#3088f0', '#f05868', '#a27dfa',
        '#259073', '#3cbeb7', '#f5ac78', '#9db7f5', '#f05868', '#a27dfa', '#259073', '#3cbeb7', '#f5ac78', '#9db7f5'];
    var fromcolorindex = 0;
    var minimum = 10000;
    for (var k = 0; k < singe_traceback_trace_information_sel.length; k++) {
        var from_box_sel = boxnumberfromcoordinates(singe_traceback_trace_information_sel[k][0][0], singe_traceback_trace_information_sel[k][0][1], boxes_sel.length);
        if (from_box_sel < minimum) {
            minimum = from_box_sel;
        }

        for (var l = 0; l < singe_traceback_trace_information_sel[k][1].length; l++) {
            var to_box_sel = boxnumberfromcoordinates(singe_traceback_trace_information_sel[k][1][l][0], singe_traceback_trace_information_sel[k][1][l][1], boxes_sel.length);
            var link_sel = new joint.dia.Link({
                source: {id: idfromboxnumber(from_box_sel, boxes_sel)},
                target: {id: idfromboxnumber(to_box_sel, boxes_sel)}
            });
            link_sel.attr({
                '.connection': {stroke: 'red', 'stroke-width': 3},
                '.marker-target': {fill: 'yellow', d: 'M 10 0 L 0 5 L 10 10 z'}
            });
            live_links_sel.push(link_sel);
            if(boxes_sel[to_box_sel] != undefined)
                boxes_sel[to_box_sel].attr('rect/fill', colorarray[fromcolorindex]);

        }

        fromcolorindex++;
    }
    graph.addCells(live_links_sel);
    if(boxes[minimum] != undefined)
        boxes_sel[minimum].attr('rect/fill', 'red');


}

/**
 * Listner function
 *
 *
 */

function listen_to_click_events(paper, graph, matrix_object, boxes, information_box, headerbasepairboxes, traceback_info_box, single_optimal_traceback_button,
                                all_optimal_traceback_button, opt, live_links_sel) {

    live_links = [];
    var single_traceback_structure_string = opt.structure;
    var singe_traceback_trace_information = opt.traces;
    var all_optimal_traceback_select_box = [];
    var type_of_box_clicked = 0;
    var givemelinksforarray = [9999, -1];
    var last_clicked_box = 9999;
    var plusonebox = [0, 0];
    var temp_error = 0; //TEMPORARY ERROR MESSAGE;

    paper.on('cell:pointerup',
        function (cellView, evt, x, y) {
            refresher(graph, boxes, sequence, live_links, plusonebox, headerbasepairboxes, information_box,
                single_optimal_traceback_button, all_optimal_traceback_button, all_optimal_traceback_select_box, live_links_sel);
            ////////////TEMPORARY ERROR MESSAGE
            if (graph.getCell(temp_error[0]) != undefined) {
                graph.getCell(temp_error[0]).remove();
            } //TEMPORARY
            //////////////////////
            live_links = [];

            //First we have to determine what is clicked on the paper object.
            type_of_box_clicked = type_of_box_clicked_function(cellView.model.id, boxes, single_optimal_traceback_button.id, all_optimal_traceback_button.id);
            //console.log("type of box clicked:", type_of_box_clicked);
            if (type_of_box_clicked == 1 && boxnumberfromid(cellView.model.id, boxes) > sequence.length) ///////////////// CLICKED THE MATRIX
            {
                //Determine what box number was clicked.
                var boxclicked = boxnumberfromid(cellView.model.id, boxes);
                //Determine the coordinates of the box clicked.
                var coordinatesclicked = coordinates_from_box_number(boxclicked, sequence.length);
                // Get the 2 element array of the coordinate values
                var coordinates_clicked_value = coordinate_values_from_box_numbers(boxclicked, sequence.length);
                //Change the color of the cliced box to red
                boxes[boxclicked].attr('rect/fill', 'red');
                // Determine the number of ways the parents can be formed
                var numberofways = 0;
                if (matrix_object.getTraces(coordinates_clicked_value[0], coordinates_clicked_value[1]) != null) {
                    numberofways = matrix_object.getTraces(coordinates_clicked_value[0], coordinates_clicked_value[1]).length;
                }
                if (numberofways == undefined) {
                    numberofways = -1;
                }
                //Get a 2 element array that have the box number and the number of way that should be displayed.
                givemelinksforarray = givemelinksfor(givemelinksforarray, boxclicked, last_clicked_box, numberofways);
                last_clicked_box = boxclicked;
                //console.log("Give me links for box number and way number "+givemelinksforarray);
                if (numberofways != 0) //more then one ways
                {
                    var number_of_parents = 0;
                    var base_pairs_array = [];
                    number_of_parents = matrix_object.getTraces(coordinates_clicked_value[0], coordinates_clicked_value[1])[givemelinksforarray[1]].parents.length;
                    base_pairs_array = matrix_object.getTraces(coordinates_clicked_value[0], coordinates_clicked_value[1])[givemelinksforarray[1]].bps;
                    var parents_boxes_values = [];
                    for (var i = 0; i < number_of_parents; i++) {
                        var point_to_boxnumber = boxnumberfromcoordinates(matrix_object.getTraces(coordinates_clicked_value[0], coordinates_clicked_value[1])[givemelinksforarray[1]].parents[i][0], matrix_object.getTraces(coordinates_clicked_value[0], coordinates_clicked_value[1])[givemelinksforarray[1]].parents[i][1], boxes.length);
                        boxes[point_to_boxnumber].attr('rect/fill', 'yellow');
                        var link = new joint.dia.Link({
                            source: {id: cellView.model.id},
                            target: {id: boxes[point_to_boxnumber].id}
                        });
                        link.attr({
                            '.connection': {stroke: 'red', 'stroke-width': 3},
                            '.marker-target': {fill: 'yellow', d: 'M 10 0 L 0 5 L 10 10 z'}
                        });
                        live_links.push(link);
                        var parent_value = matrix_object.getValue(matrix_object.getTraces(coordinates_clicked_value[0], coordinates_clicked_value[1])[givemelinksforarray[1]].parents[i][0], matrix_object.getTraces(coordinates_clicked_value[0], coordinates_clicked_value[1])[givemelinksforarray[1]].parents[i][1]);
                        parents_boxes_values.push(parent_value);
                    }
                    graph.addCells(live_links);
                    plusonebox = plusone(boxes, boxclicked, matrix_object.getValue(coordinates_clicked_value[0], coordinates_clicked_value[1]), parents_boxes_values, graph);
                    var parent_strings = " ";  // INFORMATION BOX

                    for (var iparents = 0; iparents < number_of_parents; iparents++) {
                        parent_strings = parent_strings + "(" + matrix_object.getTraces(coordinates_clicked_value[0], coordinates_clicked_value[1])[givemelinksforarray[1]].parents[iparents] + ") "
                    }
                    information_box[0].attr('text/text', 'Cell (' + coordinates_clicked_value + ") :  Recursion to " + parent_strings); // INFORMATION BOX  							}
                    information_box[0].attr('rect/fill', '#7F7F7F');
                    information_box[0].attr('text/fill', 'white');
                    //BASE PAIR INFORMATION
                    if (base_pairs_array[0] != undefined && base_pairs_array != undefined) {
                        headerbasepairboxes[base_pairs_array[0][0] - 1].attr('text/text', '(');
                        headerbasepairboxes[base_pairs_array[0][1] - 1].attr('text/text', ')');

                    }


                }
                else // No parents
                {
                    information_box[0].attr('text/text', 'Cell (' + coordinates_clicked_value + ")"); // INFORMATION BOX
                    information_box[0].attr('rect/fill', '#7F7F7F');
                    information_box[0].attr('text/fill', 'white');
                }


            }
            else if (type_of_box_clicked == 2) //clicked the Display Single Traceback Button
            {
                //Change the color of the traceback button and the color of the text to white.
                single_optimal_traceback_button.attr('rect/fill', '#7F7F7F');
                single_optimal_traceback_button.attr('text/fill', 'white');
                //put the structure on the headerbase pair boxes.
                //Update header base pair information
                for (i = 0; i < headerbasepairboxes.length; i++) {
                    headerbasepairboxes[i].attr('text/text', single_traceback_structure_string[i]);
                }
                //console.log("SINGLE TRACE __________________", JSON.stringify(singe_traceback_trace_information));
                //Now display the traceback on the matrix
                var colorarray = ['#a292bc', '#a868c0', '#40b868', '#f0d048', '#3088f0', '#f05868', '#a27dfa', '#259073', '#3cbeb7', '#f5ac78', '#9db7f5', '#f05868', '#a27dfa', '#259073', '#3cbeb7', '#f5ac78', '#9db7f5', '#a292bc', '#a868c0', '#40b868', '#f0d048', '#3088f0', '#f05868', '#a27dfa', '#259073', '#3cbeb7', '#f5ac78', '#9db7f5', '#f05868', '#a27dfa', '#259073', '#3cbeb7', '#f5ac78', '#9db7f5', '#a292bc', '#a868c0', '#40b868', '#f0d048', '#3088f0', '#f05868', '#a27dfa', '#259073', '#3cbeb7', '#f5ac78', '#9db7f5', '#f05868', '#a27dfa', '#259073', '#3cbeb7', '#f5ac78', '#9db7f5'];
                var fromcolorindex = 0;
                var minimum = 10000;
                for (var k = 0; k < singe_traceback_trace_information.length; k++) {
                    var from_box = boxnumberfromcoordinates(singe_traceback_trace_information[k][0][0], singe_traceback_trace_information[k][0][1], boxes.length);
                    if (from_box < minimum) {
                        minimum = from_box;
                    }
                    for (var l = 0; l < singe_traceback_trace_information[k][1].length; l++) {
                        var to_box = boxnumberfromcoordinates(singe_traceback_trace_information[k][1][l][0], singe_traceback_trace_information[k][1][l][1], boxes.length);
                        var link = new joint.dia.Link({
                            source: {id: idfromboxnumber(from_box, boxes)},
                            target: {id: idfromboxnumber(to_box, boxes)}
                        });
                        link.attr({
                            '.connection': {stroke: 'red', 'stroke-width': 3},
                            '.marker-target': {fill: 'yellow', d: 'M 10 0 L 0 5 L 10 10 z'}
                        });
                        live_links.push(link);
                        boxes[to_box].attr('rect/fill', colorarray[fromcolorindex]);

                    }
                    fromcolorindex++;
                }
                graph.addCells(live_links);
                if(boxes[minimum] != undefined)
                    boxes[minimum].attr('rect/fill', 'red');

            }

            else if (type_of_box_clicked == 3) //Click the Display all optimal Traceback button
            {
                all_optimal_traceback_button.attr('rect/fill', '#7F7F7F');
                all_optimal_traceback_button.attr('text/fill', 'white');

                // Create a HTML box here that will have the drop down containing the list of all the optimal tracebacks.
                // ------------------------------------------------------------------------------------------------------

                var struct = "";
                for (var i in allOpt) {
                    struct += allOpt[i].structure + ",";
                }

                var supposition = struct.split(",");
                supposition.splice(-1, 1);
                //console.log("SUPPOSDITION", supposition);

                joint.shapes.html = {};
                // joint.shapes.html.Element = joint.shapes.basic.Rect.extend({
                //                                                                 defaults: joint.util.deepSupplement({
                //                                                                 type: 'html.Element',
                //                                                                 attrs: {
                //                                                                 	position: { x: 400, y: 250 },
                //                                                                         rect: { fill: 'red' }, text: { text:'Select box' , fill: 'black' }
                //                                                                         }
                //                                                                 }, joint.shapes.basic.Rect.prototype.defaults)
                //                                                             });
                joint.shapes.html.Element = joint.shapes.basic.Rect.extend({
                    defaults: joint.util.deepSupplement({
                        type: 'html.Element',
                        attrs: {
                            rect: {stroke: 'none', 'fill-opacity': 0}
                            //text: { text:'Select a Structure to see its Traceback' , fill: 'black' }
                        }
                    }, joint.shapes.basic.Rect.prototype.defaults)
                });

                var allopt_string = "hello";
                var selectortext = '';
                //if (supposition.length >= 10){
                //    supposition.length = 10;
                //}
                for (var i = 0; i < supposition.length; i++) {
                    sele = "selectorfunction(" + i + ")";
                    ids = "list" + i
                //    selectortext = selectortext + '<option value="' + supposition[i] + '">' + supposition[i] + '</option>'
                    selectortext = selectortext + '<li id=' + ids + ' style="letter-spacing:3px ;border-bottom:5px solid #F2F3F1;cursor:pointer" onmousedown=' + sele + ' value="' + supposition[i] + '">' + supposition[i] + '</li>'
                }

                joint.shapes.html.ElementView = joint.dia.ElementView.extend({


                template: [
                        '<div id="selecterbox" class="html-element">',
                        //'<select  style="" onChange="selectorfunction(this)" name="Please Select"><option selected disabled>Select a Structure</option>' + selectortext + '</select>',
                        '<ul class="element" id="traces">' + selectortext + '</ul>',
                        '</div>'

                    ].join(''),

                    initialize: function () {

                        _.bindAll(this, 'updateBox');
                        joint.dia.ElementView.prototype.initialize.apply(this, arguments);

                        this.$box = $(_.template(this.template)());
                        // Prevent paper from handling pointerdown.
                        this.$box.find('input,select').on('mousedown click', function (evt) {
                            evt.stopPropagation();
                        });
                        // This is an example of reacting on the input change and storing the input data in the cell model.
                        this.$box.find('input').on('change', _.bind(function (evt) {
                            this.model.set('input', $(evt.target).val());
                        }, this));
                        this.$box.find('select').on('change', _.bind(function (evt) {
                            this.model.set('select', $(evt.target).val());
                        }, this));
                        this.$box.find('select').val(this.model.get('select'));
                        this.$box.find('.delete').on('click', _.bind(this.model.remove, this.model));
                        // Update the box position whenever the underlying model changes.
                        this.model.on('change', this.updateBox, this);
                        // Remove the box when the model gets removed from the graph.
                        this.model.on('remove', this.removeBox, this);

                        this.updateBox();
                    },
                    render: function () {
                        joint.dia.ElementView.prototype.render.apply(this, arguments);
                        this.paper.$el.prepend(this.$box);

                        this.updateBox();
                        return this;
                    },
                    updateBox: function () {
                        // Set the position and dimension of the box so that it covers the JointJS element.
                        var bbox = this.model.getBBox();
                        // Example of updating the HTML with a data stored in the cell model.
                        this.$box.find('label').text(this.model.get('label'));
                        this.$box.find('span').text(this.model.get('select'));
                        this.$box.css({
                            width: bbox.width,
                            height: bbox.height,
                            //left: bbox.x,
                            //top: bbox.y,
                            transform: 'rotate(' + (this.model.get('angle') || 0) + 'deg)'
                        });
                    },
                    removeBox: function (evt) {
                        this.$box.remove();
                    }


                    ///////
                });


                all_optimal_traceback_select_box[0] = new joint.shapes.html.Element({
                    position: {x: 0, y: 0},
                    size: {width: 200, height: 20},
                    label: 'Select a Structure'
                });
                // all_optimal_traceback_select_box[1] = new joint.shapes.html.Element({ position: { x: 0, y: 0 }, size: { width: 170, height: 100 }, label: 'Select a Struccture', select: 'one' });
                graph.addCells(all_optimal_traceback_select_box);


            }


        }
    )
}
