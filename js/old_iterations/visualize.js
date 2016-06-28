/**
 * @Act as the interface for the visualization library, index.html and nussinovmatrix class. Makes it easier and simple to use
 visualization functions.
 */

// Global variables
allOpt = Object.create(Traceback);
graph = new joint.dia.Graph;
matrix_boxes = [];
boxes = [];
looplength = 0;
live_links_sel = [];

function visualize()
{
//GLOBAL Variable
    looplength = $( "#looplength" ).val();
    looplength++;
    looplength--;
    delta =	 $( "#delta" ).val();
    if(delta == undefined)
        delta = 0;

    // Compute the matrix object
    var matrix_object = availableAlgorithms[formula_clicked].computeMatrix(userInput,looplength);
    // Generate Single optimal Traceback information
    var optimalTrace = Object.create(Traceback);
    optimalTrace.init(matrix_object.sequence.length);
    var opt = matrix_object.getOptimalTraceback(1, matrix_object.sequence.length, optimalTrace);

    if(delta>matrix_object.sequence.length)
        delta = matrix_object.sequence.length;

    allOpt = wuchty_2nd(matrix_object, delta, formula_clicked);

    if(opt == undefined){
        var tmp_opt = JSON.stringify(allOpt[0]);
        opt = JSON.parse(tmp_opt);
    }

    live_links_sel = [];

    // TIMER START
    var d = new Date();
    var n = d.getTime();

    // Generate the JOINTJS Graph object
    // JOINTJS Object
    graph = new joint.dia.Graph;

    // Generate the JOINTJS new paper object;
    var paper = new joint.dia.Paper({el: $('#matrix_out'),width: 800,height: 600,model: graph,interactive: false});

    // Function with arguments: visualize_matrix(matrix_object,sequence,graph,paper,x_location,y_location)
    // Defined in the rnamatrixvisualizer.js file
    matrix_boxes = visualize_matrix(matrix_object,userInput,graph,paper,0,49);  //GLOBAL VARIABLE OF MATRIX OBJECTS

    //Visualize the Information box
    // Function with arguments: visualize_information_box(graph,paper,x_cord,y_cord,width,height)
    // Defined in rnamatrixvisualizer.js file
    var endInfobox = 360;
    var widthOptbox = 200;
    var startAllopt = endInfobox + widthOptbox;
    var widthAlloptbox = 240;
    var boxHeight = 45;

    var information_box = visualize_information_box(graph,paper,0,0,endInfobox,boxHeight,'#D8D8D8');

    // TRACEBACK INFORMATION BOX.
    var traceback_info_box = visualize_traceback_info_box(graph,paper,endInfobox,0," ",0,boxHeight,'#D8D8D8');

    // Display the single optimal trace back button (visualize_traceback_info_box(graph, paper, x_cord, y_cord, text, width, height, color))
    var single_optimal_traceback_button = visualize_traceback_info_box(graph,paper,endInfobox,0,"Display an optimal traceback",widthOptbox,boxHeight,'#D8D8D8');

    // Display the all possible traceback options there are
    var all_optimal_traceback_button = visualize_traceback_info_box(graph,paper,startAllopt,0,"Display <= 10 optimal tracebacks",widthAlloptbox,boxHeight,'#D8D8D8');

    boxes = matrix_boxes[0]; //boxes are global

//Now Listen to the events.
// Defined in the rnamatrixvisualizer.js file
    if(allOpt != undefined) {
        console.log("opt is defined");

        listen_to_click_events(paper,
            graph,
            matrix_object,
            matrix_boxes[0],
            information_box,
            matrix_boxes[1],
            traceback_info_box[0],
            single_optimal_traceback_button[0],
            all_optimal_traceback_button[0],
            opt,
            live_links_sel
        );
    }
//getTraces (i,j) returns an array of object
// and .bps and .parents can only be called on objects.
}
