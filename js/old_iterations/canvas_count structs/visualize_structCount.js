/**
 * @Act as the interface for the visualization library, index.html and nussinovmatrix class. Makes it easier and simple to use
 visualization functions.
 */

// each cell is of size 40*30
var cellWidth = 40;
var cellHeight = 30;
var seqLen = null;
var canvasContext = null;
var cH = null; //canvasHeight
var cW = null; //canvasWidth


function loadCanvas(id) {
    console.log("seqLen", seqLen);
    var canvas = document.createElement('canvas');
    canvas.id = "MatrixLayer";
    canvas.width = (seqLen + 2) * 40;
    canvas.height = (seqLen + 1) * 30;
    canvas.style.zIndex = 8;
    canvas.style.border = "1px solid";

    canvasContext = canvas.getContext("2d");
    cW = canvas.width;
    cH = canvas.height;

    document.getElementById(id).appendChild(canvas);

    return canvas;
}

function write_seq() {

    var sequence = userInput;
    console.log(sequence);
    canvasContext.textBaseline = "middle";
    canvasContext.textAlign = "center";
    canvasContext.fillStyle = "#000000";
    canvasContext.font = " 18px sans-serif";

    for(var i=0; i<sequence.length; i++){
        canvasContext.fillText(sequence[i], (i+2)*cellWidth+cellWidth/2, cellHeight/2);
        canvasContext.fillText(sequence[i], cellWidth/2, (i+1)*cellHeight+cellHeight/2);
    }

}

function fill_matrix(matrix) {

    var sequence = userInput;
    console.log(sequence);
    canvasContext.textBaseline = "middle";
    canvasContext.textAlign = "center";
    canvasContext.fillStyle = "#000000";
    canvasContext.font = " 18px sans-serif";

    for (var i = 1; i <= sequence.length; i++) {
        for (var j = 1; j <= sequence.length; j++) {
            canvasContext.fillText(matrix.getValue(i, j),
                (j+1) * cellWidth + cellWidth / 2,
                (i) * cellHeight + cellHeight / 2);
        }

    }
}
function drawMatrix(matrix) {

    canvasContext.clearRect(0, 0, cW, cH);

    canvasContext.fillStyle = "#d8d8d8";
    canvasContext.fillRect(0, 0, cW, cellHeight);
    canvasContext.fillRect(0, 0, cellWidth, cH);

    // write sequence on top rox and left column
    write_seq();

    // fillmatrix
    fill_matrix(matrix);

    /* vertical lines */
    for (var x = 0; x <= cW; x += cellWidth) {
        canvasContext.moveTo(0.5 + x, 0);
        canvasContext.lineTo(0.5 + x, cH);
    }
    /* horizontal lines */
    for (var y = 0; y <= cH; y += cellHeight) {
        canvasContext.moveTo(0, 0.5 + y);
        canvasContext.lineTo(cW, 0.5 + y);
    }
    /* draw it! */
    canvasContext.strokeStyle = "#ccc";
    canvasContext.stroke();


}

function visualize() {
    var looplength = $("#looplength").val();
    looplength++;
    looplength--;
    // Compute the matrix object
    var matrix = NussinovMatrix_structuresCount.computeMatrix(userInput, looplength);
    seqLen = userInput.length;
    var canvas = loadCanvas("matrix_out");
    drawMatrix(matrix);
    //console.log("canvas width:", canvas.width);
//  Now Listen to the events.
};

