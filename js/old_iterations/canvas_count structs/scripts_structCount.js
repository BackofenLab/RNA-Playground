/**
 *
 * Created by moni on 30.05.16.
 */

// Save the recursion formulas in a Global list called formula_array.
var x;
var formula_clicked;
looplength = $("#looplength").val(); //loop

function rerendermath(idval) {
    MathJax.Hub.Queue(["Typeset", MathJax.Hub]);
}

function getFormula_structCount(){
    document.getElementById("formula").innerHTML = NussinovMatrix_structuresCount.getRecursionInLatex();
    formula_clicked = 1;
    rerendermath();
}

function validate(evt) {
    var theEvent = evt || window.event;
    var key = theEvent.keyCode || theEvent.which;
    var buKey = key;

    $("#matrix_out").remove();
    $("#output_title").remove();

    key = String.fromCharCode( key );
    var regex = /[^gcauGCAU]|\./;
    if( regex.test(key) && buKey != 8) {
        console.log(regex.test(key));
        theEvent.returnValue = false;
        if(theEvent.preventDefault) theEvent.preventDefault();
        return;
    }
    console.log(buKey);
    if($("#userInput").val().length > 0 && formula_clicked != undefined) {
        enableGO();
    }

    if($("#userInput").val().length-1  <= 0 && buKey == 8) {
        disableGO();
    }

}

function chkLength(evt){
    if($("#userInput").val().length  <= 0) {
        disableGO();
    }
}

function getGoColor() {
    var colorval = $("#GO").css("background-color");
    var parts = colorval.match(/^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/);
    delete(parts[0]);
    for (var i = 1; i <= 3; ++i) {
        parts[i] = parseInt(parts[i]).toString(16);
        if (parts[i].length == 1) parts[i] = '0' + parts[i];
    }
    color = '#' + parts.join('');
    return color;
}

function buttonCheck() {
    console.log("buttonCheck color is:", getGoColor());
    if($("#output_title").length > 0){
        return;
    }
    else if(getGoColor() == "#d6d6d6"){
        alert("Please make sure that you have entered a RNA sequence and selected a recursion formula.");
    }
    else{
        disableGO();
    }
}


function disableGO(){
    console.log("disable go");
    $("#GO").css('pointer-events', 'none');
    $("#GO").css('background-color', "#D6D6D6");
}

function enableGO(){
    console.log("enable go");
    $("#GO").css('pointer-events', 'auto');
    $("#GO").css('background-color', "grey");
}




