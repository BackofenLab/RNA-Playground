/**
 * Created by moni on 30.05.16.
 */

// Save the recursion formulas in a Global list called formula_array.

function rerendermath(idval) {
    MathJax.Hub.Queue(["Typeset", MathJax.Hub]);
}

function showRecursion() {
    console.log("val:", $("#formula").val());
    setTimeout(function () {      //introduces a little delay, so that the thing will close slowly and neatly
        document.getElementById("recursion").innerHTML = "<b>" + $("#formula").val() + "</b><br> <br>" + availableAlgorithms[$("#formula").val()].getRecursionInLatex();
        rerendermath();
        $(".animate1").empty();
    }, 450);

}

function getFormula_structCount(){
    document.getElementById("recursion").innerHTML = availableAlgorithms[$(rec_select).text()].getRecursionInLatex();
    //document.getElementById("recursion_qb").innerHTML = availableAlgorithms[$(rec_select_qb).text()].getRecursionInLatex();

    rerendermath();
}

function validate(evt) {
    var theEvent = evt || window.event;
    var key = theEvent.keyCode || theEvent.which;
    var buKey = key;

    key = String.fromCharCode( key );
    var regex = /[^gcauGCAU]|\./;
    if( regex.test(key) && buKey != 8 && buKey != 13) {
        console.log(regex.test(key));
        theEvent.returnValue = false;
        if(theEvent.preventDefault) theEvent.preventDefault();
        return;
    }
}


