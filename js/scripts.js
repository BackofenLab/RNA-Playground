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
        document.getElementById("recursion").innerHTML =  availableAlgorithms[$("#formula").val()].getRecursionInLatex();
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

function parsePSFile(sequence, table) {

    var points = "";

    for (var i in table) {
        for (var j in table[i]) {
            if (points.length > 0) {
                points += "\n";
            }
            points += format("{0} {1} {2} ubox", i, j, table[i][j].value);
        }
    }

    return  "%!PS-Adobe-3.0 EPSF-3.0\n%%Title: RNA Dot Plot\n%%Creator: PS_dot.c,v 1.38 2007/02/02 15:18:13 ivo Exp $, ViennaRNA-2.1.9\n%%CreationDate: Mon Jul  4 15:35:24 2016\n%%BoundingBox: 66 211 518 662\n%%DocumentFonts: Helvetica\n%%Pages: 1\n%%EndComments\n\n%Options: -d2 \n% \n%This file contains the square roots of the base pair probabilities in the form\n% i  j  sqrt(p(i,j)) ubox\n\n%%BeginProlog\n/DPdict 100 dict def\nDPdict begin\n/logscale false def\n/lpmin 1e-05 log def\n\n/box { %size x y box - draws box centered on x,y\n   2 index 0.5 mul sub            % x -= 0.5\n   exch 2 index 0.5 mul sub exch  % y -= 0.5\n   3 -1 roll dup rectfill\n} bind def\n\n/ubox {\n   logscale {\n      log dup add lpmin div 1 exch sub dup 0 lt { pop 0 } if\n   } if\n   3 1 roll\n   exch len exch sub 1 add box\n} bind def\n\n/lbox {\n   3 1 roll\n   len exch sub 1 add box\n} bind def\n\n/drawseq {\n% print sequence along all 4 sides\n[ [0.7 -0.3 0 ]\n  [0.7 0.7 len add 0]\n  [-0.3 len sub -0.4 -90]\n  [-0.3 len sub 0.7 len add -90]\n] {\n   gsave\n    aload pop rotate translate\n    0 1 len 1 sub {\n     dup 0 moveto\n     sequence exch 1 getinterval\n     show\n    } for\n   grestore\n  } forall\n} bind def\n\n/drawgrid{\n  0.01 setlinewidth\n  len log 0.9 sub cvi 10 exch exp  % grid spacing\n  dup 1 gt {\n     dup dup 20 div dup 2 array astore exch 40 div setdash\n  } { [0.3 0.7] 0.1 setdash } ifelse\n  0 exch len {\n     dup dup\n     0 moveto\n     len lineto \n     dup\n     len exch sub 0 exch moveto\n     len exch len exch sub lineto\n     stroke\n  } for\n  [] 0 setdash\n  0.04 setlinewidth \n  currentdict /cutpoint known {\n    cutpoint 1 sub\n    dup dup -1 moveto len 1 add lineto\n    len exch sub dup\n    -1 exch moveto len 1 add exch lineto\n    stroke\n  } if\n  0.5 neg dup translate\n} bind def\n\nend\n%%EndProlog\nDPdict begin\n%delete next line to get rid of title\n270 665 moveto /Helvetica findfont 14 scalefont setfont (dot.ps) show\n\n%GAUUUAGCUGCAUGUGCUAG\\\n/sequence { (\\\n" +
        sequence +
        "\\\n) } def\n/len { sequence length } bind def\n\n72 216 translate\n72 6 mul len 1 add div dup scale\n/Helvetica findfont 0.95 scalefont setfont\n\ndrawseq\n0.5 dup translate\n% draw diagonal\n0.04 setlinewidth\n0 len moveto len 0 lineto stroke \n\n/min { 2 copy gt { exch } if pop } bind def\n\n/utri{ % i j prob utri\n  gsave\n  1 min 2 div\n  0.85 mul 0.15 add 0.95  0.33\n  3 1 roll % prepare hsb color\n  sethsbcolor\n  % now produce the coordinates for lines\n  exch 1 sub dup len exch sub dup 4 -1 roll dup 3 1 roll dup len exch sub\n  moveto lineto lineto closepath fill\n  grestore\n} bind def\n\n%data starts here\n\n%start of quadruplex data\n\n%draw the grid\ndrawgrid\n\n%start of base pair probability data\n% i j pt ubox\n" +
        points +
        "\n\nshowpage\nend\n%%EOF\n";
}