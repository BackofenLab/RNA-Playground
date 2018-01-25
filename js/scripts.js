/**
 * Created by moni on 30.05.16.
 */


// Save the recursion formulas in a Global list called formula_array.

function rerendermath(idval) {
    MathJax.Hub.Queue(["Typeset", MathJax.Hub]);
}

//function showRecursion() {
//    console.log("val:", $("#formula").val());
//    setTimeout(function () {      //introduces a little delay, so that the thing will close slowly and neatly
//        document.getElementById("recursion").innerHTML =  availableAlgorithms[$("#formula").val()].getRecursionInLatex();
//        rerendermath();
//        $(".animate1").empty();
//    }, 450);
//
//}

function getFormula_structCount(){
    document.getElementById("recursion").innerHTML = availableAlgorithms[$(rec_select).text()].Tables[$(rec_id).text()].getRecursionInLatex();//.getRecursionInLatex();
    //document.getElementById("recursion_qb").innerHTML = availableAlgorithms[$(rec_select_qb).text()].getRecursionInLatex();

    rerendermath();
}

function validate(evt) {
    var theEvent = evt || window.event;
    var key = theEvent.keyCode || theEvent.which;
    var buKey = key;
    //console.log(key);

    key = String.fromCharCode( key );
    var regex = /[^gcauGCAU]|\./;
    if( regex.test(key) && buKey != 8 && buKey != 13 && buKey != 37 && buKey != 39 && buKey != 9) {
        //console.log(regex.test(key));
        theEvent.returnValue = false;
        if(theEvent.preventDefault) theEvent.preventDefault();
        return;
    }
}

/*
#######################--- Methods for Downloading Tables ---############################
 */
var csv = "";

function generate_tables() {
    var file = new File([csv], "table.csv", {type: "text/plain;charset=utf-8"});
    saveAs(file);
}

// Convert Matrix to CSV
function matrixToCSV(sequence, matrices) {

    var sequence_column = true;
    var res = "";
    //console.log((matrices));
    for (var idx in matrices)
    {
        if (res.length != 0) {
            res += "\n\n";
        }
        var matrix = matrices[idx];

        if (sequence_column)
            res += ",";
        // add first row: sequence
        for (var j in sequence) {
            res += "," + sequence[j];
        }
        res += "\n";
        for (var i = 1, mat_cells_len = matrix.cells.length; i < mat_cells_len; ++i) {
            if (sequence_column && i > 0)
                res += sequence[i - 1];

            // adding row
            for (var j in matrix.cells[i]) {
                if (sequence_column || j > 0)
                    res += ",";
                var mat_cell_value = matrix.cells[i][j].value;
                if (!isNaN(mat_cell_value) && mat_cell_value != null)
                    res += mat_cell_value;
            }
            res += "\n";
        }
    }
    //console.log(csv);
    csv = res;
    return res;
}

/*
 #######################--- END Download Tables ---############################
 */



/*
 #######################--- BEGIN interaction visualization ---############################
 */

function allor3dots(str) {
    return ( 10 < str.length ? ( str.substr( 0, 3 ) + "..." + str.substr( str.length - 4 ) ) : str );
}

function repres(res) {
    ret = "";
    for (var i = 0, res_len = res.length; i < res_len; ++i) {
        ret += "\t";
        for (var j in res[i]) {
            ret += res[i][j];
        }
        ret += "\n";
        //ret += "<br>"
    }
    return ret;
}


/**
 * 4d Visualization
 * ps is a set of pairs of indices from str1 and str2 respectively.
 * The indiceies represent alignment between between both characters.
 * Indicies are 1-based.
 *
 * if (a1, b1) and (a2, b2) are allignments (in ps)
 *  then a1 < a2 if and only if b2 < b1, to ensure it's a valid alignment.
 * test using console.log(repres(visualize4d("axxbxxcxxdxxexxfxxgxx", "cyydyyhyyjyykyylyymyyoyy", [[16, 16], [17, 11], [19, 7]])));
 */

function visualize4d(str1, str2, ps) {
/*
    return JSON.stringify(ps) + "\n" +
        JSON.stringify(str1) + "\n" +
        JSON.stringify(str2) + "\n" +
    "axx...xexx        x     xx" + "\n" +
    "          f    x   g" + "\n" +
    "          |    |   |" + "\n" +
    "          l    y   h" + "\n" +
    "yyoyymyy   yyky jyy   yydyyc";

*/
    if (ps.length == 0) {
        return [""];
    }
    var marked1 = new Array(str1.length + 1);
    var marked2 = new Array(str2.length + 1);

    var st1i = 0, st1f = 0, st2i = 0, st2f = 0;
    for (var i = 0; i < str1.length; ++i) {
        marked1[i] = false;
    }
    for (var i = 0; i < str2.length; ++i) {
        marked2[i] = false;
    }

    for (var i in ps) {
        marked1[ps[i][0]] = true;
        marked2[ps[i][1]] = true;
        if (st1i == 0 || ps[i][0] < st1i) {
            st1i = ps[i][0];
        }
        if (st1f == 0 || ps[i][0] > st1f) {
            st1f = ps[i][0];
        }
        if (st2i == 0 || ps[i][1] < st2i) {
            st2i = ps[i][1];
        }
        if (st2f == 0 || ps[i][1] > st2f) {
            st2f = ps[i][1];
        }
    }

    var gapBegin = Math.max(Math.min(10, st1i - 1), Math.min(10, str2.length - st2f));
    var gapEnd = Math.max(Math.min(10, str1.length - st1f), Math.min(10, st2i - 1));
    var portion = Math.max(0, Math.max(st1f - st1i, st2f - st2i) + 1);
//    var portion = Math.max(3, Math.max(st1f - st1i, st2f - st2i) + 1);
    var portionOffset = 0;

    var res = new Array(6);
    for (var i = 0; i < res.length; ++i) {
        res[i] = new String();
        for (var j = 0; j < portion + gapBegin + gapEnd + portionOffset; ++j) {
            //res[i] += '.';
            res[i][j] = ' ';
            //res[i][j] = '_';//'&nbsp;';
        }
    };
    // Fill portion
    var st_idx = gapBegin + portionOffset;
    for (var pt1 = st1i, pt2 = st2f, i = 0; (pt1 <= st1f || pt2 >= st2i) && i < portion; ++i) {

        if (marked1[pt1] && marked2[pt2]) {
            res[2][st_idx + i] = str1[pt1 - 1];
            res[4][st_idx + i] = str2[pt2 - 1];
            res[3][st_idx + i] = '|';
            ++pt1;
            --pt2;
        } else {
            if (!marked2[pt2]) {
                res[5][st_idx + i] = str2[pt2 - 1];
                --pt2;
            }
            if (!marked1[pt1]) {
                res[1][st_idx + i] = str1[pt1 - 1];
                pt1++;
            }
        }
        /*
         res[1][st_idx + i] = '_';
         res[5][st_idx + i] = '_';
         */
    }


    // Fill GapBegin Str1

    var gapB1 = allor3dots(str1.substr(0, st1i - 1)); //[0, st1i)
    var st_idx = gapBegin - gapB1.length;
    for (var i = 0; i < gapB1.length; ++i) {
        res[1][st_idx + i] = gapB1[i];
    }
    // Fill GapBegin Str2
    var st_idx = gapBegin - 1;

    var gapB2 = allor3dots(str2.substr(st2f)); // (st2f, str2.length)
    for (var i = 0; i < gapB2.length; ++i) {
        res[5][st_idx - i] = gapB2[i];
    }



    // Fill GapEnd Str1
    var st_idx = gapBegin + portion + portionOffset * 2;

    var gapE1 = allor3dots(str1.substr(st1f)); //(st1f, str1.length)
    for (var i = 0; i < gapE1.length; ++i) {
        res[1][st_idx + i] = gapE1[i];
    }
    // Fill GapEnd Str2

    var gapE2 = allor3dots(str2.substr(0, st2i - 1)); // [0, st2f)
    var st_idx = gapBegin + portion + portionOffset * 2 + gapE2.length - 1;

    for (var i = 0; i < gapE2.length; ++i) {
        res[5][st_idx - i] = gapE2[i];
    }


    return res;
}


/*
 #######################--- END interaction visualization ---############################
 */

/*
 #######################--- Methods for Generating DotPlots ---############################
 */

function mouseover(p) {
    var par = "#" + this.parentElement.parentElement.parentElement.parentElement.id; // get svg name
    var rowText = par + ' .rowtext';
    var colText = par + ' .coltext';
    var rect    = par + ' .row .dp_cell rect';
    var text    = par + ' .row .dp_cell text';

    d3.selectAll(rowText).classed("active", function(d, i) { return i == p.y; });
    d3.selectAll(colText).classed("active", function(d, i) { return i == p.x; });
    d3.selectAll(rect).classed("active", function(d, i) { return d.y == p.y && d.x == p.x; });
    d3.selectAll(text).attr("display",function(d, i) { if (d.y == p.y && d.x == p.x) return "true"; else return "none" });
}

function mouseout() {
    var par = "#" + this.parentElement.parentElement.parentElement.parentElement.id;
    var rect    = par + ' .row .dp_cell rect';

    d3.selectAll("text").classed("active", false);
    d3.selectAll(rect).classed("active", false);
    d3.selectAll(".row .dp_cell text").attr("display","none");
}

function showtext() {
    var par = "#" + this.parentElement.parentElement.parentElement.parentElement.id;
    var text    = par + ' .row .dp_cell text';
    d3.selectAll(text)
        .attr("display","true")
}

function hidetext() {
    var par = "#" + this.parentElement.parentElement.parentElement.parentElement.id;
    var text    = par + ' .row .dp_cell text';
    d3.selectAll(text)
        .attr("display","none")
}

function dotplot(sequence, table, pname, use_log_value) {
    //console.log('hi', sequence, table, pname);
    var maindic = {};

    var bpp = [];
    var mlp = [];
    var keys=["source","target","value"];

    var maxV = 0.001;

    if (use_log_value) {
        for (var i = 1; i < table.length; i++) {
            for (var j = 1; j < table[i].length; j++) {
                maxV = Math.max(maxV, table[i][j].logValue);
            }
        }
    }
    for(var i=1; i<table.length; i++){
        for(var j=1; j<table[i].length; j++){
            var a = table[i][j].i;
            var b = table[i][j].j;
            //var c = Math.pow(parseFloat(table[i][j].value), 0.5);
            var c = parseFloat(table[i][j].value);
            //console.log(table[i][j].value, table[i][j].logValue);
            if (use_log_value) {
                c = parseFloat(table[i][j].logValue);
            }
            c = parseFloat(c.toFixed(3));
            var roww = {};
            roww[keys[0]]=a;
            roww[keys[1]]=b;
            roww[keys[2]]=Math.max(c, 0.000001);
            bpp.push(roww);
            mlp.push(roww);
        }
    }

    maindic['sequence'] = sequence;
    maindic["base-pairing-probabilities"] = bpp;
    maindic["optimal-structure"] = mlp;

    var bpm=maindic;
    var cell_length = 50;
    var margin = {top: 20, right: 80, bottom: 10, left: 80},
        width = cell_length * sequence.length,
        height = cell_length * sequence.length;
    var x = d3.scale.ordinal().rangeBands([0, width]),
        z = d3.scale.linear().domain([0, 1]).clamp(true),
        c = d3.scale.category10().domain(d3.range(10));

    var force = d3.layout.force()
        .charge(-120)
        .linkDistance(30)
        .size([width, height]);

    //var svg = d3.select("#output").append("svg")
    var dev = document.createElement("div");
//    d3.select(dev).style("overflow-x", "scroll");
    var svg = d3.select(dev).append("svg")
        .attr("id", pname)
        .style("width", width + margin.left + margin.right + "px")
        .style("height", height + margin.top + margin.bottom + "px")
        .style("margin-left", -margin.left/2 + "px")
        .append("g")
        .attr("fill", "black")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var matrix = [];

    var seq = [bpm.sequence],
        isequence = bpm.sequence.split(""),
        n = isequence.length;

    var seq = bpm["sequence"];
    var isequence = bpm["sequence"];
    var n = seq.length;
    isequence=seq.split("");

    x.domain(d3.range(n));
    isequence.forEach(function(base, i) {
        base.index = i;
        matrix[i] = d3.range(n).map(function(j) { return {x: j, y: i, z: 0}; });
    });
    bpm["base-pairing-probabilities"].forEach(function(link) {
        //matrix[link.source-1][link.target-1].z += link.value;
        matrix[link.source-1][link.target-1].z = link.value;
        //console.log(link.source, link.target, link.value);
        //matrix[link.target-1][link.source-1].z += link.value;
    });

    //bpm["optimal-structure"].forEach(function(link) {
    //    matrix[link.target-1][link.source-1].z = link.value;
    //});


    var column = svg.selectAll(".column")
        .data(matrix)
        .enter().append("g")
        .attr("class", "column")
        .attr("transform", function(d, i) { return "translate(" + x(i) + ")rotate(-90)"; });

    column.append("line")
        .attr("x1", -width);

    column.append("text")
        .attr("class","coltext")
        .attr("x", 6)
        .attr("y", x.rangeBand() / 2)
        .attr("dy", ".32em")
        .attr("text-anchor", "start")
        .style("font-size", "140%")
        .text(function(d, i) { return isequence[i]; });

    var row = svg.selectAll(".row")
        .data(matrix)
        .enter().append("g")
        .attr("class", "row")
        .attr("transform", function(d, i) { return "translate(0," + x(i) + ")"; })
        .each(row);

    row.append("line")
        .attr("x2", width);

    row.append("text")
        .attr("class","rowtext")
        .attr("x", -6)
        .attr("y", x.rangeBand() / 2)
        .attr("dy", ".32em")
        .attr("text-anchor", "end")
        .style("font-size", "140%")
        .text(function(d, i) { return isequence[i]; });

    function row(row) {
        var cell = d3.select(this).selectAll(".cell")
            .data(row.filter(function(d) { return d.z; }))
            .enter().append("g")
            .attr("class", "dp_cell");

        cell.append("rect")
	        .attr("x", function(d) { return x(d.x) ; })
	        //.attr("cy", 26)
	        .attr("width", x.rangeBand())
	        .attr("height",x.rangeBand())
	        //.attr("r", function(d) { return x.rangeBand() * Math.pow(z(d.z), 0.2)/2;})
	        //.attr("ry", function(d) { return x.rangeBand() * Math.pow(z(d.z), 0.5);})
	        .attr("border", "1px solid")
	        .style("stroke", "lightgrey")
	        .style("stroke-width", 1)
	        .style("fill", "white" )
	        .style("fill-opacity", 1 )
	        //.style("fill-opacity", function(d) { return Math.pow(z(d.z), 0.5); })
	        .on("mouseover", mouseover)
	        .on("mouseout", mouseout);
        cell.append("circle")
            .attr("cx", function(d) { return x(d.x) + 26 ; })
            .attr("cy", 26)
            //.attr("width", function(d) {return x.rangeBand() * Math.pow(z(d.z), 0.5);})
            //.attr("height", function(d) {return x.rangeBand() * Math.pow(z(d.z), 0.5);})
            .attr("r", function(d) { return use_log_value ? (0.9*x.rangeBand() * Math.pow(z(d.z / (maxV + 0.01)), 0.5)/2) : (0.9*x.rangeBand() * Math.pow(z(d.z), 0.5)/2);})
            //.attr("ry", function(d) { return x.rangeBand() * Math.pow(z(d.z), 0.5);})
            //.attr("border", "1px solid red")
            //.style("stroke", "rgb(0, 165, 255)")
            .style("stroke-width", 2)
            //.style("fill-opacity", function(d) { return z(d.z); })
            .style("fill-opacity", function(d) { return Math.pow(z(d.z), 0.5); })
            .on("mouseover", mouseover)
            .on("mouseout", mouseout);
        cell.append("text")
            .attr("x", function(d) { return x(d.x)+10; })
            .attr("y", x.rangeBand()/4)
            .attr("dy", ".32em")
            .attr("text-anchor", "start")
            .attr("display", "none")
            .attr("fill", "red")
            .text(function(d) { return "i="+(d.y+1); })
            .append("svg:tspan")
            .attr("x", function(d) { return x(d.x)+10; })
            .attr("y", x.rangeBand()/4*2)
            .attr("dy", ".32em")
            .text(function(d) { return "j="+(d.x+1); })
            .append("svg:tspan")
            .attr("x", function(d) { return x(d.x)+10; })
            .attr("y", x.rangeBand()/4*3)
            .attr("dy", ".32em")
            .text(function(d) { return (d.z).toFixed(4); })
            
    }
    return dev;
}



function getHybridSequence( seq1, seq2, minLoopLength ) {
	// add first sequence
	var hybrid = seq1;
	// add minLoopLength+1 spacers
	for (var i=minLoopLength; i>-1; i--) {
		hybrid += 'X';
	}
	// add second sequence
	hybrid += seq2;
	return hybrid;
}


/*
 #######################--- END DotPlots ---############################
 */
