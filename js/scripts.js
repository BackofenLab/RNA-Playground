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
    console.log(key);

    key = String.fromCharCode( key );
    var regex = /[^gcauGCAU]|\./;
    if( regex.test(key) && buKey != 8 && buKey != 13 && buKey != 37 && buKey != 39 ) {
        console.log(regex.test(key));
        theEvent.returnValue = false;
        if(theEvent.preventDefault) theEvent.preventDefault();
        return;
    }
}

function mouseover(p) {
    var par = "#" + this.parentElement.parentElement.parentElement.parentElement.id; // get svg name
    var rowText = par + ' .rowtext';
    var colText = par + ' .coltext';
    var rect    = par + ' .row .cell rect';
    var text    = par + ' .row .cell text';

    d3.selectAll(rowText).classed("active", function(d, i) { return i == p.y; });
    d3.selectAll(colText).classed("active", function(d, i) { return i == p.x; });
    d3.selectAll(rect).classed("active", function(d, i) { return d.y == p.y && d.x == p.x; });
    d3.selectAll(text).attr("display",function(d, i) { if (d.y == p.y && d.x == p.x) return "true"; else return "none" });
}

function mouseout() {
    var par = "#" + this.parentElement.parentElement.parentElement.parentElement.id;
    var rect    = par + ' .row .cell rect';

    d3.selectAll("text").classed("active", false);
    d3.selectAll(rect).classed("active", false);
    d3.selectAll(".row .cell text").attr("display","none");
}

function showtext() {
    var par = "#" + this.parentElement.parentElement.parentElement.parentElement.id;
    var text    = par + ' .row .cell text';
    d3.selectAll(text)
        .attr("display","true")
}

function hidetext() {
    var par = "#" + this.parentElement.parentElement.parentElement.parentElement.id;
    var text    = par + ' .row .cell text';
    d3.selectAll(text)
        .attr("display","none")
}

function dotplot(sequence, table, pname) {
    console.log('hi');
    var maindic = {};


    var bpp = [];
    var mlp = [];
    var keys=["source","target","value"];

    for(var i=1; i<table.length; i++){
        for(var j=1; j<table[i].length; j++){
            var a = table[i][j].i;
            var b = table[i][j].j;
            //var c = Math.pow(parseFloat(table[i][j].value), 0.5);
            var c = parseFloat(table[i][j].value);
            c = parseFloat(c.toFixed(3));
            var roww = {};
            roww[keys[0]]=a;
            roww[keys[1]]=b;
            roww[keys[2]]=c;
            bpp.push(roww);
            mlp.push(roww);
        }
    }

    maindic['sequence'] = sequence;
    maindic["base-pairing-probabilities"] = bpp;
    maindic["optimal-structure"] = mlp;

    var bpm=maindic;

    var margin = {top: 80, right: 80, bottom: 10, left: 80},
        width = 300,
        height = 300;
    var x = d3.scale.ordinal().rangeBands([0, width]),
        z = d3.scale.linear().domain([0, 1]).clamp(true),
        c = d3.scale.category10().domain(d3.range(10));

    var force = d3.layout.force()
        .charge(-120)
        .linkDistance(30)
        .size([width, height]);

    //var svg = d3.select("#output").append("svg")
    var dev = document.createElement("div");
    var svg = d3.select(dev).append("svg")
        .attr("id", pname)
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
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
        .text(function(d, i) { return isequence[i]; });

    function row(row) {
        var cell = d3.select(this).selectAll(".cell")
            .data(row.filter(function(d) { return d.z; }))
            .enter().append("g")
            .attr("class", "cell");
        cell.append("rect")
            .attr("x", function(d) { return x(d.x); })
            .attr("width", x.rangeBand())
            .attr("height", x.rangeBand())
            //.style("fill-opacity", function(d) { return z(d.z); })
            .style("fill-opacity", function(d) { return Math.pow(z(d.z), 0.3); })
            .on("mouseover", mouseover)
            .on("mouseout", mouseout);
        cell.append("text")
            .attr("x", function(d) { return x(d.x)+15; })
            .attr("y", x.rangeBand()/2)
            .attr("dy", ".32em")
            .attr("text-anchor", "start")
            .attr("display", "none")
            .attr("color", "blue")
            //.text("hello");
            .text(function(d) { return d.z; })
    }
    return dev;
}

// Convert Matrix to CSV
function matrixToCSV(sequence, matrix) {
    var sequence_column = true;
    var res = "";

    if (sequence_column)
        res += ",";
    // add first row: sequence
    for (var j in sequence) {
        res += "," + sequence[j];
    }
    res += "\n";
    for (var i = 1; i < matrix.length; ++i) {
        if (sequence_column && i > 0)
            res += sequence[i - 1];
        
        // adding row
        for (var j in matrix[i]) {
            if (sequence_column || j > 0)
                res += ",";
            if (!isNaN(matrix[i][j].value))
                res += matrix[i][j].value;
        }
        res += "\n";
    }

    return res;
}