/**
 * @Act as the interface for the visualization library, index.html and nussinovmatrix class. Makes it easier and simple to use
 visualization functions.
 */

"use strict";

function NussinovMatrixViewModel() {
    var self = this;
    var cellWidth = 32;
    var cellHeight = 24;
    var ctx;
    var colors = ['lightseagreen', 'lightslategrey', 'lightsalmon', 'lightcoral',
                  'lightsteelblue', 'lightseagreen', 'lightslategrey', 'lightsalmon', 'lightcoral'];
    //var colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'grey', 'red', 'blue'];
    var colors2 = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'grey', 'red', 'blue'];
    var color = 0;

    self.currCell = {
        i: null,
        j: null
    };
    self.currTrace = 0;

    self.input = {
        sequence: ko.observable("GCACGACG"),
        loopLength: ko.observable(0),
        delta: ko.observable(0),
        recursion: ko.observable("nussinovUnique"),
        allowTraceback: true
    };

    self.mckaskillRecursion = ko.computed(function(){
        //console.log($(rec_select).text());
        if($(rec_select).text()=="mcKaskill") {
            console.log("mcKaskill");
            self.input.recursion("mcKaskill");
            self.input.allowTraceback = false;
        }
    });

    self.countingRecursion = ko.computed(function(){
        //console.log($(rec_select).text());
        if($(rec_select).text()=="mycounting") {
            console.log("mycounting");
            self.input.recursion("mycounting");
            self.input.allowTraceback = false;
        }
    });

    self.seqList = ko.computed(function(){
        return self.input.sequence().toUpperCase().split("");
    }, this);


    self.formula = ko.computed(function(){
        return availableAlgorithms[self.input.recursion()];
    }, this);

    self.renderer = function(matrix){
        //var res = JSON.parse(JSON.stringify(matrix));
        console.log(matrix);
        for (var i = 0; i < matrix.cells.length; ++i) {
            var p = "";
            for (var j = 0; j < matrix.cells[i].length; ++j) {
                matrix.cells[i][j].value = parseInt(matrix.cells[i][j].value);
                p += matrix.cells[i][j].value + " ";
            }
            //console.log(p);
        }
        return matrix;
    };

    self.matrix = ko.computed(function(){
        var seq = self.input.sequence().toUpperCase();
        var ll = parseInt(self.input.loopLength());
        console.log("seq len:", seq.length);
        if(seq.length==0){
            $("#matrix").hide();
            return false;
        }
        $("#matrix").show();
         if (self.input.allowTraceback) {
            $('th.cell_th, td.cell').css({
                'width': cellWidth + 'px',
                'height': cellHeight + 'px',
                'padding': '0px 0px'
            });
             ctx = $('#CanvasLayer')[0].getContext("2d");
             ctx.clearRect(0, 0, $('#CanvasLayer')[0].width, $('#CanvasLayer')[0].height);
             //ctx.stroke();
        }

        var tables = self.formula().computeMatrix(seq, ll);
        for(var i = 0; i < tables.length; ++i){
            tables[i] = self.renderer(tables[i]);
        }
        return tables;
    }, this);



    self.tracebacks = ko.computed(function(){
        var del = parseInt(self.input.delta());
        console.log("in traceback:", self.input.allowTraceback);
        if (self.input.allowTraceback) {// exclusive case for nussinov recursions
            console.log("tb allowed");
            return wuchty_2nd(self.matrix()[0], del, self.formula().Tables[0]);
        }
        return false;
    }, this);

    self.cells = ko.computed(function() {
        var tables = [];

        for(var i in self.matrix()){
            if(self.matrix()[i].cells == undefined)return;
            tables.push(self.matrix()[i].cells.slice(1));        // slice is hack to skip first row(investigate later)
        }
        return tables;

    }, this);

    // functions for visualization and interaction
    self.clickStructure = function(clicked_cell, dom) {
        var offset =  $("#matrix_body").position();
        color +=1;
        if(color >= colors.length-1) color = 0;
        console.log(clicked_cell, dom);
        $('td#structTableCells').css({'background': '#FFF'});
        $(dom.target).css({'background': colors[color]});

        var cell = JSON.stringify(clicked_cell);
        cell = JSON.parse(cell);

        $('td.cell').css("background", "white");
        drawFullTrace(offset , cell);
    };

    function drawFullTrace(location, cell){
            //console.log('loading canvas');
            // add one to cell dims because of borders
            var cH = cellHeight+1;
            var cW = cellWidth+1;

            var seqLen = self.seqList().length;
            var top = location.top - (2 * cH);
            var left = location.left + cW;
            var width = (seqLen + 1) * (cW);
            var height = (seqLen + 2) * (cH);

            $('#CanvasLayer')[0].width = width;
            $('#CanvasLayer')[0].height = height;
            $('#CanvasLayer').css({'top': top, 'left': left, 'width': width, 'height': height});
            ctx.textBaseline = "middle";
            ctx.textAlign = "center";
            ctx.font = " 18px sans-serif";
            ctx.fillStyle = "#2a6ebb";

            $(info).text("");
            $(info).hide();
            for(var i=1; i<=seqLen;i++) { // show base pairs formed at top of matrix
                ctx.fillText(cell.structure[i-1], (i) * cW + cW / 2, cH / 2 - 2);
            }

            for(var t in cell.traces){
                var child = cell.traces[t][0];
                var parents = cell.traces[t][1];

                for(var p in parents){
                    if((child[0])<parents[p][0] && (child[1]-1)==parents[p][1]){
                        ctx.font = " 10px sans-serif";
                        ctx.fillText('+1', child[1] * cW + (3*cW/4+2), (child[0]+1) * cH + (cH/4+2));
                    }
                }

                self.currCell.i = child[0];
                self.currCell.j = child[1];
                color +=1;
                if(color >= colors.length-1)color = 0;

                drawArrows(parents, cW, cH);
            }
    }

    self.clickCell = function(clicked_cell) {
        if(clicked_cell.i > clicked_cell.j) {
            return;
        }

        var offset =  $("#matrix_body").position();
        var cell = JSON.stringify(clicked_cell);
        cell = JSON.parse(cell);

        $(structTableCells).css({'background': '#FFF'});
        $('td.cell').css("background", "white");
        drawSingleTraceStep(offset , cell);
    };

    function drawSingleTraceStep(location, cell) {
        //console.log('loading canvas');
        // add one to cell dims because of borders
        var cH = cellHeight+1;
        var cW = cellWidth+1;

        var seqLen = self.seqList().length;
        var top = location.top - (2 * cH);
        var left = location.left + cW;
        var width = (seqLen + 1) * (cW);
        var height = (seqLen + 2) * (cH);

        $('#CanvasLayer')[0].width = width;
        $('#CanvasLayer')[0].height = height;
        $('#CanvasLayer').css({'top': top, 'left': left, 'width': width, 'height': height});
        ctx.textBaseline = "middle";
        ctx.textAlign = "center";
        ctx.font = " 18px sans-serif";
        ctx.fillStyle = "#2a6ebb";

        // cycle through parents pairs
        var cellI = cell.i;
        var cellJ = cell.j;
        if(self.currCell.i == cellI && self.currCell.j == cellJ ){
            if(self.currTrace == cell.traces.length-1){
                self.currTrace = 0;
            }
            else{
                self.currTrace += 1;
            }
        }
        else{
            self.currCell.i = cellI;
            self.currCell.j = cellJ;
            self.currTrace = 0;
        }

        $(info).text("Cell:[" + [cellI, cellJ] + "]  Ancestors:" + JSON.stringify(cell.traces[self.currTrace].parents));
        $(info).show();
        // show base pairs formed at top of matrix
        ctx.fillStyle = "#2a6ebb";
        for(var i=1; i<=seqLen;i++) {
                if (cell.traces[self.currTrace].bps.length != 0) {
                    ctx.font = " 10px sans-serif";
                    ctx.fillText('+1', cellJ * cW + (3*cW/4+2), (cellI+1) * cH + (cH/4+2));
                    ctx.font = " 18px sans-serif";
                    //console.log('bps', cell.traces[self.currTrace].bps[0][0], i);
                    if (i == cell.traces[self.currTrace].bps[0][0]) {

                        ctx.fillText('(', (i) * cW + cW / 2, cH / 2);
                    }
                    else if (i == cell.traces[self.currTrace].bps[0][1]) {
                        ctx.fillText(')', (i) * cW + cW / 2, cH / 2);
                    }
                }
        }
        color = 0;

        // draw arrows
        var ancestors = cell.traces[self.currTrace].parents;
        drawArrows(ancestors, cW, cH);

    }

    // http://deepliquid.com/blog/archives/98
    function drawArrows(ancestors, cW, cH){

        var left_currCell = [self.currCell.j * cW, (self.currCell.i + 1) * cH + cH/2];
        var bottom_currCell = [self.currCell.j * cW + cW/2, (self.currCell.i + 2) * cH];
        var corner_currCell = [self.currCell.j * cW, (self.currCell.i + 2) * cH];

        for (var a in ancestors){
            $('#matrix_body tr').eq(self.currCell.i-1).find('td').eq(self.currCell.j).css('background', colors[color]);
            $('#matrix_body tr').eq(ancestors[a][0]-1).find('td').eq(ancestors[a][1]).css('background', colors[color+1]);
            //console.log(self.currCell.i, self.currCell.j, ancestors[a]);
            var top = [ancestors[a][1] * cW + cW/2, (ancestors[a][0] + 1) * cH];
            var right = [(ancestors[a][1] + 1) * cW, (ancestors[a][0] + 1) * cH + cH/2];
            var corner = [(ancestors[a][1] + 1) * cW, (ancestors[a][0] + 1) * cH];

            if(self.currCell.i == ancestors[a][0]){
                drawLineArrow(left_currCell[0]+1, left_currCell[1]+0.5, right[0], right[1]+0.5);
            }
            else if(self.currCell.j == ancestors[a][1]){
                drawLineArrow(bottom_currCell[0]+0.5, bottom_currCell[1], top[0]+0.5, top[1]+1);
            }
            else if(self.currCell.i+1 < ancestors[a][0]){
                drawLineArrow(bottom_currCell[0]+0.5, bottom_currCell[1], top[0]+0.5, top[1]+1);
            }
            else{
                drawLineArrow(corner_currCell[0], corner_currCell[1], corner[0]-1, corner[1]+1);
            }
        }
    }

    function drawLineArrow(x1,y1,x2,y2) {
        var arrow = [
            [ 4, 0 ],
            [ -8, -5 ],
            [ -8, 5]
        ];

        ctx.lineWidth=2;
        ctx.strokeStyle = colors2[color];

        ctx.beginPath();
        ctx.moveTo(x1,y1);
        ctx.lineTo(x2,y2);
        ctx.stroke();
        var ang = Math.atan2(y2-y1,x2-x1);
        drawFilledPolygon(translateShape(rotateShape(arrow,ang),x2,y2));
    };

    function drawFilledPolygon(shape) {
        ctx.fillStyle = colors2[color+1];
        ctx.beginPath();
        ctx.moveTo(shape[0][0],shape[0][1]);

        for(var p in shape)
            if (p > 0) ctx.lineTo(shape[p][0],shape[p][1]);

        ctx.lineTo(shape[0][0],shape[0][1]);
        ctx.fill();
    };

    function translateShape(shape,x,y) {
        var rv = [];
        for(var p in shape)
            rv.push([ shape[p][0] + x, shape[p][1] + y ]);
        return rv;
    };

    function rotateShape(shape,ang) {
        var rv = [];
        for(var p in shape)
            rv.push(rotatePoint(ang,shape[p][0],shape[p][1]));
        return rv;
    };

    function rotatePoint(ang,x,y) {
        return [
            (x * Math.cos(ang)) - (y * Math.sin(ang)),
            (x * Math.sin(ang)) + (y * Math.cos(ang))
        ];
    };

    // custom binding to get subscript
    ko.bindingHandlers.writeSeq = {
        update: function(element, valueAccessor) {
            var value = valueAccessor();
            var valueUnwrapped = ko.unwrap(value);
            element.innerHTML = "<b>" + valueUnwrapped[0] + "</b><sub>" + valueUnwrapped[1] + "</sub>";
        }
    };

};
ko.applyBindings(new NussinovMatrixViewModel());


