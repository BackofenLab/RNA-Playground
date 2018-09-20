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

    var formulas = [];

    var maxStructures = 15;

    self.currCell = {
        i: null,
        j: null
    };
    self.currTrace = 0;
    self.fired = false;

    self.rawSeq = ko.observable("GGUCCAC");
    self.rawSeq2 = ko.observable("CCGAGG");
    self.loopLength = ko.observable($(rec_select).text() == "hybrid" ? 0 : 1);

    self.rawSequence = ko.computed(function () {
        return self.rawSeq().toUpperCase();
    });
    self.rawSequence2 = ko.computed(function () {
        var res = self.rawSeq2().toUpperCase();
        if ($(rec_select).text() == "rnaup") {
            res = res.split("").reverse().join("");
        }
        return res;
    });

    self.input = {
        loopLength: ko.computed(function () {
            return self.loopLength();
        }),
        delta: ko.observable(0),
        recursion: ko.observable("nussinovUnique"),
        allowTraceback: true,
        energy: ko.observable(-1),
        energy_normal: ko.observable(1),
        gamma: ko.observable(2.0),
        sequence: ko.computed(function () {
            var ll = self.loopLength();
            if (self.rawSeq() == undefined)
                return;
            if ($(rec_select).text() == "coFold" || $(rec_select).text() == "hybrid" || $(rec_select).text() == "rnaup") {
                if (self.rawSeq2() == undefined)
                    return;
                var linker = '';
                for (var i = 0; i <= ll; i++) linker += 'X';
                return self.rawSequence() + linker + self.rawSequence2();
            }

            return self.rawSeq().toUpperCase();
        }),

    };


    self.mcCaskillRecursion = ko.computed(function () {
        //console.log($(rec_select).text());
        if ($(rec_select).text() === "mcCaskill") {
            console.log("mcCaskill");
            self.input.recursion("mcCaskill");
            self.input.allowTraceback = false;
        }
    });

    self.countingRecursion = ko.computed(function () {
        //console.log($(rec_select).text());
        if ($(rec_select).text() == "nussinovCounting") {
            console.log("nussinovCounting");
            self.input.recursion("nussinovCounting");
            self.input.allowTraceback = false;
        }
    });

    self.maxExpAcc = ko.computed(function () {
        //console.log($(rec_select).text());
        if ($(rec_select).text() == "MaxExpAcc") {
            console.log("MaxExpAcc");
            self.input.recursion("MaxExpAcc");
            self.input.allowTraceback = true;
            cellWidth = 48;
            cellHeight = 28;
        }
    });

    self.coFold = ko.computed(function () {
        //console.log($(rec_select).text());
        if ($(rec_select).text() == "coFold") {
            console.log("coFold");
            self.input.recursion("coFold");
            self.input.allowTraceback = true;
            //cellWidth = 48;
            //cellHeight = 28;
        }
    });
    self.hybrid = ko.computed(function () {
        //console.log($(rec_select).text());
        if ($(rec_select).text() == "hybrid") {
            console.log("hybrid");
            self.input.recursion("hybrid");
            self.input.allowTraceback = true;
            //self.input.allowTraceback = false;
            //console.log("setttt");
            //cellWidth = 48;
            //cellHeight = 28;
        }
    });
    self.rnaup = ko.computed(function () {
        //console.log($(rec_select).text());
        if ($(rec_select).text() == "rnaup") {
            console.log("rnaup");
            self.input.recursion("rnaup");
            self.input.allowTraceback = true;
            //self.input.allowTraceback = false;
            //cellWidth = 48;
            //cellHeight = 28;
        }
    });


    self.seqList = ko.computed(function () {
        return self.input.sequence().toUpperCase().split("");
    }, this);

    self.formula = ko.computed(function () {
        //console.log("setting self.formula ", self.input.recursion(), availableAlgorithms);
        return availableAlgorithms[self.input.recursion()];
    }, this);

    self.matrix = ko.computed(function () {
        //var seq = self.input.sequence().toUpperCase();
        //var ll = parseInt(self.input.loopLength());
        //console.log("seq len:", self.input.sequence().length);
        if (self.input.sequence().length == 0) {
            $("#output").hide();

            return false;
        }
        $("#output").show();

        if (self.input.allowTraceback && (self.input.recursion() != "hybrid" && self.input.recursion() != "rnaup")) {

            $('th.cell_th, td.cell').css({
                'width': cellWidth + 'px',
                'height': cellHeight + 'px',
                'padding': '0px 0px'
            });
            ctx = $('#CanvasLayer')[0].getContext("2d");
            ctx.clearRect(0, 0, $('#CanvasLayer')[0].width, $('#CanvasLayer')[0].height);
            //ctx.stroke();
        }

        //console.log($("#rec_select").html());

        //console.log('input:', self.input.sequence(), self.input.loopLength(), self.input.delta(), self.input.recursion());
        //console.log('debugging', self.formula(), self.input);
        var tables = self.formula().computeMatrix(self.input);

        $("#4dVisual").text("");
        //console.log("matrix compute");
        if (self.input.recursion() === "mcCaskill") {
            $("#paired_dotplot").html(dotplot(self.input.sequence(), tables[2].cells, 'pd', false));
            $("#unpaired_dotplot").html(dotplot(self.input.sequence(), tables[3].cells, 'ud', false));

        }
        if (self.input.recursion() === "rnaup") {
            $("#dotplot_seq1").html(dotplot(self.rawSequence(), tables[1].cells, 'up1', true));
            $("#dotplot_seq2").html(dotplot(self.rawSequence2(), tables[2].cells, 'up2', true));

        }

        // add latex formulas to array
        if (!self.fired) {
            for (var i in tables) {
                formulas.push(tables[i].getRecursionInLatex());
            }
        }

        return tables;
    }, this);


    self.tracebacks = ko.computed(function () {
        var del = parseInt(self.input.delta());
        //console.log("in traceback:", self.input.allowTraceback);
        if (self.input.allowTraceback) {// exclusive case for nussinov recursions
            //console.log("tb allowed");

            if ($(rec_select).text() == "hybrid" || $(rec_select).text() == "rnaup") {
                var res = wuchty4d(self.matrix()[0], maxStructures);
                //console.log('wuchty out', res);
                return res;
            }
            if ($(rec_select).text() == "MaxExpAcc") {
                return wuchty_unlimited(self.matrix()[0], del, self.formula().Tables[0], maxStructures);
            }
            //return wuchty_2nd(self.matrix()[0], del, self.formula().Tables[0]);
            return wuchty_2nd_limited(self.matrix()[0], del, self.formula().Tables[0], maxStructures);
        }
        return false;
    }, this);

    self.cells = ko.computed(function () {
        var tables = [];

        for (var r in self.matrix()) {
            //matrixToCSV(self.input.sequence(), self.matrix()[i].cells);
            if (self.matrix()[r].cells == undefined)return;
            var oldMatrix = self.matrix()[r].cells.slice(1);  // slice is hack to skip first row(investigate later)
            var newMatrix = JSON.parse(JSON.stringify(oldMatrix));
            for (var i = 0; i < newMatrix.length; ++i) {
                for (var j = 0; j < newMatrix[i].length; ++j)
                {
                    if (self.input.recursion() === "mcCaskill" || self.input.recursion() === "MaxExpAcc") {
                        if (newMatrix[i][j].value === null || isNaN(newMatrix[i][j].value)) {
                            newMatrix[i][j].value = "";//0.0;
                        } else {
                            newMatrix[i][j].value = parseFloat(newMatrix[i][j].value).toFixed(3);
                        }
                    }
                }
            }
            tables.push(newMatrix);
        }
        //console.log(tables.simpleRepresentation());
        matrixToCSV(self.input.sequence(), self.matrix());
        return tables;

    }, this);

    self.latex = ko.computed(function () {
        //    var formulas = [];
        //    for (var i in self.matrix()) {
        //        formulas.push(self.matrix()[i].getRecursionInLatex());
        //    }
        return formulas;

    }, this);

    // functions for visualization and interaction
    self.clickStructure4d = function (clicked_cell, dom) {
        $(".col_table").scrollLeft(0);
        var offset = $("#matrix_body").position();
        color += 1;
        if (color >= colors.length - 1) color = 0;
        //console.log("CC", clicked_cell);
        $('td#structTableCells').css({'background': '#FFF'});
        $("#4dVisual").text("");
        $("#4dVisual").text(clicked_cell.rep4d);
        if (clicked_cell.value % 1 === 0) {
            $("#output_value").text(clicked_cell.value);
        } else {
            $("#output_value").text(parseFloat(clicked_cell.value).toFixed(5));
        }

        $(dom.target).css({'background': colors[color]});

        var cell = JSON.stringify(clicked_cell);
        cell = JSON.parse(cell);

        $('td.cell').css("background", "white");
        drawFullTrace(offset, cell);
    };

    // functions for visualization and interaction
    self.clickStructure = function (clicked_cell, dom) {
        $(".col_table").scrollLeft(0);
        var offset = $("#matrix_body").position();
        color += 1;
        if (color >= colors.length - 1) color = 0;
        //console.log(clicked_cell, dom);
        $('td#structTableCells').css({'background': '#FFF'});
        $(dom.target).css({'background': colors[color]});

        var cell = JSON.stringify(clicked_cell);
        cell = JSON.parse(cell);

        $('td.cell').css("background", "white");
        $('td.cell_mea').css("background", "white");

        drawFullTrace(offset, cell);
    };

    function drawFullTrace(location, cell) {
        //console.log('loading canvas');
        // add one to cell dims because of borders
        var cH = cellHeight + 1;
        var cW = cellWidth + 1;

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
        for (var i = 1; i <= seqLen; i++) { // show base pairs formed at top of matrix
            ctx.fillText(cell.structure[i - 1], (i) * cW + cW / 2, cH / 2 - 2);
        }

        for (var t in cell.traces) {
            var child = cell.traces[t][0];
            var parents = cell.traces[t][1];

            for (var p in parents) {
                if ((child[0]) < parents[p][0] && (child[1] - 1) == parents[p][1]) {
                    ctx.fillStyle = "#000";
                    ctx.font = "bold 10px arial";
                    ctx.fillText('+1', child[1] * cW + (1 * cW / 4 - 1), (child[0] + 1) * cH + (cH / 4 + 2));
                }
            }

            self.currCell.i = child[0];
            self.currCell.j = child[1];
            color += 1;
            if (color >= colors.length - 1)color = 0;

            drawArrows(parents, cW, cH);
        }
    }

    self.clickCell = function (clicked_cell) {
        if (clicked_cell.i > clicked_cell.j) {
            return;
        }
        $(".col_table").scrollLeft(0);
        var offset = $("#matrix_body").position();
        var cell = JSON.stringify(clicked_cell);
        cell = JSON.parse(cell);

        $("#structTableCells").css({'background': '#FFF'});
        $('td.cell').css("background", "white");
        $('td.cell_mea').css("background", "white");
        drawSingleTraceStep(offset, cell);
    };

    function drawSingleTraceStep(location, cell) {
        //console.log('loading canvas');
        // add one to cell dims because of borders
        var cH = cellHeight + 1;
        var cW = cellWidth + 1;

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
        //ctx.font = "18px sans-serif";
        //ctx.fillStyle = "#2a6ebb";

        // cycle through parents pairs
        var cellI = cell.i;
        var cellJ = cell.j;
        if (self.currCell.i == cellI && self.currCell.j == cellJ) {
            if (self.currTrace == cell.traces.length - 1) {
                self.currTrace = 0;
            }
            else {
                self.currTrace += 1;
            }
        }
        else {
            self.currCell.i = cellI;
            self.currCell.j = cellJ;
            self.currTrace = 0;
        }

        $(info).text("Cell:[" + [cellI, cellJ] + "]  Ancestors:" + JSON.stringify(cell.traces[self.currTrace].parents));
        $(info).show();
        // show base pairs formed at top of matrix

        for (var i = 1; i <= seqLen; i++) {
            if (cell.traces[self.currTrace].bps.length != 0) {
                ctx.fillStyle = "#000";
                ctx.font = "10px arial";
                ctx.fillText('+1', cellJ * cW + (1 * cW / 4 + 2), (cellI + 1) * cH + (cH / 4 + 2));
                ctx.fillStyle = "#2a6ebb";
                ctx.font = "18px sans-serif";
                if (self.rawSeq2 == undefined) {
                    var fillText = ['(', ')'];
                }
                else {
                    if (cell.traces[self.currTrace].bps[0][0] < self.rawSeq().length + 1 &&
                        cell.traces[self.currTrace].bps[0][1] > self.rawSeq().length + 1) {
                        var fillText = ['[', ']'];
                    }
                    else {
                        var fillText = ['(', ')'];
                    }

                }

                if (i == cell.traces[self.currTrace].bps[0][0]) {

                    ctx.fillText(fillText[0], (i) * cW + cW / 2, cH / 2);
                }
                else if (i == cell.traces[self.currTrace].bps[0][1]) {
                    ctx.fillText(fillText[1], (i) * cW + cW / 2, cH / 2);
                }
            }
        }
        color = 0;

        // draw arrows
        var ancestors = cell.traces[self.currTrace].parents;
        drawArrows(ancestors, cW, cH);

    }

    // http://deepliquid.com/blog/archives/98
    function drawArrows(ancestors, cW, cH) {

        var left_currCell = [self.currCell.j * cW, (self.currCell.i + 1) * cH + cH / 2];
        var bottom_currCell = [self.currCell.j * cW + cW / 2, (self.currCell.i + 2) * cH];
        var corner_currCell = [self.currCell.j * cW, (self.currCell.i + 2) * cH];

        for (var a in ancestors) {
            $('#matrix_body tr').eq(self.currCell.i - 1).find('td').eq(self.currCell.j).css('background', colors[color]);
            $('#matrix_body tr').eq(ancestors[a][0] - 1).find('td').eq(ancestors[a][1]).css('background', colors[color + 1]);
            //console.log(self.currCell.i, self.currCell.j, ancestors[a]);
            var top = [ancestors[a][1] * cW + cW / 2, (ancestors[a][0] + 1) * cH];
            var right = [(ancestors[a][1] + 1) * cW, (ancestors[a][0] + 1) * cH + cH / 2];
            var corner = [(ancestors[a][1] + 1) * cW, (ancestors[a][0] + 1) * cH];

            if (self.currCell.i == ancestors[a][0]) {
                drawLineArrow(left_currCell[0] + 1, left_currCell[1] + 0.5, right[0], right[1] + 0.5);
            }
            else if (self.currCell.j == ancestors[a][1]) {
                drawLineArrow(bottom_currCell[0] + 0.5, bottom_currCell[1], top[0] + 0.5, top[1] + 1);
            }
            else if (self.currCell.i + 1 < ancestors[a][0]) {
                drawLineArrow(bottom_currCell[0] + 0.5, bottom_currCell[1], top[0] + 0.5, top[1] + 1);
            }
            else {
                drawLineArrow(corner_currCell[0], corner_currCell[1], corner[0] - 1, corner[1] + 1);
            }
        }
    }

    function drawLineArrow(x1, y1, x2, y2) {
        var arrow = [
            [4, 0],
            [-8, -5],
            [-8, 5]
        ];

        ctx.lineWidth = 2;
        ctx.strokeStyle = colors2[color];

        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();
        var ang = Math.atan2(y2 - y1, x2 - x1);
        drawFilledPolygon(translateShape(rotateShape(arrow, ang), x2, y2));
    };

    function drawFilledPolygon(shape) {
        ctx.fillStyle = colors2[color + 1];
        ctx.beginPath();
        ctx.moveTo(shape[0][0], shape[0][1]);

        for (var p in shape)
            if (p > 0) ctx.lineTo(shape[p][0], shape[p][1]);

        ctx.lineTo(shape[0][0], shape[0][1]);
        ctx.fill();
    };

    function translateShape(shape, x, y) {
        var rv = [];
        for (var p in shape)
            rv.push([shape[p][0] + x, shape[p][1] + y]);
        return rv;
    };

    function rotateShape(shape, ang) {
        var rv = [];
        for (var p in shape)
            rv.push(rotatePoint(ang, shape[p][0], shape[p][1]));
        return rv;
    };

    function rotatePoint(ang, x, y) {
        return [
            (x * Math.cos(ang)) - (y * Math.sin(ang)),
            (x * Math.sin(ang)) + (y * Math.cos(ang))
        ];
    };

    // custom binding to get subscript
    ko.bindingHandlers.writeSeq = {
        update: function (element, valueAccessor) {
            var value = valueAccessor();
            var valueUnwrapped = ko.unwrap(value);
            element.innerHTML = "<b>" + valueUnwrapped[0] + "</b><sub>" + valueUnwrapped[1] + "</sub>";
        }
    };

};
ko.applyBindings(new NussinovMatrixViewModel());
