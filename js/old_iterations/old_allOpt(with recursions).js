/**
 * 
 * Created by moni on 09.05.16.
 */


/*********************************************************************************************
 All optimal Tracebacks (Start)
 **********************************************************************************************/

/**
 * Output stŕing of dots of length len
 * @param (int) len length of string to be created.
 */
dotString: function (len) { //returns string of dots of length len.
    var dstring;
    if (len <= 0) return "";
    for (var i = 0; i < len; i++) {
        if (dstring === undefined)
            dstring = ".";
        else
            dstring = dstring + ".";
    }
    return dstring;
}
,
setCharAt: function (str, index, chr) {
    //console.log("str:",str,"index:",index,'chr:',chr);

    if (index > str.length - 1) {
        return str;
    }
    return str.substr(0, index) + chr + str.substr(index + chr.length);
}
,


addBasePairs: function (structures, bpList, leftShift) {
    //console.log("addBasepairs enter:", JSON.stringify(structures));
    if (bpList != undefined) {
        for (var i in structures) {
            structures[i].structure = this.setCharAt(structures[i].structure, bpList[0] - leftShift, "(");
            structures[i].structure = this.setCharAt(structures[i].structure, bpList[1] - leftShift, ")");

        }
        ;
    }
    ;

    //console.log(structures);
    return structures;

}
,

traceallOpt: function (i, j) { //main entry function for finding all optimal structures
    //console.log("traceAllOpt:", [i, j]);
    if (this.getValue(i, j) === 0) {
        //console.log("finding ij parents:", this.getTraces(i,j));
        var optTrace = Object.create(Traceback);
        optTrace.traces = [];
        optTrace.setStructure(this.dotString(j - i + 1));
        //console.log("x1:", [optTrace]);
        //optTrace.setTraces([i,j]);
        return [optTrace];
    }
    ;
    var traces = []; //list of parent cells
    for (var x in this.getTraces(i, j)) { //get complete list of parent cells of cell (i,j)
        traces.push(this.getTraces(i, j)[x]);
    }
    ;
    //console.log("traces:", JSON.stringify(traces), "i,j:",[i,j]);
    var structures = [];
    for (var t in traces) { //iterate through list of parents
        //console.log("t in traces:", traces[t], "i,j:", [i, j]);
        //if (i === 1 && j === this.sequence.length)
        //    console.log("ij:", [i, j], "\n");
        var traceStructures = this.combineSubstructure(traces[t].parents, i, j);
        traceStructures = this.addBasePairs(traceStructures, traces[t].bps[0], i);
        //console.log("tS:", JSON.stringify(traceStructures));
        //console.log("ot", JSON.stringify(optimalTrace));
        for (var tS in traceStructures) {
            structures.push(traceStructures[tS]);
        }
        //if (i === 1 && j === this.sequence.length)
        //    console.log("i,j:", [i, j], "\n");
    }
    ;
    //console.log("x2:", JSON.stringify(structures));
    return structures;
}
,

combineSubstructure: function (parents, i, j) {
    var trace = [[i, j], parents];
    var traces = [];
    //console.log("css:", trace);

    if (parents.length === 0) {
        traces.push(trace);
        var ret = [{structure: this.dotString(j - i + 1), traces: traces}];
        //console.log("1:", JSON.stringify(ret));
        //console.log("1:");
        return ret; // return this.dotString(j - i + 1);
    }

    var allpinv = true;
    for (var p in parents) {
        allpinv = allpinv && (parents[p][0] > parents[p][1]);
    }
    if (allpinv) {
        traces.push(trace);
        var ret = [{structure: this.dotString(j - i + 1), traces: traces}];
        //console.log("2:", JSON.stringify(ret));
        //console.log("2:");
        return ret; //return this.dotString(j - i + 1);
    }

    var parentSubstructures = [];
    var numNonEmplist = 0;
    traces.push(trace);
    for (var parent in parents) {
        //console.log("tao call:", parents[parent][0], parents[parent][1]);
        var pSS = this.traceallOpt(parents[parent][0], parents[parent][1]);
        //console.log("pss1:", JSON.stringify(pSS));

        var temp = pSS;
        if (pSS[0].structure != "") {
            numNonEmplist++;
        }
        //console.log("pss3:",pSS);
        parentSubstructures.push(pSS);

    }

    //console.log("\nPSS:", JSON.stringify(parentSubstructures));
    //console.log("nel:", numNonEmplist);
    //console.log("trace:", JSON.stringify(trace))
    //console.log("traces:", JSON.stringify(traces));
    if (numNonEmplist === 1) {
        var nonEmp;
        for (nonEmp in parents) { //select one of the parents
            if (parentSubstructures[nonEmp][0].structure != "") {
                break;
            }
        }
        var pi = parents[nonEmp][0];
        var pj = parents[nonEmp][1];
        var prefix = this.dotString(pi - i);
        var suffix = this.dotString(j - pj);
        if (prefix == undefined) prefix = "";
        if (suffix == undefined) suffix = "";
        //console.log("ij:", i , j , "parent:", parents[nonEmp], " prefix:", prefix, " suffix:", suffix);
        //var struct = this.constructStruct2(prefix, suffix, parentSubstructures[nonEmp]);
        //console.log("pss bs", JSON.stringify(parentSubstructures));

        for (var nE in parentSubstructures[nonEmp]) {
            //console.log("p:",parentSubstructures[nonEmp][nE].structure);
            parentSubstructures[nonEmp][nE].structure = prefix + parentSubstructures[nonEmp][nE].structure + suffix;

            parentSubstructures[nonEmp][nE].traces.push(trace);
        }
        //console.log("pss as", JSON.stringify(parentSubstructures));
        for (var pS in parentSubstructures) { //select one of the parents
            if (parentSubstructures[pS][0].structure == "") {
                //console.log("true:",parentSubstructures[pS][0]);
                parentSubstructures.splice(pS, 1);
                break;
            }
        }
        //console.log("3:", JSON.stringify(parentSubstructures[0]));
        //console.log("3:");
        //  debugger;

        return parentSubstructures[0];//return struct;
    }
    var ds = {structure: this.dotString(j - i + 1), traces: traces};
    //console.log("\ncall addSubstructure:");
    var struct = this.addSubstructure(ds, i, parents, parentSubstructures, 0);
    //console.log("4:", JSON.stringify(struct), "\n");
    //console.log("4:");
    return struct;
}
,

addSubstructure: function (structure, leftShift, parents, parentSubstructure, parentIndex) {
    //console.log("\nass:", structure, leftShift, parentSubstructure, parentIndex);
    if (parentIndex >= parentSubstructure.length) {
        //console.log("a1:", JSON.stringify(structure), "psslen:", parentSubstructure.length, "pI:", parentIndex);
        return [structure];
    }
    ;
    //console.log("\nass:", JSON.stringify(structure), leftShift, JSON.stringify(parentSubstructure), parentIndex);
    var structures = [];
    var temp_struct = [];
    var parent = parents[parentIndex];
    for (var substructure in parentSubstructure[parentIndex]) {
        for (var i = parent[0]; i < parent[1]; i++) {
            if (parentSubstructure[parentIndex][i - parent[0]] != undefined) {
                //console.log("structure:", structure.structure, i - leftShift, parentSubstructure[parentIndex][substructure].structure);
                //console.log("traces:", JSON.stringify(structure.traces), " pss trace ", JSON.stringify(parentSubstructure[parentIndex][substructure].traces));
                //console.log("-------------------------");
                //console.log(JSON.stringify(parentSubstructure[parentIndex][substructure].traces));
                for (var tr in parentSubstructure[parentIndex][substructure].traces){
                    structure.traces.push(parentSubstructure[parentIndex][substructure].traces[tr]);
                    //console.log(parentSubstructure[parentIndex][substructure].traces[tr])
                }
                structure.structure = this.setCharAt(structure.structure, i - leftShift, parentSubstructure[parentIndex][substructure].structure);
                //structure.traces.push(parentSubstructure[parentIndex][substructure].traces);

                //console.log("-------------------------");
                //console.log("structure:", structure.structure, i - leftShift, parentSubstructure[parentIndex][substructure].structure);
                //console.log("traces:", JSON.stringify(structure.traces), " pss trace ", JSON.stringify(parentSubstructure[parentIndex][substructure].traces));
                structure.structure = this.setCharAt(structure.structure, i - leftShift, parentSubstructure[parentIndex][substructure].structure);
                i = parent[1];
            }
            ;
        }
        ;
        var temp = Object.create(structure);
        temp.structure = structure.structure;
        temp.traces = structure.traces;
        temp_struct = this.addSubstructure(temp, leftShift, parents, parentSubstructure, parentIndex + 1);
        //console.log("temp_struct:", JSON.stringify(ptemp_struct[0]));
        for (var tmp in temp_struct) {
            structures.push(temp_struct[tmp]);
        }
        //structures.push(temp_struct[0]);
        //console.log("here structures:",JSON.stringify(structures));
        for (var j = parent[0]; j < parent[1]; j++) {
            //console.log(structure.structure, j - leftShift);
            structure.structure = this.setCharAt(structure.structure, j - leftShift, ".");
            //console.log(structure.structure, j - leftShift);
            break;
        }
        ;
    }
    ;
    //console.log("a2:", JSON.stringify(structures), "psslen:", parentSubstructure.length, "pI:", parentIndex);
    return structures;
}
,
breakTb: function (tb) {
    var tbstring = [];
    for (var i in tb) {
        var divided = tb[i].split(",");
        //console.log(divided);
        for (var j in divided) {
            tbstring.push(divided[j]);
        }
    }
    return tbstring;
}
