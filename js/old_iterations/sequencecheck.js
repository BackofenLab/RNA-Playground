//////////////////////////////////////////////////////////////////////////////////////
// RNA Sequence checker                                                             //
// Institute: Department of BioInformatics, University of Freiburg, Germany         //
//////////////////////////////////////////////////////////////////////////////////////


/**
 * @Called on clicking the "GO!" button. Checks if the sequence is correct or not. Also it transforms the input sequence into upper case letters.
 */


/**
 * Introduces a loading GIF when GO button is pressed and keeps it there till the end of processing.
 */
function presequencecheck()
{
    console.log("enter preseqcheck");
    $("#GO").css('pointer-events', 'none');
    //disableGO();

    $("#output_title").remove();
    $("#loadinggif").show();
    $("#loadingtext").show();

    var gifshowingflag = 0;
    var processor = setInterval( function()
        {
            if(gifshowingflag == 1)
            {
                sequencecheck();
                clearInterval(processor);
            }
            if(gifshowingflag == 0)
            {
                $("#output_title").remove();
                $("#matrix_out").remove();


                $( "#output" ).before( "<h1 id='output_title'> Output</h1>");
                $("#output").after('<div class="row" id="matrix_out"></div>');
                gifshowingflag = 1;
            }
        }
        ,300);
    //enableGO();
    console.log("existing preseqcheck()");
}

/**
 * Checks the sequence, if it is correct then it call the visualize function.
 */
function sequencecheck()
{
    console.log("enter preseqchk");
    //Global Variable Userinput
    userInput = document.getElementById("userInput").value; //getting the sequence value entered by the user
    var sequencearray="";
    for(var i=0;i<userInput.length;i++)
    {
        sequencearray=sequencearray+userInput[i].toUpperCase(); //set the letter to UPPERCASE if its in small case.
    }
    userInput = sequencearray;
    sequence = userInput;

    visualize();
    $("#loadinggif").hide(); // REMOVING THE GIF NOW
    $("#loadingtext").hide();   //REMOVING THE LOADING TEXT
    console.log("exiting sequence check");

}