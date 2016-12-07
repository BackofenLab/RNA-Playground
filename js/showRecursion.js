/**
 * Created by moni on 07.12.16.
 */
<script>
function showRecursion() {
    console.log("val:", $("#formula").val());
    setTimeout(function () {      //introduces a little delay, so that the thing will close slowly and neatly
        document.getElementById("recursion").innerHTML =  availableAlgorithms[$("#formula").val()].getRecursionInLatex();
        rerendermath();
        $(".animate1").empty();
    }, 450);

};
</script>