<!DOCTYPE HTML>
<html lang="en" xmlns="http://www.w3.org/1999/html">
<head>
    <meta charset="UTF-8">
    <title>StructureFinder</title>

    <link rel="stylesheet" href="http://{{my_ip}}/static/w2ui/w2ui-1.4.3.min.css">
    <link rel="stylesheet" href="http://{{my_ip}}/static/bootstrap-3.3.7/css/bootstrap.min.css">
    <link rel="stylesheet" href="http://{{my_ip}}/static/bootstrap-3.3.7/css/bootstrap-theme.min.css">
    <script src="http://{{my_ip}}/static/jquery/jquery-3.2.1.min.js"></script>
    <script src="http://{{my_ip}}/static/bootstrap-3.3.7/js/bootstrap.min.js"></script>
    <script src="http://{{my_ip}}/static/jsmol/JSmol_dk.nojq.lite.js"></script>
    <script src="http://{{my_ip}}/static/w2ui/w2ui-1.4.3.min.js"></script>
    <script src="http://{{my_ip}}/static/clipboard/clipboard.min.js"></script>

<style type="text/css">

body {
    /*background-color: #ffffff;*/
    font-size: 12px;
    line-height: inherit;
}

.border-right {
    border-right: 1px solid #bcbcbc;
}

.collapsing {
    transition: 1ms;
}

#resitable1 tr td {
    line-height: 13px;
    height: 13px;
    overflow: auto;
}

#resitable2 tr td {
    line-height: 13px;
    height: 13px;
    overflow: auto;
}

.btn-group {
    padding-bottom: 4px;
    padding-top: 4px;
}

.input-group {
    padding-bottom: 2px;
    padding-top: 2px;
}

.input-group-addon {
    min-width:75px;
    text-align:left;
}

#resfile {
    font-family: "Bitstream Vera Sans Mono", Monaco, "Courier New", Courier, monospace;
}

#jsmolcolumn {
    margin-right: 15px;
    margin-left: 15px;
    width: 360px;
    height: 320px;
}

</style>


<script>
var cgifile = 'http://{{my_ip}}';

// height of cif list to 35% of the screen height:
$(document).ready(function(){
    //gets the window's height
    var b = $(window).height();
    var h = b * 0.35;
    if (h < 200) {
        h = 220;
    }
    $("#mygrid").css("height", h);
});

$(document).ready(function(){
    $('[data-toggle="cell_tooltip"]').tooltip();
});

elements = ['X',  'H',  'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
            'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',  'Ca', 'Sc', 'Ti', 'V', 'Cr',
            'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
            'Rb', 'Sr', 'Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
            'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
            'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
            'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
            'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es'];

</script>


</head>

<body>


<div class="container">
<h2>StructureFinder</h2>


<a type="button" class="btn btn-primary btn-xs" data-toggle="collapse" data-target="#adv-search"
        id="toggle_advsearch-button">Advanced Search</a>
<a type="button" class="btn btn-warning btn-xs" id="all_structures">Show All</a>
<!-- ------------  The collapsible for simple search options: -----------------  -->
<div class="form-group row" id="mainsearch">
    <div class="col-sm-6">
        <div class="input-group input-group-sm">
            <span class="input-group-addon" data-toggle="tooltip" title="Search for a Unit Cell">Unit Cell</span>
            <input type="text" class="form-control" placeholder="a b c &alpha; &beta; &gamma;" style="font-style: italic" id="smpl_cellsrch" name="cell">
                <div class="input-group-sm input-group-btn">
                <button class="btn btn-default" type="submit" id="smpl_cellsrchbutton">
                    <i class="glyphicon glyphicon-search"></i>
                </button>
                </div>
            </input>
        </div>
    </div>
    <div class="col-sm-6">
        <div class="input-group input-group-sm">
            <span class="input-group-addon" data-toggle="tooltip" title="Search for a Unit Cell">Text</span>
            <input type="text" class="form-control" placeholder="Search Text" id="smpl_textsrch" name="text">
                <div class="input-group-sm input-group-btn">
                    <button class="btn btn-default" type="submit" id="smpl_textsrchbutton">
                        <i class="glyphicon glyphicon-search"></i>
                    </button>
                </div>
            </input>
        </div>
    </div>
</div>
<!-- ---------------------    End of simple search     ----------------------    -->


<!-- The collapsible for Advanced search options: -->
<div id="adv-search" class="collapse">
    <div class="row">
        <div class="col-sm-12">
            <div class="btn-group btn-group-sm" role="group">
                <a href="#" class="badge" id="more_info_badge">info</a> &nbsp;&nbsp;&nbsp;
                <input type="checkbox" value="" id="more_results">More cell search results &nbsp;&nbsp;&nbsp;
                <input type="checkbox" value="" id="supercells">Find supercells &nbsp;&nbsp;&nbsp;
                {{!space_groups}} Find by space group
            </div>
        </div>
    </div>

    <div id="more-cell-info" class="collapse">
        <div class="row">
            <div class="col-sm-6">
                <b>regular</b><br>
                volume: &plusmn;3 %, length: 0.06&nbsp;&angst;, angle: 1.0&deg;<br>
                <br>
                <b>more results option</b><br>
                volume: &plusmn;9 &percnt;, length: 0.2&nbsp;&angst;, angle: 2.0&deg;<br>
            </div>
            <div class="col-sm-6">
                <b>Supercells</b>
                <br>
                Find also unit cells of 1, 2, 3, 4, 6, 8, 10 times the volume.
                <br>
                <b>Space group search</b>
                <br>
                Be aware that not every cif file before SHELXL-2013 has a space group number. These will not be found.
            </div>
        </div>
    </div>

    <div class="row">
        <div class="column col-xs-6">
            <div class="input-group input-group-sm advsearchfields">
                <span class="input-group-addon">Uni Cell</span>
                <input type="text" class="form-control form-sm" style="font-style: italic" placeholder="a b c &alpha; &beta; &gamma;" id="cell_adv">
            </div>
        </div>
        <div class="column col-md-6">
            <div class="input-group input-group-sm w2ui-field advsearchfields">
                <span class="input-group-addon" data-toggle="tooltip"
                      title="Search for structures that were modified between two dates">Date from</span>
                <input class="input-sm" type="text" id="date1" style="width: 95%">
                <span class="input-group-addon">to</span>
                <input class="input-sm" type="text" id="date2" style="width: 95%">
                <a type="button" class="btn btn-sm btn-standard input-group-addon"
                   data-toggle="tooltip" title="Search for structures modified during the last month."
                   id="lastmsearchlink"> From Last Month</a>
            </div>
        </div>

        <script>
            var d1 = $('input[id=date1]');
            d1.w2field('date', {
                format: 'yyyy-mm-dd',
                end: d1
            });
            $('input[id=date2]').w2field('date', {
                format: 'yyyy-mm-dd',
                start: d1
            });
        </script>
    </div>

    <div class="row">
        <div class="column col-xs-6">
            <div class="input-group input-group-sm has-success advsearchfields">
                <span class="input-group-addon" data-toggle="tooltip" title="should contain">Elements</span>
                <input type="text" class="form-control form-sm" placeholder="C H O ... (should contain)"
                       pattern="^[A-z]{1,}$" id="elements_in">
            </div>
        </div>
        <div class="column col-xs-6">
            <div class="input-group input-group-sm has-error advsearchfields">
                <span class="input-group-addon" data-toggle="tooltip" title="should not contain">Elements</span>
                <input type="text" class="form-control form-sm" placeholder="C H O ... (should not contain)"
                       pattern="^[A-z]{1,}$" id="elements_out">
            </div>
        </div>
    </div>

    <div class="row">
        <div class="column col-xs-6">
            <div class="input-group input-group-sm has-success advsearchfields">
                <span class="input-group-addon" data-toggle="tooltip" title="should contain">Text</span>
                <input type="text" class="form-control form-sm" placeholder="should contain" id="text_in">
            </div>
        </div>
        <div class="column col-xs-6">
            <div class="input-group input-group-sm has-error advsearchfields">
                <span class="input-group-addon" data-toggle="tooltip" title="should not contain">Text</span>
                <input type="text" class="form-control form-sm" placeholder="should not contain" id="text_out">
            </div>
        </div>
    </div>

    <div class="row">
        <div class="column col-xs-12">
            <a type="button" class="btn btn-sm btn-success" id="advsearch-button" style="min-width:90px"> Search </a>
        </div>
    </div>

<br>
</div>

<!-- End of collapsible for search options. -->
<div class="row">
    <div class="column col-lg-12">
        <div id="mygrid" style="height: 450px"></div>
    </div>
    <div class="column col-lg-12">
        <div class="panel panel-default">
            <div class="panel-body">
                <span class="invisible btn-group" id="cellrow" style="font-size: 14px"> </span>
                <span class="btn btn-default glyphicon glyphicon-copy invisible" id="cell_copy_btn"
                      data-toggle="tooltip" title="Copy cell to clipboard." data-clipboard-target="#hidden-cell">
                </span>
            </div>
        </div>
    </div>
</div>

<div class="row">
    <div class="column col-md-4 panel panel-default invisible" id="jsmolcolumn"
         data-toggle="tooltip" title="Asymmetric Unit">
    </div>
    <div class="column col-md-4" id="residualstable1"> </div>
    <div class="column col-md-4" id="residualstable2"> </div>
</div>


<div class="row">
    <div class="col-sm-10">
        <span id="residuals"></span>
    </div>
    <div class="col-sm-2">

    </div>
</div>

<address>
    <strong><a href="https://www.xs3.uni-freiburg.de/research/structurefinder">StructureFinder</a> by Daniel Kratzert</strong><br>
  <a href="mailto:daniel.kratzert@ac.uni-freiburg.de">daniel.kratzert@ac.uni-freiburg.de</a>
</address>

</div>  <!-- End of the main container div -->


</body>

<!-- Javascript functions start here: -->
<script>

jQuery(document).ready(function($) {

    // The main structures table:
    $('#mygrid').w2grid({
        name: 'mygrid',
        header: 'StructureFinder',
        url: cgifile+"/all",
        method: 'GET',
        show: {
            toolbar: false,
            footer: true
        },
        columns: [
            {field: 'recid',    caption: 'ID',        size: '45px', sortable: false, attr: 'align=center'},
            {field: 'filename', caption: 'filename',  size: '20%',  sortable: false, resizable: true},
            {field: 'dataname', caption: 'dataname',  size: '15%',  sortable: false, resizable: true},
            {field: 'path',     caption: 'directory', size: '65%',  sortable: false, resizable: true}
        ],
        searches: [
            {field: 'filename', caption: 'filename', type: 'text'},
            {field: 'dataname', caption: 'dataname', type: 'text'},
            {field: 'path', caption: 'directory', type: 'text'}
        ],
        //sortData: [{field: 'dataname', direction: 'ASC'}],
        onSelect:function(event) {
            showprop(event.recid);
            //console.log(event);
        }
    });

    // Do advanced search:
    var advanced_search_button = $("#advsearch-button");
    advanced_search_button.click(function(event) {
        var txt_in = document.getElementById("text_in").value;
        var txt_out = document.getElementById("text_out").value;
        var elements_in = document.getElementById("elements_in").value;
        var elements_out = document.getElementById("elements_out").value;
        var cell_adv = document.getElementById("cell_adv").value;
        var more_res = $('#more_results').is(':checked');
        var supercell = $('#supercells').is(':checked');
        var datefield1 = document.getElementById("date1").value;
        var datefield2 = document.getElementById("date2").value;
        var itnum = $("#IT_number").val().split(" ")[0];
        advanced_search(txt_in, txt_out, elements_in, elements_out, cell_adv, more_res,
                        supercell, datefield1, datefield2, itnum);
    });

    function is_elem_doubled(elements_in, elements_out) {
        var sumlist = elements_in.split(" ");
        var outlist = elements_out.split(" ");
        var ok = true;
        ok = validateSumForm(elements_in);
        for (i = 0; i < sumlist.length; i++) {
            var el = sumlist[i];
            if ($.inArray(el, outlist) >= 0) {
                // A space character is allowed:
                if (el === "") {
                    continue
                }
                ok = false;
            }
        }
        return ok;
    }

    function check_elin() {
        var elements_in = document.getElementById("elements_in").value;
        var elements_out = document.getElementById("elements_out").value;
        return is_elem_doubled(elements_in, elements_out);
    }

    function check_elex() {
        var elements_in = document.getElementById("elements_in").value;
        var elements_out = document.getElementById("elements_out").value;
        return is_elem_doubled(elements_out, elements_in)
    }

    function elements_red() {
        var elinform = $("#elements_in.form-control");
        var elexform = $("#elements_out.form-control");
        elinform.css("color", "#f35e59");
        elinform.css("font-weight", "bold");
        elexform.css("color", "#f35e59");
        elexform.css("font-weight", "bold");
    }

    function elements_regular() {
        var elinform = $("#elements_in.form-control");
        var elexform = $("#elements_out.form-control");
        elinform.css("color", "#000000");
        elinform.css("font-weight", "normal");
        elexform.css("color", "#000000");
        elexform.css("font-weight", "normal");
    }

    function validate_element_input() {
        if (!check_elin() || !check_elex()){
            elements_red();
        } else {
            elements_regular();
        }

    }

    // Validators for chemical elemets included search field:
    $("#elements_in").keyup(function() {
        validate_element_input();
    });

    // Validators for chemical elemets excluded search field:
    $("#elements_out").keyup(function() {
        validate_element_input();
    });

    function validateSumForm(sumform) {
        // Validates if sumform contains only valid chemical elements
        // Space characters are allowed
        var ok = true;
        //console.log(sumform);
        if (sumform.length === 0) {
            return true;
        }
        var sumlist = sumform.split(" ");
        //console.log(sumlist);
        for (i = 0; i < sumlist.length; i++) {
            var el = sumlist[i];
            //console.log(el);
            if ($.inArray(el, elements) === -1) {
                // A space character is allowed:
                if (el.length === 0) {continue}
                ok = false;
            }
        }
        return ok;
    }

    // Toggle search info:
    var more_info_button = $('#more_info_badge');
    more_info_button.click(function(){
        var button_text = more_info_button.text();
        $("#more-cell-info").toggle(1);
    });

    $('#all_structures').click(function () {
        w2ui['mygrid'].reload();
    });

    // Switch between advanced and simple search:
    var advbutton = $('#toggle_advsearch-button');
    advbutton.click(function(){
        var button_text = advbutton.text();
        $("#mainsearch").toggle(1);
        if (button_text.split(" ")[0] === "Advanced") {
            advbutton.html("Simple Search");
            document.getElementById("cell_adv").value = document.getElementById("smpl_cellsrch").value;
        } else {
            advbutton.html("Advanced Search");
            document.getElementById("smpl_cellsrch").value = document.getElementById("cell_adv").value;
        }
    });

    // Text search Button clicked:
    $("#smpl_textsrchbutton").click(function(event) {
        var txt = document.getElementById("smpl_textsrch").value;
        txtsearch(txt);
        //console.log(txt);
    });

    // Cell search Button clicked:
    $("#smpl_cellsrchbutton").click(function(event) {
        var cell = document.getElementById("smpl_cellsrch").value;
        cellsearch(cell);
    });

    // Enter key pressed in the simple text search field:
    $('#smpl_textsrch').keypress(function(e) {
        if (e.which === 13) {  // enter key
            var txt = document.getElementById("smpl_textsrch").value;
            txtsearch(txt);
            //console.log(txt);
        }
    });

    // Enter key pressed in the simple cell search field:
    $('#smpl_cellsrch').keypress(function(e) {
        if (e.which === 13) {  // enter key
            var cell = document.getElementById("smpl_cellsrch").value;
            cellsearch(cell);
            //console.log(cell);
        }
    });

    // Enter key pressed in one of the advanced search fields:
    $('.advsearchfields').keypress(function(e) {
        if (e.which === 13) {  // enter key
            advanced_search_button.click();
            //console.log(cell);
        }
    });

    function displayresultnum(result) {
        numresult = result.total;
        if (typeof numresult === 'undefined') numresult = 0;
        $("#cellrow").removeClass('invisible');
        $("#cell_copy_btn").addClass('invisible');
        document.getElementById("cellrow").innerHTML = "Found " + numresult + " structures";
    }

    function showprop(idstr) {
        /*
        This function uses AJAX POST calls to get the data of a structure and displays
        them below the main table.
        */
        // Get residuals table 1:
        $.post(url = cgifile, data = {id: idstr, residuals1: true}, function (result) {
            document.getElementById("residualstable1").innerHTML = result;
        });

        // Get residuals table 2:
        $.post(url = cgifile, data = {id: idstr, residuals2: true}, function (result) {
            document.getElementById("residualstable2").innerHTML = result;
        });

        // Get molecule data and display the molecule:
        var jsmolcol = $("#jsmolcolumn");
        $.post(url = cgifile+'/molecule', data = {id: idstr}, function (result) {
            Jmol._document = null;
            Jmol.getTMApplet("jmol", Info);
            jsmolcol.html(jmol._code);
            jmol.__loadModel(result);
            jsmolcol.removeClass('invisible');
            var tbl = $('#residualstable2');
            jsmolcol.css("height", tbl.height()-20);
        });

        // Get unit cell row:
        $.post(url = cgifile, data = {id: idstr, unitcell: true}, function (result) {
            $("#cellrow").removeClass('invisible');
            $("#cell_copy_btn").removeClass('invisible');
            document.getElementById("cellrow").innerHTML = result;

            var clipboard = new Clipboard('.btn');
            clipboard.on('success',
                function(e) {
                    e.clearSelection();
                }
            );
        });

        // display the big cif data table:
        $.post(url = cgifile, data = {id: idstr, all: true},
            function (result) {
                document.getElementById("residuals").innerHTML = result;
            }
        );
    }

    // some options for JSmol:
    var bgcolor = $(this.body).css("background-color");
    var Info;
    Info = {
        width: 320,
        height: 300,
        color: bgcolor,
        //color: "0xf0f0f0",
        shadeAtoms: false,
        addSelectionOptions: false,
        use: "HTML5",
        readyFunction: null,
        defaultModel: "",
        bondWidth: 3,
        zoomScaling: 5,
        pinchScaling: 5.0,
        mouseDragFactor: 0.9,
        touchDragFactor: 0.9,
        multipleBondSpacing: 0,
        spinRateX: -0.08,
        spinRateY: 0.05,
        spinFPS: 20,
        spin: false,
        infodiv: false,
        debug: false,
        j2sPath: "."
    };

    function advanced_search(text_in, text_out, elements_in, elements_out, cell_adv, more_res, supercell,
                             date1, date2, itnum) {
        var cell = cell_adv.replace(/\s+/g, ' ').trim();
        cell = cell.replace(/,/g, '.');  // replace comma with point
        if (!isValidCell(cell)) {
            cell = "";
        }
        var gridparams = {cell_search: cell, text_in: text_in, text_out: text_out, elements_in: elements_in,
                          elements_out: elements_out, more: more_res, supercell: supercell, date1: date1, date2: date2
                          ,it_num: itnum};
        //console.log(gridparams);
        var url;
        w2ui['mygrid'].request('get-records', gridparams,
            url = cgifile + "/adv_srch",
            function (result) {
                displayresultnum(result);
                //console.log(result);
            }
        );
    }

    // Search for structures of last month:
    $('#lastmsearchlink').click(function(){
        var date_now = new Date();
        var month = date_now.getUTCMonth();
        var day = date_now.getUTCDate();
        var lastmonth = new Date(date_now.getUTCFullYear(), month-1, day);
        var lastmdate = lastmonth.toISOString().split("T")[0];
        //console.log(lastmdate+ ' '+ date_now.toISOString().split("T")[0]);
        // From last month to now():
        advanced_search("", "", "", "", "", "", "", lastmdate, date_now.toISOString().split("T")[0]);
    });

    // Test if a valid unit cell is in cell:
    function isValidCell(cell) {
        var scell = cell.split(" ");
        //console.log(scell);
        if (isNumericArray(scell)) {
            return !(scell.length !== 6); // return True if 6 values
        } else {
            return false;
        }
    }

    // Test if all values in array are numeric:
    function isNumericArray(array) {
        var isal = true;
        for (var i=0; i<array.length; i++) {
            if (!$.isNumeric(array[i])) {
                isal = false;
            }
        }
        return isal;
    }

    function cellsearch(cell) {
        var more_res = $('#more_results').is(':checked');
        var supercell = $('#supercells').is(':checked');
        cell = cell.replace(/\s+/g, ' ').trim();  // replace multiple spaces with one
        cell = cell.replace(/,/g, '.');  // replace comma with point
        //console.log(cell);
        var params;
        var url;
        if (isValidCell(cell)) {
            w2ui['mygrid'].request('get-records',
                params = {cell_search: cell, more: more_res, supercell: supercell},
                url = cgifile + "/cellsrch",
                function (result) {
                    displayresultnum(result);
                    //console.log(result.total);
                    //console.log(more_res);
                }
            );
        }
    }

    function txtsearch(text) {
        var advanced = false;
        var params;
        var url;
        //console.log(text);
        w2ui['mygrid'].request('get-records',
            params = {text_search: text},
            url = cgifile + "/txtsrch",
            function (result) {
                displayresultnum(result);
                //console.log(result);
            }
        );
    }

});

</script>



</html>
