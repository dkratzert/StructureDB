<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>StructureFinder</title>
    
    <link rel="stylesheet" href="http://{{my_ip}}/static/w2ui/w2ui-1.4.3.min.css">
    <link rel="stylesheet" href="http://{{my_ip}}/static/bootstrap-3.3.7/css/bootstrap.min.css">
    <link rel="stylesheet" href="http://{{my_ip}}/static/bootstrap-3.3.7/css/bootstrap-theme.min.css">
    <script src="http://{{my_ip}}/static/jquery/jquery-3.2.1.min.js"></script>
    <script src="http://{{my_ip}}/static/bootstrap-3.3.7/js/bootstrap.min.js"></script>
    <script src="http://{{my_ip}}/static/w2ui/w2ui-1.4.3.min.js"></script>
    <script src="http://{{my_ip}}/static/clipboard/clipboard.min.js"></script>
    <script src="http://{{my_ip}}/static/strf_web.js"></script>
    <script src="http://{{my_ip}}/static/strf_web_csd.js"></script>
    <script> 
        var cgifile = 'http://{{my_ip}}';
    </script>
    
    
</head>
<body>

<div class="container"  id="dropZone">
    <h2>CellCheckCSD</h2>
    <div class="row">
        <div class="column col-sm-7">
            <div class="input-group input-group-sm advsearchfields">
                <span class="input-group-addon">Unit Cell</span>
                <input type="text" class="form-control form-sm" style="font-style: italic" placeholder="a b c &alpha; &beta; &gamma;" 
                       id="cell_csd_inp" value="{{str_id}}">
                <span class="input-group-addon input-group-btn btn-success" id="csd_search_btn">Search</span>
            </div>
        </div>
        <div class="column col-sm-2">
            <select id="centering_drop" class="dropdown dropdown-toggle btn-md" data-toggle="dropdown">
                <option value="0">Primitive (P)</option>
                <option value="1">A-centered (A)</option>
                <option value="2">B-centered (B)</option>
                <option value="3">C-centered (C)</option>
                <option value="4">Face-centered (F)</option>
                <option value="5">Body-centered (I)</option>
                <option value="6">Rhombohedral (R)</option>
            </select>
        </div>
        <div class="column col-sm-2">
            <a type="button" class="btn btn-default btn-sm" id="backtomain_button" href="http://{{my_ip}}">Back to StructureFinder</a>
        </div>
    </div>
    <div class="row"> <br></div>
    <!-- ------------- The main table: ------------ -->
    <div class="row">
        <div class="column col-sm-12">
            <div class="w2ui-grid" id="my_ccdc_grid" style="height: 450px">
                
            </div>
        </div>
    </div>    


    <div class="row">
    <div class="col-sm-12">
        <address>
            <strong><a href="https://www.xs3.uni-freiburg.de/research/structurefinder">StructureFinder</a> 
                <span id="version"></span> by Daniel Kratzert</strong><br>
          <a href="mailto:daniel.kratzert@ac.uni-freiburg.de">daniel.kratzert@ac.uni-freiburg.de</a>
        </address>
    </div>
</div>
    

</div>  <!-- End of main container -->



<script>
$(document).ready(function($){
    

    
    
});

</script>

</body>
</html>