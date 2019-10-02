% rebase('cgi_ui/views/strf_base.tpl', title='StructureFinder')

<!-- "dropZone" adds drag&drop support for the web site -->
<div class="container" id="dropZone">
    <h2>StructureFinder</h2>

    <div class="row">
        <div class="col-12">
            <a role="button" class="btn btn-success btn-sm ml-0" data-toggle="collapse" data-target="#adv-search"
               id="toggle_advsearch-button">Advanced Search</a>
            <a role="button" class="btn btn-warning btn-sm" id="all_structures">Show All</a>
            <a role="button" class="btn btn-primary btn-sm ml-3 {{!'' if host == '127.0.0.1' else 'invisible'}}"
               id="cellsearchcsd_button"
               href="http://{{my_ip}}/csd" target="_blank">CellCheckCSD</a>
        </div>
    </div>
    <!-- ------------  The collapsible for simple search options: -----------------  -->
    <div class="row" id="mainsearch">
        <div class="column col-6">
            <div class="input-group input-group-sm">
                <div class="input-group-prepend">
                    <span class="input-group-text" data-toggle="tooltip" title="Search for a Unit Cell">Unit Cell</span>
                </div>
                <input type="text" class="form-control"
                       placeholder="a  b  c  &alpha;  &beta;  &gamma;    (or drag&drop .p4p, .res, cif file)"
                       style="font-style: italic" id="smpl_cellsrch" name="cell">
                <div class="input-group-append">
                    <button role="button" class="btn btn-primary ml-0" type="submit" id="smpl_cellsrchbutton">Search
                    </button>
                </div>
            </div>
        </div>
        <div class="column col-6">
            <div class="input-group input-group-sm">
                <div class="input-group-prepend">
                    <span class="input-group-text" data-toggle="tooltip" title="Search for a Unit Cell">Text</span>
                </div>
                <input type="text" class="form-control" placeholder="Search Text" id="smpl_textsrch" name="text">
                <div class="input-group-append">
                    <button role="button" class="btn btn-primary ml-0 mr-0" type="submit" id="smpl_textsrchbutton">
                        Search
                    </button>
                </div>
            </div>
        </div>
    </div>

    <!-- ---------------------    End of simple search     ----------------------    -->


    <!------ Advanced search collapsible: -------->
    <div id="adv-search" class="collapse">
        <div class="row">
            <div class="column col-12">
                <div class="btn-group btn-group-sm" role="group">
                    <span class="badge badge-secondary mr-2" id="more_info_badge">Info</span>
                    <div class="form-check form-check-inline">
                        <input class="form-check-input" type="checkbox" id="more_results" value="">
                        <label class="form-check-label" for="more_results">More cell search results</label>
                    </div>
                    <div class="form-check form-check-inline">
                        <input class="form-check-input" type="checkbox" id="supercells" value="">
                        <label class="form-check-label" for="supercells">Find supercells</label>
                    </div>

                    <div class="d-inline">
                        %include('cgi_ui/views/spgr.tpl')
                        Find by space group
                    </div>
                </div>
            </div>
        </div>

        <div id="more-cell-info" class="collapse">

            <div class="card mb-2">
                <div class="card-body p-1">
                    <div class="row">
                        <div class="col-6">
                            <b>regular</b><br>
                            volume: &plusmn;3 %, length: 0.06&nbsp;&angst;, angle: 1.0&deg;<br>
                            <br>
                            <b>more results option</b><br>
                            volume: &plusmn;9 &percnt;, length: 0.2&nbsp;&angst;, angle: 2.0&deg;<br>
                        </div>
                        <div class="col-6">
                            <b>Supercells</b>
                            <br>
                            Find also unit cells of 1, 2, 3, 4, 6, 8, 10 times the volume.
                            <br>
                            <b>Space group search</b>
                            <br>
                            Be aware that not every cif file before SHELXL-2013 has a space group number. These will not
                            be
                            found.
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <div class="row">
            <div class="col-6">
                <div class="form-row">
                    <div class="col">
                        <div class="input-group input-group-sm mb-2 mr-sm-2">
                            <div class="input-group-prepend">
                                <span class="input-group-text input-group-sm" data-toggle="tooltip"
                                      title="">Unit Cell</span>
                            </div>
                            <input type="text" class="form-control form-control-sm" style="font-style: italic"
                                   id="cell_adv" placeholder="a b c &alpha; &beta; &gamma;">
                        </div>
                    </div>
                </div>
            </div>
            <div class="col-6">
                <div class="form-row">
                    <div class="col">
                        <div class="input-group input-group-sm w2ui-field">

                            <div class="input-group-prepend" data-toggle="tooltip">
                                <span class="input-group-text input-group-sm" data-toggle="tooltip"
                                      title="">Date from</span>
                            </div>
                            <input type="text" class="form-control form-control-sm" style="font-style: italic"
                                   id="date1" placeholder="">

                            <span class="input-group-text input-group-sm">to</span>
                            <input class="input-sm" title="Date" type="text" id="date2" style="width: 95%">
                            <a type="button" class="input-group-addon"
                               data-toggle="tooltip" title="Search for structures modified during the last month."
                               id="lastmsearchlink"> From Last Month</a>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <div class="row">
            <div class="column col-xs-6">
                <div class="input-group input-group-sm input-group-prepend has-success">
                    <span class="input-group-addon" data-toggle="tooltip" title="should contain">Elements</span>
                    <input type="text" class="form-control form-sm" placeholder="C H O ... (should contain)"
                           pattern="^[A-z]{1,}$" id="elements_in">
                    <input class="checkbox-addon" type="checkbox" aria-label="Only these elements"
                           title="Only above Elements" id="onlythese_elem"> Only above Elements
                </div>
            </div>
            <div class="column col-xs-6">
                <div class="input-group input-group-sm has-error">
                    <span class="input-group-addon" data-toggle="tooltip" title="should not contain">Elements</span>
                    <input type="text" class="form-control form-sm" placeholder="C H O ... (should not contain)"
                           pattern="^[A-z]{1,}$" id="elements_out">
                </div>
            </div>
        </div>

        <div class="row">
            <div class="column col-xs-6">
                <div class="input-group input-group-sm has-success">
                    <span class="input-group-addon" data-toggle="tooltip" title="should contain">Text</span>
                    <input type="text" class="form-control form-sm" placeholder="should contain" id="text_in">
                </div>
            </div>
            <div class="column col-xs-6">
                <div class="input-group input-group-sm has-error">
                    <span class="input-group-addon" data-toggle="tooltip" title="should not contain">Text</span>
                    <input type="text" class="form-control form-sm" placeholder="should not contain" id="text_out">
                </div>
            </div>
        </div>

        <div class="row">
            <div class="column col-xs-12">
                <a type="button" class="btn btn-sm btn-success" id="advsearch-button" style="min-width:90px">
                    Search </a>
            </div>
        </div>

        <br>
    </div>

    <!-- End of collapsible for search options. -->
    <div class="row">
        <div class="column col-12">
            <div id="mygrid" style="height: 450px"></div>
        </div>

        <br>
    </div>

    <!-- End of collapsible for search options. -->
    <div class="row">
        <div class="column col-sm-12">
            <div class="card mt-2">
                <div class="card-header">Unit Cell</div>
                <div class="card-body">
                    <span class="invisible btn-group" id="cellrow"> </span>
                    <span class="btn btn-default glyphicon glyphicon-copy invisible" id="cell_copy_btn"
                          data-toggle="tooltip" title="Copy cell to clipboard." data-clipboard-target="#hidden-cell">
                    </span>
                </div>
            </div>
        </div>
    </div>

    <div class="row mt-2">
        <div class="col-4">
            <div class="card w-100 p-0 mr-0 invisible" id="jsmolMaincolumn">
                <div class="card-body p-0 ml-0" id="jsmolcolumn">
                </div>
                <div class="card-footer">
                    <div class="form-check form-check-inline">
                        <input class="form-check-input" type="checkbox" value="" id="growCheckBox">
                        <label class="form-check-label" for="growCheckBox">Grow Structure
                        </label>
                    </div>
                </div>
            </div>
        </div>
        <div class="col-4 mx-0" id="residualstable1"></div>
        <div class="col-4 mx-0" id="residualstable2"></div>
    </div>


    <div class="row">
        <div class="col-sm-10">
            <span id="residuals"></span>
        </div>
        <div class="col-sm-2">

        </div>
    </div>

    <address>
        <strong>
            <a href="https://www.xs3.uni-freiburg.de/research/structurefinder">StructureFinder</a>
            <span id="version"></span> by Daniel Kratzert
        </strong><br>
        <a href="mailto:daniel.kratzert@ac.uni-freiburg.de">daniel.kratzert@ac.uni-freiburg.de</a><br>
        <p><a href="http://{{my_ip}}/dbfile.sqlite" download="structurefinder.sqlite" type="application/*">Download
            database file</a></p>
    </address>

</div>  <!-- End of the main container div -->
