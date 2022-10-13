//tmp
const runID = 1;
var uploadName = "GSM3860733_E10";
const nDims = 50;
var expTitle = "P14_Prime"; //TODO: change to undefined
var selectedCluster = "0" //TODO: change to undefined;
var selectedClusterDEGAsList = undefined;
var selectedClusterDEGListAsStr = undefined;

var selectedDBPhase = undefined;
var selectedDBSpecies = undefined;

var savedDEGLists = [];
var selectedSavedDEGList = undefined;

var createDEGListManuallyModal = undefined;
var saveAsDEGListModal = undefined;

$(document).ready(function() {

    $('#btn-accordion-filtering-run').click(function () {

        //show spinner
        $('#accordion-filtering-spinner').removeClass("d-none")

        const minNFeatureRNA = $('#txt-min-nfeature-rna').val();
        const maxNFeatureRNA = $('#txt-max-nfeature-rna').val();
        const percentMT = $('#txt-percent-mt').val();

        json_url = `/run_r_script/qc_metrics/${runID}/${uploadName}/${minNFeatureRNA}/${maxNFeatureRNA}/${percentMT}`
        $.get(json_url, function(response) {
        })
            .done(function(response) {
                if (response) {

                    let data = response;
                    //console.log(data);
                    const img_src = data["img_src"];

                    //remove previous response result then show new response
                    //$('#accordion-filtering-container img').attr("src", "");
                    $('#accordion-filtering-container img').attr("src", img_src);

                    //hide spinner
                    $('#accordion-filtering-spinner').addClass("d-none");

                    //show img
                    $('#accordion-filtering-container').hide();
                    $('#accordion-filtering-container').removeClass('d-none')
                    $('#accordion-filtering-container').fadeIn();


                }
            })
            .fail(function() {
                alert( "Error. Please try again later." );
            })
    });

    $('#btn-accordion-elbow-run').click(function () {

        //show spinner
        $('#accordion-elbow-spinner').removeClass("d-none")

        const minNFeatureRNA = $('#txt-min-nfeature-rna').val();
        const maxNFeatureRNA = $('#txt-max-nfeature-rna').val();
        const percentMT = $('#txt-percent-mt').val();
        const nDims = $('#txt-elbow-n-dims').val();

        json_url = `/run_r_script/elbow_plot/${runID}/${uploadName}/${minNFeatureRNA}/${maxNFeatureRNA}/${percentMT}/${nDims}`;
        $.get(json_url, function(response) {
        })
            .done(function(response) {
                if (response) {

                    let data = response;
                    //console.log(data);
                    const img_src = data["img_src"];

                    //remove previous response result then show new response
                    //$('#accordion-elbow-container img').attr("src", "");
                    $('#accordion-elbow-container img').attr("src", img_src);

                    //hide spinner
                    $('#accordion-elbow-spinner').addClass("d-none");

                    //show img
                    $('#accordion-elbow-container').hide();
                    $('#accordion-elbow-container').removeClass('d-none')
                    $('#accordion-elbow-container').fadeIn();


                }
            })
            .fail(function() {
                alert( "Error. Please try again later." );
            })
    });

     $('#btn-accordion-run-exp-run').click(function () {

        //show spinner
        $('#accordion-run-exp-spinner').removeClass("d-none")

        const minNFeatureRNA = $('#txt-min-nfeature-rna').val();
        const maxNFeatureRNA = $('#txt-max-nfeature-rna').val();
        const percentMT = $('#txt-percent-mt').val();
        //TODO: validate exp title: only letters, numbers and _ (underscore)
        expTitle = $("#txt-run-exp-title").val(); // update the global var here for next operations
        const nDims = $('#txt-run-exp-n-dims').val();
        const clusteringRes = $('#txt-run-exp-clustering-res').val();


        json_url = `/run_r_script/run_experiment/${runID}/${expTitle}/${uploadName}/${minNFeatureRNA}/${maxNFeatureRNA}/${percentMT}/${nDims}/${clusteringRes}`
        $.get(json_url, function(response) {
        })
            .done(function(response) {
                if (response) {

                    let data = response;
                    //console.log(data);
                    const umap_img_src = data["umap_img_src"];
                    const tsne_img_src = data["tsne_img_src"];
                    const data_rds_href = data["data_rds_href"];

                    //remove previous response result then show new response
                    //('#img-run-exp-umap').attr("src", "");
                    $('#img-run-exp-umap').attr("src", umap_img_src);
                    //$('#img-run-exp-tsne').attr("src", "");
                    $('#img-run-exp-tsne').attr("src", tsne_img_src);
                    //$('#lnk-run-exp-data-rds').attr("href", "");
                    $('#lnk-run-exp-data-rds').attr("href", data_rds_href);

                    //hide spinner
                    $('#accordion-run-exp-spinner').addClass("d-none");

                    //show img
                    $('#accordion-run-exp-container').hide();
                    $('#accordion-run-exp-container').removeClass('d-none')
                    $('#accordion-run-exp-container').fadeIn();


                }
            })
            .fail(function() {
                alert( "Error. Please try again later." );
            })
    });


    $("#lnk-exp-degs").click(function() {
        //show spinner
        $('#find-degs-spinner').removeClass("d-none")



        json_url = `/run_r_script/list_clusters_dea_first/${runID}/${expTitle}`
        $.get(json_url, function(response) {
        })
            .done(function(response) {
                if (response) {

                    //re-init HTML items
                    $('#ddl-find-clusters-degs ul').html("");
                    $('#selected-cluster-wrapper').addClass('d-none');
                    $('#ddl-find-clusters-degs').addClass("d-none");

                    let data = response;
                    //console.log(data);

                    let clusters = data["clusters"] // avaiable clusters
                    let clusterDEGAsList = data["cluster_0_degs_as_list"]; //dict: clusterNum  -> HTML tbody content
                    selectedClusterDEGAsList = clusterDEGAsList;

                    let clusterDEGAsStr = data["cluster_0_degs_as_str"];
                    selectedClusterDEGListAsStr = clusterDEGAsStr;

                    //fill table with DEGs for cluster 0
                    var cluster_0_degs_html = "";
                    let deg_file_rows = clusterDEGAsStr.split("\n");
                    for(var i = 0; i < deg_file_rows.length; i++) {
                        cluster_0_degs_html += "<tr class='text-center'>"
                        let tds = deg_file_rows[i].split(",");
                        for(var j = 0; j < tds.length; j++) {
                            cluster_0_degs_html += "<td>" + tds[j] + "</td>";
                        }
                        cluster_0_degs_html += "</tr>";
                    }

                    //fill table with DEGs for first cluster
                    $('#tbl-cluster-degs tbody').html(cluster_0_degs_html);

                    //hide spinner
                    $('#find-degs-spinner').addClass("d-none");


                    //fill ddl-find-clusters-degs with clusters and set cluster 0 as default
                    var ddlFindClustersDegsHTML = ""
                    for (var i = 0; i < clusters.length; i++) {
                        ddlFindClustersDegsHTML += '<li><a class="dropdown-item find-cluster-degs-item" href="javascript:;" data-cluster="' + clusters[i] + '">' + 'cluster_' + clusters[i] + '</a></li>';
                    }
                    $('#ddl-find-clusters-degs ul').html(ddlFindClustersDegsHTML);
                    $("#ddl-find-clusters-degs").val(clusters[0]);
                    $("#selected-cluster").text("cluster_" + clusters[0]);
                    $('#selected-cluster-wrapper').removeClass('d-none');
                    $('#ddl-find-clusters-degs').removeClass("d-none");

                    //fill table with DEGs for cluster 0
                    $('#tbl-cluster-degs tbody').html(cluster_0_degs_html);
                    $('#tbl-cluster-degs').removeClass("d-none");

                    //hide spinner
                    $('#find-degs-spinner').addClass("d-none");


                }
            })
            .fail(function() {
                alert( "Error. Please try again later." );
            })
    });

    //On selected cluster change: generate and load DEGs table
    $(document).on("click", ".find-cluster-degs-item", function(){


        $('#tbl-cluster-degs tbody').html("<tr class='text-center'><td>..</td><td>..</td><td>..</td><td>..</td></tr>");
        $('#find-degs-spinner').removeClass("d-none");

        selectedCluster = $(this).data("cluster");
        $("#selected-cluster").text("cluster_" + selectedCluster);

        //TODO: remove next line . it's temp
        const expTitle = 'P14_Prime';

        json_url = `/run_r_script/dea_cluster/${runID}/${expTitle}/${selectedCluster}`
        $.get(json_url, function(response) {
        })
            .done(function(response) {
                if (response) {

                    let data = response;
                    //console.log(data);

                    let clusterDEGAsList = data["cluster_degs_as_list"]; //dict: clusterNum  -> HTML tbody content
                    selectedClusterDEGAsList = clusterDEGAsList;

                    let clusterDEGAsStr = data["cluster_degs_as_str"];
                    selectedClusterDEGListAsStr = clusterDEGAsStr;

                    //fill table with DEGs for cluster
                    var cluster_degs_html = "";
                    let deg_file_rows = clusterDEGAsStr.split("\n");
                    for(var i = 0; i < deg_file_rows.length; i++) {
                        cluster_degs_html += "<tr class='text-center'>"
                        let tds = deg_file_rows[i].split(",");
                        for(var j = 0; j < tds.length; j++) {
                            cluster_degs_html += "<td>" + tds[j] + "</td>";
                        }
                        cluster_degs_html += "</tr>";
                    }
                    $('#tbl-cluster-degs tbody').html(cluster_degs_html);

                    //hide spinner
                    $('#find-degs-spinner').addClass("d-none");


                    $('#selected-cluster-wrapper').removeClass('d-none');
                    $('#ddl-find-clusters-degs').removeClass("d-none");

                    //fill table with DEGs for selected cluster
                    $('#tbl-cluster-degs').removeClass("d-none");

                    //hide spinner
                    $('#find-degs-spinner').addClass("d-none");

                }
            })
            .fail(function() {
                alert( "Error. Please try again later." );
            })

   });

    function loadSavedDEGLists() {
        json_url = `/json_lst_saved_deg_lists/${runID}`;
        $.get(json_url, function(response) {
        })
            .done(function(response) {
                if (response) {

                    let data = response;
                    //console.log(data);

                    let savedDEGLists = data["saved_deg_lists"]

                    // A- do it for ddl-navigation
                    //empty ddl-saved-deg-lists as intermediate step
                    $('.saved-deg-lists-item').remove();

                    //1-get template HTML, 2- replace placeholder, 3- add to all over HTML
                    var html = "";

                    for (var i=0; i < savedDEGLists.length; i++) {
                        html += "<li class='saved-deg-lists-item'>";
                        html += $('#saved-deg-lists-item-template').html();
                        html = html.replace("saved_deg_list_id_placeholder", savedDEGLists[i]);
                        let href = `/saved_deg_list/${runID}/${savedDEGLists[i]}`;
                        html = html.replace("saved_deg_list_href_placeholder", href);
                        html += "</li>";
                    }

                    $('#create-deg-list-manually-wrapper').before(html);

                    // B- do it for ddl-cross-with-deg-list inside accordion
                    var html = "";
                    for (var i=0; i < savedDEGLists.length; i++) {
                        html += `<li><a class="dropdown-item" href="javascript:;" data-list-id="${savedDEGLists[i]}">${savedDEGLists[i]}</a></li>`;
                    }
                    $("#ddl-cross-with-deg-list .dropdown-menu").html(html);


                }
            })
            .fail(function() {
                alert( "Error. Please try again later." );
            })
    }

    //save current DEG list as a SAVED DEG LIST
    $('#lnk-save-as-deg-list').click(function () {
        $('#txt-modal-save-deg-list-id').val("deg_lst_" + expTitle + "_cluster_" + selectedCluster);
        saveAsDEGListModal = new bootstrap.Modal(document.getElementById('modal-save-as-deg-list'), {
            keyboard: false
        });
        saveAsDEGListModal.toggle();

    });

    //bind click event to submit button
    $('#modal-save-as-deg-list-submit').click(function () {

        let degListID = $('#txt-modal-save-deg-list-id').val();

        console.log(degListID);

        if (savedDEGLists.includes(degListID) == true) {
            alert("There is already a saved DEG list with the same name. Please choose another name.")

        } else {
            json_url = `/save_deg_list/${runID}/${expTitle}/${selectedCluster}/${degListID}`;

            $.get(json_url, function (response) {
            })
                .done(function (response) {
                    if (response) {

                        let data = response;
                        alert("Saved!")
                        savedDEGLists.push(degListID);
                        loadSavedDEGLists();

                        saveAsDEGListModal.toggle();


                    }
                })
                .fail(function () {
                    alert("Error. Please try again later.");
                })
        }

    });

    //TODO: implement delete
    $('.deg_list_del').click(function () {

        let degListID = $(this).data("deg-list-id");

        if (confirm("Are you sure you want to delete?") == False) {
            alert("Canceled.")
            return;
        }

        json_url = `/workspace/run/${runID}/del_deg_list/${degListID}`;

        $.get(json_url, function (response) {
        })
            .done(function (response) {
                if (response) {

                    let data = response;

                    alert("Deleted!");

                    loadSavedDEGLists();

                }
            })
            .fail(function () {
                alert("Error. Please try again later.");
            })
    });


    //intersect DEG list from current cluster with one or more SAVED DEG LIST(S)
    $('#lnk-intersect-with').click(function () {
        var intersectClusterWithListModal = new bootstrap.Modal(document.getElementById('modal-intersect-cluster-with-list'), {
            keyboard: false
        })
        intersectClusterWithListModal.toggle();
    });

    //start experiment by subsetting a cluster. actually it will only subset and add to raw uploads
    $('#lnk-cluster-exp').click(function () {

        const minNFeatureRNA = $('#txt-min-nfeature-rna').val();
        const maxNFeatureRNA = $('#txt-max-nfeature-rna').val();
        const percentMT = $('#txt-percent-mt').val();
        const nDims = $('#txt-run-exp-n-dims').val();
        const clusteringRes = $('#txt-run-exp-clustering-res').val();

        let data_source_name = uploadName + "_" + expTitle + "_cluster_" + selectedCluster;

        //scroll to data source to set attention
        $([document.documentElement, document.body]).animate({
            scrollTop: $("#workspace-data-sources").offset().top
        }, 1000);

        //use workspace-data-source-template to add data source. Create an item... just remove the spinner when everything is fine.
        var html = $('#workspace-data-source-template').html();
        html = html.replaceAll("data_source_name_placeholder", data_source_name);

        $('#workspace-data-sources').append(html);

        json_url = `/run_r_script/subset_cluster/${runID}/${uploadName}/${expTitle}/${selectedCluster}/${minNFeatureRNA}/${maxNFeatureRNA}/${percentMT}/${nDims}/${clusteringRes}`;

        $.get(json_url, function (response) {
        })
            .done(function (response) {
                if (response) {

                    let data = response;

                    //everything is fine. just remove the spinner ;)
                    $('#workspace-data-source-new-spinner-' + data_source_name).remove();
                    alert("Data source created from cluster!")

                }
            })
            .fail(function () {
                alert("Error. Please try again later.");
            })



    });

    //attach a new data source when click on load
    $('.lnk-load-upload').click(function() {
        let data_folder_name = $(this).data("raw-upload");
        alert("Loaded " + data_folder_name);

        uploadName = data_folder_name

        $('#loaded-data').html(uploadName);
    });


    //save current DEG list as a SAVED DEG LIST
    $('#btn-create-deg-list-manually').click(function () {
        createDEGListManuallyModal = new bootstrap.Modal(document.getElementById('modal-create-deg-list-manually'), {
            keyboard: false
        });
        createDEGListManuallyModal.toggle();
    });

    $('#modal-create-deg-list-manually-submit').click(function () {

        const del_deg_list = $('#txt-modal-create-deg-list-manually-list-id').val();
        const genes = $('#txt-modal-create-deg-list-manually-genes').val();

        json_url = `/create_deg_list_manually/${runID}/${del_deg_list}/${genes}`;

        $.get(json_url, function (response) {
        })
            .done(function (response) {
                if (response) {

                    let data = response;

                    //refresh ddl list
                    loadSavedDEGLists();

                    alert("Created Successfully");

                    createDEGListManuallyModal.toggle();

                }
            })
            .fail(function () {
                alert("Error. Please try again later.");
            })

    });

    $('#ddl-cross-with-marker-db li a').click(function() {
        selectedDBSpecies = $(this).data("species");
        selectedDBPhase = $(this).data("phase");
        alert("DB selected successfully.");
    })

    $('#btn-accordion-cross-with-marker-db-run').click(function() {

        if(selectedDBSpecies === undefined || selectedDBPhase === undefined) {
            alert("Please select DB first.")
            return;
        }

        json_url = `/json_find_cluster_to_db_gene_intersection/${runID}/${selectedCluster}/${selectedDBSpecies}/${selectedDBPhase}`;

        $.get(json_url, function (response) {
        })
            .done(function (response) {
                if (response) {

                    let data = response;
                    //console.log(data)

                    const tbl_html = data["tbl_intersection"]
                    $('#accordion-cross-with-marker-db-container tbody').html(tbl_html);
                }
            })
            .fail(function () {
                alert("Error. Please try again later.");
            })


    });

    $(document).on("click", "#ddl-cross-with-deg-list .dropdown-menu li a", function() {
        selectedSavedDEGList =  $(this).data("list-id");
        alert("Selected successfully!");
    });

    $('#btn-accordion-cross-with-saved-deg-list-run').click(function() {

        if(selectedSavedDEGList === undefined) {
            alert("You have to select a DEG list from the list above.")
            return;
        }

        json_url = `/json_find_exp_to_list_gene_intersection/${runID}/${expTitle}/${selectedSavedDEGList}`;

        $.get(json_url, function (response) {
        })
            .done(function (response) {
                if (response) {

                    let data = response;
                    console.log(data)

                    const dict_intersection_result = data["dict_intersection_result"];

                    //TODO: make available for download
                    var html = "";
                    for (const [key, value] of Object.entries(dict_intersection_result)) {
                        let cluster = key;
                        let genes = value;
                        let nGenes = genes.length;

                        html += "<tr>";
                        html += `<td>${cluster}</td>`;
                        html += `<td>${nGenes}</td>`;
                        html += `<td style="white-space: pre-wrap">${genes.join(", ")}</td>`;
                        html += "</tr>";
                    }
                    $('#accordion-cross-deg-lists-container tbody').html(html);
                }
            })
            .fail(function () {
                alert("Error. Please try again later.");
            })


    });



    //bind on document ready

    //fill Saved DEG Lists ddl with content
    loadSavedDEGLists();

});
