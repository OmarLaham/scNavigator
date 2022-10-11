//tmp
const runID = 1;
const uploadName = "GSM3860733_E10";
const nDims = 50;
var expTitle = undefined;
var selectedCluster = undefined;

$(document).ready(function() {
   $(".lnk-elbow-plot").click(function() {
        let rawUpload = $(this).data("raw-upload");
        $('#modal-elbow-plot').data("raw-upload", rawUpload);

        var elbowPlotModal = new bootstrap.Modal(document.getElementById('modal-elbow-plot'), {
            keyboard: false
        })
        elbowPlotModal.toggle();
   });

    $(".lnk-run-experiment").click(function() {
        let rawUpload = $(this).data("raw-upload");
        $('#modal-run-experiment').data("raw-upload", rawUpload);

        var runExperimentModal = new bootstrap.Modal(document.getElementById('modal-run-experiment'), {
            keyboard: false
        })
        runExperimentModal.toggle();
   });


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

        //TODO: remove next line . it's temp
        const expTitle = 'P14_Prime';

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
                    let cluster_1_degs_html = data["cluster_1_degs_html"] //dict: clusterNum  -> HTML tbody content

                    //fill table with DEGs for first cluster
                    $('#tbl-cluster-degs tbody').html(cluster_1_degs_html);

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
                    $('#tbl-cluster-degs tbody').html(cluster_1_degs_html);
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

                    let cluster_degs_html = data["cluster_degs_html"] //dict: clusterNum  -> HTML tbody content

                    //fill table with DEGs for cluster 0
                    $('#tbl-cluster-degs tbody').html(cluster_degs_html);

                    //hide spinner
                    $('#find-degs-spinner').addClass("d-none");


                    $('#ddl-find-clusters-degs ul').html(ddlFindClustersDegsHTML);
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

    //save current DEG list as a SAVED DEG LIST
    $('#lnk-save-as-deg-list').click(function () {

    });

    //intersect DEG list from current cluster with one or more SAVED DEG LIST(S)
    $('#lnk-intersect-with').click(function () {

    });
});
