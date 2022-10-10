//tmp
const runID = 1;
const uploadName = "GSM3860733_E10";
const nDims = 50;

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
        const expTitle = $("#txt-run-exp-title").val();
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

});
