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

        json_url = `/run_r_script/elbow_plot/${runID}/${uploadName}/${minNFeatureRNA}/${maxNFeatureRNA}/${percentMT}/${nDims}`
        $.get(json_url, function(response) {
        })
            .done(function(response) {
                if (response) {

                    let data = response;
                    //console.log(data);
                    const img_src = data["img_src"];

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

    $('#btn-run-exp-run').click(function () {

        //show spinner
        $('#accordion-elbow-spinner').removeClass("d-none")

        const minNFeatureRNA = $('#txt-min-nfeature-rna').val();
        const maxNFeatureRNA = $('#txt-max-nfeature-rna').val();
        const percentMT = $('#txt-percent-mt').val();
        const nDims = $('#txt-elbow-n-dims').val();

        json_url = `/run_r_script/elbow_plot/${runID}/${uploadName}/${minNFeatureRNA}/${maxNFeatureRNA}/${percentMT}/${nDims}`
        $.get(json_url, function(response) {
        })
            .done(function(response) {
                if (response) {

                    let data = response;
                    //console.log(data);
                    const img_src = data["img_src"];

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

});
