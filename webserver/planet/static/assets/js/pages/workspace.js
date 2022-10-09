$(document).ready(function() {
   $(".lnk-run-experiment").click(function() {
        let rawUpload = $(this).data("raw-upload");
        $('#modal-run-experiment').data("raw-upload", rawUpload);

        var runExperimentModal = new bootstrap.Modal(document.getElementById('modal-run-experiment'), {
            keyboard: false
        })
        runExperimentModal.toggle();
   });
});
