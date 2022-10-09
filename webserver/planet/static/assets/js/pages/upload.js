$(document).ready(function() {
    $('#btn-go-to-upload').click(function() {
        $('#upload').hide();
        $('#experiment-settings').fadeOut(function() {
            $('#upload').removeClass("d-none");
            $('#upload').fadeIn();
        });
    });
});
