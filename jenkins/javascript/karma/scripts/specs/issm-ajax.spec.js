window.jasmine.DEFAULT_TIMEOUT_INTERVAL = 10000; // Change timeout for Jasmine tests
var AJAX_TIMEOUT = 5000;

function onAjaxSuccess(data, status, jqxhr) {
    console.log("Success");
}

function onAjaxError(jqxhr, status, err) {
    console.log("Unexpected error: " + err);
}

describe("test101", function() {
    it("contains test101", function(done) {
        $.ajax({
            url: 'http://localhost:8080/test101.js',
            dataType: 'script',
            success: onAjaxSuccess,
            error: onAjaxError,
            complete: function(jqxhr, status) { done(); },
            timeout: AJAX_TIMEOUT
        });
    });
});
