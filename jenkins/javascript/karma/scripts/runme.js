$(document).ready(function() {
    $.ajaxSetup({
        cache: true
    });

    var PORT=8080;
    var tests = [];

    for (var i = 101; i < 109; ++i) {
        tests.push(i);
    }

    $.each(tests, function(i, v) {
        $.getScript('http://localhost:'+PORT+'/test'+v+'.js', function(src, status) {
            console.log('='.repeat(30));
        });
    });
});
