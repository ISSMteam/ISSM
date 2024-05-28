$(document).ready(function() {
    $.ajaxSetup({
        cache: true
    });

    var PORT=80;
    var tests = []; //last test to check

    // Adding translated tests to array. Hardcoded for now.
    for (var i = 101; i < 109; ++i) {
        tests.push(i);
    }

    for (var i = 111; i < 115; ++i) {
        tests.push(i);
    }

    for (var i = 201; i < 209; ++i) {
        tests.push(i);
    }


    $.each(tests, function(i, v) {
        var btn = $('<input type="button" class="btn btn-primary" value="test'+v+'" id="test'+v+'"/>');
        $(".btn-container").append(btn);
    });

    $("[id^=test]").each(function () {
        $(this).click(function () {
            $('#debug').empty();
            var id = this.id.replace(/[^\d]/g, ''); 
            $.getScript('../../../test/NightlyRun/test'+id+'.js', function(src, status) {
                console.log('='.repeat(30));
                console.log('Status: ' + status);
                console.log('Script executed: test' + id + '.js');
            });
        });
    });

    //Show console log in div
    if (typeof console  != "undefined") 
        if (typeof console.log != 'undefined')
            console.olog = console.log;
    else
        console.olog = function() {};

    console.log = function(message) {
        console.olog(message);
        $('#debug').append('<p>' + message + '</p>');
    };

    console.error = console.debug = console.info =  console.log;
});
