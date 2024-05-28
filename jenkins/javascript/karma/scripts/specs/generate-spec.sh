#!/bin/sh

if [ -z $1 ]; then
        echo "Usage: $0 [NUMBER[-NUMBER]][,NUMBER[-NUMBER]]..."
        exit 1
fi

# Create array to store numbers of the tests
TESTNUMBERS=()

# Location of test scripts
SCRIPTDIR=$ISSM_DIR/jenkins/javascript/karma/scripts

OLDIFS=$IFS
IFS=, # Overwrite the in-field-separator to detect ranges of numbers as numbers separated by commas

# Add the test numbers to the array
for range in $1; do
    if [[ ! "$range" =~ "-" ]]; then # check if it is a range of numbers or just a single number
        if [ ! -f "$SCRIPTDIR/test$range.js" ]; then
            >&2 echo "Warning: test$num.js does not exist."
        else
            TESTNUMBERS+=($range)
        fi
    else
        SEQUENCE=($(seq -w -s ' ' $(sed "s/-/$IFS/" <<< ${range})))

        IFS=' '
        for num in ${SEQUENCE[@]}; do
            if [ ! -f "$SCRIPTDIR/test$num.js" ]; then
                >&2 echo "Warning: test$num.js does not exist."
            fi
        done

        TESTNUMBERS=("${TESTNUMBERS[@]}" "${SEQUENCE[@]}")
    fi
done

IFS=$OLDIFS # Reset the in-field-separator

# Include necessary functions and constants
cat << EOF
window.jasmine.DEFAULT_TIMEOUT_INTERVAL = 10000; // Change timeout for Jasmine tests
var AJAX_TIMEOUT = 5000;
function onAjaxSuccess(data, status, jqxhr) {
    console.log("Success");
}
function onAjaxError(jqxhr, status, err) {
    console.log("Unexpected error: " + err);
}
EOF

# Create stubs for individual tests
for num in ${TESTNUMBERS[@]}; do
cat  << EOF
describe("test$num", function() {
    it("contains test$num", function(done) {
        $.ajax({
            url: 'http://localhost:8080/test$num.js',
            dataType: 'script',
            success: onAjaxSuccess,
            error: onAjaxError,
            complete: function(jqxhr, status) { done(); },
            timeout: AJAX_TIMEOUT
        });
    });
});
EOF
done
