# Run a file server for serving the tests
PORT="${1:-8081}"
python2.7 -m SimpleHTTPServer $PORT
