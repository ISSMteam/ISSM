// Karma configuration
module.exports = function(config) {
  config.set({

    basePath: './', 

    // frameworks to use
    // available frameworks: https://npmjs.org/browse/keyword/karma-adapter
    frameworks: ['jasmine-jquery', 'jasmine'],


    // list of files / patterns to load in the browser
    files: [
      './../../../externalpackages/emscripten/src/node/4.1.1_64bit/lib/node_modules/jquery/dist/jquery.min.js',
      './../../../externalpackages/emscripten/src/node/4.1.1_64bit/lib/node_modules/mathjs/dist/math.min.js',
      'lib/Exp/Square.js',
      'lib/Par/SquareShelfConstrained.js',
      'lib/Data/SquareShelfConstrained.data.js',
      'scripts/specs/issm.spec.js',
      './../../../bin/issm-bin.js',
      './../../../bin/issm-prebin.js',
      './../../../bin/IssmModule.js'
      //'scripts/specs/temp.spec.js'
      //'scripts/specs/3.spec.js'
    ],

    // list of files to exclude
    exclude: [
    ],


    // preprocess matching files before serving them to the browser
    // available preprocessors: https://npmjs.org/browse/keyword/karma-preprocessor
    preprocessors: {
    },


    // test results reporter to use
    // possible values: 'dots', 'progress'
    // available reporters: https://npmjs.org/browse/keyword/karma-reporter
    reporters: ['dots', 'junit'],
    junitReporter: {
        outputFile: 'test-results.xml'
    },


    // web server port
    port: 9876,


    // enable / disable colors in the output (reporters and logs)
    colors: true,


    // level of logging
    // possible values: config.LOG_DISABLE || config.LOG_ERROR || config.LOG_WARN || config.LOG_INFO || config.LOG_DEBUG
    logLevel: config.LOG_INFO,


    // enable / disable watching file and executing tests whenever any file changes
    autoWatch: false,


    // start these browsers
    // available browser launchers: https://npmjs.org/browse/keyword/karma-launcher
    //browsers: ['Chrome'],
    browsers: ['Firefox'],


    // Continuous Integration mode
    // if true, Karma captures browsers, runs the tests and exits
    singleRun: true,

    // Concurrency level
    // how many browser should be started simultaneous
    concurrency: Infinity
  });
};
