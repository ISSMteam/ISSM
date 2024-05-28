/*
A collection of functions that replicate the behavior of MATLAB built-in
functions of the same, respective name.

Where possible, users are encouraged to use native and/or the most efficient 
methods in JavaScript, but we provide these functions as a way to make 
translations from the MATLAB to the JavaScript ISSM API more seamless.

NOTE:
- We cannot implement the following MATLAB built-in functions by name as their
names are reserved keywords in JavaScript,
    class

TODO:
- Implement,
    sort
    unique

Sources:
- https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Lexical_grammar#keywords
*/

/**
 * FUNCTION any - Determine if any array elements are nonzero
 *
 * Replicates behavior of MATLAB's 'any' function.
 *
 * Usage:
 *     B = any(A);
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/ref/any.html
 *
 * NOTE:
 * - Only basic functionality is implemented
 */
function any(A) {//{{{
    for (let i = 0; i < A.length; ++i) {
        if (A[i] !== 0) {
            return true;
        }
    }
    return false;
} //}}}

/**
 * FUNCTION diff - Differences and approximate derivatives
 *
 * Replicates behavior of MATLAB's 'diff' function.
 * 
 * Y = diff(X) calculates differences between adjacent elements of X.
 *
 * Usage:
 *     Y = diff(X);
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/ref/diff.html
 *
 * NOTE:
 * - Not all functionality is implemented
 */
function diff(X) {//{{{
    let diffs = [];
    for (let i = 0; i < (X.length - 1); ++i) {
        diffs[i] = X[i + 1] - X[i];
    }
    return diffs;
} //}}}

/**
 * FUNCTION disp - Display value of variable
 *
 * Replicates behavior of MATLAB's 'disp' function.
 *
 * Output is logged to console.
 *
 * Usage:
 *     disp(X);
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/ref/disp.html
 */
function disp(X) {//{{{
    console.log(X);
} //}}}

/**
 * FUNCTION error - Throw error and display message
 *
 * Replicates behavior of MATLAB's 'error' function
 *
 * Usage:
 *     error(msg);
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/ref/error.html
 *
 * NOTE:
 * - Only basic functionality is implemented
 */
function error(msg) {//{{{
    throw new Error(msg);
} //}}}

/**
 * FUNCTION find - Find indices and values of nonzero elements
 *
 * Replicates behavior of MATLAB's 'find' function
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/ref/find.html
 *
 * TODO:
 * - Implement support for multidimensional arrays
 * - Implement 'n' parameter
 * - Implement 'direction' parameter
 */
function find(X, n, direction) {//{{{
    if (typeof(X[0]) == 'number') {
        let indices = [];
        for (let i = 0; i < X.length; ++i) {
            if (X[i] != 0) {
                indices.push(i);
            }
        }
        return indices;
    } else { //TODO: If 2d array, assume return rows & cols - try to find a way to not always return rows/cols
        let rowindices = [];
        let colindices = [];
        for (let i = 0; i < X.length; ++i) {
            for (let j = 0; j < X[i].length; ++j) {
                if (X[i][j] != 0) {
                    rowindices.push(i);
                    colindices.push(j);
                }
            }
        }
        return [rowindices, colindices];
    }
} //}}}

/**
 * FUNCTION isempty - Determine whether array is empty
 *
 * Replicates behavior of MATLAB's 'isempty' function
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/ref/isempty.html
 *
 * TODO:
 * - Implement support for multidimensional arrays
 */
function isempty(A) {//{{{
    //NOTE: Expanded for clarity/debugging. Can reduce later.
    if (A === undefined) {
        return 1; //TODO: Fix this: for now, treat undefined as empty
    } else {
        return A.length == 0;
    }
} //}}}

/**
 * FUNCTION ismember - Array elements that are members of set array
 *
 * Replicates basic behavior of MATLAB's 'ismember' function with an important 
 * modification: because of the way we use ismember under MATLAB (to determine 
 * if *any* element of A is in B) and because the truthiness of an empty array 
 * in JavaScript is true, we return 0 if the A is not in B.
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/ref/double.ismember.html
 *
 * TODO:
 * - Implement support for multidimensional arrays
 */
function ismember(A, B) {//{{{
    let b = B;
    if (typeof(B) == 'number' || typeof(B) == 'string') {
        b = [B];
    }
    let indices = zeros(A.length);
    for (let i = 0; i < A.length; ++i) {
        for (let j = 0; j < b.length; ++j) {
            if (A[i] === b[j]) {
                indices[i] = 1;
            }
        }
    }
    if (indices.length) {
        return indices;
    } else {
        return 0;
    }
} //}}}

/**
 * FUNCTION isnan - Determine which array elements are NaN
 *
 * Replicates behavior of MATLAB's 'isnan' function
 *
 * Usage:
 *     TF = isnan(A)
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/ref/isnan.html
 * - https://medium.com/coding-in-simple-english/how-to-check-for-nan-in-javascript-4294e555b447
 * 
 * NOTE:
 * - Not to be confused with JavaScript built-in function isNaN
 * - Not to be confused with function Number.isNaN
 */
function isnan(A) {//{{{
    if (A.constructor !== Array) {
        error('isnan: argument must be an array')
    }
    is_nan = [];
    for (let i = 0; i < A.length; ++i) {
        if (Number.isNaN(A[i])) {
            is_nan.push(1);
        } else {
            is_nan.push(0);
        }
    }
    return is_nan;
} //}}}/*

/**
 * FUNCTION length - Length of largest array dimension
 *
 * Replicates behavior of MATLAB's 'length' function
 *
 * Usage:
 *     L = length(X)
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/ref/length.html
 *
 * TODO:
 * - Implement support for multidimensional arrays
 */
function length(A) {//{{{
    return A.length;
} //}}}/*

/**
 * FUNCTION max - Maximum elements of an array
 *
 * Replicates behavior of MATLAB's 'max' function
 *
 * Usage:
 *     M = max(A)
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/ref/max.html
 *
 * TODO:
 * - Implement support for multidimensional arrays
 * - Implement 'nanflag' parameter
 */
function max(A) {//{{{
    return Math.max.apply(null, A);
} //}}}/*

/**
 * FUNCTION min - Minimum elements of an array
 *
 * Replicates behavior of MATLAB's 'min' function
 *
 * Usage:
 *     M = min(A)
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/ref/min.html
 *
 * TODO:
 * - Implement support for multidimensional arrays
 * - Implement 'nanflag' parameter
 */
function min(A) {//{{{
    return Math.min.apply(null, A);
} //}}}/*

/**
 * FUNCTION ones - Create array of all ones
 *
 * Replicates behavior of MATLAB's 'zeros' function
 * 
 * Usage:
 *     X = ones
 *     X = ones(n)
 *     X = ones(sz1,...,szN)
 *     X = ones(sz)
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/ref/ones.html
 *
 * TODO:
 * - Create lower-level function to handle both ones and zeros functions as 
 * they are essentially the same
 * - Implement functionality for more than 2 dimensions
 */
function ones() {//{{{
    nargs = arguments.length;
    if (nargs == 0) {
        return 1;
    } else if (nargs == 1) {
        let arg = arguments[0];
        if (typeof(arg) == 'number') {
            return NewArrayFill(arg, 1);
        } else if (arg.constructor == Array) {
            return ones(...arg); // spread array of sizes
        } else {
            error('ones: functionality for greater than 2 dimensions is not currently implemented')
        }
    } else if (nargs == 2) {
        return NewArrayFill2D(arguments[0], arguments[1], 1);
    } else {
        error('ones: functionality for greater than 2 dimensions is not currently implemented')
    }
} //}}}

/**
 * FUNCTION size - Array size
 *
 * Replicates behavior of MATLAB's 'size' function
 * 
 * Usage:
 *     sz = size(A);
 *     szdim = size(A, dim);
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/ref/size.html
 * 
 * NOTE:
 * - Not all functionality is implemented
 * - In a 2D array, if length of all columns are not equal, unexpected behavior 
 * may result from subsequent processing based on the values returned from this 
 * function.
 */
function size(A, dim) {//{{{
    let nargs = arguments.length;
    if (nargs == 0) {
        error('size: at least one argument is required');
    } else if (nargs == 1) {
        if (typeof(A) == 'number') { // scalar numbers are of size [1, 1]
            return [1, 1];
        } else if (A.constructor != Array) {
            error('size: A argument must be a number or an Array');
        }
        if (typeof(A) == 'undefined' || isNaN(A)) {
            return [0];
        } else if (A.length && A[0].constructor == Array) {
            return [A.length, A[0].length];
        } else {
            return [A.length];
        }
    } else if (nargs == 2) {
        if (typeof(dim) != 'number' || dim < 0 || dim > 1) {
            error('size: dim argument must be a number between 0 and 1, inclusive');
        }
        if (dim == 0) {
            if (typeof(A) == 'number') { // scalar numbers are of size [1, 1]
                return 1;
            } else {
                return A.length;
            }
        } else {
            if (typeof(A) == 'number') { // scalar numbers are of size [1, 1]
                return 1;
            } else if (typeof(A) == 'undefined' || isNaN(A)) {
                return 0;
            } else if (A[0].constructor != Array) {
                error('size: A[0] is not an Array');
            }
            return A[0].length;
        }
    } else {
        error('size: functionality for more than 2 arguments is not currently implemented');
    }
} //}}}

/**
 * FUNCTION strcmpi - Compare strings (case insensitive)
 *
 * Replicates behavior of MATLAB's 'strcmpi' function
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/ref/strcmpi.html
 */
function strcmpi(s1, s2) {//{{{
    return s1.toLowerCase() == s2.toLowerCase();
} //}}}

/**
 * FUNCTION sum - Sum of elements of an array
 *
 * Replicates behavior of MATLAB's 'sum' function
 *
 * Usage:
 *     S = sum(A)
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/ref/sum.html
 *
 * TODO:
 * - Implement support for multidimensional arrays
 * - Implement 'nanflag' parameter
 */
function sum(A) {//{{{
    return ArraySum(A);
} //}}}

/**
 * FUNCTION zeros - Create array of all zeros
 *
 * Replicates behavior of MATLAB's 'zeros' function
 * 
 * Usage:
 *     X = zeros
 *     X = zeros(n)
 *     X = zeros(sz1,...,szN)
 *     X = zeros(sz)
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/ref/zeros.html
 *
 * TODO:
 * - Create lower-level function to handle both ones and zeros functions as 
 * they are essentially the same
 * - Implement functionality for more than 2 dimensions
 */
function zeros() {//{{{
    nargs = arguments.length;
    if (nargs == 0) {
        return 0;
    } else if (nargs == 1) {
        let arg = arguments[0];
        if (typeof(arg) == 'number') {
            return NewArrayFill(arg, 0);
        } else if (arg.constructor == Array) {
            return zeros(...arg); // spread array of sizes
        } else {
            error('zeros: functionality for greater than 2 dimensions is not currently implemented')
        }
    } else if (nargs == 2) {
        return NewArrayFill2D(arguments[0], arguments[1], 0);
    } else {
        error('zeros: functionality for greater than 2 dimensions is not currently implemented')
    }
} //}}}

/**
 * FUNCTIONS sin, cos, tan, asin, acos, atan2 - trig functions that work with radians
 *
 * Replicates behavior of MATLAB's trig functions
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/trigonometry.html
 *
 */
function sin(X) {//{{{
    let result = NewArrayFill(size, X.length);
    for (let i = 0; i < X.length; ++i) {
        result[i] = Math.sin(X[i]);
    }
    return result;
} //}}}

function cos(X) {//{{{
    let result = NewArrayFill(size, X.length);
    for (let i = 0; i < X.length; ++i) {
        result[i] = Math.cos(X[i]);
    }
    return result;
} //}}}

function tan(X) {//{{{
    let result = NewArrayFill(size, X.length);
    for (let i = 0; i < X.length; ++i) {
        result[i] = Math.tan(X[i]);
    }
    return result;
} //}}}

function asin(X) {//{{{
    let result = NewArrayFill(size, X.length);
    for (let i = 0; i < X.length; ++i) {
        result[i] = Math.asin(X[i]);
    }
    return result;
} //}}}

function acos(X) {//{{{
    let result = NewArrayFill(size, X.length);
    for (let i = 0; i < X.length; ++i) {
        result[i] = Math.acos(X[i]);
    }
    return result;
} //}}}

// TODO: Test that the arguments do not need to be reversed, as they are in MATLAB and Python
function atan2(X, Y) {//{{{
    let result = NewArrayFill(size, X.length);
    for (let i = 0; i < X.length; ++i) {
        result[i] = Math.atan2(X[i], Y[i]);
    }
    return result;
} //}}}

/**
 * FUNCTIONS sind, cosd, tand, asind, acosd, atan2d - trig functions that work with degrees
 *
 * Replicates behavior of MATLAB's trig functions
 *
 * Sources:
 * - https://www.mathworks.com/help/matlab/trigonometry.html
 *
 */
function sind(X) {//{{{
    let result = NewArrayFill(size, X.length);
    for (let i = 0; i < X.length; ++i) {
        result[i] = Math.sin(X[i] * DEG2RAD);
    }
    return result;
} //}}}

function cosd(X) {//{{{
    let result = NewArrayFill(size, X.length);
    for (let i = 0; i < X.length; ++i) {
        result[i] = Math.cos(X[i] * DEG2RAD);
    }
    return result;
} //}}}

function tand(X) {//{{{
    let result = NewArrayFill(size, X.length);
    for (let i = 0; i < X.length; ++i) {
        result[i] = Math.tan(X[i] * DEG2RAD);
    }
    return result;
} //}}}

function asind(X) {//{{{
    let result = NewArrayFill(size, X.length);
    for (let i = 0; i < X.length; ++i) {
        result[i] = RAD2DEG * Math.asin(X[i]);
    }
    return result;
} //}}}

function acosd(X) {//{{{
    let result = NewArrayFill(size, X.length);
    for (let i = 0; i < X.length; ++i) {
        result[i] = RAD2DEG * Math.acos(X[i]);
    }
    return result;
} //}}}

// TODO: Test that the arguments do not need to be reversed, as they are in MATLAB and Python
function atan2d(X, Y) {//{{{
    let result = NewArrayFill(size, X.length);
    for (let i = 0; i < X.length; ++i) {
        result[i] = RAD2DEG * Math.atan2(X[i], Y[i]);
    }
    return result;
} //}}}
