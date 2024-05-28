function ArrayMax(array){ //{{{
	//Calculate array max using for loop instead of Math.max.apply(null, array) to avoid recursive stack overflow in mobile browsers
	
	var max=-Infinity;
	
	for (var i=0;i<array.length; i++) {
		max=Math.max(max,array[i]);
	}
	
	return max;
} //}}}
function ArrayMax2D(array){ //{{{
	
	var max=-Infinity;

	for (var i=0;i<array.length;i++){
		var subarray=array[i];
		max=Math.max(max,ArrayMax(subarray));
	}

	return max;
} //}}}
function ArrayMin(array){ //{{{
	//Calculate array min using for loop instead of Math.min.apply(null, array) to avoid recursive stack overflow in mobile browsers
	
	var min=Infinity;
	
	for (var i=0;i<array.length; i++) {
		min=Math.min(min,array[i]);
	}
	
	return min;
} //}}}
function ArrayMin2D(array){ //{{{
	
	var min=Infinity;

	for (var i=0;i<array.length;i++){
		var subarray=array[i];
		min=Math.min(min,ArrayMin(subarray));
	}

	return min;
} //}}}
function ArraySum(array){ //{{{
	var sum=0;
	for(var i=0;i<array.length;i++)sum+=array[i];
	return sum;
} //}}}
function ArrayAdd(){ //{{{
    //Takes in any number of scalars or arrays, and calculates the sum. Scalars are treated as similar length arrays of the scalar.
    //Determine reference array and size
    var size, array, arg, initial;
	for (var a = 0; a < arguments.length; a++) {
        arg = arguments[a];
        if (arg instanceof Array || arg instanceof Float64Array) {
            size = arg.length;
            array = arg;
			initial = a;
            break;
        }
    }
	//check internal consistency of arrays provided!: 
	for (var a = 0; a < arguments.length; a++) {
        arg = arguments[a];
        if ((arg instanceof Array || arg instanceof Float64Array) && arg.length != size) {
            throw Error("ArrayAdd error message: arrays provided as arguments are not of the same length!");
        } else if (!(arg instanceof Array || arg instanceof Float64Array) && typeof arg != 'number') {
            throw Error("ArrayAdd error message: arguments provided are not of the type Array or Number!");
        }
	}
	//do the result:
	var result = array.slice(0);
	for (var a = 0; a < arguments.length; a++) {
		if (a != initial) {
			arg = arguments[a];
			if (arg instanceof Array || arg instanceof Float64Array) {
				for(var i = 0; i < result.length; i++){
					result[i] += arg[i];
				}
			} else if (typeof arg === 'number') {
				for(var i = 0; i < result.length; i++){
					result[i] += arg;
				}
			}
        }
	}
	return result;
} //}}}
function ArraySubtract(){ //{{{
    //Takes in any number of scalars or arrays, and calculates the subtraction. Scalars are treated as similar length arrays of the scalar.
    //Determine reference array and size
    var size, array, arg, initial;
	for (var a = 0; a < arguments.length; a++) {
        arg = arguments[a];
        if (arg instanceof Array || arg instanceof Float64Array) {
            size = arg.length;
            array = arg;
			initial = a;
            break;
        }
    }
	//check internal consistency of arrays provided!: 
	for (var a = 0; a < arguments.length; a++) {
        arg = arguments[a];
        if ((arg instanceof Array || arg instanceof Float64Array) && arg.length != size) {
            throw Error("ArrayAdd error message: arrays provided as arguments are not of the same length!");
        } else if (!(arg instanceof Array || arg instanceof Float64Array) && typeof arg != 'number') {
            throw Error("ArrayAdd error message: arguments provided are not of the type Array or Number!");
        }
	}
	//calculate the result, using the first argument to initialize:
	var result = array.slice(0);
	for (var a = 0; a < arguments.length; a++) {
		if (a !== initial) {
			arg = arguments[a];
			if (arg instanceof Array || arg instanceof Float64Array) {
				for(var i = 0; i < result.length; i++){
					result[i] -= arg[i];
				}
			} else if (typeof arg === 'number') {
				for(var i = 0; i < result.length; i++){
					result[i] -= arg;
				}
			}
        }
	}
	return result;
} //}}}
function ArraySubtract2D(){ //{{{
    //Takes in any number of scalars or arrays, and calculates the subtraction. Scalars are treated as similar length arrays of the scalar.
    //Determine reference array and size
    var size, array, arg, initial;
	for (var a = 0; a < arguments.length; a++) {
        arg = arguments[a];
        if (arg instanceof Array || arg instanceof Float64Array) {
            size = arg.length;
            array = arg;
			initial = a;
            break;
        }
    }
	//check internal consistency of arrays provided!: 
	for (var a = 0; a < arguments.length; a++) {
        arg = arguments[a];
        if ((arg instanceof Array || arg instanceof Float64Array) && arg.length != size) {
            throw Error("ArrayAdd error message: arrays provided as arguments are not of the same length!");
        } else if (!(arg instanceof Array || arg instanceof Float64Array) && typeof arg != 'number') {
            throw Error("ArrayAdd error message: arguments provided are not of the type Array or Number!");
        }
	}
	//calculate the result, using the first argument to initialize:
	var result = [];
	for (var a = 0; a < arguments.length; a++) {
		if (a !== initial) {
			arg = arguments[a];
			if (arg instanceof Array || arg instanceof Float64Array) {
				for(var i = 0; i < array.length; i++){
					result[i] = [];
					for(var j = 0; j < array[i].length; j++){
					    result[i][j] = array[i][j] - arg[i][j];
					}
				}
			} else if (typeof arg === 'number') {
				for(var i = 0; i < array.length; i++){
					result[i] = [];
					for(var j = 0; j < array[i].length; j++){
					    result[i][j] = array[i][j] - arg;
					}
				}
			}
        }
	}
	return result;
} //}}}
function ArrayMultiply(){ //{{{
    //Takes in any number of scalars or arrays, and calculates the product. Scalars are treated as similar length arrays of the scalar.
    //Determine reference array and size
    var size, array, arg, initial;
	for (var a = 0; a < arguments.length; a++) {
        arg = arguments[a];
        if (arg instanceof Array || arg instanceof Float64Array) {
            size = arg.length;
            array = arg;
			initial = a;
            break;
        }
    }
	//check internal consistency of arrays provided!: 
	for (var a = 0; a < arguments.length; a++) {
        arg = arguments[a];
        if ((arg instanceof Array || arg instanceof Float64Array) && arg.length != size) {
            throw Error("ArrayAdd error message: arrays provided as arguments are not of the same length!");
        } else if (!(arg instanceof Array || arg instanceof Float64Array) && typeof arg != 'number') {
            throw Error("ArrayAdd error message: arguments provided are not of the type Array or Number!");
        }
	}
	//do the result:
	var result = array.slice(0);
	for (var a = 0; a < arguments.length; a++) {
		if (a !== initial) {
			arg = arguments[a];
			if (arg instanceof Array || arg instanceof Float64Array) {
				for(var i = 0; i < result.length; i++){
					result[i] *= arg[i];
				}
			} else if (typeof arg === 'number') {
				for(var i = 0; i < result.length; i++){
					result[i] *= arg;
				}
			}
		}
	}
	return result;
} //}}}
function ArrayDivide(){ //{{{
    //Takes in any number of scalars or arrays, and calculates the quotient. Scalars are treated as similar length arrays of the scalar.
    //Determine reference array and size
    var size, array, arg, initial;
	for (var a = 0; a < arguments.length; a++) {
        arg = arguments[a];
        if (arg instanceof Array || arg instanceof Float64Array) {
            size = arg.length;
            array = arg;
			initial = a;
            break;
        }
    }
	//check internal consistency of arrays provided!: 
	for (var a = 0; a < arguments.length; a++) {
        arg = arguments[a];
        if ((arg instanceof Array || arg instanceof Float64Array) && arg.length != size) {
            throw Error("ArrayAdd error message: arrays provided as arguments are not of the same length!");
        } else if (!(arg instanceof Array || arg instanceof Float64Array) && typeof arg != 'number') {
            throw Error("ArrayAdd error message: arguments provided are not of the type Array or Number!");
        }
	}
	//calculate the result, using the first argument to initialize:
	var result = array.slice(0);
	for (var a = 0; a < arguments.length; a++) {
		if (a !== initial) {
			arg = arguments[a];
			if (arg instanceof Array || arg instanceof Float64Array) {
				for(var i = 0; i < result.length; i++){
					result[i] /= arg[i];
				}
			} else if (typeof arg === 'number') {
				for(var i = 0; i < result.length; i++){
					result[i] /= arg;
				}
			}
		}
	}
	return result;
} //}}}
function ArrayMean(array){ //{{{
    //Calculate the mean:
    var sum = ArraySum(array);
    return sum / array.length;

} //}}}
function ArraySTD(array){ //{{{
    //Calculate the standard deviation:
    var sum = ArraySum(array);
    var differences = 0;
    for(var i = 0; i < array.length; i++){
        differences += (array[i] - sum) ** 2;
    }
    var variance = differences / array.length;
    return Math.sqrt(variance);

} //}}}
function ArrayLerp(array1, array2, weight1, weight2, alpha){ //{{{
    //Returns linear combination of array1 and array2, based on lerp coefficients weight1, weight2, and alpha.
    //Parameters:
    //  array1  - array 
    //  array2  - array 
    //  weight1 - value associated with array1.
    //  weight2 - value associated with array2
    //  alpha   - value associated with output array
    //Output:
    //  result - linear combination of array1 and array2 based on the equation array1 * alpha + 
    
    var range = weight2 - weight1;
    var normalizedAlpha = (alpha - weight1) / range;
    var result = array1.slice();
    for(var i = 0; i < array1.length; i++){
        result[i]  = array1[i] * normalizedAlpha + array2[i] * (1 - normalizedAlpha);
    }
    return result;

} //}}}
function ArrayTranspose(array){ //{{{
    //Transposes 2d array.
    var rows = array.length;
    var cols = array[0].length;
    var result = Create2DArray(cols, rows);

    for(var i = 0; i < rows; i++) {
        for(var j = 0; j < cols; j++) {
            result[j][i]  = array[i][j];
        }
    }
    return result;

} //}}}
function ArrayXPY(){ //{{{
    if (arguments.length<2)throw Error("ArrayXPY error message: sum has to be for at least two arrays!");

	//check internal consistency of arrays provided!: 
	var firstarray=arguments[0];
	var firstsize=firstarray.length;
	
	for(var a=1;a<arguments.length;a++){
		var array=arguments[a];
		if(array.length!=firstsize)throw Error("ArrayXPY error message: arrays provided as arguments are not of the same length!");
	}

	//do the sum:
	var sum=NewArrayFill(firstsize,0);
	for(var a=0;a<arguments.length;a++){
		var array=arguments[a];
		for(var i=0;i<array.length;i++){
			sum[i]+=array[i];
		}
	}
	return sum;

} //}}}
function ArrayConcat(a, b) { //{{{
	// Make sure that both typed arrays are of the same type
	if(Object.prototype.toString.call(a) !== Object.prototype.toString.call(b)) {
		throw 'The types of the two arguments passed for parameters a and b do not match.';
	}
	var c;
	if(a instanceof Float64Array) {
		c = new a.constructor(a.length + b.length);
		c.set(a);
		c.set(b, a.length);
	}
	else {
		c = clone(a).concat(clone(b));
	}	
	return c;
} //}}}
function ArrayCol(matrix, cols) { //{{{
    var columns = [];
    if (cols instanceof Array || cols instanceof Float64Array) {
        for (var i = 0; i < matrix.length; i++){
            var col = [];
            for (var j = 0; j < cols.length; j++){
                col.push(matrix[i][cols[j]]);
            }
            columns.push(col); 
        }

    } else if (typeof cols == 'number') {
        for (var i = 0; i < matrix.length; i++){
            columns.push(matrix[i][cols]);
        }
    } else {
        throw new Error("ArrayCol error: cols must be a single integer or an array with 2 integers!");
    }
   return columns;
} //}}}
function ListToMatrix(list, elementsPerSubArray) { //{{{
	var matrix = [], i, k;

	for (i = 0, k = -1; i < list.length; i++) {
		if (i % elementsPerSubArray === 0) {
			k++;
			matrix[k] = [];
		}

		matrix[k].push(list[i]);
	}

	return matrix;
} //}}}
function MatrixToList(matrixin) { //{{{

	var matrix=matrixin;

	if (!IsArray(matrix[0])) return matrix;
	else{
		var width = matrix[0].length;
		var length = matrix.length;
		var list= new Array(width*length);

		for(var i=0;i<length;i++){
			for(var j=0;j<width;j++){
				list[i*width+j]=matrix[i][j];
			}
		}
		return list;
	}
} //}}}
function IsArray(object) { //{{{

	var type=Object.prototype.toString.call( object );
	if( type === '[object Array]' ) return 1;
	if( type === '[object Float64Array]' ) return 1;
	if( type === '[object Float32Array]' ) return 1;
	if( type === '[object Int32Array]' ) return 1;
	if( type === '[object Int16Array]' ) return 1;
	if( type === '[object Uint32Array]' ) return 1;
	if( type === '[object Uint16Array]' ) return 1;
	if( type === '[object Uint8Array]' ) return 1;
	return 0;

} //}}}
function ArrayNot(array) { //{{{

	var notarray=array.slice();
	for (var i=0;i<array.length;i++)notarray[i]=-array[i];
	return notarray;
} //}}}
function ArrayFlip(array) { //{{{
	
	var notarray=array.slice();
	for (var i=0;i<array.length;i++)notarray[i]=array[i]^1;
	return notarray;
} //}}}
function ArrayCopy(array) { //{{{

	var copy=[];
	for(var i=0;i<array.length;i++)copy[i]=array[i];
	return copy;
} //}}}
function ArrayPow(array,coefficient) { //{{{

	var powarray=array.slice();
	for (var i=0;i<array.length;i++)powarray[i]=Math.pow(array[i],coefficient);
	return powarray;
} //}}}
function ArraySqrt(array) { //{{{

	var sqrtarray=array.slice();
	for (var i=0;i<array.length;i++)sqrtarray[i]=Math.sqrt(array[i]);
	return sqrtarray;
} //}}}
function ArrayScale(array,alpha) { //{{{

	var scalearray=array.slice();
	for (var i=0;i<array.length;i++)scalearray[i]=array[i]*alpha;
	return scalearray;
} //}}}
function ArrayMag(array1,array2) { //{{{

	var arraymag=NewArrayFill(array1.length,0);
	for (var i=0;i<array1.length;i++)arraymag[i]=Math.sqrt(Math.pow(array1[i],2)+Math.pow(array2[i],2));
	return arraymag;
} //}}}
function ArrayAnyNaN(array) { //{{{

    if(IsArray(array[0])){
        for(var i=0;i<array.length;i++){
            for(var j=0;j<array[0].length;j++){
                if (isNaN(array[i][j])) return 1;
            }
        }
    }
    else{
        for(var i=0;i<array.length;i++){
            if (isNaN(array[i])) return 1;
        }
    }
    return 0;
} //}}}
function ArrayUnique(arr,rows) { //{{{
	if (arguments.length == 2){
		if (rows == 'rows') {
			//See Matlab unique function and https://stackoverflow.com/a/20339709/1905613
			let equals = (a, b) => JSON.stringify(a) === JSON.stringify(b);
			let uniques = [];
			let indexA = [];
			let indexC = [];
			let itemsFound = {};;
			for(let i = 0, l = arr.length; i < l; i++) {
				let stringified = JSON.stringify(arr[i]);
				if (typeof(itemsFound[stringified]) != 'undefined') {
					indexC.push(itemsFound[stringified]);
					continue;
				}
				uniques.push(arr[i]);
				indexA.push(i);
				itemsFound[stringified] = uniques.length-1;
				indexC.push(itemsFound[stringified]);
			}
			//assert arr == uniques[indexC,:];
			for (let i = 0; i < indexC.length; i++) {
				if (!equals(arr[i], uniques[indexC[i]])) {
					throw new Error('bad implementation');	
				}
			}
			//assert uniques == arr[indexA, :];
			for (let i = 0; i < indexA.length; i++) {
				if (!equals(uniques[i], arr[indexA[i]])) {
					throw new Error('bad implementation');	
				}
			}
			let [uniquesSorted, indexInToOut, indexOutToIn] = ArraySortWithIndices(uniques);
			//indexMapping is the index of the edge in the old array
			indexCSorted = []; //indexC.length == arr.length
			//assert uniquesSorted[i,:] = uniques[indexInToOut[i],:]
			for (let i = 0; i < indexInToOut.length; i++) {
				if (!equals(uniquesSorted[i], uniques[indexInToOut[i]])) {
					console.log(i, uniquesSorted[indexInToOut[i]], uniques[i]);
					throw new Error('bad implementation');	
				}
			}
			//assert uniques[i,:] = uniquesSorted[indexOutToIn[i],:]
			for (let i = 0; i < indexOutToIn.length; i++) {
				if (!equals(uniques[i], uniquesSorted[indexOutToIn[i]])) {
					console.log(i, uniques[indexOutToIn[i]], uniquesSorted[i]);
					throw new Error('bad implementation');	
				}
			}
			//GOAL: assert arr[i,:] == uniquesSorted[indexCSorted[i], :]
			//GIVEN: assert arr[i,:] == uniques[indexC[i],:];
			//GIVEN: assert uniquesSorted[i,:] = uniques[indexInToOut[i],:]
			//GIVEN: assert uniques[i,:] = uniquesSorted[indexOutToIn[i],:]
			//assert uniques[indexC[i],:] == uniquesSorted[indexOutToIn[indexC[i]],:]
			//assert uniquesSorted[indexCSorted[i],:]; == uniquesSorted[indexOutToIn[indexC[i]],:];
			for (let i = 0; i < arr.length; i++) {
				indexCSorted[i] = indexOutToIn[indexC[i]];
			}
			for (let i = 0; i < indexC.length; i++) {
				if (!equals(arr[i], uniquesSorted[indexCSorted[i]])) {
					console.log(i, arr[i], uniquesSorted[indexCSorted[i]]);
					throw new Error('bad implementation');	
				}
			}

			indexASorted = []; //indexA.length == uniques.length
			//GOAL: uniquesSorted[i, :] == arr[indexASorted[i], :]
			//GIVEN: assert arr[i,:] == uniques[indexC[i],:];
			//GIVEN: assert arr[indexA[i],:] == uniques[i,:];
			//GIVEN: assert uniques[indexInToOut[i],:] == uniquesSorted[i,:]
			//GIVEN: assert uniques[i,:] = uniquesSorted[indexOutToIn[i],:]
			//assert uniquesSorted[i,:] == uniques[indexInToOut[i],:]
			//assert uniques[indexInToOut[i],:] == arr[indexA[indexInToOut[i]],:];
			//assert indexA[indexInToOut] == indexASorted
			//indexASorted == indexA[indexMapping[i]]
			for (let i = 0; i < indexA.length; i++) {
				indexASorted[i] = indexA[indexInToOut[i]];
			}
			//assert uniques == arr[indexA, :];
			for (let i = 0; i < indexASorted.length; i++) {
				if (!equals(uniquesSorted[i], arr[indexASorted[i]])) {
					throw new Error('bad implementation');	
				}
			}
			console.log('Good uniques');
			return [uniquesSorted, indexASorted, indexCSorted];
		} else {
			throw new Error('ArrayUnique non "rows" not supported');	
		}
	} else {
		return arr.reverse().filter(function (e, i, arr) {
				return arr.indexOf(e, i+1) === -1;
		}).reverse();
	}
} //}}}
function ArraySortWithIndices(toSort, sortingFunction) { //{{{
	//returns the sorted and index such that toSort[index[i]] == sorted[i]
	let toSortCopy = [];
	for (var i = 0; i < toSort.length; i++) {
	    toSortCopy[i] = [toSort[i], i];
	}
	if (typeof(sortingFunction) == 'undefined') {
		let numeric2DFunction = function(a, b) {
			if (a[0][0] == b[0][0]) {
                return a[0][1] - b[0][1];
			} else {
			    return a[0][0] - b[0][0];
			}
		};
		sortingFunction = numeric2DFunction;
	}
	toSortCopy.sort(sortingFunction);
	let indicesInToOut = [];
	let indicesOutToIn = [];
	let sorted = [];
	for (var j = 0; j < toSortCopy.length; j++) {
	    indicesInToOut[j] = toSortCopy[j][1];
	    indicesOutToIn[toSortCopy[j][1]] = j;
	    sorted[j] = toSortCopy[j][0];
	}
	return [sorted, indicesInToOut, indicesOutToIn];
} //}}}
function ArraySort(array,dim) { //{{{
	let numericFunction = function(a, b) {
	    return a - b;
	};
	let numeric2DFunction = function(a, b) {
	    return a[0] - b[0];
	};
	if (arguments.length == 2){
		if (dim == 1) {
			array.sort(numeric2DFunction);
		} else if (dim == 2) {
			for (let i = 0; i < array.length; i++) {
				array[i].sort(numericFunction);
			}
		} else {
			throw new Error('ArraySort dim > 2 not yet supported')
		}
		return array;
	} else {
		return array.sort(numericFunction);
	}
} //}}}
function ArrayRange(lower, upper) { //{{{

    var range = upper - lower + 1;
    return Array.apply(null, Array(range)).map(function (val, ind) {return ind + lower;});

} //}}}
function ArrayAny(array) { //{{{
    //Emulates Matlab 'any' function
    
    for(var i=0;i<array.length;i++){
        if (array[i]!=0)return 1;
    }
    return 0;
} //}}}
function ArrayAnyEqual(array,value) { //{{{
	
	if(!isNaN(value)){
		for(var i=0;i<array.length;i++){
			if (array[i]==value)return 1;
		}
	}
	else{
		for(var i=0;i<array.length;i++){
			if (isNaN(array[i]))return 1;
		}
	}
	return 0;
} //}}}
function ArrayAnyBelowOrEqual(array,value) { //{{{

	for(var i=0;i<array.length;i++){
		if (array[i]<=value)return 1;
	}
	return 0;
} //}}}
function ArrayAnyBelowStrict(array,value) { //{{{

	for(var i=0;i<array.length;i++){
		if (array[i]<value)return 1;
	}
	return 0;
} //}}}
function ArrayAnyAboveOrEqual(array,value) { //{{{

	for(var i=0;i<array.length;i++){
		if (array[i]>=value)return 1;
	}
	return 0;
} //}}}
function ArrayAnyAboveStrict(array,value) { //{{{

	for(var i=0;i<array.length;i++){
		if (array[i]>value)return 1;
	}
	return 0;
} //}}}
function ArrayAnd(array1,array2) { //{{{
	var array = new Array(array1.length);
	for (var i=0;i<array1.length;i++) {
		array[i]=array1[i] & array2[i];
	}
	return array;
} //}}}
function ArrayOr(array1,array2) { //{{{
	var array = new Array(array1.length);
	for (var i=0;i<array1.length;i++) {
		array[i]=array1[i] | array2[i];
	}
	return array;
} //}}}
function ArrayEqual(array1,array2) { //{{{
	var array = new Array(array1.length);

	if (typeof(array1[0]) == 'number') {
		if (typeof(array2) == 'number') {
			for(var i=0;i<array1.length;i++){
				array[i] = array1[i] == array2;
			}
		} else {
			for(var i=0;i<array1.length;i++){
				array[i] = array1[i] == array2[i];
			}
		}
	} else { //provide support for 2d arrays
		if (typeof(array2) == 'number') {
			for(var i=0;i<array1.length;i++){
				array[i] = new Array(array1[i].length);
				for(var j=0;j<array1[i].length;j++){
					array[i][j] = array1[i][j] == array2;
				}
			}
		} else {
			for(var i=0;i<array1.length;i++){
				array[i] = new Array(array1[i].length);
				for(var j=0;j<array1[i].length;j++){
					array[i][j] = array1[i][j] == array2[i][j];
				}
			}
		}
	}
	return array;
} //}}}
function ArrayLessThan(array1,array2) { //{{{
	var array = new Array(array1.length);

	if (typeof(array2) == 'number') {
		for(var i=0;i<array1.length;i++){
			array[i] = array1[i] < array2;
		}
	} else {
		for(var i=0;i<array1.length;i++){
			array[i] = array1[i] < array2[i];
		}
	}
	return array;
} //}}}
function ArrayGreaterThan(array1,array2) { //{{{
	var array = new Array(array1.length);
	if (typeof(array2) == 'number') {
		for(var i=0;i<array1.length;i++){
			array[i] = array1[i] > array2;
		}
	} else {
		for(var i=0;i<array1.length;i++){
			array[i] = array1[i] > array2[i];
		}
	}
	return array;
} //}}}
function ArrayLessEqualThan(array1,array2) { //{{{
	var array = new Array(array1.length);
	if (typeof(array2) == 'number') {
		for(var i=0;i<array1.length;i++){
			array[i] = array1[i] <= array2;
		}
	} else {
		for(var i=0;i<array1.length;i++){
			array[i] = array1[i] <= array2[i];
		}
	}
	return array;
} //}}}
function ArrayGreaterEqualThan(array1,array2) { //{{{
	var array = new Array(array1.length);
	if (typeof(array2) == 'number') {
		for(var i=0;i<array1.length;i++){
			array[i] = array1[i] >= array2;
		}
	} else {
		for(var i=0;i<array1.length;i++){
			array[i] = array1[i] >= array2[i];
		}
	}
	return array;
} //}}}
function ArrayIsMember(array1,array2) { //{{{

	var array=NewArrayFill(array1.length,0);
	for (var i=0;i<array1.length;i++){
		for(var j=0;j<array2.length;j++){
			if (array1[i] == array2[j]){
				array[i]=1;
				break;
			}
		}
	}
	return array;
} //}}}
function NewArrayFill(size,value) { //{{{
	var array = new Array(size);
	
	for (var i = 0; i < size; i++) {
		array[i] = value;
	}
	
	return array;
} //}}}
function NewArrayFill2D(rows,cols,value) { //{{{
	var arr=new Array(rows); 

	for (var i=0;i<rows;i++) {
		arr[i] = NewArrayFill(cols,value);
	}

	return arr;
} //}}}
function NewArrayFillIncrement(start,size,increment) { //{{{
	var array=new Array(size); 

	for(var i=0;i<size;i++){
		array[i]=start+i*increment;
	}

	return array;
} //}}}
function ones(size) { //{{{
	return NewArrayFill(size,1);
} //}}}
function zeros(size) { //{{{
	return NewArrayFill(size,0);
} //}}}
function ArrayFind(array,value) { //{{{
	var count=0;
	var indices=[];

	for (var i=0;i<array.length;i++){
		if(array[i]===value){
			indices.push(count);
		}
		count++;
	}
	return indices;
} //}}}
function ArrayFind2D(array,value) { //{{{
	var count=0;
	var indices=[];

	for (var i=0;i<array.length;i++){
		for (var j=0;j<array[i].length;j++){
			if(array[i][j]===value){
				indices.push(count);
			}
			count++;
		}
	}
	return indices;
} //}}}
function ArrayFindNot(array,value) { //{{{
	var count=0;
	var indices=[];

	for (var i=0;i<array.length;i++){
		if(array[i]!==value){
			indices.push(count);
		}
		count++;
	}
	return indices;
} //}}}
function ArrayFindNot2D(array,value) { //{{{
	var count=0;
	var indices=[];

	for (var i=0;i<array.length;i++){
		for (var j=0;j<array[i].length;j++){
			if(array[i][j]!==value){
				indices.push(count);
			}
			count++;
		}
	}
	return indices;
} //}}}
function ArrayIndex(array1,array2,value) { //{{{
	//Change behavior between get (if no value is provided) to set (if value to set is provided)
	if (arguments.length == 2){
		let data = []
		if (typeof(array2[0]) == 'number') {
			for (let i=0;i<array2.length;i++){
				data.push(array1[array2[i]]);
			}
		} else {
			//2d index array
			for (let i=0;i<array2.length;i++){
				let data2 = [];
				for (let j=0;j<array2[i].length;j++){
				    data2.push(array1[array2[i][j]]);
				}
				data.push(data2);
			}
		}
		return data;
	} else {
		for (var i=0;i<array2.length;i++){
			array1[array2[i]]=value;
		}
		return array1;
	}
} //}}}
function Create2DArray(rows,cols) { //{{{
	var arr=new Array(rows); 

	for (var i=0;i<rows;i++) {
		arr[i] = new Array(cols);
	}

	return arr;
} //}}}
function MapIsEmpty(map) { //{{{
	for (var key in map){
		if(map.hasOwnProperty(key)){
			return false;
		}
	}
	return true;
} //}}}
function clone(obj) {//{{{
	
	var copy;

	// Handle the 3 simple types, and null or undefined
	if (null == obj || "object" != typeof obj) return obj;

	// Handle Date
	if (obj instanceof Date) {
		copy = new Date();
		copy.setTime(obj.getTime());
		return copy;
	}

	// Handle Array
	if (obj instanceof Array || arg instanceof Float64Array) {
		copy = [];
		for (var i = 0, len = obj.length; i < len; i++) {
			copy[i] = clone(obj[i]);
		}
		return copy;
	}

	// Handle Object
	if (obj instanceof Object) {
		copy = {};
		for (var attr in obj) {
			if (obj.hasOwnProperty(attr)) copy[attr] = clone(obj[attr]);
		}
		return copy;
	}

	throw new Error("Unable to copy obj! Its type isn't supported.");
} //}}}
function FloatFix(pointer,size) {//{{{

	var buffer=new Float64Array(size);
	for(var i=0;i<size;i++)buffer[i]=pointer[i];
	return buffer;


} //}}}
function NullFix(pointer,value) {//{{{

	if(pointer==null)return value;
	else{
		//check that the pointer values are not null: 
		if(IsArray(pointer)){
			if(IsArray(pointer[0])){
				for(var i=0;i<pointer.length;i++){
					for(var j=0;j<pointer[0].length;j++){
						if(pointer[i][j]==null)pointer[i][j]=value;
					}
				}	
			}
			else{
				for(var i=0;i<pointer.length;i++){
					if(pointer[i]==null)pointer[i]=value;
				}
			}
		}
		return pointer;
	}

} //}}}
function typedArraySliceSupport() { //{{{
	//TypedArray compatibility for Safari/IE
	if (typeof Int8Array !== 'undefined') {
		if (!Int8Array.prototype.fill) { Int8Array.prototype.fill = Array.prototype.fill; }
		if (!Int8Array.prototype.slice) { Int8Array.prototype.slice = Array.prototype.slice; }
	}
	if (typeof Uint8Array !== 'undefined') {
		if (!Uint8Array.prototype.fill) { Uint8Array.prototype.fill = Array.prototype.fill; }
		if (!Uint8Array.prototype.slice) { Uint8Array.prototype.slice = Array.prototype.slice; }
	}
	if (typeof Uint8ClampedArray !== 'undefined') {
		if (!Uint8ClampedArray.prototype.fill) { Uint8ClampedArray.prototype.fill = Array.prototype.fill; }
		if (!Uint8ClampedArray.prototype.slice) { Uint8ClampedArray.prototype.slice = Array.prototype.slice; }
	}
	if (typeof Int16Array !== 'undefined') {
		if (!Int16Array.prototype.fill) { Int16Array.prototype.fill = Array.prototype.fill; }
		if (!Int16Array.prototype.slice) { Int16Array.prototype.slice = Array.prototype.slice; }
	}
	if (typeof Uint16Array !== 'undefined') {
		if (!Uint16Array.prototype.fill) { Uint16Array.prototype.fill = Array.prototype.fill; }
		if (!Uint16Array.prototype.slice) { Uint16Array.prototype.slice = Array.prototype.slice; }
	}
	if (typeof Int32Array !== 'undefined') {
		if (!Int32Array.prototype.fill) { Int32Array.prototype.fill = Array.prototype.fill; }
		if (!Int32Array.prototype.slice) { Int32Array.prototype.slice = Array.prototype.slice; }
	}
	if (typeof Uint32Array !== 'undefined') {
		if (!Uint32Array.prototype.fill) { Uint32Array.prototype.fill = Array.prototype.fill; }
		if (!Uint32Array.prototype.slice) { Uint32Array.prototype.slice = Array.prototype.slice; }
	}
	if (typeof Float32Array !== 'undefined') {
		if (!Float32Array.prototype.fill) { Float32Array.prototype.fill = Array.prototype.fill; }
		if (!Float32Array.prototype.slice) { Float32Array.prototype.slice = Array.prototype.slice; }
	}
	if (typeof Float64Array !== 'undefined') {
		if (!Float64Array.prototype.fill) { Float64Array.prototype.fill = Array.prototype.fill; }
		if (!Float64Array.prototype.slice) { Float64Array.prototype.slice = Array.prototype.slice; }
	}
	if (typeof TypedArray !== 'undefined') {
		if (!TypedArray.prototype.fill) { TypedArray.prototype.fill = Array.prototype.fill; }
		if (!TypedArray.prototype.slice) { TypedArray.prototype.slice = Array.prototype.slice; }
	}
} //}}}
