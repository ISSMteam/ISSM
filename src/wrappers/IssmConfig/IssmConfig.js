function IssmConfig(string){
/*IssmConfig 
	   usage: var config = IssmConfig('_HAVE_PETSC_');
*/

	//output
	var pvalue= Module._malloc(4); 
	var psvalue= Module._malloc(4); 

	//Declare IssmConfig module: 
	IssmConfigModule = Module.cwrap('IssmConfigModule','number',['number','string','string']);
	
	//Call IssmConfig module: 
	IssmConfigModule(pvalue, psvalue, string);
	
	/*Dynamic copying from heap: {{{*/
	var value = Module.getValue(pvalue, 'double');
	/*}}}*/

	/*Free resources: */
	Module._free(pvalue); 

	return value;
}
