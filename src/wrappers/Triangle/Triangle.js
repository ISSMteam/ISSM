function Triangle(md,domain,rifts, area){
/*Triangle 
	   usage: var array = Triangle(domain,rifts,area);
	      where: array is made of [index,x,y,segments,segmentmarkers]
		  and index,x,y defines a triangulation, segments is an array made 
	      of exterior segments to the mesh domain outline, segmentmarkers is an array 
		  flagging each segment, domain a js array defining the domain outline  (sames for 
		  rifts) and area is the maximum area desired for any element of the resulting mesh.

		  Ok, for now, we are not dealing with rifts. Also, the domain is made of only one 
		  profile, this to avoid passing a double** pointer to js. 
*/

	//Dynamic allocations: {{{
	//Retrieve domain arrays, and allocate on Module heap: 
	
	//input
	var dx=new Float64Array(domain['x']); var nx=dx.length * dx.BYTES_PER_ELEMENT;
	var dxPtr= Module._malloc(nx); var domainxHeap = new Uint8Array(Module.HEAPU8.buffer,dxPtr,nx);
	domainxHeap.set(new Uint8Array(dx.buffer)); var domainx=domainxHeap.byteOffset;

	var dy=new Float64Array(domain['y']); var ny=dy.length * dy.BYTES_PER_ELEMENT;
	var dyPtr = Module._malloc(ny); var domainyHeap = new Uint8Array(Module.HEAPU8.buffer,dyPtr,ny);
	domainyHeap.set(new Uint8Array(dy.buffer)); var domainy=domainyHeap.byteOffset;
	
	//output
	var nel,indexlinear,index,nods,x,y;
	var pnel= Module._malloc(4); 
	var pindex= Module._malloc(4); 
	var pnods= Module._malloc(4); 
	var px= Module._malloc(4); 
	var py= Module._malloc(4); 
	var psegments= Module._malloc(4); 
	var psegmentmarkers= Module._malloc(4); 
	var pnsegs= Module._malloc(4); 
	//}}}

	//Declare Triangle module: 
	TriangleModule = Module.cwrap('TriangleModule','number',['number','number','number','number','number','number','number','number','number','number','number','number']);
	
	//Call Triangle module: 
	TriangleModule(pindex,px,py,pnel,pnods,psegments,psegmentmarkers,pnsegs, domainx,domainy,dx.length,area);
	
	/*Dynamic copying from heap: {{{*/
	//recover mesh: 
	nel = Module.getValue(pnel, 'i32');
	var indexptr = Module.getValue(pindex,'i32');
	indexlinear = Module.HEAPF64.slice(indexptr /8, indexptr/8 + nel*3);
	index = ListToMatrix(indexlinear,3);

	nods = Module.getValue(pnods, 'i32');
	var xptr = Module.getValue(px,'i32');
	var yptr = Module.getValue(py,'i32');
	x = Module.HEAPF64.slice(xptr /8, xptr/8 + nods);
	y = Module.HEAPF64.slice(yptr /8, yptr/8 + nods);
	
	nsegs = Module.getValue(pnsegs, 'i32');
	var segmentsptr = Module.getValue(psegments,'i32');
	segmentslinear = Module.HEAPF64.slice(segmentsptr /8, segmentsptr/8 + nsegs*3);
	segments = ListToMatrix(segmentslinear,3);
	
	var segmentmarkersptr = Module.getValue(psegmentmarkers,'i32');
	segmentmarkers = Module.HEAPF64.slice(segmentmarkersptr /8, segmentmarkersptr/8 + nsegs);
	/*}}}*/

	var return_array=[index,x,y,segments,segmentmarkers];

	/*Free resources: */
	Module._free(pindex); 
	Module._free(indexlinear); 
	Module._free(px); 
	Module._free(x); 
	Module._free(py); 
	Module._free(y); 
	Module._free(pnel); 
	Module._free(pnods); 
	Module._free(psegments); 
	Module._free(psegmentmarkers); 
	Module._free(pnsegs); 

	return return_array;
}
