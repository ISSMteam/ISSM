function plot_contour(md, datain, options, canvas) { //{{{
	//PLOT_MESH - Function for plotting wireframe contour.
	//
	//   Usage:
	//      plot_contour(md, options, canvas);
	//
	//   See also: PLOTMODEL, PLOT_MANAGER

	//{{{
	//Process data and model
	//let meshresults = processmesh(md, [], options);
	//let [x, y, z, index, is2d, isplanet, vertices, scale] = scaleMesh(md, meshresults, options);
	let [x, y, z, index, is2d, isplanet]=processmesh(md,[],options);
	options.removefield('log',0);;
	let [data, datatype]=processdata(md,datain,options);
	if (vesl.helpers.isEmptyOrUndefined(data)) error('data provided is empty');
	
	//check is2d
	if (!is2d && !isplanet) {
		error('plot_contour error message: contour not supported for 3d meshes, project on a layer');
	}

	//first, process data: must be on nodes
	if (datatype==1) {
		//elements -> take average
		//data=averaging(md,data,0);
		error('plot_contour error message: contour not supported for element data yet');
	} else if (datatype==2) {
		//nodes -> do nothing
	} else if (datatype==3) {
		//quiver -> take norm
		//data=sqrt(sum(datain.*datain,2)); //(original)
		//data=ArraySqrt(ArraySum(ArrayMultiply(datain, datain),2)); //js version
		error('plot_contour error message: contour not supported for quiver data yet');
	} else {
		error('datatype not supported yet');
	}

	//prepare colors
	if (options.exist('contouronly')) {
		//remove the previous plots
		//cla
		error('contouronly not supported yet');
	}
	let color=options.getfieldvalue('contourcolor','yellow');
	let linewidth=options.getfieldvalue('linewidth',1);

	//get contours levels
	let contourlevels=options.getfieldvalue('contourlevels');
	let levels;
	if (typeof(contourlevels) == 'number') {
		//levels=round_ice(linspace(max(data),min(data),contourlevels),2);
		error('numeric contourlevels not supported yet - must provide levels as an array');
	} else {
		//levels=sort(unique(levels),'descend');
		levels=ArraySort(ArrayUnique(contourlevels)).reverse();
	}
	let numlevels=levels.length;

	//initialization of some variables
	let numberofelements=index.length; //same as size(index,1)
	let elementslist=NewArrayFillIncrement(0,numberofelements,1); //1:numberofelements;
	let c=[];
	let h=[];;

	//get unique edges in mesh
	//1: list of edges
	index = ArraySubtract2D(index, 1);
	//edges=[index[:,[1,2]); index(:,[2,3]); index(:,[3,1])];
	let edges=ArrayConcat(ArrayConcat(ArrayCol(index,[0,1]), ArrayCol(index,[1,2])), ArrayCol(index,[2,0]));
	//2: find unique edges
	//[edges,I,J]=unique(sort(edges,2),'rows');
	[edges,I,J]=ArrayUnique(ArraySort(edges,2),'rows');
	//3: unique edge numbers
	let vec=J;
	//4: unique edges numbers in each triangle (2 triangles sharing the same edge will have
	//   the same edge number)
	let edges_tria=ArrayTranspose([ArrayIndex(vec,elementslist), ArrayIndex(vec,ArrayAdd(elementslist,numberofelements)), ArrayIndex(vec,ArrayAdd(elementslist,2*numberofelements))]);

	//segments [nodes1 nodes2]
	let Seg1=ArrayCol(index,[0,1]);
	let Seg2=ArrayCol(index,[1,2]);
	let Seg3=ArrayCol(index,[2,0]);
	
	//segment numbers [1;4;6;...]
	let Seg1_num=ArrayCol(edges_tria,0);
	let Seg2_num=ArrayCol(edges_tria,1);
	let Seg3_num=ArrayCol(edges_tria,2);

	//value of data on each tips of the segments
	let Data1=ArrayIndex(data,Seg1);
	let Data2=ArrayIndex(data,Seg2);
	let Data3=ArrayIndex(data,Seg3);

	//get the ranges for each segment
	let Range1=ArraySort(Data1,2);
	let Range2=ArraySort(Data2,2);
	let Range3=ArraySort(Data3,2);
	
	let hx = [];
	let hy = [];
	let hz = [];
	for (let i=0; i<numlevels; i++) {
		let level=levels[i];

		//find the segments that contain this value
		let pos1=ArrayAnd(ArrayLessThan(ArrayCol(Range1,0),level), ArrayGreaterThan(ArrayCol(Range1,1),level)); //pos1=(Range1(:,1)<level & Range1(:,2)>level);
		let pos2=ArrayAnd(ArrayLessThan(ArrayCol(Range2,0),level), ArrayGreaterThan(ArrayCol(Range2,1),level)); //pos2=(Range2(:,1)<level & Range2(:,2)>level);
		let pos3=ArrayAnd(ArrayLessThan(ArrayCol(Range3,0),level), ArrayGreaterThan(ArrayCol(Range3,1),level)); //pos3=(Range3(:,1)<level & Range3(:,2)>level);

		//get elements
		let poselem12=ArrayAnd(pos1, pos2);
		let poselem13=ArrayAnd(pos1, pos3);
		let poselem23=ArrayAnd(pos2, pos3);
		let poselem=find(ArrayOr(ArrayOr(poselem12, poselem13), poselem23));
		let numelems=length(poselem);

		//if no element has been flagged, skip to the next level
		if (numelems==0) {
			continue;
		}

		//go through the elements and build the coordinates for each segment (1 by element)
		let x1=zeros(numelems,1);
		let x2=zeros(numelems,1);
		let y1=zeros(numelems,1);
		let y2=zeros(numelems,1);
		let z1=zeros(numelems,1);
		let z2=zeros(numelems,1);

		let edge_l=zeros(numelems,2);

		for (let j=0; j < numelems; j++) {

			let weight1=(level-Data1[poselem[j]][0])/(Data1[poselem[j]][1]-Data1[poselem[j]][0]);
			let weight2=(level-Data2[poselem[j]][0])/(Data2[poselem[j]][1]-Data2[poselem[j]][0]);
			let weight3=(level-Data3[poselem[j]][0])/(Data3[poselem[j]][1]-Data3[poselem[j]][0]);

			if (poselem12[poselem[j]]) {
				x1[j]=x[Seg1[poselem[j]][0]]+weight1*(x[Seg1[poselem[j]][1]]-x[Seg1[poselem[j]][0]]);
				x2[j]=x[Seg2[poselem[j]][0]]+weight2*(x[Seg2[poselem[j]][1]]-x[Seg2[poselem[j]][0]]);
				y1[j]=y[Seg1[poselem[j]][0]]+weight1*(y[Seg1[poselem[j]][1]]-y[Seg1[poselem[j]][0]]);
				y2[j]=y[Seg2[poselem[j]][0]]+weight2*(y[Seg2[poselem[j]][1]]-y[Seg2[poselem[j]][0]]);
				z1[j]=z[Seg1[poselem[j]][0]]+weight1*(z[Seg1[poselem[j]][1]]-z[Seg1[poselem[j]][0]]);
				z2[j]=z[Seg2[poselem[j]][0]]+weight2*(z[Seg2[poselem[j]][1]]-z[Seg2[poselem[j]][0]]);
				edge_l[j][0]=Seg1_num[poselem[j]];
				edge_l[j][1]=Seg2_num[poselem[j]];

			} else if (poselem13[poselem[j]]) {
                x1[j]=x[Seg1[poselem[j]][0]]+weight1*(x[Seg1[poselem[j]][1]]-x[Seg1[poselem[j]][0]]);
				x2[j]=x[Seg3[poselem[j]][0]]+weight3*(x[Seg3[poselem[j]][1]]-x[Seg3[poselem[j]][0]]);
				y1[j]=y[Seg1[poselem[j]][0]]+weight1*(y[Seg1[poselem[j]][1]]-y[Seg1[poselem[j]][0]]);
				y2[j]=y[Seg3[poselem[j]][0]]+weight3*(y[Seg3[poselem[j]][1]]-y[Seg3[poselem[j]][0]]);
				z1[j]=z[Seg1[poselem[j]][0]]+weight1*(z[Seg1[poselem[j]][1]]-z[Seg1[poselem[j]][0]]);
				z2[j]=z[Seg3[poselem[j]][0]]+weight3*(z[Seg3[poselem[j]][1]]-z[Seg3[poselem[j]][0]]);
				edge_l[j][0]=Seg1_num[poselem[j]];
				edge_l[j][1]=Seg3_num[poselem[j]];

			} else if (poselem23[poselem[j]]) {
                x1[j]=x[Seg2[poselem[j]][0]]+weight2*(x[Seg2[poselem[j]][1]]-x[Seg2[poselem[j]][0]]);
				x2[j]=x[Seg3[poselem[j]][0]]+weight3*(x[Seg3[poselem[j]][1]]-x[Seg3[poselem[j]][0]]);
				y1[j]=y[Seg2[poselem[j]][0]]+weight2*(y[Seg2[poselem[j]][1]]-y[Seg2[poselem[j]][0]]);
				y2[j]=y[Seg3[poselem[j]][0]]+weight3*(y[Seg3[poselem[j]][1]]-y[Seg3[poselem[j]][0]]);
				z1[j]=z[Seg2[poselem[j]][0]]+weight2*(z[Seg2[poselem[j]][1]]-z[Seg2[poselem[j]][0]]);
				z2[j]=z[Seg3[poselem[j]][0]]+weight3*(z[Seg3[poselem[j]][1]]-z[Seg3[poselem[j]][0]]);
				edge_l[j][0]=Seg2_num[poselem[j]];
				edge_l[j][1]=Seg3_num[poselem[j]];
			} else {
				//it shoud not go here
			}
			let test = new THREE.Vector3(x1[j], y1[j], z1[j]);
			//console.log('helo', test.length());
		}
		//now that we have the segments, we must try to connect them...

		//loop over the subcontours
		while (!isempty(edge_l)) {
			//take the right edge of the second segment and connect it to the next segments if any
			let e1=edge_l[0][0];   let e2=edge_l[0][1];
			let xc=[x1[0],x2[0]]; let yc=[y1[0],y2[0]]; let zc=[z1[0],z2[0]]; 

			//erase the lines corresponding to this edge
			//edge_l(1,:)=[];
			//x1(1)=[]; x2(1)=[];
			//y1(1)=[]; y2(1)=[];
			//z1(1)=[]; z2(1)=[]
			edge_l.splice(0,1);
			x1.splice(0,1); x2.splice(0,1);
			y1.splice(0,1); y2.splice(0,1);
			z1.splice(0,1); z2.splice(0,1);
			let [ro1,co1]=find(ArrayEqual(edge_l,e1));

			while (!isempty(ro1)) {
				if (co1==0) {
					xc=[x2[ro1]].concat(xc); yc=[y2[ro1]].concat(yc); zc=[z2[ro1]].concat(zc);

					//next edge:
					e1=edge_l[ro1][1];

				} else {
					xc=[x1[ro1]].concat(xc); yc=[y1[ro1]].concat(yc); zc=[z1[ro1]].concat(zc);

					//next edge:
					e1=edge_l[ro1][0];
				}

				//erase the lines of this
				edge_l.splice(ro1,1);
				x1.splice(ro1,1); x2.splice(ro1,1);
				y1.splice(ro1,1); y2.splice(ro1,1);
				z1.splice(ro1,1); z2.splice(ro1,1);

				//next connection
				[ro1,co1]=find(ArrayEqual(edge_l,e1));
			}

			//same thing the other way (to the right)
			let [ro2,co2]=find(ArrayEqual(edge_l,e2));

			while (!isempty(ro2)) {

				if (co2==0) {
					xc=[x2[ro2]].concat(xc); yc=[y2[ro2]].concat(yc); zc=[z2[ro2]].concat(zc);

					//next edge:
					e2=edge_l[ro2][1];

				} else {
					xc=[x1[ro2]].concat(xc); yc=[y1[ro2]].concat(yc); zc=[z1[ro2]].concat(zc);

					//next edge:
					e2=edge_l[ro2][0];
				}

				//erase the lines of this
				edge_l.splice(ro2,1);
				x1.splice(ro2,1); x2.splice(ro2,1);
				y1.splice(ro2,1); y2.splice(ro2,1);
				z1.splice(ro2,1); z2.splice(ro2,1);

				//next connection
				[ro2,co2]=find(ArrayEqual(edge_l,e2));
			}

			//we now have one subcontour ready to be plotted
			if (options.getfieldvalue('contouronly',0)) {
				if (isplanet) {
					hx = ArrayConcat(hx, ArrayConcat(xc, [NaN]));
					hy = ArrayConcat(hy, ArrayConcat(yc, [NaN]));
					hz = ArrayConcat(hz, ArrayConcat(zc, [NaN]));
					//h=[h;patch('Xdata',[xc;NaN],'Ydata',[yc;NaN],'Zdata',[zc;NaN],'facecolor','none','linewidth',linewidth)];
				} else {
					hx = ArrayConcat(hx, ArrayConcat(xc, [NaN]));
					hy = ArrayConcat(hy, ArrayConcat(yc, [NaN]));
					hz = ArrayConcat(hz, ArrayConcat(zc, [NaN]));
					//h=[h;patch('Xdata',[xc;NaN],'Ydata',[yc;NaN],'Zdata',zc,'Cdata',zc,'facecolor','none','edgecolor','flat','linewidth',linewidth)];
				}
				//hold on      
			} else {
				//dist = 5000;
				dist = 0;
				if (isplanet) {
					if ((max(xc)-min(xc)+max(yc)-min(yc)+max(zc)-min(zc))<dist) { continue; }
					hx = ArrayConcat(hx, ArrayConcat(xc, [NaN]));
					hy = ArrayConcat(hy, ArrayConcat(yc, [NaN]));
					hz = ArrayConcat(hz, ArrayConcat(zc, [NaN]));
					//h=[h;patch('Xdata',[xc;NaN],'Ydata',[yc;NaN],'Zdata',[zc;NaN],'facecolor','none','edgecolor',color,'linewidth',linewidth)];
					//c = horzcat([level, xc'; length(xc), yc'; length(xc), zc']);
				} else {
					if ((max(xc)-min(xc)+max(yc)-min(yc))<dist) { continue; }
					hx = ArrayConcat(hx, ArrayConcat(xc, [NaN]));
					hy = ArrayConcat(hy, ArrayConcat(yc, [NaN]));
					hz = ArrayConcat(hz, ArrayConcat(zc, [NaN]));
					//h=patch('Xdata',[xc;NaN],'Ydata',[yc;NaN],'facecolor','none','edgecolor',color,'linewidth',linewidth);
					//c = horzcat([level, xc'; length(xc), yc']);
				}
				//clabel(c,h,'FontSize',10,'labelspacing',20000,'color',color);
				//hold on
			}
			//break;
		}
	}

	//Compute gl variables:
	let scale = options.getfieldvalue('contourscale', 1.0005);
	let node = new Node(
		'canvas', canvas,
		'options', options,
		'md', md,
		'name', 'contours',
		'replaceNaN', false,
		'shaderName', 'line_strip',
		'lineWidth', options.getfieldvalue('linewidth', 1),
		'scale', [scale, scale, scale]
	);

	node.patch('Vertices', [hx, hy, hz], 'FaceColor', 'none', 'EdgeColor', color);
	//node.patch('Faces', elements, 'Vertices', vertices, 'FaceColor', 'none', 'EdgeColor', edgecolor);
	//}}}
} //}}}
