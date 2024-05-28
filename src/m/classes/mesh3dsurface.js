//MESH3DSURFACE class definition
//
//   Usage:
//      mesh3dsurface= new mesh3dsurface();

function mesh3dsurface () {
	//methods 
		this.setdefaultparameters = function (){ //{{{

			//the connectivity is the averaged number of nodes linked to a
			//given node through an edge. This connectivity is used to initially
			//allocate memory to the stiffness matrix. A value of 16 seems to
			//give a good memory/time ration. This value can be checked in
			//trunk/test/Miscellaneous/runme.m
			this.average_vertex_connectivity=25;
		}
		// }}}
		this.disp = function () { //{{{
			console.log(sprintf("   2D tria Mesh (3D surface):")); 

			console.log(sprintf("\n      Elements and vertices:"));
			fielddisplay(this,"numberofelements","number of elements");
			fielddisplay(this,"numberofvertices","number of vertices");
			fielddisplay(this,"elements","vertex indices of the mesh elements");
			fielddisplay(this,"x","vertices x coordinate [m]");
			fielddisplay(this,"y","vertices y coordinate [m]");
			fielddisplay(this,"z","vertices z coordinate [m]");
			fielddisplay(this,"edges","edges of the 2d mesh (vertex1 vertex2 element1 element2)");
			fielddisplay(this,"numberofedges","number of edges of the 2d mesh");

			console.log(sprintf("\n      Properties:"));
			fielddisplay(this,"vertexonboundary","vertices on the boundary of the domain flag list");
			fielddisplay(this,"segments","edges on domain boundary (vertex1 vertex2 element)");
			fielddisplay(this,"segmentmarkers","number associated to each segment");
			fielddisplay(this,"vertexconnectivity","list of elements connected to vertex_i");
			fielddisplay(this,"elementconnectivity","list of elements adjacent to element_i");
			fielddisplay(this,"average_vertex_connectivity","average number of vertices connected to one vertex");

			console.log(sprintf("\n      Extracted model:"));
			fielddisplay(this,"extractedvertices","vertices extracted from the model");
			fielddisplay(this,"extractedelements","elements extracted from the model");

			console.log(sprintf("\n      Projection:"));
			fielddisplay(this,"lat","vertices latitude [degrees]");
			fielddisplay(this,"long","vertices longitude [degrees]");
			fielddisplay(this,"r","vertices radius [m]");
		} //}}}
		this.classname = function () { //{{{
			return "mesh3dsurface";
		} //}}}
		this.domaintype=function (){ // {{{
			return '3Dsurface';
		} // }}}
		this.dimension = function () { //{{{
			return 2;
		} //}}}
		this.elementtype = function() {//{{{
			return 'Tria';
		} // }}}
		this.checkconsistency = function(md,solution,analyses){ //{{{

			checkfield(md,'fieldname','mesh.x','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
			checkfield(md,'fieldname','mesh.y','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
			checkfield(md,'fieldname','mesh.z','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
			checkfield(md,'fieldname','mesh.lat','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
			checkfield(md,'fieldname','mesh.long','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
			checkfield(md,'fieldname','mesh.r','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
			checkfield(md,'fieldname','mesh.elements','NaN',1,'Inf',1,'>',0,'values',NewArrayFillIncrement(md.mesh.numberofvertices,1,1));
			checkfield(md,'fieldname','mesh.elements','size',[md.mesh.numberofelements, 3]);
			if(ArrayAnyEqual(ArrayIsMember(NewArrayFillIncrement(md.mesh.numberofvertices,1,1),ArraySort(ArrayUnique(MatrixToList(md.mesh.elements)))),0)){
				md.checkmessage('orphan nodes have been found. Check the mesh outline');
			}
			checkfield(md,'fieldname','mesh.numberofelements','>',0);
			checkfield(md,'fieldname','mesh.numberofvertices','>',0);
			checkfield(md,'fieldname','mesh.average_vertex_connectivity','>=',9,'message',"'mesh.average_vertex_connectivity' should be at least 9 in 2d");
			checkfield(md,'fieldname','mesh.segments','NaN',1,'Inf',1,'>',0,'size',[NaN, 3]);

		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{
			WriteData(fid,prefix,'name','md.mesh.domain_type','data','Domain' + this.domaintype(),'format','String');
			WriteData(fid,prefix,'name','md.mesh.domain_dimension','data',this.dimension(),'format','Integer');
			WriteData(fid,prefix,'name','md.mesh.elementtype','data',this.elementtype(),'format','String');
			WriteData(fid,prefix,'object',this,'fieldname','x','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','y','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','z','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','lat','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','long','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','r','format','DoubleMat','mattype',1);
			//WriteData(fid,prefix,'object',this,'fieldname','z','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','elements','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',this,'fieldname','numberofelements','format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','numberofvertices','format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','average_vertex_connectivity','format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','vertexonboundary','format','DoubleMat','mattype',1);
			//WriteData(fid,prefix,'object',this,'fieldname','segments','format','DoubleMat','mattype',3);
		}//}}}
		this.fix=function() { //{{{
			//Transform objects into Float64Arrays:
			this.x=FloatFix(this.x,this.numberofvertices); 
			this.y=FloatFix(this.y,this.numberofvertices); 
			this.z=FloatFix(this.y,this.numberofvertices); 
			this.r=FloatFix(this.y,this.numberofvertices); 
			this.edges=NullFix(this.edges,NaN);
			this.vertexonboundary=FloatFix(this.vertexonboundary,this.numberofvertices); 
			this.segmentmarkers=FloatFix(this.segmentmarkers,this.segments.length);
			this.extractedvertices=NullFix(this.extractedvertices,NaN);
			this.extractedelements=NullFix(this.extractedelements,NaN);
			this.lat=NullFix(this.lat,NaN);
			this.long=NullFix(this.long,NaN);
		}//}}}

	//properties 
	// {{{
		this.x                           = NaN;
		this.y                           = NaN;
		this.z                           = NaN;
		this.elements                    = NaN;
		this.numberofelements            = 0;
		this.numberofvertices            = 0;
		this.numberofedges               = 0;

		this.lat                         = NaN;
		this.long                        = NaN;
		this.r                           = NaN;

		this.vertexonboundary            = NaN;
		this.edges                       = NaN;
		this.segments                    = NaN;
		this.segmentmarkers              = NaN;
		this.vertexconnectivity          = NaN;
		this.elementconnectivity         = NaN;
		this.average_vertex_connectivity = 0;

		this.extractedvertices           = NaN;
		this.extractedelements           = NaN;

		this.setdefaultparameters();
		//}}}
}
