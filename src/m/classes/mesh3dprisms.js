//MESH3DPRISMS class definition
//
//   Usage:
//      mesh=mesh3dprisms();

function mesh3dprisms() {
		this.setdefaultparameters = function() { // {{{)

			//the connectivity is the averaged number of nodes linked to a
			//given node through an edge. This connectivity is used to initially
			//allocate memory to the stiffness matrix. A value of 16 seems to
			//give a good memory/time ration. This value can be checked in
			//trunk/test/Miscellaneous/runme.m
			this.average_vertex_connectivity=25;
		} // }}}
		this.checkconsistency = function(md,solution,analyses) { //{{{[)

			checkfield(md,'fieldname','mesh.x','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
			checkfield(md,'fieldname','mesh.y','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
			checkfield(md,'fieldname','mesh.z','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
			checkfield(md,'fieldname','mesh.elements','NaN',1,'Inf',1,'>',0,'values',NewArrayFillIncrement(md.mesh.numberofvertices,1,1));
			checkfield(md,'fieldname','mesh.elements','size',[md.mesh.numberofelements, 6]);

            if(ArrayAnyEqual(ArrayIsMember(NewArrayFillIncrement(md.mesh.numberofvertices,1,1),ArraySort(ArrayUnique(MatrixToList(md.mesh.elements)))),0)){
				//checkmessage(md,'orphan nodes have been found. Check the mesh outline'); @TODO
                md.checkmessage('orphan nodes have been found. Check the mesh outline');
			}
			checkfield(md,'fieldname','mesh.numberoflayers','>=',0);
			checkfield(md,'fieldname','mesh.numberofelements','>',0);
			checkfield(md,'fieldname','mesh.numberofvertices','>',0);
			checkfield(md,'fieldname','mesh.vertexonbase','size',[md.mesh.numberofvertices, 1],'values',[0, 1]);
			checkfield(md,'fieldname','mesh.vertexonsurface','size',[md.mesh.numberofvertices, 1],'values',[0, 1]);
			checkfield(md,'fieldname','mesh.z','>=',md.geometry.base-Math.pow(10, -10),'message','\'mesh.z\' lower than bedrock');
			checkfield(md,'fieldname','mesh.z','<=',md.geometry.surface+Math.pow(10, -10),'message','\'mesh.z\' higher than surface elevation');
			checkfield(md,'fieldname','mesh.average_vertex_connectivity','>=',24,'message','\'mesh.average_vertex_connectivity\' should be at least 24 in 3d');
			if(this.scale_factor.length>1) checkfield(md,'fieldname','mesh.scale_factor','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
		} // }}}
		this.disp = function()  { // {{{
			console.log(sprintf('   3D prism Mesh:')); 

			console.log(sprintf('\n      Elements and vertices of the original 2d mesh:'));
			fielddisplay(this,'numberofelements2d','number of elements');
			fielddisplay(this,'numberofvertices2d','number of vertices');
			fielddisplay(this,'elements2d','vertex indices of the mesh elements');
			fielddisplay(this,'x2d','vertices x coordinate [m]');
			fielddisplay(this,'y2d','vertices y coordinate [m]');

            console.log(sprintf('\n      Elements and vertices of the extruded 3d mesh:'));
			fielddisplay(this,'numberofelements','number of elements');
			fielddisplay(this,'numberofvertices','number of vertices');
			fielddisplay(this,'elements','vertex indices of the mesh elements');
			fielddisplay(this,'x','vertices x coordinate [m]');
			fielddisplay(this,'y','vertices y coordinate [m]');
			fielddisplay(this,'z','vertices z coordinate [m]');

			console.log(sprintf('\n      Properties:'));
			fielddisplay(this,'numberoflayers','number of extrusion layers');
			fielddisplay(this,'vertexonbase','lower vertices flags list');
			fielddisplay(this,'vertexonsurface','upper vertices flags list');
			fielddisplay(this,'uppervertex','upper vertex list (NaN for vertex on the upper surface)');
			fielddisplay(this,'upperelements','upper element list (NaN for element on the upper layer)');
			fielddisplay(this,'lowervertex','lower vertex list (NaN for vertex on the lower surface)');
			fielddisplay(this,'lowerelements','lower element list (NaN for element on the lower layer');
			fielddisplay(this,'vertexonboundary','vertices on the boundary of the domain flag list');

			fielddisplay(this,'vertexconnectivity','list of elements connected to vertex_i');
			fielddisplay(this,'elementconnectivity','list of elements adjacent to element_i');
			fielddisplay(this,'average_vertex_connectivity','average number of vertices connected to one vertex');

			console.log(sprintf('\n      Extracted model:'));
			fielddisplay(this,'extractedvertices','vertices extracted from the model');
			fielddisplay(this,'extractedelements','elements extracted from the model');

			console.log(sprintf('\n      Projection:'));
			fielddisplay(this,'lat','vertices latitude [degrees]');
			fielddisplay(this,'long','vertices longitude [degrees]');
			fielddisplay(this,'epsg','EPSG code (ex: 3413 for UPS Greenland, 3031 for UPS Antarctica)');
			fielddisplay(self,"scale_factor","Projection correction for volume, area, etc. computation)");
		} // }}}
		this.marshall = function(md,prefix,fid) { // {{{
			WriteData(fid,prefix,'name','md.mesh.domain_type','data','Domain' + this.domaintype(),'format','String');
			WriteData(fid,prefix,'name','md.mesh.domain_dimension','data',this.dimension(),'format','Integer');
			WriteData(fid,prefix,'name','md.mesh.elementtype','data',this.elementtype(),'format','String');
			WriteData(fid,prefix,'object',this,'class','mesh','fieldname','x','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'class','mesh','fieldname','y','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'class','mesh','fieldname','z','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'class','mesh','fieldname','elements','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',this,'class','mesh','fieldname','numberoflayers','format','Integer');
			WriteData(fid,prefix,'object',this,'class','mesh','fieldname','numberofelements','format','Integer');
			WriteData(fid,prefix,'object',this,'class','mesh','fieldname','numberofvertices','format','Integer');
			WriteData(fid,prefix,'object',this,'class','mesh','fieldname','vertexonbase','format','BooleanMat','mattype',1);
			WriteData(fid,prefix,'object',this,'class','mesh','fieldname','vertexonsurface','format','BooleanMat','mattype',1);
			WriteData(fid,prefix,'object',this,'class','mesh','fieldname','lowerelements','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',this,'class','mesh','fieldname','upperelements','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',this,'class','mesh','fieldname','average_vertex_connectivity','format','Integer');
			WriteData(fid,prefix,'object',this,'class','mesh','fieldname','elements2d','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',this,'class','mesh','fieldname','numberofvertices2d','format','Integer');
			WriteData(fid,prefix,'object',this,'class','mesh','fieldname','numberofelements2d','format','Integer');
			WriteData(fid,prefix,'object',this,'class','mesh','fieldname','scale_factor','format','DoubleMat','mattype',1);
		} // }}}
        this.fix=function() { //{{{
            //Transform objects into Float64Arrays:
            this.x=FloatFix(this.x,this.numberofvertices); 
            this.y=FloatFix(this.y,this.numberofvertices); 
            this.z=FloatFix(this.y,this.numberofvertices); 
            this.r=FloatFix(this.y,this.numberofvertices); 
            this.edges=NaNFix(this.edges,NaN);
            this.vertexonboundary=FloatFix(this.vertexonboundary,this.numberofvertices); 
            this.segmentmarkers=FloatFix(this.segmentmarkers,this.segments.length);
            this.extractedvertices=NaNFix(this.extractedvertices,NaN);
            this.extractedelements=NaNFix(this.extractedelements,NaN);
            this.lat=NaNFix(this.lat,NaN);
            this.long=NaNFix(this.long,NaN);
				this.scale_factor=NullFix(this.scale_factor,NaN);
        }//}}}
		this.domaintype = function() { // {{{)
			return '3D';
		} // }}}
		this.dimension = function() { // {{{)
			return 3;
		} // }}}
		this.elementtype = function() { // {{{)
			return 'Penta';
		} // }}}
        this.classname = function () { //{{{
            return "mesh3dprisms";
        } //}}}

	//properties 
	// {{{
        this.x                           = NaN;
        this.y                           = NaN;
        this.z                           = NaN;
        this.elements                    = NaN;
        this.numberoflayers              = 0;
        this.numberofelements            = 0;
        this.numberofvertices            = 0;

        this.lat                         = NaN;
        this.long                        = NaN;
        this.epsg                        = 0;
		this.scale_factor                = NaN;

        this.vertexonbase                = NaN;
        this.vertexonsurface             = NaN;
        this.lowerelements               = NaN;
        this.lowervertex                 = NaN;
        this.upperelements               = NaN;
        this.uppervertex                 = NaN;
        this.vertexonboundary            = NaN;

        this.vertexconnectivity          = NaN;
        this.elementconnectivity         = NaN;
        this.average_vertex_connectivity = 0;

        this.x2d                         = NaN;
        this.y2d                         = NaN;
        this.elements2d                  = NaN;
        this.numberofvertices2d          = 0;
        this.numberofelements2d          = 0;

        this.extractedvertices           = NaN;
        this.extractedelements           = NaN;
	//}}}
}
