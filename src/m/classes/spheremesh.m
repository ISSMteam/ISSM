%SPHEREMESH class definition
%
%   Usage:
%      spheremesh=spheremesh();

classdef spheremesh
	properties (SetAccess=public) 
		x                           = NaN;
		y                           = NaN;
		z                           = NaN;
		r                           = NaN;
		theta                       = NaN;
		phi                         = NaN
		elements                    = NaN
		numberoflayers              = 0;
		numberofelements            = 0;
		numberofvertices            = 0;

		vertexconnectivity          = NaN
		elementconnectivity         = NaN
		average_vertex_connectivity = 0;
	end
	methods
		function self = spheremesh(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%the connectivity is the avergaded number of nodes linked to a
			%given node through an edge. This connectivity is used to initially
			%allocate memory to the stiffness matrix. A value of 16 seems to
			%give a good memory/time ration. This value can be checked in
			%trunk/test/Miscellaneous/runme.m
			self.average_vertex_connectivity=25;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','spheremesh.x','NaN',1,'Inf',1,'size',[md.spheremesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','spheremesh.y','NaN',1,'Inf',1,'size',[md.spheremesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','spheremesh.z','NaN',1,'Inf',1,'size',[md.spheremesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','spheremesh.r','NaN',1,'Inf',1,'size',[md.spheremesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','spheremesh.theta','NaN',1,'Inf',1,'size',[md.spheremesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','spheremesh.phi','NaN',1,'Inf',1,'size',[md.spheremesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','spheremesh.elements','NaN',1,'Inf',1,'>',0,'values',1:md.spheremesh.numberofvertices);
			md = checkfield(md,'fieldname','spheremesh.elements','size',[md.spheremesh.numberofelements 3]);
			if any(~ismember(1:md.spheremesh.numberofvertices,sort(unique(md.spheremesh.elements(:)))));
				md = checkmessage(md,'orphan nodes have been found. Check the spheremesh outline');
			end
			md = checkfield(md,'fieldname','spheremesh.numberoflayers','>=',0);
			md = checkfield(md,'fieldname','spheremesh.numberofelements','>',0);
			md = checkfield(md,'fieldname','spheremesh.numberofvertices','>',0);
			md = checkfield(md,'fieldname','spheremesh.elementconnectivity','size',[md.spheremesh.numberofelements 3],'NaN',1,'Inf',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Mesh:')); 

			disp(sprintf('\n      Elements and vertices:'));
			fielddisplay(self,'numberofelements','number of elements');
			fielddisplay(self,'numberofvertices','number of vertices');
			fielddisplay(self,'elements','vertex indices of the mesh elements');
			fielddisplay(self,'x','vertices x coordinate [m]');
			fielddisplay(self,'y','vertices y coordinate [m]');
			fielddisplay(self,'z','vertices z coordinate [m]');
			fielddisplay(self,'r','vertices r coordinate [m]');
			fielddisplay(self,'theta','vertices theta coordinate [degrees]');
			fielddisplay(self,'phi','vertices phi coordinate [degrees]');

			disp(sprintf('\n      Properties:'));
			fielddisplay(self,'numberoflayers','number of extrusion layers');

			fielddisplay(self,'vertexconnectivity','list of elements connected to vertex_i');
			fielddisplay(self,'elementconnectivity','list of elements adjacent to element_i');
			fielddisplay(self,'average_vertex_connectivity','average number of vertices connected to one vertex');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'fieldname','x','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','y','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','z','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','r','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','theta','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','phi','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','elements','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'fieldname','numberoflayers','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','numberofelements','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','numberofvertices','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','elementconnectivity','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'fieldname','average_vertex_connectivity','format','Integer');
		end % }}}
	end
end
