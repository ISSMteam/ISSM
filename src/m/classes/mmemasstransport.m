%MMEMASSTRANSPORT class definition
%
%   Usage:
%      masstransport=mmemasstransport();
%  A little bit of info here. thickness is an ensemble set of element based time series (cell array 
%  of matrice size nel x ntimes. 
%  We supply indices into this ensemble: ex: [1 5 10] and a partition vector (values 0 to nindices-1). 
%  This means, that when we resolve this MME, we'll have: 
%  thickness = thickness{0} modified by: 
%		 - thickness{1} on elements such that their partition value is 0
%		 - thickness{5} on elements such that their partition value is 1
%		 - thickness{10} on elements such that their partition value is 2
%  The goal is to allow to index into different Mmes. Here, partition value 0 might be for GlacierMIP2 
%  elements, partition value 1 for ISMIP6 Antarctica and value 2 for ISMIP6 Greenland. 

classdef mmemasstransport 
	properties (SetAccess=public)
		thickness; %multi model ensemble time series of ice thickness evolution
		ids; %indices into the ensemble 
		partition; %partition vector 
		requested_outputs      = {};
	end
	methods
		function self = mmemasstransport(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
			self.ids=[];
			self.thickness={};
			self.partition=[];
			%default output
			self.requested_outputs={'default'};
		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'DeltaIceThickness'};
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('MmemasstransportAnalysis',analyses) |  (strcmp(solution,'TransientSolution') & md.transient.ismmemasstransport==0), return; end

			nids=length(self.ids); 
			nt=length(self.thickness);
			if size(self.ids,1) < size(self.ids,2),
				error('mmemasstransport checkconsistency error message: ids should be a column vector');
			end
			md = checkfield(md,'field',self.partition,'fieldname','partition','NaN',1,'Inf',1,'>=',-1,'<=',nids-1,'size',[md.mesh.numberofelements,1]);
			md = checkfield(md,'field',self.ids,'fieldname','ids','NaN',1,'Inf',1,'>=',1,'<=',nt);
			for i=1:nt,
				md = checkfield(md,'field',self.thickness{i},'fieldname',['self.thickness{' num2str(i) '}'],'NaN',0,'Inf',1,'timeseries',1,'timeserieslength',md.mesh.numberofelements+1); %don't check for NaN, as it's a legit value here for vertices that do not belong to the current Mme. 
			end
			md = checkfield(md,'fieldname','mmemasstransport.requested_outputs','stringrow',1);

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   mmemasstransport: '));
			fielddisplay(self,'thickness','ensemble set of thickness changes');
			fielddisplay(self,'ids','indices  into the ensemble set');
			fielddisplay(self,'partition','partition vector');
			fielddisplay(self,'requested_outputs','additional outputs requested');
		end % }}}
		function s=nummodels(self) % {{{
			s=numel(self.thickness);
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			
			WriteData(fid,prefix,'object',self,'fieldname','ids','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'fieldname','partition','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','thickness','format','MatArray','timeseries',1,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);

			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.mmemasstransport.requested_outputs','format','StringArray');

		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
			error('mmemasstransport error message: not implemented yet');
		end % }}}
		function self = extrude(self,md) % {{{
		end % }}}
	end
end
