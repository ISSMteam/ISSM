%LOVENUMBERS numbers class definition
%
%   Usage:
%      lovenumbers=lovenumbers() 
%      lovenumbers=lovenumbers('maxdeg',10000,'referenceframe','CF'); 
% 
%   choose numbers of degrees required (1000 by default) and reference frame (between CF and CM; CM by default) 

classdef lovenumbers
	properties (SetAccess=public) 

		%loading love numbers:
		h             = []; %provided by PREM model
		k             = []; %idem
		l             = []; %idem
		
		%tidal love numbers for computing rotational feedback:
		th            = [];
		tk            = [];
		tl            = [];
		tk2secular    = 0; %deg 2 secular number.
		pmtf_colinear = [];
		pmtf_ortho    = [];

		%time/frequency for visco-elastic love numbers
		timefreq      = [];
		istime        = 1;

	end
	methods
		function self = lovenumbers(varargin) % {{{
			options=pairoptions(varargin{:});
			maxdeg=getfieldvalue(options,'maxdeg',1000);
			referenceframe=getfieldvalue(options,'referenceframe','CM');
			self=setdefaultparameters(self,maxdeg,referenceframe);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   lovenumbers parameters:'));

			fielddisplay(self,'h','load Love number for radial displacement');
			fielddisplay(self,'k','load Love number for gravitational potential perturbation');
			fielddisplay(self,'l','load Love number for horizontal displacements');

			fielddisplay(self,'th','tidal load Love number (deg 2)');
			fielddisplay(self,'tk','tidal load Love number (deg 2)');
			fielddisplay(self,'tl','tidal load Love number (deg 2)');
			fielddisplay(self,'tk2secular','secular fluid Love number');
			fielddisplay(self,'pmtf_colinear','Colinear component of the Polar Motion Transfer Function (e.g. x-motion due to x-component perturbation of the inertia tensor)');
			fielddisplay(self,'pmtf_ortho','Orthogonal component of the Polar Motion Transfer Function (couples x and y components, only used for Chandler Wobble)');


			fielddisplay(self,'istime','time (default: 1) or frequency love numbers (0)');
			fielddisplay(self,'timefreq','time/frequency vector (yr or 1/yr)');

		end % }}}
		function self = setdefaultparameters(self,maxdeg,referenceframe) % {{{

			%initialize love numbers:
			self.h=getlovenumbers('type','loadingverticaldisplacement','referenceframe',referenceframe,'maxdeg',maxdeg);
			self.k=getlovenumbers('type','loadinggravitationalpotential','referenceframe',referenceframe,'maxdeg',maxdeg);
			self.l=getlovenumbers('type','loadinghorizontaldisplacement','referenceframe',referenceframe,'maxdeg',maxdeg);
			self.th=getlovenumbers('type','tidalverticaldisplacement','referenceframe',referenceframe,'maxdeg',maxdeg);
			self.tk=getlovenumbers('type','tidalgravitationalpotential','referenceframe',referenceframe,'maxdeg',maxdeg);
			self.tl=getlovenumbers('type','tidalhorizontaldisplacement','referenceframe',referenceframe,'maxdeg',maxdeg);

			%secular fluid love number: 
			self.tk2secular=0.942;

			self.pmtf_colinear=0.0;
			self.pmtf_ortho=0.0;
			if maxdeg>=2
				self.pmtf_colinear=(1.0+self.k(3,:))/(1.0-self.tk(3,:)/self.tk2secular); %valid only for elastic regime, not viscous. Also neglects chandler wobble
				self.pmtf_ortho=0.0;
			end
			%time: 
			self.istime=1; %temporal love numbers by default
			self.timefreq=0; %elastic case by default.

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ismember('SealevelchangeAnalysis',analyses) | (strcmp(solution,'TransientSolution') & md.transient.isslc==0), 
				return; 
			end

			md = checkfield(md,'fieldname','solidearth.lovenumbers.h','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','solidearth.lovenumbers.k','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','solidearth.lovenumbers.l','NaN',1,'Inf',1);
			
			
			md = checkfield(md,'fieldname','solidearth.lovenumbers.th','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','solidearth.lovenumbers.tk','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','solidearth.lovenumbers.tl','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','solidearth.lovenumbers.tk2secular','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','solidearth.lovenumbers.pmtf_colinear','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','solidearth.lovenumbers.pmtf_ortho','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','solidearth.lovenumbers.timefreq','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','solidearth.lovenumbers.istime','NaN',1,'Inf',1,'values',[0 1]);

			%check that love numbers are provided at the same level of accuracy: 
			if (size(self.h,1)~=size(self.k,1) | size(self.h,1)~=size(self.l,1)),
				error('lovenumbers error message: love numbers should be provided at the same level of accuracy');
			end

			ntf=length(self.timefreq);
			if( size(self.h,2) ~= ntf | size(self.k,2) ~= ntf | size(self.l,2) ~= ntf | size(self.th,2) ~= ntf | size(self.tk,2) ~= ntf | size(self.tl,2) ~= ntf | size(self.pmtf_colinear,2) ~= ntf | size(self.pmtf_ortho,2) ~= ntf),
				error('lovenumbers error message: love numbers should have as many time/frequency steps as the time/frequency vector');
			end

			if self.istime && self.timefreq(1)~=0
				error('temporal love numbers must start with elastic response, i.e timefreq(1)=0');
			end

		end % }}}
		function list=defaultoutputs(self,md) % {{{
			list = {};
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			WriteData(fid,prefix,'object',self,'fieldname','h','name','md.solidearth.lovenumbers.h','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','k','name','md.solidearth.lovenumbers.k','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','l','name','md.solidearth.lovenumbers.l','format','DoubleMat','mattype',1);

			WriteData(fid,prefix,'object',self,'fieldname','th','name','md.solidearth.lovenumbers.th','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','tk','name','md.solidearth.lovenumbers.tk','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','tl','name','md.solidearth.lovenumbers.tl','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','pmtf_colinear','name','md.solidearth.lovenumbers.pmtf_colinear','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','pmtf_ortho','name','md.solidearth.lovenumbers.pmtf_ortho','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'data',self.tk2secular,'fieldname','lovenumbers.tk2secular','format','Double');

			if self.istime,
				scale=md.constants.yts; 
			else
				scale=1.0/md.constants.yts;
			end
			WriteData(fid,prefix,'object',self,'fieldname','istime','name','md.solidearth.lovenumbers.istime','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','timefreq','name','md.solidearth.lovenumbers.timefreq','format','DoubleMat','mattype',1,'scale',scale);

		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
			writejs1Darray(fid,[modelname '.lovenumbers.h'],self.h);
			writejs1Darray(fid,[modelname '.lovenumbers.k'],self.k);
			writejs1Darray(fid,[modelname '.lovenumbers.l'],self.l);
			writejs1Darray(fid,[modelname '.lovenumbers.istime'],self.istime);
			writejs1Darray(fid,[modelname '.lovenumbers.time'],self.time);
		end % }}}
		function self = extrude(self,md) % {{{
		end % }}}
	end
end
