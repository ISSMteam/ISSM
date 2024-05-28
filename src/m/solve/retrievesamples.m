function varargout=retrievesamples(varargin)

	options=pairoptions(varargin{:});

	directory=getfieldvalue(options,'directory');
	name=getfieldvalue(options,'name');
	nsamples=getfieldvalue(options,'nsamples');
	step=getfieldvalue(options,'step');
	fields=getfieldvalue(options,'fields');

	if ~isa(fields,'cell'),
		error('retrievesamples error message: ''fields'' should be a cell array of field names');
	end
	
	nout=length(fields);
	for n=1:nout,
		field=fields{n};
							
		[sample,fpos]=loadresultfromdisk(sprintf('%s/%s.outbin.%i',directory,name,1),step,field);
		nv=length(sample);
		
		%initialize: 
		samples=zeros(nv,nsamples);
		samples(:,1)=sample;
	
		for i=2:nsamples,
			if mod(i,10)==0, disp(i/nsamples*100); end
			samples(:,i)=loadresultfromdisk(sprintf('%s/%s.outbin.%i',directory,name,i),step,field,fpos);
		end
		varargout{n}=samples;
	end
