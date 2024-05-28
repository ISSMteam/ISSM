function [variable fpos]=loadresultfromdisk(filename,step,name,varargin)
%LOADRESULTFROMDISK - load specific result of solution sequence from disk file "filename"
%
%   Usage:
%      variable=loadresultfromdisk(filename,step,name);

	%Open file
	fid=fopen(filename,'rb');
	if(fid==-1),
		error(['loadresultsfromdisk error message: could not open ',filename,' for binary reading']);
	end

	if nargin==4,
		%Put the pointer on the right position in the file:
		fpos=varargin{1};
		fseek(fid,fpos,'bof');
	end

	while 1,

		%read field
		fpos=ftell(fid);
		[length,count]=fread(fid,1,'int');

		if count==0,
			break;
			%we are done, break;
		else
			fieldname=fread(fid,length,'char');
			fieldname=fieldname(1:end-1)';
			fieldname=char(fieldname);
			rtime=fread(fid,1,'double');
			rstep=fread(fid,1,'int');

			%check on field: 
			if ((step==rstep) & (strcmpi(name,fieldname))),
				%ok, go read the result really: 
				type=fread(fid,1,'int');
				M=fread(fid,1,'int');
				if type==1,
					field=fread(fid,M,'double');
				elseif type==2,
					field=fread(fid,M,'char');
					field=char(field(1:end-1)');
				elseif type==3,
					N=fread(fid,1,'int');
					field=fread(fid,[N M],'double')';
				elseif type==4,
					N=fread(fid,1,'int');
					field=fread(fid,[N M],'int')';
				elseif type==5,
					N=fread(fid,1,'int');
					fieldr=fread(fid,[N M],'double')';
					fieldi=fread(fid,[N M],'double')';
					field=complex(fieldr,fieldi);
				else
					error(['cannot read data of type ' num2str(type) ]);
				end
				variable=field;
				break;
			else
				%just skim to next results.
				type=fread(fid,1,'int');
				M=fread(fid,1,'int');
				if type==1,
					fseek(fid,8*M,'cof');
				elseif type==2,
					fseek(fid,M,'cof');
				elseif type==3,
					N=fread(fid,1,'int');
					fseek(fid,M*N*8,'cof');
				elseif type==4,
					N=fread(fid,1,'int');
					fseek(fid,M*N*4,cof);
				else
					error(['cannot read data of type ' num2str(type) ]);
				end
			end
		end
	end
	fclose(fid);
