function self=read_netCDF(filename)

	% Different types in the netcdf standard are:
	%   2 for char
	%   4 for integer
	%   6 for doubles	

	ncid=netcdf.open(filename,'NC_NOWRITE');
	groupIDs=netcdf.inqGrps(ncid);%retrieve group IDs
	self=model;
	%loop on groups
	for i=1:length(groupIDs)
		whichclass = netcdf.getAtt(groupIDs(i),netcdf.getConstant('NC_GLOBAL'),'classtype');
		groupName = netcdf.inqGrpName(groupIDs(i));		
		%results needs a special treatment as it is a structure
		if strcmp(whichclass,'results'),
			subgroupIDs=netcdf.inqGrps(groupIDs(i));%retrieve group IDs
			%define the model structure
			self=setfield(self,groupName,struct);
			for j=1:length(subgroupIDs)
				subclass = netcdf.getAtt(subgroupIDs(j),netcdf.getConstant('NC_GLOBAL'),'classtype');
				self.results=setfield(self.results,subclass,struct);
				[ndims nvar natts]=netcdf.inq(subgroupIDs(j));
				varIDs=netcdf.inqVarIDs(subgroupIDs(j));
				%first loop on group atributes
				for k=1:natts,
					attname = netcdf.inqAttName(subgroupIDs(j),netcdf.getConstant('NC_GLOBAL'),k-1);
					[xtype,attlen] = netcdf.inqAtt(subgroupIDs(j),netcdf.getConstant('NC_GLOBAL'),attname);
					disp(sprintf('In %s, Treating attribute %s of type %i',subclass,attname,xtype));
					%classtype have done is job, no need to keep it any more
					if ~strcmp(attname,'classtype'),
						attval=netcdf.getAtt(subgroupIDs(j),netcdf.getConstant('NC_GLOBAL'),attname);
						if strcmp(attval,'False'),
							self.(groupName).(subclass).(attname)=false;
						elseif strcmp(attval,'True')
							self.(groupName).(subclass).(attname)=true;
						else
							self.(groupName).(subclass).(attname)=attval;
						end
					end
				end
				%now loop on variable in group
				count=0;
				for k=1:length(varIDs),
					[varname, xtype, varDimIDs, varAtts] =netcdf.inqVar(subgroupIDs(j),varIDs(k));
					disp(sprintf('In %s, Treating variable %s of type %i',whichclass,varname,xtype));
					%time dimension seems to be last in our construction
					for l=1:length(varDimIDs),
						[dimname, dimlen] = netcdf.inqDim(ncid,varDimIDs(l));
						count(l)=[dimlen];
					end
					startpoint=zeros(size(varDimIDs));
					timestep=count(end);
					count(end)=1;
					for l=1:timestep,
						data=netcdf.getVar(subgroupIDs(j),varIDs(k),startpoint,count);
						self.(groupName).(subclass)(l).(varname)=data;
						startpoint(end)=startpoint(end)+1;
						self.(groupName).(subclass)(l).('errlog')='';
						self.(groupName).(subclass)(l).('outlog')='';
						self.(groupName).(subclass)(l).('SolutionType')=subclass;
					end
					count=0;
				end
			end
		else,
			%define the model structure
			self.(groupName)=eval(whichclass);
			varIDs=netcdf.inqVarIDs(groupIDs(i));
			[ndims nvar natts]=netcdf.inq(groupIDs(i));
			%first loop on group atributes
			for j=1:natts,
				attname = netcdf.inqAttName(groupIDs(i),netcdf.getConstant('NC_GLOBAL'),j-1);
				[xtype,attlen] = netcdf.inqAtt(groupIDs(i),netcdf.getConstant('NC_GLOBAL'),attname);
				disp(sprintf('In %s, Treating attribute %s of type %i',whichclass,attname,xtype));
				%classtype have done is job, no need to keep it any more
				if ~strcmp(attname,'classtype'),
					attval=netcdf.getAtt(groupIDs(i),netcdf.getConstant('NC_GLOBAL'),attname);
					if strcmp(attval,'False'),
						self.(groupName).(attname)=false;
					elseif strcmp(attval,'True')
						self.(groupName).(attname)=true;
					else
						self.(groupName).(attname)=attval;
					end
				end
			end
			%now loop on variable in group
			for j=1:length(varIDs),
				[varname, xtype, varDimIDs, varAtts] =netcdf.inqVar(groupIDs(i),varIDs(j));
				disp(sprintf('In %s, Treating variable %s of type %i',whichclass,varname,xtype));			
				%if the value is a single string, we need to transpose it (cross check with python file is necessary)
				if xtype==2
					varval=netcdf.getVar(groupIDs(i),varIDs(j))';
					varval=cellstr(varval)';
					if strcmp(varval{1},'emptystruct'),
						self.(groupName).(varname)=struct;
					elseif strcmp(varval{1},'emptycell'),
						self.(groupName).(varname)=cell(0,0);
					else
						self.(groupName).(varname)=varval;
					end
				else
					self.(groupName).(varname)=netcdf.getVar(groupIDs(i),varIDs(j));
				end
			end
		end
	end
	netcdf.close(ncid)
end
