function plot_vstime(md,options,nlines,ncols,i)

%plot mesh
subplot(nlines,ncols,i); 

%getting the variable structure
datastruct  = getfieldvalue(options,'Input');
structnames = strsplit_strict(datastruct,'.');

%getting the position of the to be ploted point
location=getfieldvalue(options,'position');
if ~isnumeric(location) | numel(location)~=2,
	error('location provided not supported (should be [x y])');
end
xpoint=location(1);
ypoint=location(2);

%gathering a few time related parameters
FirstTime=md.timestepping.start_time;
LastTime=md.timestepping.final_time;
Dt=md.timestepping.time_step;
OutputDt=md.settings.output_frequency*Dt;
%Constructiion of the result structure(one step above variable in
%the structure tree)
solstruct = md;
for i=2:numel(structnames)-1,
	solstruct=solstruct.(char(structnames(i)));
end

timesteps   = numel(solstruct);
%now build a table with the variable needed in the end plot
Value(1:timesteps)=NaN*ones(timesteps,1);
time(1:timesteps)=NaN*ones(timesteps,1);

for t=1:timesteps,
	
	data=solstruct(t).(char(structnames(numel(structnames))));

	[TimeData datatype]=processdata(md,data,options);
	clear data
	Value(t)=PointValues(md,TimeData,xpoint,ypoint);
	time(t)=FirstTime+(t-1)*OutputDt+Dt;
end

plot(time,Value);
clear Value time

options=addfielddefault(options,'title',md.miscellaneous.name);
options=addfielddefault(options,'ylabel',char(structnames(numel(structnames))));
options=addfielddefault(options,'xlabel','Time');
options=addfielddefault(options,'colorbar',0);
options=addfielddefault(options,'axis','auto');
options=addfielddefault(options,'xlim',[FirstTime LastTime]);
options=addfielddefault(options,'view',2);
applyoptions(md,[],options);
