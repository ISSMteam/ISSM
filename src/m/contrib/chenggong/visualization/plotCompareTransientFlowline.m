function plotCompareTransientFlowline(velList, titleList, md, flowline, TStart, TEnd, xl, yl, ca)
%plotCompareTransientFlowline - to compare the solutions along a given
%  flowine for the given time periods
%
%   velList:	a cell of transient velocities on the model mesh
%   titleList: names of the velocity data
%   md:			ISSM model
%   flowline:	a flowline data structure contains flowline.x and flowline.y. 
%					If possible, it should also contains flowline.Xmain, which is 
%					the distance along the flowline from upstream to the terminus position.
%   TStart:		a list of start time points
%   TEnd:		a list of end time points
%   ca:			caixs value in plot
%   xl:			xlim value
%   yl:			ylim value
%

N = length(velList);
vel_flowline = cell(N,1);
if ~isfield(flowline, 'Xmain')
    flowline.Xmain = cumsum([0;sqrt((flowline.x(2:end)- flowline.x(1:end-1)).^2 + (flowline.y(2:end)- flowline.y(1:end-1)).^2)]')/1000;
end

for i = 1:N
    vel_flowline{i} = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,velList{i},flowline.x,flowline.y);
end

for p = 1:N
    subplot(N,1,p)
    hold on
    for i = 1:length(TStart)
        overlay = [vel_flowline{p}(:,i), vel_flowline{p}(:,i)]';
        imAlpha = ones(size(overlay));
        imAlpha(isnan(overlay)) = 0;
        imageNonUni(flowline.Xmain, [TStart(i), TEnd(i)], overlay);
    end
    title(titleList{p}, 'Interpreter', 'latex');
    xlim(xl)
    xlabel('x(km)')
    ylabel('year')
    ylim(yl);
    caxis(ca)
    colormap('jet')
    h=colorbar;
    title(h,'m/a')
end
