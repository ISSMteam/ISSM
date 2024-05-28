function showregion(md,insetpos)
%SHOWREGION - show region on plot
%
%   Usage:
%      showregion(md,insetpos);

%get inset relative position (x,y,width,height)
%insetpos=getfieldvalue(options,'insetpos',[0.02 0.70 0.18 0.18]);

%get current plos position
cplotpos=get(gca,'pos');
%compute inset position
PosInset=[cplotpos(1)+insetpos(1)*cplotpos(3),cplotpos(2)+insetpos(2)*cplotpos(4), insetpos(3)*cplotpos(3), insetpos(4)*cplotpos(4)];
axes('pos',PosInset);
axis equal off
%box off
if md.mesh.epsg==3413,
	A=expread(['/u/astrid-r1b/ModelData/Exp/Greenland.exp']);
elseif md.mesh.epsg==3031,
	A=expread(['/u/astrid-r1b/ModelData/Exp/Antarctica.exp']);
else
	error('md.mesh.epsg not supported yet');
end

Ax=[min(A.x) max(A.x)];
Ay=[min(A.y) max(A.y)];

mdx=[min(md.mesh.x) max(md.mesh.x)];
mdy=[min(md.mesh.y) max(md.mesh.y)];

line(A.x,A.y,'color','b');
patch([Ax(1)  Ax(2)  Ax(2)  Ax(1) Ax(1)],[Ay(1)  Ay(1)  Ay(2)  Ay(2) Ay(1)],[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1,'FaceLighting','none')
patch( [mdx(1) mdx(2) mdx(2) mdx(1)],[mdy(1) mdy(1) mdy(2) mdy(2)],ones(4,1),'EdgeColor',[0 0 0],'FaceColor','r','FaceAlpha',0.5)
colorbar('off');
