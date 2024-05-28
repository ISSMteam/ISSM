%Test Name: EISMINTMassConservation
%This test is a test from the EISMINT for Ice shelves Vincent Rommelaere 1996.
printingflag=false;

results={};

for stabilization=1:3;
	%The goal is to test the masstransport model
	md=bamg(model(),'domain','../Exp/SquareEISMINT.exp','hmax',3000.);
	md=setmask(md,'all','');
	md=parameterize(md,'../Par/SquareEISMINT.par');
	md.smb.mass_balance(:)=0.;
	md=setflowequation(md,'SSA','all');
	md.cluster=generic('name',oshostname(),'np',8);

	disp('      initial velocity');
	md.initialization.vx=zeros(md.mesh.numberofvertices,1);
	md.initialization.vy=-400.*ones(md.mesh.numberofvertices,1);

	%Stabilization
	if stabilization==2,
		md.masstransport.stabilization=0;
	else
		md.masstransport.stabilization=stabilization;
	end

	%spc thickness
	pos=find(md.mesh.y>199999.9);
	times=0:1:500;
	md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices+1,length(times));
	md.masstransport.spcthickness(end,:)=times;
	md.masstransport.spcthickness(pos,:)=repmat(500.+100.*sin(2.*pi*times/200.),length(pos),1);
	if stabilization==3,
		pos=find(isnan(md.masstransport.spcthickness)); md.masstransport.spcthickness(pos)=500.; %No NaN for DG
	end

	%solve
	md.transient.isstressbalance=0;
	md.settings.output_frequency=500; %keep only last step
	md.verbose=verbose();
	md=solve(md,'Transient');
	results{stabilization}=(md.results.TransientSolution(end).Thickness);
end

%plot results
[elements,x,y,z,s,h1]=SectionValues(md,results{1},'../Exp/CrossLineEISMINT.exp',100.);
[elements,x,y,z,s,h2]=SectionValues(md,results{2},'../Exp/CrossLineEISMINT.exp',100.);
[elements,x,y,z,s,h3]=SectionValues(md,results{3},'../Exp/CrossLineEISMINT.exp',100.);
[elements,x,y,z,s,hth]=SectionValues(md, 500+100*sin(2*pi/200*(500-md.mesh.y/400)),'../Exp/CrossLineEISMINT.exp',100.);
plot(s,h1,'r',s,h2,'b',s,h3,'g',s,hth,'k')
legend('Art. diff.','No Art. diff.','D.G.','Theoretical')
if printingflag,
	set(gcf,'Color','w')
	export_fig([issmdir() '/website/doc_pdf/validation/Images/EISMINT/IceShelf/eismintmasscthickness.pdf']);
end

%Fields and tolerances to track changes
field_names     ={ ...
	'ThicknessArtDiff','ThicknessNoArtDiff','ThicknessDG' ...
};
field_tolerances={...
	1e-13, 1e-13, 1e-13...
};
field_values={
	results{1}, ...
	results{2}, ...
	results{3}, ...
};
