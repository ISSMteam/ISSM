function plotgemb(results,fieldname,varargin);

	options=pairoptions(varargin{:});

	zerolevel= getfieldvalue(options,'zerolevel',0);
	numlevels= getfieldvalue(options,'numlevels',-1);
	element= getfieldvalue(options,'element',1);
	maxstep=getfieldvalue(options,'maxstep',length(results));

	%number of results: 
	nt=length(results); 

	for i=1:nt, 

		z0=zerolevel;

		%retrieve time, and delta around time: 
		time=results(i).time;
		if i<nt,
			deltat=results(i+1).time-time;
		else
			deltat=time-results(i-1).time;			
		end

		%figure out number of levels: 
		dz=results(i).SmbDz(element,:);
		if numlevels==-1,
			nlevels=length(dz); 
		else
			nlevels=numlevels;
		end
        dz=flipud(dz(1:nlevels)') ;

		%retrieve values: 
		field=results(i).(fieldname);
		T=flipud(field(element,1:nlevels)');

		%build vertical values: 
		nz=length(dz); 

		xA=(time-deltat/2)*ones(nz,1);
		xB=(time+deltat/2)*ones(nz,1);
		xC=(time+deltat/2)*ones(nz,1);
		xD=(time-deltat/2)*ones(nz,1);

		zA=zeros(nz,1);
		zB=zeros(nz,1);
		zC=zeros(nz,1);
		zD=zeros(nz,1);
		
        
		for j=1:nz, 
			zA(j)=z0;
			zB(j)=z0;
			zC(j)=z0+dz(j);
			zD(j)=z0+dz(j);
			z0=z0+dz(j);
        end       
    

		patch([xA,xB,xC,xD]',[zA,zB,zC,zD]',[T,T,T,T]','EdgeColor','none');

		if i>=maxstep,
			break;
		end
    end
    
