%%
%  create kml linestrings for the partition edges.
%
%  [kfold]=kml_part_edges(md,params)
%
%  where the required input is:
%    md            (model, model class object)
%
%  the optional input is:
%    params        (string/numeric, parameter names and values)
%
%  and the optional input is:
%    latsgn        (numeric, +1/-1 for north/south latitude)
%    data          (numeric, element or nodal results data)
%    alt           (numeric, altitude for polygons, default 10000)
%    cmin          (numeric, minimum of color map)
%    cmax          (numeric, maximum of color map)
%    cmap          (char or numeric, colormap definition)
%    prtplt        (char, 'off'/'no' for partition segment plot)
%
%  and the required output is:
%    kfold         (kml_folder, folder of linestring placemarks)
%
function [kfold]=kml_part_edges(varargin)

if ~nargin
    help kml_part_edges
    return
end

%%  process input data

iarg=1;
if (nargin >= 1)
    md=varargin{1};
end
if ~exist('md','var') || isempty(md) || ~isa(md,'model')
    error(['Model ''' inputname(iarg) ''' is unrecognized class ''' class(md) '''.']);
end

%  parameters

iarg=iarg+1;
while (iarg <= nargin-1)
    if ischar(varargin{iarg})
        if ~isempty(strmatch(varargin{iarg},...
                {'latsgn','data','alt',...
                 'cmin','cmax','cmap','prtplt'},...
                'exact'))
            eval([varargin{iarg} '=varargin{iarg+1};']);
            disp([varargin{iarg} '=' any2str(varargin{iarg+1},20) ';']);
        else
            warning([varargin{iarg} '=' any2str(varargin{iarg+1},20) ' is not recognized.']);
        end
    else
        error(['''' any2str(varargin{iarg}) ''' is not a parameter name.']);
    end
    if strcmpi(varargin{iarg},'data')
        cdata=inputname(iarg+1);
    end
    iarg=iarg+2;
end

if isempty(md.mesh.lat)  || ((numel(md.mesh.lat) == 1)  && isnan(md.mesh.lat)) || ...
   isempty(md.mesh.long) || ((numel(md.mesh.long) == 1) && isnan(md.mesh.long))
    if     ~exist('latsgn','var')
        error(['Missing ''latsgn'' parameter to calculate missing lat/long data.']);
    elseif (abs(latsgn) ~= 1)
        error(['Incorrect latsgn=' num2str(latsgn) ' parameter to calculate missing lat/long data.']);
    else
        display('Converting x/y data to lat/long data.');
        [md.mesh.lat,md.mesh.long]=xy2ll(md.x,md.y,latsgn);
    end
end

if exist('data','var') && ~isempty(data)
    if     (numel(data)==md.mesh.numberofelements)
        edata=data;
    elseif (numel(data)==md.mesh.numberofvertices)
        ndata=data;
        display('Averaging nodal data to element data.');
        edata=zeros(1,md.mesh.numberofelements);
        for i=1:size(md.mesh.elements,1)
            for j=1:size(md.mesh.elements,2)
                edata(i)=edata(i)+ndata(md.mesh.elements(i,j));
            end
            edata(i)=edata(i)/size(md.mesh.elements,2);
        end
    else
        error(['Data has incorrect number of ' num2str(numel(data)) ' values.']);
    end
end

%  colormap command operates on a figure, so create an invisible one
%  (could also directly call colormaps, e.g. jet(64), but risky)

hfig=figure('Visible','off');
if exist('cmap','var')
    colormap(cmap)
end
cmap=colormap;
close(hfig)

if exist('edata','var')
    if ~exist('cmin','var')
        cmin=min(min(edata));
    end
    if ~exist('cmax','var')
        cmax=max(max(edata));
    end
end

if ~exist('alt','var')
    alt=10000;
end

%%  write folder for partition edges

if (~exist('prtplt','var') || strncmpi(prtplt,'on' ,2) || strncmpi(prtplt,'y',1)) && ...
    md.qmu.numberofpartitions
    kfold=kml_folder();
    kfold.name      ='Partition Edges';
    kfold.visibility=1;
    kfold.descript  =sprintf('Partitions=%d, Nodes=%d',...
        md.qmu.numberofpartitions,md.mesh.numberofvertices);
    kfold.feature   ={repmat(kml_placemark(),1,md.qmu.numberofpartitions)};

%  write each partition as a linestring multigeometry placemark

    disp(['Writing ' num2str(md.qmu.numberofpartitions) ' partitions as KML linestrings.']);
    epart=md.qmu.partition(md.mesh.elements)+1;
    if exist('ndata','var') || exist('edata','var')
        pdata=zeros(1,md.qmu.numberofpartitions);
        pdata(:)=NaN;
    end

%  loop over each partition

    for k=1:md.qmu.numberofpartitions
%        disp(['partition k=' int2str(k)])

%  for each partition, find all the included elements and determine the
%  perimeter (including those shared by another partition)

        [icol,irow]=find(epart'==k);
        if isempty(irow)
            continue;
        end
        irow=unique(irow);
        elemp=md.mesh.elements(irow,:);
        epartp=epart(irow,:);
        nodeconp=kmlnodeconnectivity(elemp,md.mesh.numberofvertices);
        [edgeadjp]=edgeadjacency(elemp,nodeconp);
        [edgeper,elemper,iloop]=edgeperimeter(elemp,nodeconp,edgeadjp);
        iloop(end+1)=size(edgeper,1)+1;

%  determine the data to be used for the colors (if any)

        if exist('ndata','var')
            pdata(k)=ndata(find(md.qmu.partition+1==k,1));
        elseif exist('edata','var')
            for i=1:size(epartp,1)
                if isempty(find(epart(i,:)~=k,1))
                    pdata(k)=edata(irow(i));
                    break
                end
            end
            if isnan(pdata(k))
                warning('Data for Partition %d is not defined.\n',k)
            end
        end

%  set up the placemark with multigeometry

        kplace=kml_placemark();
        if (length(iloop)-1 > 1)
            kplace.name      =sprintf('Partition %d (%d loops)',k,length(iloop)-1);
        else
            kplace.name      =sprintf('Partition %d',k);
        end
        kplace.visibility=1;
        if exist('pdata','var')
            kplace.descript  =sprintf('Partition data: %g',pdata(k));
            imap = fix((pdata(k)-cmin)/(cmax-cmin)*size(cmap,1))+1;
            if     (imap >= 1) && (imap <= size(cmap,1))
                kplace.styleurl  =sprintf('#MatlabColor%d',imap);
            elseif (pdata(k) == cmax)
                kplace.styleurl  =sprintf('#MatlabColor%d',size(cmap,1));
            else
                kplace.styleurl  =sprintf('#BlackLineEmptyPoly');
            end
        else
            kplace.styleurl  =sprintf('#BlackLineRandomPoly');
        end

        kmgeom=kml_multigeometry();
        kmgeom.geometry  ={repmat(kml_linestring(),1,length(iloop)-1)};

%  loop over each loop of the perimeter for the given partition

        for i=1:length(iloop)-1
            kline=kml_linestring();
            kline.extrude   =1;
            kline.tessellate=1;
            kline.altmode   ='relativeToGround';
            kline.coords    =zeros(0,3);

            elast=0;
            nlast=0;
            slast=0;
            lat=[];
            long=[];

%  loop over the element edges on the loop of the partition

            j=iloop(i);
            while (j < iloop(i+1))
%  find which side of element is referenced in perimeter list
                for l=1:size(elemp,2)
                    if ((elemp(elemper(j),l)          == edgeper(j,1)) && ...
                        (elemp(elemper(j),mod(l,3)+1) == edgeper(j,2))) || ...
                       ((elemp(elemper(j),l)          == edgeper(j,2)) && ...
                        (elemp(elemper(j),mod(l,3)+1) == edgeper(j,1)))
                        jedge=l;
                        break
                    end
                end

%  check if element side connects nodes in partition
                if (epartp(elemper(j),jedge)          == k) && ...
                   (epartp(elemper(j),mod(jedge,3)+1) == k)
%  write out specified element side
%                    disp(['segment j=' int2str(j) ' unshared edge ' int2str(edgeper(j,1)) ' to ' int2str(edgeper(j,2)) ' on side ' int2str(jedge) ' from element ' int2str(elemper(j)) ' written.'])
%  if first edge, write out first node
                    if ~elast
                        kline.coords(end+1,:)=[md.mesh.long(edgeper(j,1)) md.mesh.lat(edgeper(j,1)) alt];
                    end
                    kline.coords(end+1,:)=[md.mesh.long(edgeper(j,2)) md.mesh.lat(edgeper(j,2)) alt];
                    elast=elemper(j);
                    nlast=edgeper(j,2);
                    slast=0;
                    j=j+1;

%  element not entirely within partition, so figure out boundary
                else
%                    disp(['segment j=' int2str(j) ' from element ' int2str(elemper(j)) ' shared by other partitions.'])
                    ielem=elemper(j);

%  follow partition boundary through elements not wholly in partition
%  (may include elements not in perimeter list)

                    while 1
%  if first edge, figure out direction from perimeter edge direction
                        if ~nlast && ~slast
                            nlast=find(elemp(ielem,:)==edgeper(j,1));
                            nnext=find(elemp(ielem,:)==edgeper(j,2));
                            if     (nlast+nnext == 3)
                                slast=1;
                            elseif (nlast+nnext == 5)
                                slast=2;
                            elseif (nlast+nnext == 4)
                                slast=3;
                            end
                            if     (nnext+(6-nlast-nnext) == 3)
                                snext=1;
                            elseif (nnext+(6-nlast-nnext) == 5)
                                snext=2;
                            elseif (nnext+(6-nlast-nnext) == 4)
                                snext=3;
                            end

%  find how many nodes of current element are in current partition
%  (1 or 2, not 3) and handle each permutation separately
                            ipart=find(epartp(ielem,:)==k);
%  two nodes are in current partition, so cut off other node
                            if (length(ipart) == 2)
                                switch 6-sum(ipart)
                                    case nlast
                                        slast=6-slast-snext;
                                        nlast=0;
                                    case nnext
                                        if (epartp(ielem,nnext) == k)
                                            nlast=nnext;
                                        end
                                    otherwise
                                        slast=6-slast-snext;
                                        nlast=0;
                                end
%  one node is in current partition
                            else
%  all different, so cut through centroid
                                if (epartp(ielem,1) ~= epartp(ielem,2)) && ...
                                   (epartp(ielem,2) ~= epartp(ielem,3)) && ...
                                   (epartp(ielem,3) ~= epartp(ielem,1))
                                    switch ipart
                                        case {nlast,nnext}
                                            if (epartp(ielem,nnext) == k)
                                                nlast=nnext;
                                            end
                                        otherwise
                                            slast=6-slast-snext;
                                            nlast=0;
                                    end
%  other two are in the same partition, so cut them off
                                else
                                    switch ipart
                                        case nlast
                                            if (epartp(ielem,nnext) == k)
                                                nlast=nnext;
                                            end
                                        case nnext
                                            slast=snext;
                                            nlast=0;
                                        otherwise
                                            slast=6-slast-snext;
                                            nlast=0;
                                    end
                                end
                            end

%  last edge exited last element at node
                            if nlast
%  write out first node of first side for half-edge to midpoint
%                                disp(['segment j=' int2str(j) ' unshared half edge from node ' int2str(elemp(ielem,nlast)) ' (node ' int2str(nlast) ') on side ' int2str(slast) ' from element ' int2str(ielem) ' written.'])
                                kline.coords(end+1,:)=[md.mesh.long(elemp(ielem,nlast)) ...
                                                       md.mesh.lat(elemp(ielem,nlast)) alt];
                            end
                            nlast=0;

%  write out midpoint of first side
                            kline.coords(end+1,:)=[(md.mesh.long(elemp(ielem,slast))...
                                                   +md.mesh.long(elemp(ielem,mod(slast,3)+1)))/2. ...
                                                   (md.mesh.lat(elemp(ielem,slast))...
                                                   +md.mesh.lat(elemp(ielem,mod(slast,3)+1)))/2. alt];
                        end

%  last edge exited last element at node
                        if nlast
                            if elast
%  find where last node on previous element occurs on current element
                                nlast=find(elemp(ielem,:)==nlast,1);
                            end
%  half-edge occurs on unshared side from current node (unique unless mesh
%  is only attached at node)
                            switch nlast
                                case 1
                                    if ~edgeadjp(ielem,1)
                                        nnext=2;
                                        slast=1;
                                    else
                                        nnext=3;
                                        slast=3;
                                    end
                                case 2
                                    if ~edgeadjp(ielem,2)
                                        nnext=3;
                                        slast=2;
                                    else
                                        nnext=1;
                                        slast=1;
                                    end
                                case 3
                                    if ~edgeadjp(ielem,3)
                                        nnext=1;
                                        slast=3;
                                    else
                                        nnext=2;
                                        slast=2;
                                    end
                            end
%  write out half-edge from current node to midpoint of unshared side
%                            disp(['segment j=' int2str(j) ' unshared half edge from node ' int2str(elemp(ielem,nlast)) ' (node ' int2str(nlast) ') on side ' int2str(slast) ' from element ' int2str(ielem) ' written.'])
                            kline.coords(end+1,:)=[(md.mesh.long(elemp(ielem,nlast))...
                                                   +md.mesh.long(elemp(ielem,nnext)))/2. ...
                                                   (md.mesh.lat(elemp(ielem,nlast))...
                                                   +md.mesh.lat(elemp(ielem,nnext)))/2. alt];
                            nlast=0;

%  last edge exited last element at midpoint of side
                        elseif slast
                            if elast
%  find where last side on previous element occurs on current element
                                slast=find(edgeadjp(ielem,:)==elast,1);
                            end
                        end

%  find how many nodes of current element are in current partition
%  (1 or 2, not 3) and handle each permutation separately
                        ipart=find(epartp(ielem,:)==k);
                        if (length(ipart) == 2)
%  two nodes are in current partition, so cut off other node
                            switch 6-sum(ipart)
                                case 1
                                    snext=3+1-slast;
                                case 2
                                    snext=1+2-slast;
                                case 3
                                    snext=2+3-slast;
                            end
                        else
                            if (epartp(ielem,1) ~= epartp(ielem,2)) && ...
                               (epartp(ielem,2) ~= epartp(ielem,3)) && ...
                               (epartp(ielem,3) ~= epartp(ielem,1))
%  all different, so cut through centroid
%                                disp(['element ielem=' int2str(ielem) ' centroid written.'])
                                kline.coords(end+1,:)=[sum(md.mesh.long(elemp(ielem,:)))/3. ...
                                                       sum(md.mesh.lat(elemp(ielem,:)))/3. alt];
                            end
%  one node is in current partition, so cut off other two nodes
                            switch ipart
                                case 1
                                    snext=3+1-slast;
                                case 2
                                    snext=1+2-slast;
                                case 3
                                    snext=2+3-slast;
                            end
                        end
%  write out midpoint of opposite side
%                        disp(['segment j=' int2str(j) ' internal edge from side ' int2str(slast) ' to side ' int2str(snext) ' from element ' int2str(ielem) ' written.'])
                        kline.coords(end+1,:)=[(md.mesh.long(elemp(ielem,snext))...
                                               +md.mesh.long(elemp(ielem,mod(snext,3)+1)))/2. ...
                                               (md.mesh.lat(elemp(ielem,snext))...
                                               +md.mesh.lat(elemp(ielem,mod(snext,3)+1)))/2. alt];
                        elast=ielem;
                        nlast=0;
                        slast=snext;
%  find adjacent element to opposite side
                        ielem=edgeadjp(elast,slast);
%  if opposite side is unshared, find it in edge perimeter list
                        if ~ielem
                            jlast=find(elemper(j:end)==elast)+j-1;
                            j=0;
                            for l=1:length(jlast)
                                if ((elemp(elast,slast)          == edgeper(jlast(l),1)) && ...
                                    (elemp(elast,mod(slast,3)+1) == edgeper(jlast(l),2))) || ...
                                   ((elemp(elast,slast)          == edgeper(jlast(l),2)) && ...
                                    (elemp(elast,mod(slast,3)+1) == edgeper(jlast(l),1)))
                                    j=jlast(l);
                                    break
                                end
                            end
                            if ~j
                                j=iloop(i+1)-1;
                            end
%  write out half-edge from midpoint of unshared side to node
                            if (epartp(elast,slast) == k)
                                nnext=slast;
                            else
                                nnext=mod(slast,3)+1;
                            end
%                            disp(['segment j=' int2str(j) ' unshared half edge on side ' int2str(slast) ' to node ' int2str(elemp(elast,nnext)) ' (node ' int2str(nnext) ') from element ' int2str(elast) ' written.'])
                            kline.coords(end+1,:)=[md.mesh.long(elemp(elast,nnext)) ...
                                                   md.mesh.lat(elemp(elast,nnext)) alt];
                            break
%  if not unshared, advance perimeter list and watch for end
                        else
                            if (elast == elemper(j))
                                if (j+1 < iloop(i+1)) && ...
                                   ~isempty(find(elemper(j+1:end)~=elast,1))
                                    j=j+find(elemper(j+1:end)~=elast,1);
                                else
                                    break
                                end
                            end
                        end
                    end
                    j=j+1;
                end
            end

            kmgeom.geometry{1}(i)=kline;
            clear kline
        end

        kplace.geometry=kmgeom;
        kfold.feature{1}(k)=kplace;
        clear kmgeom kplace
    end
end

end
