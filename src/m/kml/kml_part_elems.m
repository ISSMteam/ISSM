%%
%  create kml polygons for the partition elements.
%
%  [kfold]=kml_part_elems(md,params)
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
%    kfold         (kml_folder, folder of polygon placemarks)
%
function [kfold]=kml_part_elems(varargin)

if ~nargin
    help kml_part_elems
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

%  write folder for partition elements

if (~exist('prtplt','var') || strncmpi(prtplt,'on' ,2) || strncmpi(prtplt,'y',1)) && ...
    md.qmu.numberofpartitions
    kfold=kml_folder();
    kfold.name      ='Partition Elements';
    kfold.visibility=1;
    kfold.descript  =sprintf('Partitions=%d, Nodes=%d\n',...
        md.qmu.numberofpartitions,md.mesh.numberofvertices);
    kfold.feature   ={repmat(kml_placemark(),1,md.qmu.numberofpartitions)};

%  write each partition loop as a polygon multigeometry placemark

    disp(['Writing ' num2str(md.qmu.numberofpartitions) ' partitions as KML polygons.']);
    epart=md.qmu.partition(md.mesh.elements)+1;
    if exist('ndata','var') || exist('edata','var')
        pdata=zeros(1,md.qmu.numberofpartitions);
        pdata(:)=NaN;
    end

%  loop over each partition

    for k=1:md.qmu.numberofpartitions

%  for each partition, find all the included elements

        [icol,irow]=find(epart'==k);
        if isempty(irow)
            continue;
        end
        irow=unique(irow);
        elem=md.mesh.elements(irow,:);

%  determine the data to be used for the colors (if any)

        if exist('ndata','var')
            pdata(k)=ndata(find(md.qmu.partition+1==k,1));
        elseif exist('edata','var')
            for i=1:size(epart,1)
                if isempty(find(epart(i,:)~=k,1))
                    pdata(k)=edata(i);
                    break
                end
            end
            if isnan(pdata(k))
                warning('Data for Partition %d is not defined.\n',k)
            end
        end

%  set up the placemark with multigeometry

        kplace=kml_placemark();
        kplace.name      =sprintf('Partition %d (%d elements)',k,size(elem,1));
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
        kmgeom.geometry  ={repmat(kml_polygon(),1,size(elem,1))};

%  loop over each element for the given partition

        for i=1:size(elem,1)
            kpoly=kml_polygon();
            kpoly.extrude   =1;
            kpoly.altmode   ='relativeToGround';

            kring=kml_linearring();
            kring.coords    =zeros(size(elem,2)+1,3);

%  loop over the element nodes

            for j=1:size(elem,2)
                kring.coords(j,:)=[md.mesh.long(elem(i,j)) md.mesh.lat(elem(i,j)) alt];
            end
            kring.coords(end,:)=kring.coords(1,:);

            kpoly.outer=kring;
            kmgeom.geometry{1}(i)=kpoly;
            clear kring kpoly
        end

        kplace.geometry=kmgeom;
        kfold.feature{1}(k)=kplace;
        clear kmgeom kplace
    end
end

end
