%%
%  create kml polygons for the element mesh.
%
%  [kfold]=kml_mesh_elem(md,params)
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
%
%  and the required output is:
%    kfold         (kml_folder, folder of polygon placemarks)
%
function [kfold]=kml_mesh_elem(varargin)

if ~nargin
    help kml_mesh_elem
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
                 'cmin','cmax','cmap'},...
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

%%  write folder for mesh

kfold=kml_folder();
if exist('cdata','var') && ~isempty(cdata)
    kfold.name      =sprintf('Data: %s',cdata);
else
    kfold.name      =sprintf('Mesh');
end
kfold.visibility=1;
kfold.descript  =sprintf('Elements=%d, Nodes=%d',...
    md.mesh.numberofelements,md.mesh.numberofvertices);
% see matlab_oop, "initializing a handle object array"
%kfold.feature   ={repmat(kml_placemark(),1,size(md.mesh.elements,1))};
kfeat(size(md.mesh.elements,1))=kml_placemark();
kfold.feature={kfeat};

%  write each element as a polygon placemark

disp(['Writing ' num2str(size(md.mesh.elements,1)) ' tria elements as KML polygons.']);
for i=1:size(md.mesh.elements,1)
    kplace=kml_placemark();
    kplace.name      =sprintf('Element %d',i);
    kplace.visibility=1;
    if exist('edata','var')
%        kplace.descript  =sprintf('Element data: %g',edata(i));
        kplace.descript  =sprintf('campaign{\n  deformation 1 %g quad_pol ascending right asap;\n}',edata(i));
        imap = fix((edata(i)-cmin)/(cmax-cmin)*size(cmap,1))+1;
        if     (imap >= 1) && (imap <= size(cmap,1))
            kplace.styleurl  =sprintf('#MatlabColor%d',imap);
        elseif (edata(i) == cmax)
            kplace.styleurl  =sprintf('#MatlabColor%d',size(cmap,1));
        else
            kplace.styleurl  =sprintf('#BlackLineEmptyPoly');
        end
    else
        kplace.styleurl  =sprintf('#BlackLineRandomPoly');
    end

    kpoly=kml_polygon();
    kpoly.extrude   =1;
    kpoly.altmode   ='relativeToGround';

    kring=kml_linearring();
    kring.coords    =zeros(size(md.mesh.elements,2)+1,3);

    for j=1:size(md.mesh.elements,2)
        kring.coords(j,:)=[md.mesh.long(md.mesh.elements(i,j)) md.mesh.lat(md.mesh.elements(i,j)) alt];
    end
    kring.coords(end,:)=kring.coords(1,:);

    kpoly.outer=kring;
    kplace.geometry=kpoly;
    kfold.feature{1}(i)=kplace;
    clear kring kpoly kplace

    if ~mod(i,1000)
        disp(['  ' num2str(i) ' tria elements written.']);
    end
end

end
