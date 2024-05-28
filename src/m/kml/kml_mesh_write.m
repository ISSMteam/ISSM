%%
%  write a kml file of the mesh from the model.
%
%  []=kml_mesh_write(filek,md,params)
%
%  where the required input is:
%    filek         (char, name of .kml file)
%    md            (model, model class object)
%
%  the optional input is:
%    params        (string/numeric, parameter names and values)
%
%  and the optional input is:
%    latsgn        (numeric, +1/-1 for north/south latitude)
%    data          (numeric, element or nodal results data)
%    alt           (numeric, altitude for polygons, default 10000)
%    lwidth        (numeric, line width in pixels, default 1)
%    popac         (numeric, polygon opacity, default 0.50)
%    cmin          (numeric, minimum of color map)
%    cmax          (numeric, maximum of color map)
%    cmap          (char or numeric, colormap definition)
%    prtplt        (char, 'off'/'no' for partition segment plot)
%
function []=kml_mesh_write(varargin)

if ~nargin
    help kml_mesh_write
    return
end

%%  process input data

iarg=1;
if (nargin >= 1)
    filek=varargin{1};
end

iarg=iarg+1;
if (nargin >= 2)
    md=varargin{2};
end
if ~exist('md','var') || isempty(md) || ~isa(md,'model')
    error(['Model ''' inputname(iarg) ''' is unrecognized class ''' class(md) '''.']);
end

%  parameters

iarg=iarg+1;
while (iarg <= nargin-1)
    if ischar(varargin{iarg})
        if ~isempty(strmatch(varargin{iarg},...
                {'latsgn','data','alt','lwidth','popac',...
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

%%  construct kml document

kdoc=kml_document();
kdoc.name      =sprintf('ISSM Mesh: %s',md.miscellaneous.name);
kdoc.open      =1;
ifirst=true;
for i=1:numel(md.miscellaneous.notes)
    if ~isempty(md.miscellaneous.notes{i})
        if ~ifirst
            kdoc.descript  =[kdoc.descript sprintf('\n')];
        end
        ifirst=false;
        kdoc.descript  =[kdoc.descript sprintf('%s',md.miscellaneous.notes{i})];
    end
end
clear ifirst
kdoc.style     ={repmat(kml_style(),0,0)};
kdoc.feature   ={repmat(kml_folder(),0,0)};

%  write style templates for defaults and for each color of the matlab
%  colormap (note that matlab colormap format is rgb, where each varies
%  from 0 to 1, whereas the kml color format is aabbggrr, where each
%  varies from 00 to ff.)

if ~exist('lwidth','var')
    lwidth=1;
end
if ~exist('popac','var')
    popac=0.50;
end

klsty=kml_linestyle();
klsty.color     ='ff000000';
klsty.colormode ='normal';
klsty.width     =lwidth;
kpsty=kml_polystyle();
kpsty.color     =sprintf('%02xffffff',round(popac*255));
kpsty.colormode ='random';
kstyle=kml_style();
kstyle.id        =sprintf('BlackLineRandomPoly');
kstyle.line      =klsty;
kstyle.poly      =kpsty;
kdoc.style{1}(end+1)=kstyle;
clear kstyle kpsty klsty

klsty=kml_linestyle();
klsty.color     ='ff000000';
klsty.colormode ='normal';
klsty.width     =lwidth;
kpsty=kml_polystyle();
kpsty.color     =sprintf('00ffffff');
kpsty.colormode ='random';
kstyle=kml_style();
kstyle.id        =sprintf('BlackLineEmptyPoly');
kstyle.line      =klsty;
kstyle.poly      =kpsty;
kdoc.style{1}(end+1)=kstyle;
clear kstyle kpsty klsty

klsty=kml_linestyle();
klsty.color     ='ff0000ff';
klsty.colormode ='normal';
klsty.width     =lwidth;
kpsty=kml_polystyle();
kpsty.color     =sprintf('%02x0000ff',round(popac*255));
kpsty.colormode ='random';
kstyle=kml_style();
kstyle.id        =sprintf('RedLineRedPoly');
kstyle.line      =klsty;
kstyle.poly      =kpsty;
kdoc.style{1}(end+1)=kstyle;
clear kstyle kpsty klsty

%  colormap command operates on a figure, so create an invisible one
%  (could also directly call colormaps, e.g. jet(64), but risky)

if exist('edata','var')
    hfig=figure('Visible','off');
    if exist('cmap','var')
        colormap(cmap)
    end
    cmap=colormap;
    close(hfig)

    disp(['Writing ' num2str(size(cmap,1)) ' Matlab colors as KML style templates.']);
    for i=1:size(cmap,1)
        klsty=kml_linestyle();
        klsty.color     ='ff000000';
        klsty.colormode ='normal';
        klsty.width     =lwidth;
        kpsty=kml_polystyle();
        kpsty.color     =sprintf('%02x%02x%02x%02x',round(popac*255),...
            round(cmap(i,3)*255),round(cmap(i,2)*255),round(cmap(i,1)*255));
        kpsty.colormode ='normal';
        kstyle=kml_style();
        kstyle.id        =sprintf('MatlabColor%d',i);
        kstyle.line      =klsty;
        kstyle.poly      =kpsty;
        kdoc.style{1}(end+1)=kstyle;
        clear kstyle kpsty klsty
    end
end

%  write folder for mesh

kdoc.feature{1}(end+1)=kml_mesh_elem(md,varargin{3:end});

%  write folder for partition segments

if (~exist('prtplt','var') || strncmpi(prtplt,'on' ,2) || strncmpi(prtplt,'y',1)) && ...
    md.qmu.numberofpartitions
    kdoc.feature{1}(end+1)=kml_part_flagedges(md,varargin{3:end});
end

%  write folder for unshared edges

if (~exist('prtplt','var') || strncmpi(prtplt,'on' ,2) || strncmpi(prtplt,'y',1)) && ...
    md.qmu.numberofpartitions
    kdoc.feature{1}(end+1)=kml_unsh_edges(md,varargin{3:end});
end

%  write folder for partition elements

if (~exist('prtplt','var') || strncmpi(prtplt,'on' ,2) || strncmpi(prtplt,'y',1)) && ...
    md.qmu.numberofpartitions
    kdoc.feature{1}(end+1)=kml_part_elems(md,varargin{3:end});
end

%  write folder for partition edges

if (~exist('prtplt','var') || strncmpi(prtplt,'on' ,2) || strncmpi(prtplt,'y',1)) && ...
    md.qmu.numberofpartitions
    kdoc.feature{1}(end+1)=kml_part_edges(md,varargin{3:end});
end

%  write folder for partitions

if (~exist('prtplt','var') || strncmpi(prtplt,'on' ,2) || strncmpi(prtplt,'y',1)) && ...
    md.qmu.numberofpartitions
    kdoc.feature{1}(end+1)=kml_partitions(md,varargin{3:end});
end

%%  write kml file

kml_file_write(kdoc,filek);
% kml_file_swrite(kdoc,filek);
delete(kdoc);

end
