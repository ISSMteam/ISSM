%%
%  create kml linestrings for the flagged partition edges.
%
%  [kfold]=kml_part_flagedges(md,params)
%
%  where the required input is:
%    md            (model, model class object)
%
%  the optional input is:
%    params        (string/numeric, parameter names and values)
%
%  and the optional input is:
%    latsgn        (numeric, +1/-1 for north/south latitude)
%    alt           (numeric, altitude for polygons, default 10000)
%    prtplt        (char, 'off'/'no' for partition edge plot)
%
%  and the required output is:
%    kfold         (kml_folder, folder of linestring placemarks)
%
function [kfold]=kml_part_flagedges(varargin)

if ~nargin
    help kml_part_flagedges
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
                {'latsgn','alt','prtplt'},...
                'exact'))
            eval([varargin{iarg} '=varargin{iarg+1};']);
            disp([varargin{iarg} '=' any2str(varargin{iarg+1},20) ';']);
        else
            warning([varargin{iarg} '=' any2str(varargin{iarg+1},20) ' is not recognized.']);
        end
    else
        error(['''' any2str(varargin{iarg}) ''' is not a parameter name.']);
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

if ~exist('alt','var')
    alt=10000;
end

%%  write folder for partition segments

if (~exist('prtplt','var') || strncmpi(prtplt,'on' ,2) || strncmpi(prtplt,'y',1)) && ...
    md.qmu.numberofpartitions
    [latseg,lonseg]=flagedges(md.mesh.elements,md.mesh.lat,md.mesh.long,md.qmu.partition);
    kfold=kml_folder();
    kfold.name      ='Partition Segments';
    kfold.visibility=1;
    kfold.descript  =sprintf('Partitions=%d, Segments=%d',...
        md.qmu.numberofpartitions,size(latseg,1));
    kfold.feature   ={repmat(kml_placemark(),1,size(latseg,1))};

%  write each segment as a linestring placemark

    disp(['Writing ' num2str(size(latseg,1)) ' partition segments as KML linestrings.']);
    for i=1:size(latseg,1)
        kplace=kml_placemark();
        kplace.name      =sprintf('Segment %d',i);
        kplace.visibility=1;
        kplace.styleurl  ='#RedLineRedPoly';

        kline=kml_linestring();
        kline.extrude   =1;
        kline.tessellate=1;
        kline.altmode   ='relativeToGround';
        kline.coords    =zeros(2,3);

        for j=1:2
            kline.coords(j,:)=[lonseg(i,j) latseg(i,j) alt];
        end

        kplace.geometry=kline;
        kfold.feature{1}(i)=kplace;
        clear kline kplace
    end
end

end
