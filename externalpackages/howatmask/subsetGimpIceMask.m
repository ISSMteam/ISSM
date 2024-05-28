function [M,Mx,My] = subsetGimpIceMask(x0,x1,y0,y1)
% subsetGimpIceMask samples mask tiles for a given coordinate range
%
%   [M,Mx,My] = subsetGimpIceMask(x0,x1,y0,y1) returns the mask in array M
%   with map vectors Mx and My within the bounding box specified by
%   coordinate ranges x0 < x1 and y0 < y1. The function finds tiles
%   overlappimg the box, samples the tiles and mosaics the samples into one
%   array.
%
% 	Example: Specify the following coordinate range and map the result.
%   >> x0 = 260000; x1 = 377900;  y0 = -2609400; y1 = -2543000;
%   >> [M,Mx,My] = subsetGimpIceMask(x0,x1,y0,y1);
%   >> imagesc(Mx,My,M); axis equal xy; colormap gray;
%
%   Ian Howat, Ohio State University, ihowat@gmail.com
%   The Greenland Ice Mapping Project is funded by NASA.
%   $Revision: 0 $  $Date: 09-Nov-2011 14:48:24$
%

%make coordinate index
nt = 6;
ps = 15;
xn0 = -640000 ; %xn1 = 857550;
yn0 =  -3313350; %yn1 = -665450;

c = 16640;
r = 29420;

x = xn0:ps:xn0 + c.*ps.*nt;
y = yn0:ps:yn0 + r.*ps.*nt;

X = cell(1,nt);
Y = cell(nt,1);
i=1;
for i=1:nt;
    
    X{i} = [x((i-1)*c+1),x(i*c)];
    Y{i} = [y((i-1)*r+1),y(i*r)];
    
end

X = repmat(X,[nt,1]);
Y = repmat(Y,[1,nt]);
N = cell(size(X));


i=1;
for i=1:numel(X)
    N{i} = [[X{i}(1),Y{i}(2)];...
        [X{i}(2),Y{i}(2)];...
        [X{i}(2),Y{i}(1)];...
        [X{i}(1),Y{i}(1)]];
end


%% find overlapping tiles
p = [[x0,y1];[x1,y1];[x1,y0];[x0,y0]];

n = zeros(size(N));

i=1;
for i=1:numel(N)
    n(i) = any(inpolygon(p(:,1),p(:,2),N{i}(:,1),N{i}(:,2)));
end
i=1;
for i=1:size(n,1)
    n(i,find(n(i,:),1,'first'):find(n(i,:),1,'last')) = 1;
end
i=1;
for i=1:size(n,2)
    n(find(n(:,i),1,'first'):find(n(:,i),1,'last'),i) = 1;
end

[i,j] = find(n);
row = (min(i):max(i))';
col = min(j):max(j);
col = repmat(col,[length(row),1])-1;
row = repmat(row,[1,size(col,2)])-1;

% make output cells
M = cell(size(row));
Mx = M;
My = M;

% read each overlapping file and populate output cell
j=1;
for j=1:size(col,2)
    i=1;
    for i=1:size(row,1);
        m = enviread(['GimpIceMask_15m_tile',num2str(col(i,j)),'_',...
            num2str(row(i,j))],'map_subset',[x0 x1 y0 y1]);
        M{end-i+1,j} = m.z;
        Mx{end-i+1,j} = m.x;
        My{end-i+1,j} = m.y';
        
        
    end
end

% put together cells
M = cell2mat(M);
Mx =cell2mat(Mx(1,:));
My =cell2mat(My(:,1));


function I=enviread(varargin)

file=varargin{1};
hdrfile=[deblank(file),'.hdr'];
info=read_envihdr(hdrfile);
sub = [1, info.samples, 1, info.lines];

sub = varargin{3};
subx = (sub(1:2)-info.map_info.mapx)./info.map_info.dx;
suby = (info.map_info.mapy - sub(3:4))./info.map_info.dy;
subx = round(subx);
suby = round(suby);

subx(subx < 1) = 1;
suby(suby < 1) = 1;
subx(subx > info.samples) = info.samples;
suby(suby > info.lines)   = info.lines;

sub  = [subx,suby];

sub(1:2) = sort(sub(1:2));
sub(3:4) = sort(sub(3:4));

%% Make geo-location vectors
if isfield(info.map_info,'mapx') && isfield(info.map_info,'mapy')
    xi = info.map_info.image_coords(1);
    yi = info.map_info.image_coords(2);
    xm = info.map_info.mapx;
    ym = info.map_info.mapy;
    %adjust points to corner (1.5,1.5)
    if yi > 1.5
        ym =  ym + ((yi*info.map_info.dy)-info.map_info.dy);
    end
    if xi > 1.5
        xm = xm - ((xi*info.map_info.dy)-info.map_info.dx);
    end
    
    I.x = xm + ((0:info.samples-1).*info.map_info.dx);
    I.y = ym - ((0:info.lines-1).*  info.map_info.dy);
    
end

I.x = I.x(sub(1):sub(2));
I.y = I.y(sub(3):sub(4));

%% Set binary format parameters
switch info.byte_order
    case {0}
        machine = 'ieee-le';
    case {1}
        machine = 'ieee-be';
    otherwise
        machine = 'n';
end

format = 'uint8';

tmp=zeros(sub(4)-sub(3)+1,sub(2)-sub(1)+1,info.bands,format);
fid=fopen(file,'r');

offset1=(sub(3)-1)*info.samples;
fseek(fid,offset1,'bof');
for b=1:info.bands
    for i=sub(3):sub(4)
        t=fread(fid,info.samples,format);
        tmp(i-sub(3)+1,:,b)=t(sub(1):sub(2));
    end
    offset2 = info.samples*info.lines*b+offset1;
    fseek(fid,offset2,'bof');
end

fclose(fid);

I.z=tmp;
I.info =info;

%% sub function
function info = read_envihdr(hdrfile)
% READ_ENVIHDR read and return ENVI image file header information.
%   INFO = READ_ENVIHDR('HDR_FILE') reads the ASCII ENVI-generated image
%   header file and returns all the information in a structure of
%   parameters.
%
%   Example:
%   >> info = read_envihdr('my_envi_image.hdr')
%   info =
%          description: [1x101 char]
%              samples: 658
%                lines: 749
%                bands: 3
%        header_offset: 0
%            file_type: 'ENVI Standard'
%            data_type: 4
%           interleave: 'bsq'
%          sensor_type: 'Unknown'
%           byte_order: 0
%             map_info: [1x1 struct]
%      projection_info: [1x102 char]
%     wavelength_units: 'Unknown'
%           pixel_size: [1x1 struct]
%           band_names: [1x154 char]
%
%   NOTE: This function is used by ENVIREAD to import data.
% Ian M. Howat, Applied Physics Lab, University of Washington
% ihowat@apl.washington.edu
% Version 1: 19-Jul-2007 00:50:57
fid = fopen(hdrfile);
while fid;
    line = fgetl(fid);
    if line == -1
        break
    else
        eqsn = findstr(line,'=');
        if ~isempty(eqsn)
            param = strtrim(line(1:eqsn-1));
            param(findstr(param,' ')) = '_';
            value = strtrim(line(eqsn+1:end));
            if isempty(str2num(value))
                if ~isempty(findstr(value,'{')) && isempty(findstr(value,'}'))
                    while isempty(findstr(value,'}'))
                        line = fgetl(fid);
                        value = [value,strtrim(line)];
                    end
                end
                eval(['info.',param,' = ''',value,''';'])
            else
                eval(['info.',param,' = ',value,';'])
            end
        end
    end
end
fclose(fid);

if isfield(info,'map_info')
    line = info.map_info;
    line(line == '{' | line == '}') = [];
    line = strtrim(split(line,','));
    info.map_info = [];
    info.map_info.projection = line{1};
    info.map_info.image_coords = [str2num(line{2}),str2num(line{3})];
    info.map_info.mapx = str2num(line{4});
    info.map_info.mapy = str2num(line{5});
    info.map_info.dx  = str2num(line{6});
    info.map_info.dy  = str2num(line{7});
    if length(line) == 9
        info.map_info.datum  = line{8};
        info.map_info.units  = line{9}(7:end);
    elseif length(line) == 11
        info.map_info.zone  = str2num(line{8});
        info.map_info.hemi  = line{9};
        info.map_info.datum  = line{10};
        info.map_info.units  = line{11}(7:end);
    end
end

if isfield(info,'pixel_size')
    line = info.pixel_size;
    line(line == '{' | line == '}') = [];
    line = strtrim(split(line,','));
    info.pixel_size = [];
    info.pixel_size.x = str2num(line{1});
    info.pixel_size.y = str2num(line{2});
    info.pixel_size.units = line{3}(7:end);
end

%%
function A = split(s,d)
%This function by Gerald Dalley (dalleyg@mit.edu), 2004
A = {};
while (length(s) > 0)
    [t,s] = strtok(s,d);
    A = {A{:}, t};
end




