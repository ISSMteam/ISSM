function I=enviread(varargin)
%enviread: Read binary image files using ENVI header information
%
%I=enviread('filename') where the header file is named filename.hdr
%and exists in the same directory. Otherwise use 
%
%I=enviread('filename','hdrfilename')
%
%The output structure I contains fields I.x, I.y, I.z and I.info
%containing the x-coordinate vector, y-coordinate vector,
%images data and header info, respectively. I.z will be in whatever
%number format (double, int, etc.) as in the envi file.
%
%Original version by Ian Howat, Ohio State Universtiy, ihowat@gmail.com
%Thanks to Yushin Ahn and Ray Jung
%
%% Current version 3 - adding subsetting feature to read only parts of
%% files.

file=varargin{1};
hdrfile=[deblank(file),'.hdr'];
info=read_envihdr(hdrfile);
sub = [1, info.samples, 1, info.lines];
if nargin==2
    if varargin{2}=='date'
        [pathstr, name, ext, versn] = fileparts(file);
        info.hdrname=hdrfile;
        info.ipath=pathstr;
        info.iname=name;
        info.ipname=file;
        info.year =name(1:4);
        info.month=name(5:6);
        info.day  =name(7:8);
        info.date =name(1:8);
        info.sfname=name(1:14);
        info.datenum=datenum(str2double(info.year),...
            str2double(info.month),str2double(info.day));
    end
end

% subset
if nargin == 3;
    if strcmp(varargin{2},'pixel_subset');
        sub = varargin{3};
    elseif strcmp(varargin{2},'map_subset');
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
    end
end

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
    %I.x = xm + (subE(1):subE(2)-1).*info.map_info.dx;
    %I.y = ym - (subE(3):subE(4)-1).*info.map_info.dy;
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
switch info.data_type
    case {1}
        format = 'uint8';
    case {2}
        format= 'int16';
    case{3}
        format= 'int32';
    case {4}
        format= 'single';
    case {5}
        format= 'double';
    case {6}
        disp('>> Sorry, Complex (2x32 bits)data currently not supported');
        disp('>> Importing as double-precision instead');
        format= 'double';
    case {9}
        error('Sorry, double-precision complex (2x64 bits) data currently not supported');
    case {12}
        format= 'uint16';
    case {13}
        format= 'uint32';
    case {14}
        format= 'int64';
    case {15}
        format= 'uint64';
    otherwise
        error(['File type number: ',num2str(dtype),' not supported']);
end

%% file read
% Version 2 code by Yushin Ahn - replaces resize calls with loops (except
% for BIP formats) to work on big arrays.

        %tmp=zeros(info.lines, info.samples,info.bands,format);
        tmp=zeros(sub(4)-sub(3)+1,sub(2)-sub(1)+1,info.bands,format);
        fid=fopen(file,'r');

switch lower(info.interleave)

    case {'bsq'}
        % Format:
        % [Band 1]       
        % R1: C1, C2, C3, ...
        % R2: C1, C2, C3, ...
        %  ...
        % RN: C1, C2, C3, ...
        %
        % [Band 2]
        %  ...
        % [Band N]

%% Old Script w/out subsetting
%         for b=1:info.bands
%             for i=1:info.lines
%                 t=fread(fid,info.samples,format);
%                 tmp(i,:,b)=t;    
%             end
%         end
        
%% New Subsetting script - IMH
        offset1=(sub(3)-1)*info.samples;
        fseek(fid,offset1,'bof');
        for b=1:info.bands
            for i=sub(3):sub(4)
                t=fread(fid,info.samples,format);
                %if i >=  sub(3) && i <= sub(4)
                    tmp(i-sub(3)+1,:,b)=t(sub(1):sub(2));
                %end
            end
            offset2 = info.samples*info.lines*b+offset1;
            fseek(fid,offset2,'bof');
        end

    case {'bil'}
        % Format:        
        % [Row 1]      
        % B1: C1, C2, C3, ...
        % B2: C1, C2, C3, ...
        %
        %  ...
        % [Row N]

%         for i=1:info.lines
%             for b=1:info.bands
%                  t=fread(fid,info.samples,format);
%                 tmp(i,:,b)=t;       
%             end
%         end
      
%% New Subsetting script - IMH
        for i=1:sub(4) 
            for b=1:info.bands
                t=fread(fid,info.samples,format);
                if i >=  sub(3) && i <= sub(4) 
                    tmp(i-sub(3)+1,:,b)=t(sub(1):sub(2));    
                end
            end
        end



    case {'bip'}
    
        % Row 1
        % C1: B1 B2 B3, ...
        % C2: B1 B2 B3, ...
        % ...
        % Row N
        %This section authored by Ray Jung, APL-Johns Hopkins
        Z = fread(fid,info.samples*info.lines*info.bands,format,0,machine);  
        Z = reshape(Z, [info.bands, info.samples, info.lines]);

        for k=1:info.bands
            tmp(:,:,k) = squeeze(Z(k,:,:))';
        end     
end
fclose(fid);

% tmp2=zeros(subE(4)-subE(3)+1,subE(2)-subE(1)+1,info.bands,format);
% [a,b]=size(tmp);
% for i=1:info.bands
%     tmp2(xoff+1:a,yoff+1:b,1)=tmp(1:a,1:b,1);
% end
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
