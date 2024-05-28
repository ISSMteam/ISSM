%
%  write a kml file of the kml objects.
%
%  [fid]=kml_file_write(kobj,filek,indent)
%
%  where the required input is:
%    kobj          (kml_object, kml object to be written)
%     or
%    kobj          (cell array, array of kml objects)
%
%  the optional input is:
%    filek         (char, name of .kml file)
%     or
%    filek         (numeric, file ID of already-open file,
%                            noting 1=stdout and 2=stderr)
%    indent        (char, indention string)
%
%  and the optional output is:
%    fid           (numeric, file ID of still-open file)
%
function [fid]=kml_file_write(kobj,filek,indent)

if ~nargin
    help kml_file_write
    return
end

%%  process input data

if ~iscell(kobj)
    kobj={kobj};
end

fid=0;
if ~exist('filek' ,'var') || (~ischar(filek) && ~isnumeric(filek))
    filek='';
elseif ischar(filek)
    if     strcmpi(filek,'stdout')
        fid=1;
    elseif strcmpi(filek,'stderr')
        fid=2;
    end
elseif isnumeric(filek)
    fid=filek;
    if     (fid == 1)
        filek='stdout';
    elseif (fid == 2)
        filek='stderr';
    else
        filek='';
    end
end

if ~exist('indent','var') || ~ischar(indent)
    indent='  ';
end

%%  write kml file

%  open file and write header data (if necessary)

if ~fid
    if isempty(filek)
        filek=input('kml file to write?  ','s');
    end
    [pathstr,name,ext,versn] = fileparts(filek);
    if isempty(ext)
        ext='.kml';
    end
    filek=fullfile(pathstr,[name ext versn]);

    display(sprintf('Opening kml file ''%s''.',filek));
    fid=fopen(sprintf('%s',filek),'w');
    if (fid < 0)
        error('File ''%s'' could not be opened.',filek);
    end

    fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
    fprintf(fid,'<kml xmlns="http://www.opengis.net/kml/2.2">\n');
end

%  write kml objects

if ~isempty(filek)
    display(sprintf('Writing to kml file ''%s'':',filek));
else
    display(sprintf('Writing to kml file id=%d:',fid));
end
for i=1:numel(kobj)
    if isa(kobj{i},'kml_object')
        display(sprintf('  Writing object %d of class ''%s'' and size %s.',...
            i,class(kobj{i}),string_size(kobj{i})));
        kml_write(kobj{i},fid,indent);
    else
        if ~isempty(inputname(1))
            warning('Object ''%s{%d}'' is a ''%s'' class object, not ''%s''.',...
                inputname(1),i,class(kobj{i}),'kml_object');
        else
            warning('Object {%d} is a ''%s'' class object, not ''%s''.',...
                             i,class(kobj{i}),'kml_object');
        end
    end
end

%  write trailer data and close file (if necessary)

if ~nargout && (fid >= 3)
    fprintf(fid,'</kml>\n');

    if (fclose(fid) < 0)
        if ~isempty(filek)
            error('File ''%s'' could not be closed.',filek);
        else
            error('File id=%d could not be closed.',fid);
        end
    else
        if ~isempty(filek)
            disp(['End of file ''' filek ''' successfully written.']);
        else
            disp(['End of file successfully written.']);
        end
    end
end

end
