%
%  definition for the kml_linestring sub (derived) class.
%
%  [kml]=kml_linestring(varargin)
%
%  where the optional varargin and defaults are:
%    id            (char, linestring id, '')
%    extrude       (logical, extrusion, false)
%    tessellate    (logical, tessellation, false)
%    altmode       (char, altitude mode, 'clampToGround')
%    coords        (numeric, long/lat/alt (n x 3), empty)
%
%  note that zero arguments constructs a default instance; one
%  argument of the class copies the instance; and two or more
%  arguments constructs a new instance from the arguments.
%
classdef kml_linestring < kml_geometry
    properties
        extrude   =false;
        tessellate=false;
        altmode   ='clampToGround';
        coords    =zeros(0,3);
    end

    methods
        function [kml]=kml_linestring(varargin)

            kml=kml@kml_geometry(varargin{:});

            switch nargin

%  create a default object

                case 0

%  copy the object or create the object from the input

                otherwise
                    if (nargin == 1) && isa(varargin{1},class(kml))
                        kml=varargin{1};

                    else
                        fnames=fieldnames(kml_linestring());

                        for i=length(fieldnames(kml_geometry()))+1:min(nargin,length(fnames))
                            if isa(varargin{i},class(kml.(fnames{i})))
                                if ~isempty(varargin{i})
                                    kml.(fnames{i})=varargin{i};
                                end
                            else
                                if ~isempty(inputname(i))
                                    warning('Argument ''%s'' for property ''%s'' is a ''%s'' class object, not ''%s''.',...
                                        inputname(i),fnames{i},class(varargin{i}),class(kml.(fnames{i})));
                                else
                                    warning('Argument %d for property ''%s'' is a ''%s'' class object, not ''%s''.',...
                                        i           ,fnames{i},class(varargin{i}),class(kml.(fnames{i})));
                                end
                            end
                        end
                    end

            end

        end

%  display the object

        function []=disp(kml)

            for i=1:numel(kml)
                disp(sprintf('class ''%s'' object ''%s%s'' = \n',...
                    class(kml),inputname(1),string_dim(kml,i)));
                disp@kml_geometry(kml(i));
                disp(sprintf('       extrude: %g'      ,kml(i).extrude));
                disp(sprintf('    tessellate: %g'      ,kml(i).tessellate));
                disp(sprintf('       altmode: ''%s'''  ,kml(i).altmode));
                disp(sprintf('        coords: %s %s\n' ,string_size(kml(i).coords),...
                             class(kml(i).coords)));
            end

        end

%  return the fieldnames of the object

        function [fnames]=fieldnames(kml)

%  fieldnames for a sub (derived) class list those before super (base)

            fnames=fieldnames(kml_geometry());
            fnames={fnames{:} ...
                    'extrude' ...
                    'tessellate' ...
                    'altmode' ...
                    'coords' ...
                   }';

        end

%  set the properties of the object

        function [kml]=setprops(kml,varargin)

            kmlref=feval(class(kml));
            fnames=fieldnames(kmlref);

%  loop through each parameter in the input list (comparing to the reference
%  object in case property types have been changed)

            for i=1:2:length(varargin)
                if ismember(varargin{i},fnames) && (i+1 <= length(varargin))
                    if isa(varargin{i+1},class(kmlref.(varargin{i})))
                        kml.(varargin{i})=varargin{i+1};
                    else
                        if ~isempty(inputname(i+1))
                            warning('Argument ''%s'' for property ''%s'' is a ''%s'' class object, not ''%s''.',...
                                inputname(i+2),varargin{i},class(varargin{i+1}),class(kmlref.(varargin{i})));
                        else
                            warning('Argument %d for property ''%s'' is a ''%s'' class object, not ''%s''.',...
                                i+2           ,varargin{i},class(varargin{i+1}),class(kmlref.(varargin{i})));
                        end
                    end
                else
                    warning('Property ''%s'' for class ''%s'' does not exist.',...
                        varargin{i},class(kmlref));
                end
            end

        end

%  write the object

        function []=kml_write(kml,fid,indent)

            if ~exist('fid','var') || isempty(fid)
                fid=1;
            end
            if ~exist('indent','var') || isempty(indent)
                indent='';
            end

%  loop over the linestrings

            for i=1:numel(kml)
                kmli=kml(i);
                if ~isempty(kmli.id)
                    fprintf(fid,'%s<LineString id="%s">\n',indent,kmli.id);
                else
                    fprintf(fid,'%s<LineString>\n',indent);
                end
                kml_write@kml_geometry(kmli,fid,indent);
                fprintf(fid,'%s  <extrude>%d</extrude>\n',indent,kmli.extrude);
                fprintf(fid,'%s  <tessellate>%d</tessellate>\n',indent,kmli.tessellate);
                fprintf(fid,'%s  <altitudeMode>%s</altitudeMode>\n',indent,kmli.altmode);
                fprintf(fid,'%s  <coordinates>\n',indent);

%  loop over the coordinates for each linestring

                for j=1:size(kmli.coords,1)
                    fprintf(fid,'%s    %0.16g,%0.16g,%0.16g\n',indent,kmli.coords(j,:));
                end

                fprintf(fid,'%s  </coordinates>\n',indent);
                fprintf(fid,'%s</LineString>\n',indent);
            end

        end

%  string write the object

        function [sbuf]=kml_swrite(kml,sbuf,indent)

            if ~exist('sbuf','var') || isempty(sbuf)
                sbuf=string_buf;
            end
            if ~exist('indent','var') || isempty(indent)
                indent='';
            end

%  loop over the linestrings

            for i=1:numel(kml)
                kmli=kml(i);
                if ~isempty(kmli.id)
                    sbuf=add(sbuf,sprintf('%s<LineString id="%s">\n',indent,kmli.id));
                else
                    sbuf=add(sbuf,sprintf('%s<LineString>\n',indent));
                end
                sbuf=kml_swrite@kml_geometry(kmli,sbuf,indent);
                sbuf=add(sbuf,sprintf('%s  <extrude>%d</extrude>\n',indent,kmli.extrude));
                sbuf=add(sbuf,sprintf('%s  <tessellate>%d</tessellate>\n',indent,kmli.tessellate));
                sbuf=add(sbuf,sprintf('%s  <altitudeMode>%s</altitudeMode>\n',indent,kmli.altmode));
                sbuf=add(sbuf,sprintf('%s  <coordinates>\n',indent));

%  loop over the coordinates for each linestring

                for j=1:size(kmli.coords,1)
                    sbuf=add(sbuf,sprintf('%s    %0.16g,%0.16g,%0.16g\n',indent,kmli.coords(j,:)));
                end

                sbuf=add(sbuf,sprintf('%s  </coordinates>\n',indent));
                sbuf=add(sbuf,sprintf('%s</LineString>\n',indent));
            end

        end

    end

end
