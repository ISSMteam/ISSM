%
%  definition for the kml_polygon sub (derived) class.
%
%  [kml]=kml_polygon(varargin)
%
%  where the optional varargin and defaults are:
%    id            (char, polygon id, '')
%    extrude       (logical, extrusion, false)
%    tessellate    (logical, tessellation, false)
%    altmode       (char, altitude mode, 'clampToGround')
%    outer         (kml_linearring, outer boundary)
%    inner         (kml_linearring, inner boundaries)
%
%  note that zero arguments constructs a default instance; one
%  argument of the class copies the instance; and two or more
%  arguments constructs a new instance from the arguments.
%
classdef kml_polygon < kml_geometry
    properties
        extrude   =false;
        tessellate=false;
        altmode   ='clampToGround';
        outer     =kml_linearring.empty();
        inner     =kml_linearring.empty();
    end

    methods
        function [kml]=kml_polygon(varargin)

            kml=kml@kml_geometry(varargin{:});

            switch nargin

%  create a default object

                case 0

%  copy the object or create the object from the input

                otherwise
                    if (nargin == 1) && isa(varargin{1},class(kml))
                        kml=varargin{1};

                    else
                        fnames=fieldnames(kml_polygon());

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
                disp(sprintf('         outer: %s %s'   ,string_size(kml(i).outer),...
                             class(kml(i).outer)));
                disp(sprintf('         inner: %s %s\n' ,string_size(kml(i).inner),...
                             class(kml(i).inner)));
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
                    'outer' ...
                    'inner' ...
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

%  loop over the polygons

            for i=1:numel(kml)
                kmli=kml(i);
                if ~isempty(kmli.id)
                    fprintf(fid,'%s<Polygon id="%s">\n',indent,kmli.id);
                else
                    fprintf(fid,'%s<Polygon>\n',indent);
                end
                kml_write@kml_geometry(kmli,fid,indent);
                fprintf(fid,'%s  <extrude>%d</extrude>\n',indent,kmli.extrude);
                fprintf(fid,'%s  <tessellate>%d</tessellate>\n',indent,kmli.tessellate);
                fprintf(fid,'%s  <altitudeMode>%s</altitudeMode>\n',indent,kmli.altmode);
                fprintf(fid,'%s  <outerBoundaryIs>\n',indent);
                if isa(kmli.outer,'kml_linearring')
                    kml_write(kmli.outer,fid,[indent '    ']);
                else
                    warning('kml(%d).outer is a ''%s'' class object, not ''%s''.',...
                        i,class(kmli.outer),'kml_linearring');
                end
                fprintf(fid,'%s  </outerBoundaryIs>\n',indent);

%  loop over any inner boundaries for each polygon

                if isa(kmli.inner,'kml_linearring')
                    for j=1:numel(kmli.inner)
                        fprintf(fid,'%s  <innerBoundaryIs>\n',indent);
                        kml_write(kmli.inner(j),fid,[indent '    ']);
                        fprintf(fid,'%s  </innerBoundaryIs>\n',indent);
                    end
                else
                    warning('kml(%d).inner is a ''%s'' class object, not ''%s''.',...
                        i,class(kmli.inner),'kml_linearring');
                end

                fprintf(fid,'%s</Polygon>\n',indent);
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

%  loop over the polygons

            for i=1:numel(kml)
                kmli=kml(i);
                if ~isempty(kmli.id)
                    sbuf=add(sbuf,sprintf('%s<Polygon id="%s">\n',indent,kmli.id));
                else
                    sbuf=add(sbuf,sprintf('%s<Polygon>\n',indent));
                end
                sbuf=kml_swrite@kml_geometry(kmli,sbuf,indent);
                sbuf=add(sbuf,sprintf('%s  <extrude>%d</extrude>\n',indent,kmli.extrude));
                sbuf=add(sbuf,sprintf('%s  <tessellate>%d</tessellate>\n',indent,kmli.tessellate));
                sbuf=add(sbuf,sprintf('%s  <altitudeMode>%s</altitudeMode>\n',indent,kmli.altmode));
                sbuf=add(sbuf,sprintf('%s  <outerBoundaryIs>\n',indent));
                if isa(kmli.outer,'kml_linearring')
                    sbuf=kml_swrite(kmli.outer,sbuf,[indent '    ']);
                else
                    warning('kml(%d).outer is a ''%s'' class object, not ''%s''.',...
                        i,class(kmli.outer),'kml_linearring');
                end
                sbuf=add(sbuf,sprintf('%s  </outerBoundaryIs>\n',indent));

%  loop over any inner boundaries for each polygon

                if isa(kmli.inner,'kml_linearring')
                    for j=1:numel(kmli.inner)
                        sbuf=add(sbuf,sprintf('%s  <innerBoundaryIs>\n',indent));
                        sbuf=kml_swrite(kmli.inner(j),sbuf,[indent '    ']);
                        sbuf=add(sbuf,sprintf('%s  </innerBoundaryIs>\n',indent));
                    end
                else
                    warning('kml(%d).inner is a ''%s'' class object, not ''%s''.',...
                        i,class(kmli.inner),'kml_linearring');
                end

                sbuf=add(sbuf,sprintf('%s</Polygon>\n',indent));
            end

        end

%  delete the object

        function []=delete(kml)

%  loop over the polygons

            for i=numel(kml):-1:1
                kmli=kml(i);
                if isa(kmli.outer,'kml_linearring')
                    delete(kmli.outer);
                else
                    warning('kml(%d).outer is a ''%s'' class object, not ''%s''.',...
                        i,class(kmli.outer),'kml_linearring');
                end
                kmli.outer     =kml_linearring.empty();

%  loop over any inner boundaries for each polygon

                if isa(kmli.inner,'kml_linearring')
                    for j=numel(kmli.inner):-1:1
                        delete(kmli.inner(j));
                    end
                else
                    warning('kml(%d).inner is a ''%s'' class object, not ''%s''.',...
                        i,class(kmli.inner),'kml_linearring');
                end
                kmli.inner     =kml_linearring.empty();

            end

        end

    end

end
