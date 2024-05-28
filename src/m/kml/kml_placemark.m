%
%  definition for the kml_placemark sub (derived) class.
%
%  [kml]=kml_placemark(varargin)
%
%  where the optional varargin and defaults are:
%    id            (char, placemark id, '')
%    name          (char, name, '')
%    visibility    (logical, visibility, true)
%    open          (logical, open, false)
%    snippet       (char, snippet, '')
%    descript      (char, description, '')
%    styleurl      (char, style url, '')
%    style         (cell array, styles)
%    geometry      (kml_geometry, placemark geometry)
%
%  note that zero arguments constructs a default instance; one
%  argument of the class copies the instance; and two or more
%  arguments constructs a new instance from the arguments.
%
classdef kml_placemark < kml_feature
    properties
        geometry  =kml_geometry.empty();
    end

    methods
        function [kml]=kml_placemark(varargin)

            kml=kml@kml_feature(varargin{:});

            switch nargin

%  create a default object

                case 0

%  copy the object or create the object from the input

                otherwise
                    if (nargin == 1) && isa(varargin{1},class(kml))
                        kml=varargin{1};

                    else
                        fnames=fieldnames(kml_placemark());

                        for i=length(fieldnames(kml_feature()))+1:min(nargin,length(fnames))
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
                disp@kml_feature(kml(i));
                disp(sprintf('      geometry: %s %s\n' ,string_size(kml(i).geometry),...
                             class(kml(i).geometry)));
            end

        end

%  return the fieldnames of the object

        function [fnames]=fieldnames(kml)

%  fieldnames for a sub (derived) class list those before super (base)

            fnames=fieldnames(kml_feature());
            fnames={fnames{:} ...
                    'geometry' ...
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

%  loop over the placemarks

            for i=1:numel(kml)
                kmli=kml(i);
                if ~isempty(kmli.id)
                    fprintf(fid,'%s<Placemark id="%s">\n',indent,kmli.id);
                else
                    fprintf(fid,'%s<Placemark>\n',indent);
                end
                kml_write@kml_feature(kmli,fid,indent);

%  loop over the geometry elements for each placemark

                for j=1:min(1,numel(kmli.geometry))
                    kmlij=kmli.geometry(j);
                    if ~isempty(kmlij)
                        if isa(kmlij,'kml_geometry')
                            kml_write(kmlij,fid,[indent '  ']);
                        else
                            warning('kml(%d).geometry(%d) is a ''%s'' class object, not ''%s''.',...
                                i,j,class(kmlij),'kml_geometry');
                        end
                    end
                end

                fprintf(fid,'%s</Placemark>\n',indent);
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

%  loop over the placemarks

            for i=1:numel(kml)
                kmli=kml(i);
                if ~isempty(kmli.id)
                    sbuf=add(sbuf,sprintf('%s<Placemark id="%s">\n',indent,kmli.id));
                else
                    sbuf=add(sbuf,sprintf('%s<Placemark>\n',indent));
                end
                sbuf=kml_swrite@kml_feature(kmli,sbuf,indent);

%  loop over the geometry elements for each placemark

                for j=1:min(1,numel(kmli.geometry))
                    kmlij=kmli.geometry(j);
                    if ~isempty(kmlij)
                        if isa(kmlij,'kml_geometry')
                            sbuf=kml_swrite(kmlij,sbuf,[indent '  ']);
                        else
                            warning('kml(%d).geometry(%d) is a ''%s'' class object, not ''%s''.',...
                                i,j,class(kmlij),'kml_geometry');
                        end
                    end
                end

                sbuf=add(sbuf,sprintf('%s</Placemark>\n',indent));
            end

        end

%  delete the object

        function []=delete(kml)

%  loop over the placemarks

            for i=numel(kml):-1:1
                kmli=kml(i);
                delete@kml_feature(kmli);

%  loop over the geometry elements for each placemark

                for j=min(1,numel(kmli.geometry)):-1:1
                    kmlij=kmli.geometry(j);
                    if ~isempty(kmlij)
                        if isa(kmlij,'kml_geometry')
                            delete(kmlij);
                        else
                            warning('kml(%d).geometry(%d) is a ''%s'' class object, not ''%s''.',...
                                i,j,class(kmlij),'kml_geometry');
                        end
                    end
                end
                kmli.geometry  =kml_geometry.empty();

            end

        end

    end

end
