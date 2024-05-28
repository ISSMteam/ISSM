%
%  definition for the kml_polystyle sub (derived) class.
%
%  [kml]=kml_polystyle(varargin)
%
%  where the optional varargin and defaults are:
%    id            (char, polystyle id, '')
%    color         (char, opacity/color in hex aabbggrr, 'ffffffff')
%    colormode     (char, color mode , 'normal')
%    fill          (logical, polygon fill, true)
%    outline       (logical, polygon outline, true)
%
%  note that zero arguments constructs a default instance; one
%  argument of the class copies the instance; and two or more
%  arguments constructs a new instance from the arguments.
%
classdef kml_polystyle < kml_colorstyle
    properties
        fill      =true;
        outline   =true;
    end

    methods
        function [kml]=kml_polystyle(varargin)

            kml=kml@kml_colorstyle(varargin{:});

            switch nargin

%  create a default object

                case 0

%  copy the object or create the object from the input

                otherwise
                    if (nargin == 1) && isa(varargin{1},class(kml))
                        kml=varargin{1};

                    else
                        fnames=fieldnames(kml_polystyle());

                        for i=length(fieldnames(kml_colorstyle()))+1:min(nargin,length(fnames))
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
                disp@kml_colorstyle(kml(i));
                disp(sprintf('          fill: %d'      ,kml(i).fill));
                disp(sprintf('       outline: %d\n'    ,kml(i).outline));
            end

        end

%  return the fieldnames of the object

        function [fnames]=fieldnames(kml)

%  fieldnames for a sub (derived) class list those before super (base)

            fnames=fieldnames(kml_colorstyle());
            fnames={fnames{:} ...
                    'fill' ...
                    'outline' ...
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

%  loop over the polystyles

            for i=1:numel(kml)
                kmli=kml(i);
                if ~isempty(kmli.id)
                    fprintf(fid,'%s<PolyStyle id="%s">\n',indent,kmli.id);
                else
                    fprintf(fid,'%s<PolyStyle>\n',indent);
                end
                kml_write@kml_colorstyle(kmli,fid,indent);
                fprintf(fid,'%s  <fill>%d</fill>\n',indent,kmli.fill);
                fprintf(fid,'%s  <outline>%d</outline>\n',indent,kmli.outline);
                fprintf(fid,'%s</PolyStyle>\n',indent);
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

%  loop over the polystyles

            for i=1:numel(kml)
                kmli=kml(i);
                if ~isempty(kmli.id)
                    sbuf=add(sbuf,sprintf('%s<PolyStyle id="%s">\n',indent,kmli.id));
                else
                    sbuf=add(sbuf,sprintf('%s<PolyStyle>\n',indent));
                end
                sbuf=kml_swrite@kml_colorstyle(kmli,sbuf,indent);
                sbuf=add(sbuf,sprintf('%s  <fill>%d</fill>\n',indent,kmli.fill));
                sbuf=add(sbuf,sprintf('%s  <outline>%d</outline>\n',indent,kmli.outline));
                sbuf=add(sbuf,sprintf('%s</PolyStyle>\n',indent));
            end

        end

    end

end
