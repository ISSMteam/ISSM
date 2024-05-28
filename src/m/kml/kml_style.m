%
%  definition for the kml_style sub (derived) class.
%
%  [kml]=kml_style(varargin)
%
%  where the optional varargin and defaults are:
%    id            (char, style id, '')
%    icon          (char, icon style, '')
%    label         (char, label style, '')
%    line          (char, line style, '')
%    poly          (char, poly style, '')
%    balloon       (char, balloon style, '')
%    list          (char, list style, '')
%
%  note that zero arguments constructs a default instance; one
%  argument of the class copies the instance; and two or more
%  arguments constructs a new instance from the arguments.
%
classdef kml_style < kml_styleselector
    properties
%         icon      =kml_iconstyle.empty();
%         label     =kml_labelstyle.empty();
        icon      =[];
        label     =[];
        line      =kml_linestyle.empty();
        poly      =kml_polystyle.empty();
%         balloon   =kml_balloonstyle.empty();
%         list      =kml_liststyle.empty();
        balloon   =[];
        list      =[];
    end

    methods
        function [kml]=kml_style(varargin)

            kml=kml@kml_styleselector(varargin{:});

            switch nargin

%  create a default object

                case 0

%  copy the object or create the object from the input

                otherwise
                    if (nargin == 1) && isa(varargin{1},class(kml))
                        kml=varargin{1};

                    else
                        fnames=fieldnames(kml_style());

                        for i=length(fieldnames(kml_styleselector()))+1:min(nargin,length(fnames))
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
                disp@kml_styleselector(kml(i));
                disp(sprintf('          icon: %s %s'   ,string_size(kml(i).icon),...
                             class(kml(i).icon)));
                disp(sprintf('         label: %s %s'   ,string_size(kml(i).label),...
                             class(kml(i).label)));
                disp(sprintf('          line: %s %s'   ,string_size(kml(i).line),...
                             class(kml(i).line)));
                disp(sprintf('          poly: %s %s'   ,string_size(kml(i).poly),...
                             class(kml(i).poly)));
                disp(sprintf('       balloon: %s %s'   ,string_size(kml(i).balloon),...
                             class(kml(i).balloon)));
                disp(sprintf('          list: %s %s\n' ,string_size(kml(i).list),...
                             class(kml(i).list)));
            end

        end

%  return the fieldnames of the object

        function [fnames]=fieldnames(kml)

%  fieldnames for a sub (derived) class list those before super (base)

            fnames=fieldnames(kml_styleselector());
            fnames={fnames{:} ...
                    'icon' ...
                    'label' ...
                    'line' ...
                    'poly' ...
                    'balloon' ...
                    'list' ...
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

%  loop over the styles

            for i=1:numel(kml)
                kmli=kml(i);
                if ~isempty(kmli.id)
                    fprintf(fid,'%s<Style id="%s">\n',indent,kmli.id);
                else
                    fprintf(fid,'%s<Style>\n',indent);
                end
                kml_write@kml_styleselector(kmli,fid,indent);
%                 if isa(kmli.icon,'kml_iconstyle')
%                     kml_write(kmli.icon,fid,[indent '  ']);
%                 else
%                     warning('kml(%d).icon is a ''%s'' class object, not ''%s''.',...
%                         i,class(kmli.icon),'kml_iconstyle');
%                 end
%                 if isa(kmli.label,'kml_labelstyle')
%                     kml_write(kmli.label,fid,[indent '  ']);
%                 else
%                     warning('kml(%d).label is a ''%s'' class object, not ''%s''.',...
%                         i,class(kmli.label),'kml_labelstyle');
%                 end
                if isa(kmli.line,'kml_linestyle')
                    kml_write(kmli.line,fid,[indent '  ']);
                else
                    warning('kml(%d).line is a ''%s'' class object, not ''%s''.',...
                        i,class(kmli.line),'kml_linestyle');
                end
                if isa(kmli.poly,'kml_polystyle')
                    kml_write(kmli.poly,fid,[indent '  ']);
                else
                    warning('kml(%d).poly is a ''%s'' class object, not ''%s''.',...
                        i,class(kmli.poly),'kml_polystyle');
                end
%                 if isa(kmli.balloon,'kml_balloonstyle')
%                     kml_write(kmli.balloon,fid,[indent '  ']);
%                 else
%                     warning('kml(%d).balloon is a ''%s'' class object, not ''%s''.',...
%                         i,class(kmli.balloon),'kml_balloonstyle');
%                 end
%                 if isa(kmli.list,'kml_liststyle')
%                     kml_write(kmli.list,fid,[indent '  ']);
%                 else
%                     warning('kml(%d).list is a ''%s'' class object, not ''%s''.',...
%                         i,class(kmli.list),'kml_liststyle');
%                 end
                fprintf(fid,'%s</Style>\n',indent);
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

%  loop over the styles

            for i=1:numel(kml)
                kmli=kml(i);
                if ~isempty(kmli.id)
                    sbuf=add(sbuf,sprintf('%s<Style id="%s">\n',indent,kmli.id));
                else
                    sbuf=add(sbuf,sprintf('%s<Style>\n',indent));
                end
                sbuf=kml_swrite@kml_styleselector(kmli,sbuf,indent);
%                 if isa(kmli.icon,'kml_iconstyle')
%                     sbuf=kml_swrite(kmli.icon,sbuf,[indent '  ']);
%                 else
%                     warning('kml(%d).icon is a ''%s'' class object, not ''%s''.',...
%                         i,class(kmli.icon),'kml_iconstyle');
%                 end
%                 if isa(kmli.label,'kml_labelstyle')
%                     sbuf=kml_swrite(kmli.label,sbuf,[indent '  ']);
%                 else
%                     warning('kml(%d).label is a ''%s'' class object, not ''%s''.',...
%                         i,class(kmli.label),'kml_labelstyle');
%                 end
                if isa(kmli.line,'kml_linestyle')
                    sbuf=kml_swrite(kmli.line,sbuf,[indent '  ']);
                else
                    warning('kml(%d).line is a ''%s'' class object, not ''%s''.',...
                        i,class(kmli.line),'kml_linestyle');
                end
                if isa(kmli.poly,'kml_polystyle')
                    sbuf=kml_swrite(kmli.poly,sbuf,[indent '  ']);
                else
                    warning('kml(%d).poly is a ''%s'' class object, not ''%s''.',...
                        i,class(kmli.poly),'kml_polystyle');
                end
%                 if isa(kmli.balloon,'kml_balloonstyle')
%                     sbuf=kml_swrite(kmli.balloon,sbuf,[indent '  ']);
%                 else
%                     warning('kml(%d).balloon is a ''%s'' class object, not ''%s''.',...
%                         i,class(kmli.balloon),'kml_balloonstyle');
%                 end
%                 if isa(kmli.list,'kml_liststyle')
%                     sbuf=kml_swrite(kmli.list,sbuf,[indent '  ']);
%                 else
%                     warning('kml(%d).list is a ''%s'' class object, not ''%s''.',...
%                         i,class(kmli.list),'kml_liststyle');
%                 end
                sbuf=add(sbuf,sprintf('%s</Style>\n',indent));
            end

        end

%  delete the object

        function []=delete(kml)

%  loop over the styles

            for i=numel(kml):-1:1
                kmli=kml(i);
%                 if isa(kmli.icon,'kml_iconstyle')
%                     delete(kmli.icon);
%                 else
%                     warning('kml(%d).icon is a ''%s'' class object, not ''%s''.',...
%                         i,class(kmli.icon),'kml_iconstyle');
%                 end
%                 kmli.icon      =kml_iconstyle.empty();
%                 if isa(kmli.label,'kml_labelstyle')
%                     delete(kmli.label);
%                 else
%                     warning('kml(%d).label is a ''%s'' class object, not ''%s''.',...
%                         i,class(kmli.label),'kml_labelstyle');
%                 end
%                 kmli.label     =kml_labelstyle.empty();
                if isa(kmli.line,'kml_linestyle')
                    delete(kmli.line);
                else
                    warning('kml(%d).line is a ''%s'' class object, not ''%s''.',...
                        i,class(kmli.line),'kml_linestyle');
                end
                kmli.line      =kml_linestyle.empty();
                if isa(kmli.poly,'kml_polystyle')
                    delete(kmli.poly);
                else
                    warning('kml(%d).poly is a ''%s'' class object, not ''%s''.',...
                        i,class(kmli.poly),'kml_polystyle');
                end
                kmli.poly      =kml_polystyle.empty();
%                 if isa(kmli.balloon,'kml_balloonstyle')
%                     delete(kmli.balloon);
%                 else
%                     warning('kml(%d).balloon is a ''%s'' class object, not ''%s''.',...
%                         i,class(kmli.balloon),'kml_balloonstyle');
%                 end
%                 kmli.balloon   =kml_balloonstyle.empty();
%                 if isa(kmli.list,'kml_liststyle')
%                     delete(kmli.list);
%                 else
%                     warning('kml(%d).list is a ''%s'' class object, not ''%s''.',...
%                         i,class(kmli.list),'kml_liststyle');
%                 end
%                 kmli.list      =kml_liststyle.empty();
            end

        end

    end

end
