%
%  definition for the kml_feature super (base) and sub (derived) abstract class.
%
%  [kml]=kml_feature(varargin)
%
%  where the optional varargin and defaults are:
%    id            (char, feature id, '')
%    name          (char, name, '')
%    visibility    (logical, visibility, true)
%    open          (logical, open, false)
%    snippet       (char, snippet, '')
%    descript      (char, description, '')
%    styleurl      (char, style url, '')
%    style         (cell array, styles)
%
%  note that zero arguments constructs a default instance; one
%  argument of the class copies the instance; and two or more
%  arguments constructs a new instance from the arguments.
%
classdef kml_feature < kml_object
    properties
        name      ='';
        visibility=true;
        open      =false;
        snippet   ='';
        descript  ='';
        styleurl  ='';
        style     ={};
    end

    methods
        function [kml]=kml_feature(varargin)

            kml=kml@kml_object(varargin{:});

            switch nargin

%  create a default object

                case 0

%  copy the object or create the object from the input

                otherwise
                    if (nargin == 1) && isa(varargin{1},class(kml))
                        kml=varargin{1};

                    else
                        fnames=fieldnames(kml_feature());

                        for i=length(fieldnames(kml_object()))+1:min(nargin,length(fnames))
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
                if strcmp(class(kml),'kml_feature')
                    disp(sprintf('class ''%s'' object ''%s%s'' = \n',...
                        class(kml),inputname(1),string_dim(kml,i)));
                end
                disp@kml_object(kml(i));
                disp(sprintf('          name: ''%s'''  ,kml(i).name));
                disp(sprintf('    visibility: %g'      ,kml(i).visibility));
                disp(sprintf('          open: %g'      ,kml(i).open));
                disp(sprintf('       snippet: ''%s'''  ,kml(i).snippet));
                disp(sprintf('      descript: ''%s'''  ,kml(i).descript));
                disp(sprintf('      styleurl: ''%s'''  ,kml(i).styleurl));
                if strcmp(class(kml),'kml_feature')
                    disp(sprintf('         style: %s %s\n' ,string_size(kml(i).style),...
                                 class(kml(i).style)));
                else
                    disp(sprintf('         style: %s %s'   ,string_size(kml(i).style),...
                                 class(kml(i).style)));
                end
            end

        end

%  return the fieldnames of the object

        function [fnames]=fieldnames(kml)

%  fieldnames for a sub (derived) class list those before super (base)

            fnames=fieldnames(kml_object());
            fnames={fnames{:} ...
                    'name' ...
                    'visibility' ...
                    'open' ...
                    'snippet' ...
                    'descript' ...
                    'styleurl' ...
                    'style' ...
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

%  loop over the features

            for i=1:numel(kml)
                kmli=kml(i);
                if strcmp(class(kml),'kml_feature')
                    if ~isempty(kmli.id)
                        fprintf(fid,'%s<!Feature id="%s">\n',indent,kmli.id);
                    else
                        fprintf(fid,'%s<!Feature>\n',indent);
                    end
                end
                kml_write@kml_object(kmli,fid,indent);
                if ~isempty(kmli.name)
                    fprintf(fid,'%s  <name>%s</name>\n',indent,kmli.name);
                end
                fprintf(fid,'%s  <visibility>%d</visibility>\n',indent,kmli.visibility);
                fprintf(fid,'%s  <open>%d</open>\n',indent,kmli.open);
                if ~isempty(kmli.snippet)
                    fprintf(fid,'%s  <Snippet maxLines="2">%s</Snippet>\n',indent,kmli.snippet);
                end
                if ~isempty(kmli.descript)
                    fprintf(fid,'%s  <description>%s</description>\n',indent,kmli.descript);
                end
                if ~isempty(kmli.styleurl)
                    fprintf(fid,'%s  <styleUrl>%s</styleUrl>\n',indent,kmli.styleurl);
                end

%  loop over the styles for each feature

                for j=1:numel(kmli.style)
                    if ~isempty(kmli.style{j})
                        if isa(kmli.style{j},'kml_styleselector')
                            kml_write(kmli.style{j},fid,[indent '  ']);
                        else
                            warning('kml(%d).style{%d} is a ''%s'' class object, not ''%s''.',...
                                i,j,class(kmli.style{j}),'kml_styleselector');
                        end
                    end
                end

                if strcmp(class(kml),'kml_feature')
                    fprintf(fid,'%s<!/Feature>\n',indent);
                end
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

%  loop over the features

            for i=1:numel(kml)
                kmli=kml(i);
                if strcmp(class(kml),'kml_feature')
                    if ~isempty(kmli.id)
                        sbuf=add(sbuf,sprintf('%s<!Feature id="%s">\n',indent,kmli.id));
                    else
                        sbuf=add(sbuf,sprintf('%s<!Feature>\n',indent));
                    end
                end
                sbuf=kml_swrite@kml_object(kmli,sbuf,indent);
                if ~isempty(kmli.name)
                    sbuf=add(sbuf,sprintf('%s  <name>%s</name>\n',indent,kmli.name));
                end
                sbuf=add(sbuf,sprintf('%s  <visibility>%d</visibility>\n',indent,kmli.visibility));
                sbuf=add(sbuf,sprintf('%s  <open>%d</open>\n',indent,kmli.open));
                if ~isempty(kmli.snippet)
                    sbuf=add(sbuf,sprintf('%s  <Snippet maxLines="2">%s</Snippet>\n',indent,kmli.snippet));
                end
                if ~isempty(kmli.descript)
                    sbuf=add(sbuf,sprintf('%s  <description>%s</description>\n',indent,kmli.descript));
                end
                if ~isempty(kmli.styleurl)
                    sbuf=add(sbuf,sprintf('%s  <styleUrl>%s</styleUrl>\n',indent,kmli.styleurl));
                end

%  loop over the styles for each feature

                for j=1:numel(kmli.style)
                    if ~isempty(kmli.style{j})
                        if isa(kmli.style{j},'kml_styleselector')
                            sbuf=kml_swrite(kmli.style{j},sbuf,[indent '  ']);
                        else
                            warning('kml(%d).style{%d} is a ''%s'' class object, not ''%s''.',...
                                i,j,class(kmli.style{j}),'kml_styleselector');
                        end
                    end
                end

                if strcmp(class(kml),'kml_feature')
                    sbuf=add(sbuf,sprintf('%s<!/Feature>\n',indent));
                end
            end

        end

%  delete the object

        function []=delete(kml)

%  loop over the features

            for i=numel(kml):-1:1
                kmli=kml(i);

%  loop over the styles for each feature

                for j=numel(kmli.style):-1:1
                    if ~isempty(kmli.style{j})
                        if isa(kmli.style{j},'kml_styleselector')
                            delete(kmli.style{j});
                        else
                            warning('kml(%d).style{%d} is a ''%s'' class object, not ''%s''.',...
                                i,j,class(kmli.style{j}),'kml_styleselector');
                        end
                    end
                end
                kmli.style     ={};

            end

        end

    end

end
