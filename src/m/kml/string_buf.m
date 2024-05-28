%
%  definition for the string_buf class.
%
%  [sbuf]=string_buf(varargin)
%
%  where the optional varargin and defaults are:
%    init          (numeric, initial size)
%    inc           (numeric, incremental size)
%    max           (numeric, maximum size)
%
%  and the protected properties are:
%    string        (char, string buffer)
%    size          (numeric, current size of buffer)
%    len           (numeric, current length of string in buffer)
%
%  note that zero arguments constructs a default instance; one
%  argument of the class copies the instance; and two or more
%  arguments constructs a new instance from the arguments.
%
classdef string_buf < handle
    properties
        init      =10000000;
        inc       =1000000;
        max       =100000000;
    end
%     properties (SetAccess = private, GetAccess = private)
    properties (SetAccess = private)
        string    ='';
        size      =0;
        len       =0;
    end

    methods
        function [sbuf]=string_buf(varargin)

            switch nargin

%  create a default object

                case 0

%  copy the object or create the object from the input

                otherwise
                    if (nargin == 1) && isa(varargin{1},class(string_buf))
                        sbuf=varargin{1};

                    else
                        fnames=fieldnames(string_buf());

                        for i=1:min(nargin,length(fnames))
                            if isa(varargin{i},class(sbuf.(fnames{i})))
                                if ~isempty(varargin{i})
                                    sbuf.(fnames{i})=varargin{i};
                                end
                            else
                                if ~isempty(inputname(i))
                                    warning('Argument ''%s'' for property ''%s'' is a ''%s'' class object, not ''%s''.',...
                                        inputname(i),fnames{i},class(varargin{i}),class(sbuf.(fnames{i})));
                                else
                                    warning('Argument %d for property ''%s'' is a ''%s'' class object, not ''%s''.',...
                                        i           ,fnames{i},class(varargin{i}),class(sbuf.(fnames{i})));
                                end
                            end
                        end
                    end

            end

            sbuf.string=blanks(sbuf.init);
            sbuf.size  =sbuf.init;

        end

%  display the object

        function []=disp(sbuf)

            for i=1:numel(sbuf)
                disp(sprintf('class ''%s'' object ''%s%s'' = \n',...
                    class(sbuf),inputname(1),string_dim(sbuf,i)));
                disp(sprintf('          init: %d'      ,sbuf(i).init));
                disp(sprintf('           inc: %d'      ,sbuf(i).inc));
                disp(sprintf('           max: %d'      ,sbuf(i).max));
                disp(sprintf('        string: %s'      ,any2str(sbuf(i).string,40)));
                disp(sprintf('          size: %d'      ,sbuf(i).size));
                disp(sprintf('           len: %d\n'    ,sbuf(i).len));
            end

        end

%  return the fieldnames of the object

        function [fnames]=fieldnames(sbuf)

            fnames={'init' ...
                    'inc' ...
                    'max' ...
                   }';

        end

%  set the properties of the object

        function [sbuf]=set(sbuf,varargin)

            sbufref=feval(class(sbuf));
            fnames=fieldnames(sbufref);

%  loop through each parameter in the input list (comparing to the reference
%  object in case property types have been changed)

            for i=1:2:length(varargin)
                if ismember(varargin{i},fnames) && (i+1 <= length(varargin))
                    if isa(varargin{i+1},class(sbufref.(varargin{i})))
                        sbuf.(varargin{i})=varargin{i+1};
                    else
                        if ~isempty(inputname(i+1))
                            warning('Argument ''%s'' for property ''%s'' is a ''%s'' class object, not ''%s''.',...
                                inputname(i+2),varargin{i},class(varargin{i+1}),class(sbufref.(varargin{i})));
                        else
                            warning('Argument %d for property ''%s'' is a ''%s'' class object, not ''%s''.',...
                                i+2           ,varargin{i},class(varargin{i+1}),class(sbufref.(varargin{i})));
                        end
                    end
                else
                    warning('Property ''%s'' for class ''%s'' does not exist.',...
                        varargin{i},class(sbufref));
                end
            end

        end

%  add a string to the object

        function [sbuf]=add(sbuf,str)

            if ~ischar(str)
                if ~isempty(inputname(2))
                    warning('Argument ''%s'' for string is a ''%s'' class object, not ''%s''.',...
                        inputname(2),class(str),'char');
                else
                    warning('Argument %d for string is a ''%s'' class object, not ''%s''.',...
                        2           ,class(str),'char');
                end
            end

%  check the buffer size and increase as necessary

            slen=length(str);
            while (sbuf.len+slen > sbuf.size)
                if (sbuf.size+sbuf.inc <= sbuf.max)
                    sbuf.string=[sbuf.string blanks(sbuf.inc)];
                    sbuf.size  =sbuf.size+sbuf.inc;
                else
                    error('String buffer length of %d would exceed maximum of %d.',...
                        sbuf.size+sbuf.inc,sbuf.max);
                end
            end

%  copy the string into the buffer

            sbuf.string(sbuf.len+1:sbuf.len+slen)=str;
            sbuf.len=sbuf.len+slen;

        end

%  return the string from the object

        function [str]=str(sbuf)

           str=sbuf.string(1:sbuf.len);

        end

%  reset the object

        function [sbuf]=reset(sbuf)

            string    ='';
            size      =0;
            len       =0;

        end

    end

end
