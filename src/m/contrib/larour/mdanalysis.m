function varargout = mdanalysis(varargin)
% MDANALYSIS MATLAB code for mdanalysis.fig {{{
%      MDANALYSIS, by itself, creates a new MDANALYSIS or raises the existing
%      singleton*.
%
%      H = MDANALYSIS returns the handle to a new MDANALYSIS or the handle to
%      the existing singleton*.
%
%      MDANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MDANALYSIS.M with the given input arguments.
%
%      MDANALYSIS('Property','Value',...) creates a new MDANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mdanalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mdanalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%variables: }}}
%global variables:  %{{{
global md  modelname
if strcmpi(class(varargin{1}),'model') | strcmpi(class(varargin{1}),'sealevelmodel'), 
	md=varargin{1};
	modelname=inputname(1);
end

global logvalue
global solutiontype
global comparison 
global reload
global earth
global whores
%}}}
%Initialization code{{{
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mdanalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @mdanalysis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mdanalysis is made visible.
function mdanalysis_OpeningFcn(hObject, eventdata, handles, varargin)
global md comparison logvalue solutiontype
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mdanalysis (see VARARGIN)

% Choose default command line output for mdanalysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using mdanalysis.
if strcmp(get(hObject,'Visible'),'off')
	comparison=0;
	logvalue=0;
	reload=0;
	plotm();
end

% UIWAIT makes mdanalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mdanalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%}}}
%Interp faceting {{{
% --- Executes on button press in Interp.
function Interp_Callback(hObject, eventdata, handles)
% hObject    handle to Interp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotm();
% Hint: get(hObject,'Value') returns toggle state of Interp

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Interp.
function Interp_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Interp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% }}}
%Faceted facting {{{
% --- Executes on button press in Faceted.
function Faceted_Callback(hObject, eventdata, handles)
% hObject    handle to Faceted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotm();
% Hint: get(hObject,'Value') returns toggle state of Faceted

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Faceted.
function Faceted_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Faceted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%}}}
%Flat faceting {{{
% --- Executes on button press in Flat.
function Flat_Callback(hObject, eventdata, handles)
% hObject    handle to Flat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotm();
% Hint: get(hObject,'Value') returns toggle state of Flat

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Flat.
function Flat_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Flat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%}}}
%Step {{{
% --- Executes on selection change in Step.
function Step_Callback(hObject, eventdata, handles)
global grounded ice  logvalue
% hObject    handle to Step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotm();

% Hints: contents = cellstr(get(hObject,'String')) returns Step contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Step

% --- Executes during object creation, after setting all properties.
function Step_CreateFcn(hObject, eventdata, handles)
global md grounded ice  logvalue 
% hObject    handle to Step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if issealevel(), if isearth, results=md.earth.results; else results=md.icecaps{1}.results; end; else results=md.results; end;
if isfield(results,'TransientSolution'),
	strings=cell(length(results.TransientSolution),1);
	for i=1:length(results.TransientSolution),
		strings{i}=sprintf('%4.2f     (%i)',results.TransientSolution(i).time,i);
	end
	if length(strings)>40,
		%too many strings, reduce!
		lt=length(strings);
		modlt=floor(lt/40);
		strings=strings(1:modlt:end);
	end
	set(hObject,'String',strings);
	set(hObject,'Value',1);
else
	set(hObject,'String',{});
	set(hObject,'Value',1);
end
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%}}}
%Log {{{
% --- Executes on button press in Log.
function Log_Callback(hObject, eventdata, handles)
global grounded ice  logvalue
% hObject    handle to Log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotm();
% Hint: get(hObject,'Value') returns toggle state of Log

% --- Executes on button press in Lock.
function Lock_Callback(hObject, eventdata, handles)
% hObject    handle to Lock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotm();
% Hint: get(hObject,'Value') returns toggle state of Lock
% }}}
%Collock: not sure whether that one is working anymore: {{{
% --- Executes on button press in Collock.
function Collock_Callback(hObject, eventdata, handles)
% hObject    handle to Collock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotm();
% Hint: get(hObject,'Value') returns toggle state of Collock
% }}}
%Field {{{
% --- Executes on selection change in Field.
function Field_Callback(hObject, eventdata, handles)
global md  
% hObject    handle to Field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

strings=get(hObject,'String'); value=get(hObject,'Value');
string=strings(value); 

if strcmpi(string,'IceVolume') | strcmpi(string,'IceVolumeAboveFloatation'),
	plotv();
elseif strcmpi(string,'SealevelRSLEustatic'),
	%display directly: 
	displayscalar();
else
	plotm();
end



% Hints: contents = cellstr(get(hObject,'String')) returns Field contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Field

% --- Executes during object creation, after setting all properties.
function Field_CreateFcn(hObject, eventdata, handles)
global md  
% hObject    handle to Field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if issealevel(), if isearth, results=md.earth.results; else results=md.icecaps{1}.results; end; else results=md.results; end;
if isfield(results,'TransientSolution'),
	fields=listfields(results,'TransientSolution');
	set(hObject,'String',fields);
	for i=1:length(fields), 
		if ~issealevel & strcmpi(fields{i},'Vel'),
			set(hObject,'Value',i);
			break;
		end
		if issealevel & isearth & strcmpi(fields{i},'Sealevel'),
			set(hObject,'Value',i);
			break;
		end
	end
else
	set(hObject,'String',{});
	set(hObject,'Value',0);
end

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%}}}
%Model Loading {{{
% --- Executes on selection change in ModelChoice.
function ModelChoice_Callback(hObject, eventdata, handles)
	% hObject    handle to ModelChoice (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)

	% Hints: contents = cellstr(get(hObject,'String')) returns ModelChoice contents as cell array
	%        contents{get(hObject,'Value')} returns selected item from ModelChoice
	global modelname  md whores

	%figure out which model is being loaded: 
	strings=cellstr(get(hObject,'String'));
	counter=get(hObject,'Value');
	string=strings{counter};
	counter=findstr(string,'(');
	newmodelname=string(1:counter-2); 

	%load model from base workspace:
	md= evalin('base', newmodelname);

	%replot: 
	h = findobj('Tag','Field');
	Field_Callback(h,eventdata,handles);


	%eval([newmodelname '=md;']);
	%close gcbf;
	%eval(['mdanalysis(' newmodelname ');']); 

% --- Executes during object creation, after setting all properties.
function ModelChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ModelChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global modelname


%figure out how many models are in the workspace!
whores=evalin('base','whos');
flags=zeros(length(whores),1);
for i=1:length(whores),
	if strcmpi(whores(i).class,'model') | strcmpi(whores(i).class,'sealevelmodel'),
		flags(i)=1;
	end
end
pos=find(flags);
whores=whores(pos);

strings={};
for i=1:length(whores),
	strings{end+1}=[whores(i).name ' (' whores(i).class ')'];
end

%match with varargin{1}'s name: 
counter=-1;
for i=1:length(strings),
	if strcmpi(modelname,whores(i).name),
		counter=i;
		break;
	end
end
if counter==-1, 
	error('could not find input model name matching base workspace names!');
end

set(hObject,'String',strings);
set(hObject,'Value',counter);


% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end 

%}}}
%Cmin {{{
function Cmin_Callback(hObject, eventdata, handles)
% hObject    handle to Cmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotm();
% Hints: get(hObject,'String') returns contents of Cmin as text
%        str2double(get(hObject,'String')) returns contents of Cmin as a double


% --- Executes during object creation, after setting all properties.
function Cmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String',num2str(NaN));
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%}}}
%Cmax {{{
function Cmax_Callback(hObject, eventdata, handles)
% hObject    handle to Cmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cmax as text
%        str2double(get(hObject,'String')) returns contents of Cmax as a double
plotm();

% --- Executes during object creation, after setting all properties.
function Cmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String',num2str(NaN));
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%}}}
%Xmin {{{
function xmin_Callback(hObject, eventdata, handles)
% hObject    handle to xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xmin as text
%        str2double(get(hObject,'String')) returns contents of xmin as a double

% --- Executes during object creation, after setting all properties.
function xmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%}}}
%Xmax {{{
function xmax_Callback(hObject, eventdata, handles)
% hObject    handle to xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of xmax as text
%        str2double(get(hObject,'String')) returns contents of xmax as a double

% --- Executes during object creation, after setting all properties.
function xmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% }}}
%Ymin {{{
function ymin_Callback(hObject, eventdata, handles)
% hObject    handle to ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of ymin as text
%        str2double(get(hObject,'String')) returns contents of ymin as a double


% --- Executes during object creation, after setting all properties.
function ymin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%}}}
%Ymax {{{
function ymax_Callback(hObject, eventdata, handles)
% hObject    handle to ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of ymax as text
%        str2double(get(hObject,'String')) returns contents of ymax as a double


% --- Executes during object creation, after setting all properties.
function ymax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%}}}
%Time lable {{{
function Time_Callback(hObject, eventdata, handles)
% hObject    handle to Time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Time as text
%        str2double(get(hObject,'String')) returns contents of Time as a double


% --- Executes during object creation, after setting all properties.
function Time_CreateFcn(hObject, eventdata, handles)
global md 
% hObject    handle to Time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if issealevel, if isearth, results=md.earth.results; else results=md.icecaps{1}.results; end; else results=md.results; end;
if isfield(results,'TransientSolution'),
	set(hObject,'String',sprintf('%4.2f',results.TransientSolution(1).time));
else
	set(hObject,'String',sprintf('%4.2f',0));
end
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Time.
function Time_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% }}}
%Mask: {{{
% --- Executes on selection change in Masque.
function Masque_Callback(hObject, eventdata, handles)
% hObject    handle to Masque (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotm();
% Hints: contents = cellstr(get(hObject,'String')) returns Masque contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Masque

% --- Executes during object creation, after setting all properties.
function Masque_CreateFcn(hObject, eventdata, handles)
global md 
% hObject    handle to Masque (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if issealevel, 
	if isearth, 
		mv=md.earth.mask; 
	else 
		mv=md.icecaps{1}.mask; 
	end 
else 
	mv=md.mask; 
end
strings=fieldnames(mv);
strings={'all',strings{:}};
set(hObject,'String',strings);
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end 
%}}}
%DiffStep {{{
% --- Executes on selection change in DiffStep.
function DiffStep_Callback(hObject, eventdata, handles)
% hObject    handle to DiffStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotm();

% Hints: contents = cellstr(get(hObject,'String')) returns DiffStep contents as cell array
%        contents{get(hObject,'Value')} returns selected item from DiffStep

% --- Executes during object creation, after setting all properties.
function DiffStep_CreateFcn(hObject, eventdata, handles)
global md ;
% hObject    handle to DiffStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if issealevel, if isearth, results=md.earth.results; else results=md.icecaps{1}.results; end; else results=md.results; end;
if isfield(results,'TransientSolution'),
	strings=cell(length(results.TransientSolution),1);
	for i=1:length(results.TransientSolution),
		strings{i}=sprintf('%4.2f     (%i)',results.TransientSolution(i).time,i);
	end
	if length(strings)>20,
		%too many strings, reduce!
		lt=length(strings);
		modlt=floor(lt/20);
		strings=strings(1:modlt:end);
	end
	set(hObject,'String',strings);
else
	set(hObject,'String',{});
end

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over DiffStep.
function DiffStep_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to DiffStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% }}}
%Diff button {{{
% --- Executes on button press in Diff.
function Diff_Callback(hObject, eventdata, handles)
% hObject    handle to Diff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotm();
% Hint: get(hObject,'Value') returns toggle state of Diff

% --- Executes during object creation, after setting all properties.
function Diff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Diff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%}}}
%Model Fields:  {{{
% --- Executes on selection change in ModelFields.
function ModelFields_Callback(hObject, eventdata, handles)
% hObject    handle to ModelFields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotm();
% Hints: contents = cellstr(get(hObject,'String')) returns ModelFields contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ModelFields

% --- Executes during object creation, after setting all properties.
function ModelFields_CreateFcn(hObject, eventdata, handles)
global md 
% hObject    handle to ModelFields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% fields {{{ 
fields={'geometry.thickness',...
		 'geometry.surface',...
		 'geometry.base',...
		 'geometry.bed',...
		 'inversion.vx_obs',...
		 'inversion.vy_obs',...
		 'inversion.vel_obs',...
		 'friction.coefficient',...
		 'materials.rheology_B',...
		 'mask.ice_levelset',...
		 'mask.groundedice_levelset',...
		 'slr.deltathickness'...
		 };

 if issealevel(),
	 if isearth,
		 mask=md.earth.mask;
	 else
		 mask=md.icecaps{1}.mask;
	 end
 else
	 mask=md.mask;
 end

 if  strcmpi(class(mask),'maskpsl'),
	 fields{end+1}='mask.ocean_levelset';
	 fields{end+1}='mask.land_levelset';
 end


% }}}
set(hObject,'String',fields);
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% }}}
%Mf button: {{{
% --- Executes on button press in Mf.
function Mf_CreateFcn(varargin)
function Mf_Callback(hObject, eventdata, handles)
% hObject    handle to Mf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotm();
% Hint: get(hObject,'Value') returns toggle state of Mf
%}}}
%Solution menu {{{
% --- Executes on selection change in popupmenu10.
function popupmenu10_Callback(hObject, eventdata, handles)
global solutiontype
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%figure out what results we are plotting here: 
value=get(hObject,'Value');
strings=get(hObject,'String');
solutiontype=strings(value);

plotm();
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu10

% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
global md solutiontype 
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if issealevel,
	if isearth,
		resultfields=fieldnames(md.earth.results);
	else
		resultfields=fieldnames(md.icecaps{1}.results);
	end
else
	resultfields=fieldnames(md.results);
end
if isempty(resultfields),
	resultfields={'None'};
end

%default solution: 
for i=1:length(resultfields),
	solutiontype=resultfields{i};
	if strcmpi(solutiontype,'StressbalanceSolution'),
		i0=i;
		break;
	elseif strcmpi(solutiontype,'StressbalanceSolution'),
		i0=i;
		break;
	else
		i0=1;
	end
end
set(hObject,'String',resultfields);
set(hObject,'Value',i0);

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end 
% }}}
%Icecaps menu: {{{
% --- Executes on selection change in popupmenu11.
function popupmenu11_Callback(hObject, eventdata, handles)
global md
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
strings=get(hObject,'String');
value=get(hObject,'Value');
continent=strings{value};

%update basin list: 
basins=md.basinsfromcontinent(continent);
strings={'All',basins{:}};
obj=findobj('Tag', 'popupmenu12');
set(obj,'String',strings');
set(obj,'Value',1);

plotm();
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu11 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu11


% --- Executes during object creation, after setting all properties.
function popupmenu11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global md range;

if issealevel(),
	strings=md.continents();
	set(hObject,'String',strings);
	set(hObject,'Value',1);
else
	set(hObject,'String',{});
	set(hObject,'Value',1);
end


% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end 
% }}}
%Basins menu{{{
% --- Executes on selection change in popupmenu12.
function popupmenu12_Callback(hObject, eventdata, handles)
global md 
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plotm();
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu12 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu12

% --- Executes during object creation, after setting all properties.
function popupmenu12_CreateFcn(hObject, eventdata, handles)
global md 
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if issealevel(),

	%need the continent name:  default first one.
	continents=md.continents();
	basins=md.basinsfromcontinent(continents{1});
	strings={'All',basins{:}};
	set(hObject,'String',strings);
	set(hObject,'Value',1);
else
	set(hObject,'String',{});
	set(hObject,'Value',1);
end
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% }}}
%Comparison: {{{

% --- Executes on button press in radiobutton17.
function radiobutton17_Callback(hObject, eventdata, handles)
global comparison
% hObject    handle to radiobutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
comparison=get(hObject,'Value');
plotm();
% Hint: get(hObject,'Value') returns toggle state of radiobutton17
% }}}
%Earth button: {{{

% --- Executes on button press in radiobutton18.
function radiobutton18_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global md earth 
earth=get(hObject,'value');

%Change Field: 
if issealevel, if isearth, results=md.earth.results; else results=md.icecaps{1}.results; end; else results=md.results; end;
fieldHandle=findobj('Tag', 'Field');
if isfield(results,'TransientSolution'),
	fieldss=listfields(results,'TransientSolution');
	set(fieldHandle,'String',fieldss);
	set(fieldHandle,'Value',1);
else 
	set(fieldHandle,'String',{});
	set(fieldHandle,'Value',1);
end

%Change Masque: 
if issealevel, 
	if isearth, 
		mv=md.earth.mask; 
	else 
		mv=md.icecaps{1}.mask; 
	end 
else 
	mv=md.mask; 
end
fieldHandle=findobj('Tag', 'Masque');
strings=fieldnames(mv);
strings={'all',strings{:}};
set(fieldHandle,'String',strings);

plotm();
% Hint: get(hObject,'Value') returns toggle state of radiobutton18
% --- Executes during object creation, after setting all properties.
function radiobutton18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global earth 
earth=1;
set(hObject,'Value',earth);
%}}}
%Reload model:  {{{
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global reload
reload=1;
plotm();
%}}}
%Scalar {{{
function Scalar_Callback(hObject, eventdata, handles)
% hObject    handle to Scalar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Scalar as text
%        str2double(get(hObject,'String')) returns contents of Scalar as a double

% --- Executes during object creation, after setting all properties.
function Scalar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Scalar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

	scalarHandle=findobj('Tag', 'Scalar'); 
	set(scalarHandle,'String',sprintf('%4.2f',NaN));

	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end 
% }}}
function plotv() % {{{
	global md   range

	%field:  % {{{ 
	fieldHandle=findobj('Tag', 'Field'); 
	fieldstrings=get(fieldHandle,'String');
	fieldvalue=get(fieldHandle,'Value');
	field=fieldstrings{fieldvalue};
	%}}}
		%lock limits:  % {{{
		lockHandle=findobj('Tag', 'Lock'); lockvalue=get(lockHandle,'Value');
		if lockvalue,
			xl=xlim; yl=ylim;
		end
		 %}}}
	
	cla, reset(gca);
	hold on ;
	if issealevel,
		if isearth,
			%do nothing. no ice volume.
		else
			dt=md.earth.timestepping.time_step*md.earth.settings.output_frequency;
			for j=1:length(range),
				i=range(j);
				vols=resultstomatrix(md.icecaps{i},'TransientSolution',field);
				if j==1,
					volst=vols;
				else
					volst(1,:)=volst(1,:)+vols(1,:);
				end
			end
			s1=subplot(2,1,1);
			set(s1,'Position',[0.3300    0.6100    0.4750    0.3]);
			plot(volst(2,:),(volst(1,:)-volst(1,1))*md.icecaps{1}.materials.rho_ice/1e12);
			xlabel('time (yr)'); ylabel('Mass (Gt)'); colorbar off;
			if lockvalue, xlim(xl), ylim(yl); end
			
			s2=subplot(2,1,2);
			set(s2,'Position',[0.3300    0.1100    0.4750    0.3]);
			dv=diff(volst(1,:))/dt;
			plot(volst(2,1:end-1),dv*md.icecaps{1}.materials.rho_ice/1e12);
			xlabel('time (yr)'); ylabel('Mass (Gt)'); colorbar off;
			if lockvalue, xlim(xl), ylim(yl); end
		end
	else 
		dt=md.timestepping.time_step*md.settings.output_frequency;
		vols=resultstomatrix(md,'TransientSolution',field);
		dv=diff(vols(1,:))/dt;
		
		s1=subplot(2,1,1);
		set(s1,'Position',[0.3300    0.6100    0.4750    0.3]);
		xlabel('time (yr)'); ylabel('Mass (Gt)'); colorbar off;
		plot(vols(2,:),(vols(1,:)-vols(1,1))*md.materials.rho_ice/1e12);
		if lockvalue, xlim(xl), ylim(yl); end

		s2=subplot(2,1,2);
		set(s2,'Position',[0.3300    0.1100    0.4750    0.3]);
		xlabel('time (yr)'); ylabel('Mass (Gt)'); colorbar off;
		plot(vols(2,1:end-1),dv*md.materials.rho_ice/1e12);
		if lockvalue, xlim(xl), ylim(yl); end
	end
%}}}
function displayscalar() % {{{
	global md   range

	%field:  % {{{ 
	fieldHandle=findobj('Tag', 'Field'); 
	fieldstrings=get(fieldHandle,'String');
	fieldvalue=get(fieldHandle,'Value');
	field=fieldstrings{fieldvalue};
	%}}}
	%counter:  % {{{ 
	stepHandle=findobj('Tag', 'Step'); 
	stepstrings=get(stepHandle,'String');
	stepvalue=get(stepHandle,'Value');
	if ~isempty(stepstrings),
		stepstring=stepstrings{stepvalue};
		%grab second integer: 
		A=sscanf(stepstring,'%g (%i)'); counter=A(2);
	else
		counter=1;
	end
	%}}}
			

	if issealevel,
		if isearth,
		if strcmpi(field,'SealevelRSLEustatic'),
			dt=md.earth.timestepping.time_step*md.earth.slr.geodetic_run_frequency;
			eus=md.earth.results.TransientSolution(counter-1).(field)/dt*1000; %in mm/yr
			scalarHandle=findobj('Tag', 'Scalar'); 
			set(scalarHandle,'String',sprintf('%4.5f mm/yr',eus));
		end
		else
			error('not supported yet!');
		end
	end
%}}}
function plotm() % {{{

	%retrieve field: 
	hObject=findobj('Tag', 'Field'); 
	strings=get(hObject,'String'); value=get(hObject,'Value');

	string='';
	if value, 
		if ~isempty(strings),
			string=strings(value); 
		end
	end

	if strcmpi(string,'IceVolume') | strcmpi(string,'IceVolumeAboveFloatation'),
		plotv();
	elseif strcmpi(string,'SealevelRSLEustatic'),
		displayscalar();
	else
		if issealevel(),
			plotsl();
		else 
			plotmd();
		end
	end
%}}}
	function plotmd() % {{{
		global logvalue md solutiontype comparison 

		%counter:  % {{{ 
		stepHandle=findobj('Tag', 'Step'); 
		stepstrings=get(stepHandle,'String');
		stepvalue=get(stepHandle,'Value');
		if ~isempty(stepstrings),
			stepstring=stepstrings{stepvalue};
			%grab second integer: 
			A=sscanf(stepstring,'%g (%i)'); counter=A(2);
		else
			counter=1;
		end
		%}}}
		%log:  % {{{ 
		logHandle=findobj('Tag', 'Log'); logvalue=get(logHandle,'Value');
		%}}}
		%shading:  %{{{
		interpHandle=findobj('Tag', 'Interp'); interpvalue=get(interpHandle,'Value');
		flatHandle=findobj('Tag', 'Flat'); flatvalue=get(flatHandle,'Value');
		facetedHandle=findobj('Tag', 'Faceted'); facetedvalue=get(facetedHandle,'Value');
		if interpvalue,
			shadingv='interp';
		elseif flatvalue,
			shadingv='flat';
		elseif facetedvalue,
			shadingv='faceted';
		else 
			shadingv='interp';
		end
		%}}}
		%mask:  %{{{
		maskfieldHandle=findobj('Tag', 'Masque'); 
		maskfieldstrings=get(maskfieldHandle,'String');
		maskfieldvalue=get(maskfieldHandle,'Value');
		maskfield=maskfieldstrings(maskfieldvalue); 
		maskfield=maskfield{:};
		if strcmpi(maskfield,'all'),
			maskv=ones(md.mesh.numberofvertices,1);
		else
			maskv=md.mask.(maskfield);
		end
		if strcmpi(maskfield,'all'),
			%do nothing; 
		elseif strcmpi(maskfield,'groundedice_levelset'),
			maskv=maskv>=0;
		elseif strcmpi(maskfield,'ice_levelset'),
			maskv=maskv<=0;
		elseif strcmpi(maskfield,'ocean_levelset'),
			maskv=maskv==1;
		elseif strcmpi(maskfield,'land_levelset'),
			maskv=maskv==1;
		elseif strcmpi(maskfield,'glacier_levelset'),
			if isnan(maskv),
				maskv=ones(md.mesh.numberofvertices,1);
			else
				maskv=maskv==1;
			end
		end
		maskvalue=0;
		%}}}
		%lock limits:  % {{{
		lockHandle=findobj('Tag', 'Lock'); lockvalue=get(lockHandle,'Value');
		if lockvalue,
			xl=xlim; yl=ylim;
		else
			xl=[min(md.mesh.x) max(md.mesh.x)];
			yl=[min(md.mesh.y) max(md.mesh.y)];
		end %}}}
		%time: {{{	
		timeHandle=findobj('Tag','Time');
		if strcmpi(solutiontype,'TransientSolution'),
			set(timeHandle,'String',sprintf('%4.2f',md.results.TransientSolution(counter).time));
		else
			set(timeHandle,'String',sprintf('%4.2f',0));
		end
		%}}}
		%diffcounter:  % {{{ 
		stepHandle=findobj('Tag', 'DiffStep'); 
		stepstrings=get(stepHandle,'String');
		stepvalue=get(stepHandle,'Value');
		if ~isempty(stepstrings),
			stepstring=stepstrings{stepvalue};
			%grab second integer: 
			A=sscanf(stepstring,'%g (%i)'); diffcounter=A(2);
		else
			diffcounter=1;
		end
		diffHandle=findobj('Tag', 'Diff'); 
		diff=get(diffHandle,'Value');
		%}}}
		%field:  % {{{ 
		fieldHandle=findobj('Tag', 'Field'); 
		fieldstrings=get(fieldHandle,'String');
		fieldvalue=get(fieldHandle,'Value');
		if ~isempty(fieldstrings),
			fieldv=fieldstrings{fieldvalue};
		else
			fieldv=NaN;
		end
		mfHandle=findobj('Tag', 'Mf'); 
		mf=get(mfHandle,'Value');
		if mf | strcmpi(solutiontype,'None'),
			fieldHandle=findobj('Tag', 'ModelFields'); 
			fieldstrings=get(fieldHandle,'String');
			fieldvalue=get(fieldHandle,'Value');
			fieldv=fieldstrings{fieldvalue};
			eval(['field=md.' fieldv ';']);
			if isnan(field),
				field='mesh';
			end
		else
			if strcmpi(solutiontype,'TransientSolution'),
				field=md.results.TransientSolution(counter).(fieldv);
				dfield=md.results.TransientSolution(diffcounter).(fieldv);
			elseif strcmpi(solutiontype,'StressbalanceSolution'),
				field=md.results.StressbalanceSolution.(fieldv);
				dfield=NaN;
			else 
				error('unknown solution type!');
			end
		end
		if comparison, 
			if strcmpi(solutiontype,'None'),
				error('cannot compare a solution field to model as no solution was run or selected!');
			end
			%figure out second field:
			if strcmpi(fieldv,'Vel'),
				field2=md.inversion.vel_obs;
			end
		end

		%}}}
		%color limits:  % {{{
		cminHandle=findobj('Tag','Cmin'); cmin=str2num(get(cminHandle,'String'));
		cmaxHandle=findobj('Tag','Cmax'); cmax=str2num(get(cmaxHandle,'String'));

		if isnan(cmin)
			if logvalue,
				cmin=.1;
			else
				pos=find(maskv);
				cmin=min(field(pos));
			end
		end
		if isnan(cmax)
			pos=find(maskv);
			cmax=max(field(pos));
		end
		colaxis=[cmin,cmax];
		%}}}

		cla;
		if logvalue,
			plotmodel(md,'data',field,'figurestatement','off','clf','off','mask',maskv,'caxis',colaxis,'log',10,'xlim',xl,'ylim',yl,'shading',shadingv);
		else
			if ~diff,
				if comparison,
					plotmodel(md,'data',field,'data',field2,'figurestatement','off','clf','off','mask#all',maskv,'caxis#all',colaxis,'xlim#all',xl,'ylim#all',yl,'shading#all',shadingv,'nlines',2,'ncols',1);
				else
					if ischar(field) & strcmpi(field,'mesh'),
						plotmodel(md,'data','mesh','figurestatement','off','clf','off','xlim',xl,'ylim',yl);
					else
						plotmodel(md,'data',field,'figurestatement','off','clf','off','mask',maskv,'caxis',colaxis,'xlim',xl,'ylim',yl,'shading',shadingv); 
					end
				end
			else
				plotmodel(md,'data',dfield-field,'figurestatement','off','clf','off','mask',maskv,'caxis',colaxis,'xlim',xl,'ylim',yl,'shading',shadingv);
			end

		end
		set(gca,'Position',[0.2300    0.1100    0.7750    0.8150]);
	%}}}
function plotsl() % {{{
	global logvalue md solutiontype comparison 

	%counter:  % {{{ 
	stepHandle=findobj('Tag', 'Step'); 
	stepstrings=get(stepHandle,'String');
	stepvalue=get(stepHandle,'Value');
	if ~isempty(stepstrings),
		stepstring=stepstrings{stepvalue};
		%grab second integer: 
		A=sscanf(stepstring,'%g (%i)'); counter=A(2);
	else
		counter=1;
	end
	%}}}
	%log:  % {{{ 
	logHandle=findobj('Tag', 'Log'); logvalue=get(logHandle,'Value');
	%}}}
	%shading:  %{{{
	interpHandle=findobj('Tag', 'Interp'); interpvalue=get(interpHandle,'Value');
	flatHandle=findobj('Tag', 'Flat'); flatvalue=get(flatHandle,'Value');
	facetedHandle=findobj('Tag', 'Faceted'); facetedvalue=get(facetedHandle,'Value');
	if interpvalue,
		shadingv='interp';
	elseif flatvalue,
		shadingv='flat';
	elseif facetedvalue,
		shadingv='faceted';
	else 
		shadingv='interp';
	end
	%}}}
	%time: {{{	
	timeHandle=findobj('Tag','Time');
	if strcmpi(solutiontype,'TransientSolution'),
		if isearth,
			set(timeHandle,'String',sprintf('%4.2f',md.earth.results.TransientSolution(counter).time));
		else
			set(timeHandle,'String',sprintf('%4.2f',md.icecaps{1}.results.TransientSolution(counter).time));
		end
	else
		set(timeHandle,'String',sprintf('%4.2f',0));
	end
	%}}}
	%diffcounter:  % {{{ 
	stepHandle=findobj('Tag', 'DiffStep'); 
	stepstrings=get(stepHandle,'String');
	stepvalue=get(stepHandle,'Value');
	if ~isempty(stepstrings),
		stepstring=stepstrings{stepvalue};
		%grab second integer: 
		A=sscanf(stepstring,'%g (%i)'); diffcounter=A(2);
	else
		diffcounter=1;
	end
	diffHandle=findobj('Tag', 'Diff'); 
	diff=get(diffHandle,'Value');
	%}}}

	if isearth, 
		%Earth plotting: {{{
		%field:  % {{{ 
		fieldHandle=findobj('Tag', 'Field'); 
		fieldstrings=get(fieldHandle,'String');
		fieldvalue=get(fieldHandle,'Value');
		if ~isempty(fieldstrings),
			fieldv=fieldstrings{fieldvalue};
		else
			fieldv=NaN;
		end
		colorbartitle={fieldv,''};
		if strcmpi(fieldv,'SealevelRSLRate'),
			colorbartitle={'SealevelRSLRate (mm/yr)',''};
		end

		mfHandle=findobj('Tag', 'Mf'); 
		mf=get(mfHandle,'Value');
		if mf | strcmpi(solutiontype,'None'),
			fieldHandle=findobj('Tag', 'ModelFields'); 
			fieldstrings=get(fieldHandle,'String');
			fieldvalue=get(fieldHandle,'Value');
			fieldv=fieldstrings{fieldvalue};
			eval(['field=md.earth.' fieldv ';']);
			if isnan(field),
				field='mesh';
			end
		else
			if strcmpi(solutiontype,'TransientSolution'),
				field=md.earth.results.TransientSolution(counter).(fieldv);
				dfield=md.earth.results.TransientSolution(diffcounter).(fieldv);
			elseif strcmpi(solutiontype,'StressbalanceSolution'),
				field=md.earth.results.StressbalanceSolution.(fieldv);
				dfield=NaN;
			else 
				error('unknown solution type!');
			end
		end
		if comparison, 
			if strcmpi(solutiontype,'None'),
				error('cannot compare a solution field to model as no solution was run or selected!');
			end
			%figure out second field:
			if strcmpi(fieldv,'Vel'), field2=md.earth.inversion.vel_obs; end
			if strcmpi(fieldv,'Base'), field2=md.earth.geometry.base; end
		end
		%}}}
		%mask:  %{{{
		maskfieldHandle=findobj('Tag', 'Masque'); 
		maskfieldstrings=get(maskfieldHandle,'String');
		maskfieldvalue=get(maskfieldHandle,'Value');
		maskfield=maskfieldstrings(maskfieldvalue); 
		maskfield=maskfield{:};
		if strcmpi(maskfield,'all'),
			maskv=ones(md.earth.mesh.numberofvertices,1);
		else
			maskv=md.earth.mask.(maskfield);
		end
		if strcmpi(maskfield,'all'),
			%do nothing; 
		elseif strcmpi(maskfield,'groundedice_levelset'),
			maskv=maskv>=0;
		elseif strcmpi(maskfield,'ice_levelset'),
			maskv=maskv<=0;
		elseif strcmpi(maskfield,'ocean_levelset'),
			maskv=maskv==1;
		elseif strcmpi(maskfield,'land_levelset'),
			maskv=maskv==1;
		elseif strcmpi(maskfield,'glacier_levelset'),
			if isnan(maskv),
				maskv=ones(md.earth.mesh.numberofvertices,1);
			else
				maskv=maskv==1;
			end
		end
		maskvalue=0;
		%}}}
		%color limits:  % {{{
		if ~ischar(field),
			cminHandle=findobj('Tag','Cmin'); cmin=str2num(get(cminHandle,'String'));
			cmaxHandle=findobj('Tag','Cmax'); cmax=str2num(get(cmaxHandle,'String'));

			if isnan(cmin)
				if logvalue,
					cmin=.1;
				else
					pos=find(maskv);
					cmin=min(field(pos));
				end
			end
			if isnan(cmax)
				pos=find(maskv);
				cmax=max(field(pos));
			end

			if cmin==cmax,
				colaxis=[cmin-eps,cmax+eps];
			else
				colaxis=[cmin,cmax];
			end
		end

		%}}}
		%lock limits:  % {{{
		lockHandle=findobj('Tag', 'Lock'); lockvalue=get(lockHandle,'Value');
		if lockvalue,
			xl=xlim; yl=ylim; zl=zlim;
			[az,el]=view; v=[az,el];
		else
			xl=[min(md.earth.mesh.x) max(md.earth.mesh.x)];
			yl=[min(md.earth.mesh.y) max(md.earth.mesh.y)];
			zl=[min(md.earth.mesh.z) max(md.earth.mesh.z)];
			v=[57,50];
		end %}}}

		%some quirks: 
		cla;
		if logvalue,
			plotmodel(md.earth,'data',field,'figurestatement','off','clf','off','mask',maskv,'caxis',colaxis,'log',10,'xlim',xl,'ylim',yl,'zlim',zl,'shading',shadingv,'colorbar','on','view',v,'colorbartitle',colorbartitle);
		else
			if ~diff, 
				if comparison,
					plotmodel(md.earth,'data',field,'data',field2,'figurestatement','off','clf','off','mask#all',maskv,'caxis#all',colaxis,'xlim#all',xl,'ylim#all',yl,'zlim#all',zl,'shading#all',shadingv,'nlines',2,'ncols',1,'view#all',v,'colorbartitle',colorbartitle);
				else
					if ischar(field) & strcmpi(field,'mesh'),
						plotmodel(md.earth,'data','mesh','figurestatement','off','clf','off','xlim',xl,'ylim',yl,'zlim',zl,'view#all',v,'colorbartitle',colorbartitle);
					else
						plotmodel(md.earth,'data',field,'figurestatement','off','clf','off','mask',maskv,'maskvalue',maskvalue,'caxis',colaxis,'xlim',xl,'ylim',yl,'zlim',zl,'shading',shadingv,'axes','equal','colorbarpos',[.95 .5 .01 .15],'view#all',v,'colorbartitle',colorbartitle);
					end
				end
			else
				plotmodel(md.earth,'data',dfield-field,'figurestatement','off','clf','off','mask',maskv,'caxis',colaxis,'xlim',xl,'ylim',yl,'zlim',zl,'shading',shadingv,'view#all',v,'colorbartitle',colorbartitle);
			end
		end

		%coastlines and graticule: {{{
		load coastlines
		[x,y,z]=AboveGround(coastlat,coastlon,md.earth.mesh.r(1),1000);
		hold on; p=plot3(x,y,z,'k-'); set(p,'Color','k');

		%graticule:
		grat=graticule(30,40,1);
		[x,y,z]=AboveGround(grat.lat,grat.long,md.earth.mesh.r(1),1000);
		hold on, p=plot3(x,y,z,'k.-','MarkerSize',.01); set(p,'Color','k');
		%}}}
		%}}}
	else 
		%IceCaps plotting:  {{{

		%counter:  % {{{ 
		stepHandle=findobj('Tag', 'Step'); 
		stepstrings=get(stepHandle,'String');
		stepvalue=get(stepHandle,'Value');
		if ~isempty(stepstrings),
			stepstring=stepstrings{stepvalue};
			%grab second integer: 
			A=sscanf(stepstring,'%g (%i)'); counter=A(2);
		else
			counter=1;
		end
		%}}}
		%log:  % {{{ 
		logHandle=findobj('Tag', 'Log'); logvalue=get(logHandle,'Value');
		%}}}
		%shading:  %{{{
		interpHandle=findobj('Tag', 'Interp'); interpvalue=get(interpHandle,'Value');
		flatHandle=findobj('Tag', 'Flat'); flatvalue=get(flatHandle,'Value');
		facetedHandle=findobj('Tag', 'Faceted'); facetedvalue=get(facetedHandle,'Value');
		if interpvalue,
			shadingv='interp';
		elseif flatvalue,
			shadingv='flat';
		elseif facetedvalue,
			shadingv='faceted';
		else 
			shadingv='interp';
		end
		%}}}
		%time: {{{	
		timeHandle=findobj('Tag','Time');
		if strcmpi(solutiontype,'TransientSolution'),
			set(timeHandle,'String',sprintf('%4.2f',md.icecaps{1}.results.TransientSolution(counter).time));
		else
			set(timeHandle,'String',sprintf('%4.2f',0));
		end
		%}}}
		%diffcounter:  % {{{ 
		stepHandle=findobj('Tag', 'DiffStep'); 
		stepstrings=get(stepHandle,'String');
		stepvalue=get(stepHandle,'Value');
		if ~isempty(stepstrings),
			stepstring=stepstrings{stepvalue};
			%grab second integer: 
			A=sscanf(stepstring,'%g (%i)'); diffcounter=A(2);
		else
			diffcounter=1;
		end
		diffHandle=findobj('Tag', 'Diff'); 
		diff=get(diffHandle,'Value');
		%}}}
		%range:  % {{{
		basins=getbasins();
		if strcmpi(basins,'All'),
			range=md.basinindx('continent',{getcontinent()});
		else
			range=md.basinindx('basin',basins);
		end

		%}}}
		cla;
		for j=1:length(range), 
			i=range(j);
			%mask:  %{{{
			maskfieldHandle=findobj('Tag', 'Masque'); 
			maskfieldstrings=get(maskfieldHandle,'String');
			maskfieldvalue=get(maskfieldHandle,'Value');
			maskfield=maskfieldstrings(maskfieldvalue); 
			maskfield=maskfield{:};
			if strcmpi(maskfield,'all'),
				maskv=ones(md.icecaps{i}.mesh.numberofvertices,1);
			else
				maskv=md.icecaps{i}.mask.(maskfield);
			end
			if strcmpi(maskfield,'all'),
				%do nothing; 
			elseif strcmpi(maskfield,'groundedice_levelset'),
				maskv=maskv>=0;
			elseif strcmpi(maskfield,'ice_levelset'),
				maskv=maskv<=0;
			elseif strcmpi(maskfield,'ocean_levelset'),
				maskv=maskv==1;
			elseif strcmpi(maskfield,'land_levelset'),
				maskv=maskv==1;
			elseif strcmpi(maskfield,'glacier_levelset'),
				if isnan(maskv),
					maskv=ones(md.icecaps{i}.mesh.numberofvertices,1);
				else
					maskv=maskv==1;
				end
			end
			maskvalue=0;
			%}}}
			%lock limits:  % {{{
			lockHandle=findobj('Tag', 'Lock'); lockvalue=get(lockHandle,'Value');
			if lockvalue,
				xl=xlim; yl=ylim;
			else
				if j==1,
					xl=[min(md.icecaps{i}.mesh.x) max(md.icecaps{i}.mesh.x)];
				else
					xl=[min([md.icecaps{i}.mesh.x;xl(1)]) max([md.icecaps{i}.mesh.x;xl(2)])];
				end
				if j==1,
					yl=[min(md.icecaps{i}.mesh.y) max(md.icecaps{i}.mesh.y)];
				else
					yl=[min([md.icecaps{i}.mesh.y; yl(1)]) max([md.icecaps{i}.mesh.y; yl(2)])];
				end
			end %}}}
			%field:  % {{{ 
			fieldHandle=findobj('Tag', 'Field'); 
			fieldstrings=get(fieldHandle,'String');
			fieldvalue=get(fieldHandle,'Value');
			if ~isempty(fieldstrings),
				fieldv=fieldstrings{fieldvalue};
			else
				fieldv=NaN;
			end
			mfHandle=findobj('Tag', 'Mf'); 
			mf=get(mfHandle,'Value');
			if mf | strcmpi(solutiontype,'None'),
				fieldHandle=findobj('Tag', 'ModelFields'); 
				fieldstrings=get(fieldHandle,'String');
				fieldvalue=get(fieldHandle,'Value');
				fieldv=fieldstrings{fieldvalue};
				eval(['field=md.icecaps{' num2str(i) '}.' fieldv ';']);
				if isnan(field),
					field='mesh';
				end
			else
				if strcmpi(solutiontype,'TransientSolution'),
					field=md.icecaps{i}.results.TransientSolution(counter).(fieldv);
					dfield=md.icecaps{i}.results.TransientSolution(diffcounter).(fieldv);
				elseif strcmpi(solutiontype,'StressbalanceSolution'),
					field=md.icecaps{i}.results.StressbalanceSolution.(fieldv);
					dfield=NaN;
				else 
					error('unknown solution type!');
				end
			end
			if comparison, 
				if strcmpi(solutiontype,'None'),
					error('cannot compare a solution field to model as no solution was run or selected!');
				end
				%figure out second field:
				if strcmpi(fieldv,'Vel'), field2=md.icecaps{i}.inversion.vel_obs; end
				if strcmpi(fieldv,'Base'), field2=md.icecaps{i}.geometry.base; end
				if strcmpi(fieldv,'Bed'), field2=md.icecaps{i}.geometry.bed; end
				if strcmpi(fieldv,'Thickness'), field2=md.icecaps{i}.geometry.thickness; end
				if strcmpi(fieldv,'Surface'), field2=md.icecaps{i}.geometry.surface; end
			end

			%}}}
			%color limits:  % {{{
			cminHandle=findobj('Tag','Cmin'); cminfield=str2num(get(cminHandle,'String'));
			cmaxHandle=findobj('Tag','Cmax'); cmaxfield=str2num(get(cmaxHandle,'String'));

			if isnan(cminfield)
				if logvalue,
					cmin=.1;
				else
					pos=find(maskv);
					if j==1,
						cmin=min(field(pos));
					else
						cmin=min(min(field(pos)),cmin);
					end
				end
			else
				cmin=cminfield;
			end
			if isnan(cmaxfield)
				pos=find(maskv);
				if j==1
					cmax=max(field(pos));
				else
					cmax=max(max(field(pos)),cmax);
				end
			else
				cmax=cmaxfield;
			end
			if cmin==cmax,
				cmin=cmin-1e-10; cmax=cmax+1e-10;
			end
			colaxis=[cmin,cmax];
			%}}}
			%contour levels? {{{
			contourlevels=0;
			if strncmpi(fieldv,'Mask',4),
				contourexp=[tempname '.exp'];
				isoline(md.icecaps{i},field,'output',contourexp);
				contourlevels=1;

				if diff,
					contourdiffexp=[tempname '.exp'];
					isoline(md.icecaps{i},dfield,'output',contourdiffexp);
				end
			end
			%}}}
			%display : {{{ 
			if logvalue,
				plotmodel(md.icecaps{i},'data',field,'figurestatement','off','clf','off','mask',maskv,'caxis',colaxis,'log',10,'xlim',xl,'ylim',yl,'shading',shadingv);
			else
				if ~diff,
					if comparison,
						plotmodel(md.icecaps{i},'data',field,'data',field2,'figurestatement','off','clf','off','mask#all',maskv,'caxis#all',colaxis,'xlim#all',xl,'ylim#all',yl,'shading#all',shadingv,'nlines',2,'ncols',1);
					else
						if ischar(field) & strcmpi(field,'mesh'),
							plotmodel(md.icecaps{i},'data','mesh','figurestatement','off','clf','off','xlim',xl,'ylim',yl);
						else
							plotmodel(md.icecaps{i},'data',field,'figurestatement','off','clf','off','mask',maskv,'caxis',colaxis,'xlim',xl,'ylim',yl,'shading',shadingv); 
							if contourlevels,
								contour=expread(contourexp);
								hold on; 
								for k=1:length(contour),
									plot(contour(k).x,contour(k).y,'r-');
								end
							end
						end
					end
				else
					plotmodel(md.icecaps{i},'data',dfield-field,'figurestatement','off','clf','off','mask',maskv,'caxis',colaxis,'xlim',xl,'ylim',yl,'shading',shadingv);
					if contourlevels,
						contour=expread(contourexp);
						hold on; 
						for k=1:length(contour),
							plot(contour(k).x,contour(k).y,'k-');
						end
						dcontour=expread(contourdiffexp);
						hold on; 
						for k=1:length(dcontour),
							plot(dcontour(k).x,dcontour(k).y,'r-');
						end
					end
				end
			end 
			% }}} 

		end %end for

		if ~comparison,
			set(gca,'Position',[0.2300    0.1100    0.7750    0.8150]); 
		end % }}}
	end 
% }}}
function result=isearth() % {{{
	global earth;
	result=earth;
%}}}
function result=issealevel() % {{{
	global md;
	if strcmpi(class(md),'sealevelmodel'),
		result=1;
	else 
		result=0;
	end
%}}}
function  continent=getcontinent() % {{{

	continentHandle=findobj('Tag', 'popupmenu11');
	strings=get(continentHandle,'String');
	value=get(continentHandle,'Value');
	continent=strings{value};

%}}}
function  basins=getbasins() % {{{

	basinHandle=findobj('Tag', 'popupmenu12'); 
	strings=get(basinHandle,'String');
	value=get(basinHandle,'Value');
	basins=strings{value};
% }}}
