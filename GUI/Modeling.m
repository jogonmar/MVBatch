function varargout = Modeling(varargin)
% MODELING M-file for Modeling.fig
%      MODELING, by itself, creates a new MODELING or raises the existing
%      singleton*.
%
%      H = MODELING returns the handle to a new MODELING or the handle to
%      the existing singleton*.
%
%      MODELING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODELING.M with the given input arguments.
%
%      MODELING('Property','Value',...) creates a new MODELING or raises the
%      existing singleton*.  Starting from the left, property value pairs
%      are
%      applied to the GUI before mpf_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Modeling_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Modeling

% Last Modified by GUIDE v2.5 08-Nov-2016 16:33:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Modeling_OpeningFcn, ...
                   'gui_OutputFcn',  @Modeling_OutputFcn, ...
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



% --- Executes just before Modeling is made visible.
function Modeling_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Modeling (see VARARGIN)

% Choose default command line output for Modeling

handles.output = hObject;

handles.PCAh = 0;

if length(varargin)>0,   
    handles.ParentsWindow=varargin{1};
    handles.ParentFigure = guidata(handles.ParentsWindow);

    handles.data.x = handles.ParentFigure.s_alignment.alg_batches(:,2:end,:);

    if ~isempty(handles.ParentFigure.s_calibration)
        models = handles.ParentFigure.s_calibration;
        handles.data.man_mp_group = models.man_mp_group;
        handles.data.mp_group = models.mp_group;
        handles.data.mp_group2 = models.mp_group2;
       
        set(handles.pushbuttonMod,'Enable','on');
        set(handles.pushbuttonApply,'Enable','on');
        set(handles.popupmenuMod,'Enable','on');
        set(handles.textMod,'Enable','on');
        set(handles.popupmenuPhase,'Enable','on');
        set(handles.textPhase,'Enable','on');
        set(handles.radiobutton1,'Enable','on');
        
        contents = get(handles.popupmenuCM,'String');
        contents = strvcat(contents(1:3,:));
        set(handles.popupmenuCM,'String',strvcat(contents,' Model Total Covariance'));
        contents = get(handles.popupmenuCM,'String');
        set(handles.popupmenuCM,'String',strvcat(contents,' Model Dynamic Partial Covariance'));
        contents = get(handles.popupmenuCM,'String');
        set(handles.popupmenuCM,'String',strvcat(contents,' Model Partial Covariance'));
        contents = get(handles.popupmenuCM,'String');
        set(handles.popupmenuCM,'String',strvcat(contents,' Residuals Total Covariance'));
        contents = get(handles.popupmenuCM,'String');
        set(handles.popupmenuCM,'String',strvcat(contents,' Residuals Dynamic Partial Covariance'));
        contents = get(handles.popupmenuCM,'String');
        set(handles.popupmenuCM,'String',strvcat(contents,' Residuals Partial Covariance'));

        set(handles.popupmenuMod,'String','');
        for i=1:length(handles.data.man_mp_group),
            contents = get(handles.popupmenuMod,'String');
            set(handles.popupmenuMod,'String',strvcat(contents,sprintf(' Manual MPPCA Model %d',i)));
        end
        
        for i=1:length(handles.data.mp_group),
            contents = get(handles.popupmenuMod,'String');
            set(handles.popupmenuMod,'String',strvcat(contents,sprintf(' MPPCA Model %d',i)));
        end

        set(handles.popupmenuPhase,'Value',1);
        set(handles.popupmenuPhase,'String','');
        if length(handles.data.man_mp_group)<1,
            phases=handles.data.mp_group{1}.phases;
        else
            phases=handles.data.man_mp_group{1}.phases;
        end
        for i=1:size(phases,1),
            contents = get(handles.popupmenuPhase,'String');
            set(handles.popupmenuPhase,'String',strvcat(contents,sprintf(' Phase %d: from %d to %d, %d PC(s) and %d LMV(s)',i,phases(i,4),phases(i,5),phases(i,2),phases(i,3))));
        end
   
        set(handles.text_Tm,'Enable','on');
        set(handles.edit_Tm,'Enable','on');
        set(handles.popupmenuM,'Enable','on');
        set(handles.pushbutton_SS,'Enable','on');

        if  length(handles.data.mp_group2)>0,
            for i=1:length(handles.data.mp_group2),
                contents = get(handles.popupmenuMod,'String');
                set(handles.popupmenuMod,'String',strvcat(contents,sprintf(' Merged MPPCA Model %d',i)));
            end

            if length(handles.data.mp_group2)>2,
                [P,TABLE,STATS] = anova2(handles.data.anova,1,'off');
                h=figure;
                multcompare(STATS,'ctype','lsd');
            end
        end
    else      
        handles.data.man_mp_group = {};
        handles.data.mp_group = {};
        handles.data.mp_group2 = {};
    end
 
    handles.data.text = [];
    handles.data.prep = 2;
    handles.data.preptxt = 'TCTS';

    s = size(handles.data.x);  % Fix the selection of blocks in the cross-validation
    cross = cross_parameters; 
    if cross.order.input==false,
        cross.order.input=true;
        cross.order.cols=rand(1,s(2)*s(1));
        cross.order.rows=rand(1,s(3)); 
    end
    cross.leave_m = 'ckf';
    handles.data.cross = cross;
   
    if ~isempty(find(isnan(handles.data.x))), 
        handles.data.x = missTSR3D(handles.data.x,3,s(1)-1,2);
    end
 
    % Update handles structure
    guidata(hObject, handles);
    
    % Center GUI
    set(gcf,'Units', 'pixels' );
    %get your display size
    screenSize = get(0, 'ScreenSize');
    %calculate the center of the display
    position = get( gcf,'Position' );
    position(1) = (screenSize(3)-position(3))/2;
    position(2) = (screenSize(4)-position(4))/2;
    %center the window
    set( gcf,'Position', position );
   
else
    file = LoadMenuItem_Callback(hObject, eventdata, handles);
    if ~file, error('No input data'); end;
end


% UIWAIT makes Modeling wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Modeling_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function edit_vars_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vars as text
%        str2double(get(hObject,'String')) returns contents of edit_vars as a double

val = get(hObject,'String');
if isequal(val,'all'),
    handles.cm.vars = 1:size(handles.data.x,2);
else
    handles.cm.vars = str2num(get(hObject,'String'));
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_vars_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_init_Callback(hObject, eventdata, handles)
% hObject    handle to edit_init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_init as text
%        str2double(get(hObject,'String')) returns contents of edit_init as a double

init =  str2num(get(hObject,'String'));
handles.cm.init = init;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_init_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fint_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fint as text
%        str2double(get(hObject,'String')) returns contents of edit_fint as a double

val = get(hObject,'String');
if length(val)>2 && isequal(val(end-2:end),'end'),
    val = [val(1:end-3) num2str(size(handles.data.x,1))];
end

handles.cm.fint = str2num(val);

guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function edit_fint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% function edit_init2_Callback(hObject, eventdata, handles)
% % hObject    handle to edit_init2 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edit_init2 as text
% %        str2double(get(hObject,'String')) returns contents of edit_init2 as a double
% 
% init =  str2num(get(hObject,'String'));
% handles.mv.init = init;
% guidata(hObject,handles);
% 
% radiobutton2_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_init2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_init2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% function edit_fint2_Callback(hObject, eventdata, handles)
% % hObject    handle to edit_fint2 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edit_fint2 as text
% %        str2double(get(hObject,'String')) returns contents of edit_fint2 as a double
% 
% val = get(hObject,'String');
% if length(val)>2 && isequal(val(end-2:end),'end'),
%     val = [val(1:end-3) num2str(size(handles.data.x,1))];
% end
% 
% handles.mv.fint = str2num(val);
% guidata(hObject,handles);
% 
% radiobutton2_Callback(hObject, eventdata, handles)




% --- Executes during object creation, after setting all properties.
function edit_fint2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fint2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function edit_LMVs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_LMVs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_LMVs as text
%        str2double(get(hObject,'String')) returns contents of edit_LMVs as a double

lags =  str2num(get(hObject,'String'));
handles.data.lags = lags;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_LMVs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LMVs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

edit_LMVs_Callback(hObject, eventdata, handles)


function edit_T_Callback(hObject, eventdata, handles)
% hObject    handle to edit_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_T as text
%        str2double(get(hObject,'String')) returns contents of edit_T as a double

Ts =  str2num(get(hObject,'String'));
handles.data.Ts = Ts;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_T_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

edit_T_Callback(hObject, eventdata, handles)



function edit_k_Callback(hObject, eventdata, handles)
% hObject    handle to edit_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_k as text
%        str2double(get(hObject,'String')) returns contents of edit_k as a double

ks =  str2num(get(hObject,'String'));
handles.data.ks = ks;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_k_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

edit_k_Callback(hObject, eventdata, handles)



function edit_minL_Callback(hObject, eventdata, handles)
% hObject    handle to edit_minL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_minL as text
%        str2double(get(hObject,'String')) returns contents of edit_minL as a double

minLs =  str2num(get(hObject,'String'));
handles.data.minLs = minLs;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_minL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

edit_minL_Callback(hObject, eventdata, handles)



function edit_a_Callback(hObject, eventdata, handles)
% hObject    handle to edit_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_a as text
%        str2double(get(hObject,'String')) returns contents of edit_a as a double

as =  str2num(get(hObject,'String'));
handles.data.as = as;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

edit_a_Callback(hObject, eventdata, handles)


% --- Executes on button press in pushbutton_FS.
function pushbutton_FS_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_FS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.text = [];

i=1;
for i_lags=1:length(handles.data.lags),
    for i_Ts=1:length(handles.data.Ts),
        for i_ks=1:length(handles.data.ks),
            for i_minLs=1:length(handles.data.minLs),
                for i_as=1:length(handles.data.as),
                    text = sprintf('Generating model: LMVs = %d, T = %g, k = %d, minL = %d, a = %d, %s, %s.',handles.data.lags(i_lags),handles.data.Ts(i_Ts),handles.data.ks(i_ks),handles.data.minLs(i_minLs),handles.data.as(i_as),handles.data.preptxt,handles.data.cross.leave_m);
                    handles.data.text = cprintMV(handles.console,text,handles.data.text,0);  
                    [mp_group{i},handles.data.text] = mppca2_s(handles.data.x,handles.data.lags(i_lags),handles.data.Ts(i_Ts),true,handles.data.ks(i_ks),handles.data.minLs(i_minLs),handles.data.as(i_as),handles.data.prep,handles.data.cross,handles.console,handles.data.text);
                    i = 1+i;
                end
            end
        end
    end
end

handles.data.text = cprintMV(handles.console,sprintf('%d models generated.',length(mp_group)),handles.data.text,0,2);

handles.data.mp_group = mp_group;

set(handles.pushbuttonMod,'Enable','on');
set(handles.pushbuttonApply,'Enable','on');
set(handles.popupmenuMod,'Enable','on');
set(handles.textMod,'Enable','on');
set(handles.popupmenuPhase,'Enable','on');
set(handles.textPhase,'Enable','on');
set(handles.radiobutton1,'Enable','on');

contents = get(handles.popupmenuCM,'String');
if size(contents,1)<=3,
    contents = strvcat(contents);
    set(handles.popupmenuCM,'String',strvcat(contents,' Model Total Covariance'));
    contents = get(handles.popupmenuCM,'String');
    set(handles.popupmenuCM,'String',strvcat(contents,' Model Dynamic Partial Covariance'));
    contents = get(handles.popupmenuCM,'String');
    set(handles.popupmenuCM,'String',strvcat(contents,' Model Partial Covariance'));
    contents = get(handles.popupmenuCM,'String');
    set(handles.popupmenuCM,'String',strvcat(contents,' Residuals Total Covariance'));
    contents = get(handles.popupmenuCM,'String');
    set(handles.popupmenuCM,'String',strvcat(contents,' Residuals Dynamic Partial Covariance'));
    contents = get(handles.popupmenuCM,'String');
    set(handles.popupmenuCM,'String',strvcat(contents,' Residuals Partial Covariance'));
end
    
set(handles.popupmenuMod,'Value',1);
contents = get(handles.popupmenuMod,'String');
set(handles.popupmenuMod,'String',contents(1:length(handles.data.man_mp_group),:));
for i=1:length(mp_group),
    contents = get(handles.popupmenuMod,'String');
    set(handles.popupmenuMod,'String',strvcat(contents,sprintf(' MPPCA Model %d',i)));
end

popupmenuMod_Callback(handles.popupmenuMod, eventdata, handles)

set(handles.text_Tm,'Enable','on');
set(handles.edit_Tm,'Enable','on');
set(handles.popupmenuM,'Enable','on');
set(handles.pushbutton_SS,'Enable','on');
    
% if length(mp_group)>1,
%     set(handles.text_Tm,'Enable','on');
%     set(handles.edit_Tm,'Enable','on');
%     set(handles.popupmenuM,'Enable','on');
%     set(handles.pushbutton_SS,'Enable','on');
% else
%     set(handles.text_Tm,'Enable','off');
%     set(handles.edit_Tm,'Enable','off');
%     set(handles.popupmenuM,'Enable','off');
%     set(handles.pushbutton_SS,'Enable','off');
% end

guidata(hObject,handles);


% --- Executes on button press in pushbutton_SS.
function pushbutton_SS_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_SS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mp_group2={};
s = size(handles.data.x);

handles.data.text = cprintMV(handles.console,' Merging. Please, be patient.',handles.data.text,0);

for i_Tms=1:length(handles.data.Tms),
    handles.data.text = cprintMV(handles.console,sprintf('Merging model %d...',i_Tms),handles.data.text);
    mp_group2{i_Tms} = fastmix_MP(reduce_group(handles.data.mp_group),handles.data.Tms(i_Tms),handles.data.criterium,true,handles.data.cross);
    handles.data.text = cprintMV(handles.console,'Done.',handles.data.text,2);
end
 
mp_group2 = reduce_group(mp_group2);
anova=[];
cols=[];
for o=1:length(mp_group2),
    anova=[anova squeeze(sum(sum(mp_group2{o}.pem.^2,2),1))];
    cols=[cols mp_group2{o}.arg.Tm];
end
rows=(1:s(3))';

handles.data.text = cprintMV(handles.console,sprintf('%d different models generated.',length(mp_group2)),handles.data.text);

contents = get(handles.popupmenuMod,'String');
set(handles.popupmenuMod,'Value',1);
set(handles.popupmenuMod,'String',contents(1:(length(handles.data.mp_group)+length(handles.data.man_mp_group)),:));
for i=1:length(mp_group2),
    contents = get(handles.popupmenuMod,'String');
    set(handles.popupmenuMod,'String',strvcat(contents,sprintf(' Merged MPPCA Model %d',i)));
end

handles.data.mp_group2 = mp_group2;
handles.data.anova = anova;
handles.data.rows = rows;
handles.data.cols = cols;
guidata(hObject,handles);

if length(mp_group2)>1,
    [P,TABLE,STATS] = anova2(handles.data.anova,1,'off');
    figure
    multcompare(STATS,'ctype','lsd');
end


function console_Callback(hObject, eventdata, handles)
% hObject    handle to console (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of console as text
%        str2double(get(hObject,'String')) returns contents of console as a double


% --- Executes during object creation, after setting all properties.
function console_CreateFcn(hObject, eventdata, handles)
% hObject    handle to console (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_Tm_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Tms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Tms as text
%        str2double(get(hObject,'String')) returns contents of edit_Tms as a double

Tms =  str2num(get(hObject,'String'));
handles.data.Tms = Tms;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_Tm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Tms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

edit_Tm_Callback(hObject, eventdata, handles)


% --- Executes on selection change in popupmenuM.
function popupmenuM_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenuM contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuM

contents = get(hObject,'String');
txt=contents{get(hObject,'Value')};
 
switch txt,
    case ' Model parsimony',
        criterium = 'parsimony';
    case ' Cov. matrix parsimony',
        criterium = 'parsimony2';
    case ' Minimize LMVs',
        criterium = 'LMV';
    case ' Minimize phases',
        criterium = 'phases';
    case ' Minimize LVs',
        criterium = 'pcs';
end
        
handles.data.criterium = criterium;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popupmenuM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

popupmenuM_Callback(hObject, eventdata, handles)



% --- Executes on selection change in popupmenuPhase.
function popupmenuPhase_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenuPhase contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuPhase

radiobutton1_Callback(hObject, eventdata, handles);

    

% --- Executes during object creation, after setting all properties.
function popupmenuPhase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbuttonCM.
function pushbuttonCM_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

edit_vars_Callback(handles.edit_vars, eventdata, handles)
handles=guidata(hObject);
edit_init_Callback(handles.edit_init, eventdata, handles)
handles=guidata(hObject);
edit_fint_Callback(handles.edit_fint, eventdata, handles)
handles=guidata(hObject);
popupmenuCM_Callback(handles.popupmenuCM, eventdata, handles)
handles=guidata(hObject);

data = preprocess3D(handles.data.x(handles.cm.init:handles.cm.fint,handles.cm.vars,:),handles.data.prep);
switch handles.cm.type,
    case {1,2,3}
         map = cov_map(data,[],0.9,handles.cm.type,0.1,1,handles.console);
                                                                
    case {4,5,6}
         model = get(handles.popupmenuMod,'value');
         if model <= length(handles.data.man_mp_group),             
             model = handles.data.man_mp_group{model};
         elseif model  <= length(handles.data.man_mp_group)+length(handles.data.mp_group),
             model = handles.data.mp_group{model-length(handles.data.man_mp_group)};
         else
             model = handles.data.mp_group2{model-length(handles.data.man_mp_group)-length(handles.data.mp_group)};
         end
         [t,e] = evalMP_s(model);
         map = cov_map(t(handles.cm.init:handles.cm.fint,handles.cm.vars,:),data,0.9,handles.cm.type-3,0.1,1,handles.console);                                                                
                                                               
    case {7,8,9}
         model = get(handles.popupmenuMod,'value');
         if model <= length(handles.data.man_mp_group),             
             model = handles.data.man_mp_group{model};
         elseif model  <= length(handles.data.man_mp_group)+length(handles.data.mp_group),
             model = handles.data.mp_group{model-length(handles.data.man_mp_group)};
         else
             model = handles.data.mp_group2{model-length(handles.data.man_mp_group)-length(handles.data.mp_group)};
         end
         [t,e] = evalMP_s(model);
         map = cov_map(e(handles.cm.init:handles.cm.fint,handles.cm.vars,:),data,0.9,handles.cm.type-6,0.1,1,handles.console);                                                                

end

% --- Executes on selection change in popupmenuCM.
function popupmenuCM_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuCM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenuCM contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuCM

type=get(hObject,'Value');
 
handles.cm.type = type;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popupmenuCM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuCM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popupmenuMod.
function popupmenuMod_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenuMod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        popupmenuMod

model = get(hObject,'Value');
if model <= length(handles.data.man_mp_group),             
    model = handles.data.man_mp_group{model};
elseif model  <= length(handles.data.man_mp_group)+length(handles.data.mp_group),
     model = handles.data.mp_group{model-length(handles.data.man_mp_group)};
else
    model = handles.data.mp_group2{model-length(handles.data.man_mp_group)-length(handles.data.mp_group)};
end

set(handles.popupmenuPhase,'String','');
phases=model.phases;
for i=1:size(phases,1),
    contents = get(handles.popupmenuPhase,'String');
    set(handles.popupmenuPhase,'String',strvcat(contents,sprintf(' Phase %d: from %d to %d',i,phases(i,4),phases(i,5))));
end



% --- Executes during object creation, after setting all properties.
function popupmenuMod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String',' ');


% --- Executes on button press in pushbuttonMod.
function pushbuttonMod_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = get(handles.popupmenuMod,'String');
model = get(handles.popupmenuMod,'Value');
line = deblank(contents(model,:));

cprintMV(handles.console,sprintf(cat(2,'Displaying ',line,': ')),' ',0);

if model <= length(handles.data.man_mp_group),          
    model = handles.data.man_mp_group{model};
    cprintMV(handles.console,sprintf('PRESS = %g',model.cumpress),' ',2);
elseif model <= length(handles.data.man_mp_group)+length(handles.data.mp_group),
    model = handles.data.mp_group{model-length(handles.data.man_mp_group)};
    cprintMV(handles.console,sprintf('PRESS = %g',model.cumpress),' ',2);
    cprintMV(handles.console,' ');
    cprintMV(handles.console,sprintf('LMVs = %d, T = %g, k = %d, minL = %d, a = %d, %s, %s.',model.arg.lag,model.arg.T,model.arg.gamma,model.arg.minsize,model.arg.n,handles.data.preptxt,handles.data.cross.leave_m));                       
else
    model = handles.data.mp_group2{model-length(handles.data.man_mp_group)-length(handles.data.mp_group)};
    cprintMV(handles.console,sprintf('PRESS = %g',model.cumpress),' ',2);
    cprintMV(handles.console,' ');
    cprintMV(handles.console,sprintf('Tm = %g, %s, %s',model.arg.Tm,handles.data.preptxt,handles.data.cross.leave_m)); 
end
cprintMV(handles.console,' ');

sp = size(model.phases);
for i=1:sp(1),
    cprintMV(handles.console,sprintf('From %d to %d, %d PC(s) and %d LMV(s), PRESS = %g',model.phases(i,4),model.phases(i,5),model.phases(i,2),model.phases(i,3),model.phases(i,1)));
end

% --- Executes on button press in pushbuttonDyn.
function pushbuttonDyn_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDyn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data=handles.data;
if get(handles.radiobutton1,'Value'),    
   
    model = get(handles.popupmenuMod,'value');
    if model <= length(handles.data.man_mp_group),             
        model = handles.data.man_mp_group{model};
    elseif model  <= length(handles.data.man_mp_group)+length(handles.data.mp_group),
        model = handles.data.mp_group{model-length(handles.data.man_mp_group)};
    else
        model = handles.data.mp_group2{model-length(handles.data.man_mp_group)-length(handles.data.mp_group)};
    end
    
    phase = get(handles.popupmenuPhase,'value');
    phase = model.phases(phase,:);
    
    if phase(4)>phase(5),
        errordlg('Initial sampling time is posterior to final sampling time.');
    else
        lags = min(10,floor((phase(5)-phase(4))/3));
        if lags > 0,
            x = preprocess3D(data.x,data.prep);
            mv_parreg(x(phase(4):phase(5),:,:),lags,get(handles.radiobutton5,'Value')+1,0.1,1,get(handles.radiobutton3,'Value'),handles.console);
        else
            errordlg('Insufficient phase length for time-series analysis.')
        end
    end
else
    init =  str2num(get(handles.edit_init2,'String'));
    handles.mv.init = init;
    val = get(handles.edit_fint2,'String');
    
    if length(val)>2 && isequal(val(end-2:end),'end'),
        val = [val(1:end-3) num2str(size(handles.data.x,1))];
    end
    
    handles.mv.fint = str2num(val);
    
    mv=handles.mv;

    if mv.init>mv.fint,
        errordlg('Initial sampling time is posterior to final sampling time.')
    else
        lags = min(10,floor((mv.fint-mv.init)/3));
        if lags > 0,
            x = preprocess3D(data.x,data.prep);
            mv_parreg(x(mv.init:mv.fint,:,:),lags,get(handles.radiobutton5,'Value')+1,0.1,1,get(handles.radiobutton3,'Value'),handles.console); 
        else
            errordlg('Insufficient phase length for time-series analysis.')
        end
    end
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function WorkspaceMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to WorkspaceMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

assignin('base','dataMPF',handles.data);


% --------------------------------------------------------------------
function SaveMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, pathname] = uiputfile('*.mat','Save Analysis');
if ~isequal(file, 0)
    dataMPF=handles.data;
    eval(['save ' fullfile(pathname, file) ' dataMPF']);
end


% --------------------------------------------------------------------
function file = LoadMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, pathname] = uigetfile('*.mat');
if ~isequal(file, 0)
    eval(['load ' fullfile(pathname, file) ' dataMPF']);

    handles.data = dataMPF;

    set(handles.edit_LMVs,'String',num2str(handles.data.lags,'%d '));
    set(handles.edit_a,'String',num2str(handles.data.as,'%d '));
    set(handles.edit_k,'String',num2str(handles.data.ks,'%d '));
    set(handles.edit_T,'String',num2str(handles.data.Ts,'%g '));
    set(handles.edit_minL,'String',num2str(handles.data.minLs,'%d '));
    set(handles.edit_Tm,'String',num2str(handles.data.Tms,'%d '));

    set(handles.NPMenuItem,'Checked','off')
    set(handles.TCMenuItem,'Checked','off')
    set(handles.TCTSMenuItem,'Checked','off')
    set(handles.TCVSMenuItem,'Checked','off')
    set(handles.VCMenuItem,'Checked','off')
    set(handles.VCVSMenuItem,'Checked','off')
    switch lower(handles.data.prep),
        case 0,
            set(handles.NPMenuItem,'Checked','on')
        case 1,
            set(handles.TCMenuItem,'Checked','on')
        case 2,
            set(handles.TCTSMenuItem,'Checked','on')
        case 3,
            set(handles.TCVSMenuItem,'Checked','on')
        case 4,
            set(handles.VCMenuItem,'Checked','on')
        case 5,
            set(handles.VCVSMenuItem,'Checked','on')
    end

    set(handles.LnBOMenuItem,'Checked','off')
    set(handles.LnSOMenuItem,'Checked','off')
    set(handles.LnSO2MenuItem,'Checked','off')
    set(handles.CCLnSOMenuItem,'Checked','off')
    switch lower(handles.data.cross.leave_m),
        case 'rkf',
            set(handles.LnBOMenuItem,'Checked','on')
        case 'skf',
            set(handles.LnSOMenuItem,'Checked','on')
        case 'iskf',
            set(handles.LnSO2MenuItem,'Checked','on')
        case 'cskf',
            set(handles.CCLnSOMenuItem,'Checked','on')
    end

    if  length(handles.data.mp_group)<1 && length(handles.data.man_mp_group)<1,
        set(handles.textMod,'Enable','off');
        set(handles.popupmenuPhase,'Enable','off');
        set(handles.textPhase,'Enable','off'); 
        set(handles.pushbuttonMod,'Enable','off');
        set(handles.pushbuttonApply,'Enable','off');
        set(handles.popupmenuMod,'Enable','off');
        set(handles.popupmenuMod,'String',' ');
        set(handles.popupmenuMod,'Value',1);
        set(handles.text_Tm,'Enable','off');
        set(handles.edit_Tm,'Enable','off');
        set(handles.popupmenuM,'Enable','off');
        set(handles.pushbutton_SS,'Enable','off');
        set(handles.radiobutton1,'Enable','off');

        contents = get(handles.popupmenuCM,'String');
        contents = strvcat(contents);
        set(handles.popupmenuCM,'String',contents(1:2,:));
    else
        set(handles.edit_Tm,'String',num2str(handles.data.Tms,'%g '));
        set(handles.pushbuttonMod,'Enable','on');
        set(handles.pushbuttonApply,'Enable','on');
        set(handles.popupmenuMod,'Enable','on');
        set(handles.textMod,'Enable','on');
        set(handles.popupmenuPhase,'Enable','on');
        set(handles.textPhase,'Enable','on');
        set(handles.radiobutton1,'Enable','on');
        
        contents = get(handles.popupmenuCM,'String');
        contents = strvcat(contents(1:3,:));
        set(handles.popupmenuCM,'String',strvcat(contents,' Model Total Covariance'));
        contents = get(handles.popupmenuCM,'String');
        set(handles.popupmenuCM,'String',strvcat(contents,' Model Dynamic Partial Covariance'));
        contents = get(handles.popupmenuCM,'String');
        set(handles.popupmenuCM,'String',strvcat(contents,' Model Partial Covariance'));
        contents = get(handles.popupmenuCM,'String');
        set(handles.popupmenuCM,'String',strvcat(contents,' Residuals Total Covariance'));
        contents = get(handles.popupmenuCM,'String');
        set(handles.popupmenuCM,'String',strvcat(contents,' Residuals Dynamic Partial Covariance'));
        contents = get(handles.popupmenuCM,'String');
        set(handles.popupmenuCM,'String',strvcat(contents,' Residuals Partial Covariance'));

        set(handles.popupmenuMod,'String','');
        for i=1:length(handles.data.man_mp_group),
            contents = get(handles.popupmenuMod,'String');
            set(handles.popupmenuMod,'String',strvcat(contents,sprintf(' Manual MPPCA Model %d',i)));
        end
        
        for i=1:length(handles.data.mp_group),
            contents = get(handles.popupmenuMod,'String');
            set(handles.popupmenuMod,'String',strvcat(contents,sprintf(' MPPCA Model %d',i)));
        end

        set(handles.popupmenuPhase,'Value',1);
        set(handles.popupmenuPhase,'String','');
        if length(handles.data.man_mp_group)<1,
            phases=handles.data.mp_group{1}.phases;
        else
            phases=handles.data.man_mp_group{1}.phases;
        end
        for i=1:size(phases,1),
            contents = get(handles.popupmenuPhase,'String');
            set(handles.popupmenuPhase,'String',strvcat(contents,sprintf(' Phase %d: from %d to %d, %d PC(s) and %d LMV(s)',i,phases(i,4),phases(i,5),phases(i,2),phases(i,3))));
        end
   
        set(handles.text_Tm,'Enable','on');
        set(handles.edit_Tm,'Enable','on');
        set(handles.popupmenuM,'Enable','on');
        set(handles.pushbutton_SS,'Enable','on');

        if  length(handles.data.mp_group2)>0,
            for i=1:length(handles.data.mp_group2),
                contents = get(handles.popupmenuMod,'String');
                set(handles.popupmenuMod,'String',strvcat(contents,sprintf(' Merged MPPCA Model %d',i)));
            end

            if length(handles.data.mp_group2)>2,
                [P,TABLE,STATS] = anova2(handles.data.anova,1,'off');
                h=figure;
                multcompare(STATS,'ctype','lsd');
            end
        end
    end
    guidata(hObject,handles);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

file = uiputfile('*.txt','Print Analysis Information');
if ~isequal(file, 0)
    f = fopen(file,'wt');
    for i=1:size(handles.data.text,1),
        fprintf(f,'%s \n',handles.data.text(i,:));
    end
    fclose(f);
end
open(file)



% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1

set(handles.radiobutton2,'Value',0)
set(handles.radiobutton1,'Value',1);


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2

if isequal('on',get(handles.radiobutton1,'Enable')),
    set(handles.radiobutton1,'Value',0)
end

set(handles.radiobutton2,'Value',1);


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3

set(handles.radiobutton4,'Value',0)
set(handles.radiobutton3,'Value',1);

% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4

set(handles.radiobutton3,'Value',0)
set(handles.radiobutton4,'Value',1);


% --------------------------------------------------------------------
function PreprocessingMenu_Callback(hObject, eventdata, handles)
% hObject    handle to PreprocessingMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function NPMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to NPMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.prep = 0;
handles.data.preptxt = 'NP';
guidata(hObject,handles);

set(handles.NPMenuItem,'Checked','on')
set(handles.TCMenuItem,'Checked','off')
set(handles.TCTSMenuItem,'Checked','off')
set(handles.TCVSMenuItem,'Checked','off')
set(handles.VCMenuItem,'Checked','off')
set(handles.VCVSMenuItem,'Checked','off')


% --------------------------------------------------------------------
function TCMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to TCMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.prep = 1;
handles.data.preptxt = 'TC';
guidata(hObject,handles);

set(handles.NPMenuItem,'Checked','off')
set(handles.TCMenuItem,'Checked','on')
set(handles.TCTSMenuItem,'Checked','off')
set(handles.TCVSMenuItem,'Checked','off')
set(handles.VCMenuItem,'Checked','off')
set(handles.VCVSMenuItem,'Checked','off')


% --------------------------------------------------------------------
function TCTSMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to TCTSMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.prep = 2;
handles.data.preptxt = 'TCTS';
guidata(hObject,handles);

set(handles.NPMenuItem,'Checked','off')
set(handles.TCMenuItem,'Checked','off')
set(handles.TCTSMenuItem,'Checked','on')
set(handles.TCVSMenuItem,'Checked','off')
set(handles.VCMenuItem,'Checked','off')
set(handles.VCVSMenuItem,'Checked','off')


% --------------------------------------------------------------------
function TCVSMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to TCVSMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.prep = 3;
handles.data.preptxt = 'TCVS';
guidata(hObject,handles);

set(handles.NPMenuItem,'Checked','off')
set(handles.TCMenuItem,'Checked','off')
set(handles.TCTSMenuItem,'Checked','off')
set(handles.TCVSMenuItem,'Checked','on')
set(handles.VCMenuItem,'Checked','off')
set(handles.VCVSMenuItem,'Checked','off')


% --------------------------------------------------------------------
function VCMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to VCMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.prep = 4;
handles.data.preptxt = 'VC';
guidata(hObject,handles);

set(handles.NPMenuItem,'Checked','off')
set(handles.TCMenuItem,'Checked','off')
set(handles.TCTSMenuItem,'Checked','off')
set(handles.TCVSMenuItem,'Checked','off')
set(handles.VCMenuItem,'Checked','on')
set(handles.VCVSMenuItem,'Checked','off')


% --------------------------------------------------------------------
function VCVSMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to VCVSMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.prep = 5;
handles.data.preptxt = 'VCVs';
guidata(hObject,handles);

set(handles.NPMenuItem,'Checked','off')
set(handles.TCMenuItem,'Checked','off')
set(handles.TCTSMenuItem,'Checked','off')
set(handles.TCVSMenuItem,'Checked','off')
set(handles.VCMenuItem,'Checked','off')
set(handles.VCVSMenuItem,'Checked','on')








% --------------------------------------------------------------------
function CVMenu_Callback(hObject, eventdata, handles)
% hObject    handle to CVMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % --- Executes on selection change in popupmenuC.
% function popupmenuC_Callback(hObject, eventdata, handles)
% % hObject    handle to popupmenuC (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: contents = get(hObject,'String') returns popupmenuC contents as cell array
% %        contents{get(hObject,'Value')} returns selected item from popupmenuC
% contents = get(hObject,'String');
% txt=contents{get(hObject,'Value')};
%  
% handles.data.cross.leave_m = txt(2:end);
% guidata(hObject,handles);
% 
% 
% % --- Executes during object creation, after setting all properties.
% function popupmenuC_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to popupmenuC (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: popupmenu controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% 
% set(hObject,'Value',2);
% popupmenuC_Callback(hObject, eventdata, handles)



% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5





function edit_man_init_Callback(hObject, eventdata, handles)
% hObject    handle to edit_man_init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_man_init as text
%        str2double(get(hObject,'String')) returns contents of edit_man_init as a double

init =  str2num(get(hObject,'String'));
handles.man.init = init;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_man_init_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_man_init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

edit_man_init_Callback(hObject, eventdata, handles)


function fint = edit_man_fint_Callback(hObject, eventdata, handles)
% hObject    handle to edit_man_fint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_man_fint as text
%        str2double(get(hObject,'String')) returns contents of edit_man_fint as a double

val = get(hObject,'String');
if length(val)>2 && isequal(val(end-2:end),'end'),
    val = [val(1:end-3) num2str(size(handles.data.x,1))];
end

handles.man.fint = str2num(val);
guidata(hObject,handles);

fint = handles.man.fint;

% --- Executes during object creation, after setting all properties.
function edit_man_fint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_man_fint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_man_LMVs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_man_LMVs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_man_LMVs as text
%        str2double(get(hObject,'String')) returns contents of edit_man_LMVs as a double

lmvs =  str2num(get(hObject,'String'));
handles.man.lmvs = lmvs;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_man_LMVs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_man_LMVs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

edit_man_LMVs_Callback(hObject, eventdata, handles)

    

function edit_man_PCs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_man_PCs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_man_PCs as text
%        str2double(get(hObject,'String')) returns contents of edit_man_PCs as a double

pcs =  str2num(get(hObject,'String'));
handles.man.pcs = pcs;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_man_PCs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_man_PCs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

edit_man_PCs_Callback(hObject, eventdata, handles)


% --- Executes on button press in pushbutton_man.
function pushbutton_man_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_man (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.man.fint = edit_man_fint_Callback(handles.edit_man_fint, eventdata, handles);

% Checking arguments
ok=true;
if find(isnan(handles.man.pcs)) | find(handles.man.pcs<1), 
    cprintMV(handles.console,'Incorrect number of PCs.'); ok=false; end;
if find(isnan(handles.man.lmvs)) | find(handles.man.lmvs<0), 
    cprintMV(handles.console,'Incorrect number of LMVs.'); ok=false; end;
if find(isnan(handles.man.init)) | find(handles.man.init<0) | find(handles.man.init>size(handles.data.x,1)), 
    cprintMV(handles.console,'Incorrect initial time point.'); ok=false; end;
if find(isnan(handles.man.fint)) | find(handles.man.fint<0) | find(handles.man.fint>size(handles.data.x,1)), 
    cprintMV(handles.console,'Incorrect final time point.'); ok=false; end;
n_phases = length(handles.man.pcs);
if n_phases ~= length(handles.man.lmvs), 
    cprintMV(handles.console,'Non consistent number of parameters.'); ok=false; end;
if n_phases ~= length(handles.man.init),
    cprintMV(handles.console,'Non consistent number of parameters.'); ok=false; end;
if n_phases ~= length(handles.man.fint),
    cprintMV(handles.console,'Non consistent number of parameters.'); ok=false; end;
for i=1:n_phases,
    if handles.man.init(i)>handles.man.fint(i),
    cprintMV(handles.console,sprintf('Incorrect initial/final time in phase %d.',i)); ok=false; end;
    if isinf(handles.man.lmvs(i)),
        handles.man.lmvs(i)=(handles.man.fint(i)-handles.man.init(i)); end;
    if handles.man.lmvs(i)>(handles.man.fint(i)-handles.man.init(i)),
    cprintMV(handles.console,sprintf('Incorrect number of LMVs in phase %d.',i)); ok=false; end;
    if handles.man.pcs(i)>rank(unfold(handles.data.x(handles.man.init(i):handles.man.fint(i),:,:),handles.man.lmvs(i))),
    cprintMV(handles.console,sprintf('Number of PCs in phase %d above the rank.',i)); ok=false; end;
end
    
if ok,
    pcs=ones(n_phases,1);
    pcs(:)=handles.man.pcs;
    lmvs=ones(n_phases,1);
    lmvs(:)=handles.man.lmvs;
    init=ones(n_phases,1);
    init(:)=handles.man.init;
    fint=ones(n_phases,1);
    fint(:)=handles.man.fint;
    
    d=init(2:end)-fint(1:end-1);
    if ~isempty(d) && ~isempty(find(d~=1)), ok=false; end;
    d=init(2:end)-lmvs(2:end);
    if ~isempty(find(d<1)), ok=false; end;
end

%Main code
if ok,
    arg=struct('xini',handles.data.x,'cross',handles.data.cross,'prep',handles.data.prep); 
    mp_model=struct('type','Manual','arg',arg);
    
    mp_model.phases = [ones(n_phases,1),pcs,lmvs,init,fint];

    mp_model = crossvalMP_s(mp_model,handles.data.cross);

    handles.data.man_mp_group{length(handles.data.man_mp_group)+1} = mp_model;
    
    contents = get(handles.popupmenuMod,'String');
    if isequal(contents,' '), contents=''; end;
    set(handles.popupmenuMod,'Value',1);
    set(handles.popupmenuMod,'String',strvcat(contents(1:length(handles.data.man_mp_group)-1,:),sprintf(' Manual MPPCA Model %d',length(handles.data.man_mp_group)),contents(length(handles.data.man_mp_group):end,:)));

    popupmenuMod_Callback(handles.popupmenuMod, eventdata, handles)

    
    set(handles.pushbuttonMod,'Enable','on');
    set(handles.pushbuttonApply,'Enable','on');
    set(handles.popupmenuMod,'Enable','on');
    set(handles.textMod,'Enable','on');
    set(handles.popupmenuPhase,'Enable','on');
    set(handles.textPhase,'Enable','on');
    set(handles.radiobutton1,'Enable','on');
 
    contents = get(handles.popupmenuMod,'String');
    line = deblank(contents(length(handles.data.man_mp_group),:));
    cprintMV(handles.console,sprintf(cat(2,'Displaying ',line,': ')),' ',0);
    cprintMV(handles.console,sprintf('PRESS = %g',mp_model.cumpress),' ',2);
    cprintMV(handles.console,' ');
    
    sp = size(mp_model.phases);
    for i=1:sp(1),
        cprintMV(handles.console,sprintf('From %d to %d, %d PC(s) and %d LMV(s), PRESS = %g',mp_model.phases(i,4),mp_model.phases(i,5),mp_model.phases(i,2),mp_model.phases(i,3),mp_model.phases(i,1)));
    end

    contents = get(handles.popupmenuCM,'String');
    if size(contents,1)<=3,
        contents = strvcat(contents);
        set(handles.popupmenuCM,'String',strvcat(contents,' Model Total Covariance'));
        contents = get(handles.popupmenuCM,'String');
        set(handles.popupmenuCM,'String',strvcat(contents,' Model Dynamic Partial Covariance'));
        contents = get(handles.popupmenuCM,'String');
        set(handles.popupmenuCM,'String',strvcat(contents,' Model Partial Covariance'));
        contents = get(handles.popupmenuCM,'String');
        set(handles.popupmenuCM,'String',strvcat(contents,' Residuals Total Covariance'));
        contents = get(handles.popupmenuCM,'String');
        set(handles.popupmenuCM,'String',strvcat(contents,' Residuals Dynamic Partial Covariance'));
        contents = get(handles.popupmenuCM,'String');
        set(handles.popupmenuCM,'String',strvcat(contents,' Residuals Partial Covariance'));
    end

    guidata(hObject,handles);
else
    cprintMV(handles.console,'Check the input for errors.');  
end





% --- Executes on button press in pushbutton21.
function pushbuttonPCs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.man.fint = edit_man_fint_Callback(handles.edit_man_fint, eventdata, handles);

% Checking arguments
ok=true;
if find(isnan(handles.man.pcs)) | find(handles.man.pcs<1), 
    cprintMV(handles.console,'Incorrect number of PCs.'); ok=false; end;
if find(isnan(handles.man.lmvs)) | find(handles.man.lmvs<0), 
    cprintMV(handles.console,'Incorrect number of LMVs.'); ok=false; end;
if find(isnan(handles.man.init)) | find(handles.man.init<0) | find(handles.man.init>size(handles.data.x,1)), 
    cprintMV(handles.console,'Incorrect initial time point.'); ok=false; end;
if find(isnan(handles.man.fint)) | find(handles.man.fint<0) | find(handles.man.fint>size(handles.data.x,1)), 
    cprintMV(handles.console,'Incorrect final time point.'); ok=false; end;
n_phases = length(handles.man.pcs);
if n_phases ~= length(handles.man.lmvs), 
    cprintMV(handles.console,'Non consistent number of parameters.'); ok=false; end;
if n_phases ~= length(handles.man.init),
    cprintMV(handles.console,'Non consistent number of parameters.'); ok=false; end;
if n_phases ~= length(handles.man.fint),
    cprintMV(handles.console,'Non consistent number of parameters.'); ok=false; end;
for i=1:n_phases,
    if handles.man.init(i)>handles.man.fint(i),
    cprintMV(handles.console,sprintf('Incorrect initial/final time in phase %d.',i)); ok=false; end;
    if isinf(handles.man.lmvs(i)),
        handles.man.lmvs(i)=(handles.man.fint(i)-handles.man.init(i)); end;
    if handles.man.lmvs(i)>(handles.man.fint(i)-handles.man.init(i)),
    cprintMV(handles.console,sprintf('Incorrect number of LMVs in phase %d.',i)); ok=false; end;
    if handles.man.pcs(i)>rank(unfold(handles.data.x(handles.man.init(i):handles.man.fint(i),:,:),handles.man.lmvs(i))),
    cprintMV(handles.console,sprintf('Number of PCs in phase %d above the rank.',i)); ok=false; end;
end

if ok,
    pcs=handles.man.pcs;
    lmvs=handles.man.lmvs;
    init=handles.man.init;
    fint=handles.man.fint;
    
    d=init(2:end)-fint(1:end-1);
    if ~isempty(d) && ~isempty(find(d~=1)), ok=false; end;
    d=init(2:end)-lmvs(2:end);
    if ~isempty(find(d<1)), ok=false; end;
end

%Main code
if ok,
    arg=struct('xini',handles.data.x,'cross',handles.data.cross,'prep',handles.data.prep); 
    
    [ccs,av,st] = preprocess3D(handles.data.x,handles.data.prep);
    
    txt = cprintMV(handles.console,'Computing, please wait...',[],0);
    for i=1:length(pcs),
        c_2D = unfold(ccs(init(i):fint(i),:,:),lmvs(i));
        % ----- MEDA Toolbox ----- %
        [x_var,cumpress] = var_pca(c_2D,[],0,0); 
        plot_vec([x_var cumpress/cumpress(1)],0:length(x_var)-1,[],{'#PCs',sprintf('Phase %d',i)},[],0,{'% Res. Var','ckf PRESS'});
        legend('show');
        % ----- MEDA Toolbox ----- %
    end
    
    cprintMV(handles.console,'Done.',txt,1);

else
    cprintMV(handles.console,'Check the input for errors.');
end


% --- Executes on button press in pushbuttonExp.
function pushbuttonExp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonExp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data=handles.data;
x = preprocess3D(data.x,data.prep);
if get(handles.radiobutton1,'Value'),    
  
    model = get(handles.popupmenuMod,'value');
    if model <= length(handles.data.man_mp_group),             
        model = handles.data.man_mp_group{model};
    elseif model  <= length(handles.data.man_mp_group)+length(handles.data.mp_group),
        model = handles.data.mp_group{model-length(handles.data.man_mp_group)};
    else
        model = handles.data.mp_group2{model-length(handles.data.man_mp_group)-length(handles.data.mp_group)};
    end
    
    phase = get(handles.popupmenuPhase,'value');
    phase = model.phases(phase,:);
    
    if phase(4)>phase(5),
        errordlg('Initial sampling time is posterior to final sampling time.');
    end
    
    xu = unfold(x(phase(4):phase(5),:,:),phase(3));
    if handles.PCAh~=0 & ishandle(handles.PCAh), 
        close(handles.PCAh); 
    end
    handles.PCAh = PCA(xu,1:phase(2),0);    
else
    
    init =  str2num(get(handles.edit_init2,'String'));
    handles.mv.init = init;
    val = get(handles.edit_fint2,'String');
    
    if length(val)>2 && isequal(val(end-2:end),'end'),
        val = [val(1:end-3) num2str(size(handles.data.x,1))];
    end
    
    handles.mv.fint = str2num(val);
    mv=handles.mv;
    
    if mv.init>mv.fint,
        errordlg('Initial sampling time is posterior to final sampling time.');
    end
    xu = unfold(x(mv.init:mv.fint,:,:),0);
    if handles.PCAh~=0 & ishandle(handles.PCAh), 
        close(handles.PCAh); 
    end
    handles.PCAh = PCA(xu,[],0);  
end

guidata(hObject,handles);



% --------------------------------------------------------------------
function LnBOMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to LnBOMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.cross.leave_m = 'rkf';
guidata(hObject,handles);

set(handles.LnBOMenuItem,'Checked','on')
set(handles.LnSOMenuItem,'Checked','off')
set(handles.ckfMenuItem,'Checked','off')

% --------------------------------------------------------------------
function LnSOMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to LnSOMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.cross.leave_m = 'ekf';
guidata(hObject,handles);

set(handles.LnBOMenuItem,'Checked','off')
set(handles.LnSOMenuItem,'Checked','on')
set(handles.ckfMenuItem,'Checked','off')

% --------------------------------------------------------------------
function ckfMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to ckfMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.cross.leave_m = 'ckf';
guidata(hObject,handles);

set(handles.LnBOMenuItem,'Checked','off')
set(handles.LnSOMenuItem,'Checked','off')
set(handles.ckfMenuItem,'Checked','on')


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over console.
function console_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to console (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonClose.
function pushbuttonClose_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% axes(handles.ParentFigure.main_window)
% image(handles.ParentFigure.images{5});
% axis off;
% axis image;
% handles.ParentFigure.track(4) = 1;
% handles.ParentFigure.track(5) = 0;
% guidata(handles.ParentsWindow, handles.ParentFigure);
delete(handles.figure1)

% --- Executes on button press in pushbuttonApply.
function pushbuttonApply_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonApply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

models.man_mp_group = handles.data.man_mp_group;
models.mp_group = handles.data.mp_group;
models.mp_group2 = handles.data.mp_group2;
 
handles.ParentFigure.s_calibration = models;
handles.ParentFigure.s_calibration.x = handles.data.x;
guidata(handles.ParentsWindow,handles.ParentFigure)

% Update the main layout by enabling the monitoring GUI
axes(handles.ParentFigure.main_window)
image(handles.ParentFigure.images{5});
axis off;
axis image;
% Update the tracking array of the bilinear modeling cycling (only updated
% when the user presses the button "apply"
handles.ParentFigure.track(:) = 1;
handles.ParentFigure.track(5) = 0;

% Update handles structure
guidata(handles.ParentsWindow, handles.ParentFigure);
guidata(hObject, handles);

% Allow user to proceed with the design of the monitoring scheme
set(handles.ParentFigure.pbMonitoring,'Enable','on');

% Delete the current user interface
delete(handles.figure1);
