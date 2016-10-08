function varargout = Monitoring(varargin)
% MONITORING M-file for Monitoring.fig
%      MONITORING, by itself, creates a new MONITORING or raises the existing
%      singleton*.
%
%      H = MONITORING returns the handle to a new MONITORING or the handle to
%      the existing singleton*.
%
%      MONITORING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MONITORING.M with the given input arguments.
%
%      MONITORING('Property','Value',...) creates a new MONITORING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before monitor_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Monitoring_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Monitoring

% Last Modified by GUIDE v2.5 25-Sep-2014 15:34:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Monitoring_OpeningFcn, ...
                   'gui_OutputFcn',  @Monitoring_OutputFcn, ...
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


% --- Executes just before Monitoring is made visible.
function Monitoring_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Monitoring (see VARARGIN)

% Choose default command line output for Monitoring
handles.output = hObject;

if length(varargin)<1, error('Error in the number of arguments.'); end;

handles.ParentsWindow=varargin{1};
handles.ParentFigure = guidata(handles.ParentsWindow);

handles.calibration = varargin{2};
handles.alignment = varargin{3};

% Update handles structure
guidata(hObject, handles);

set(handles.popupmenuMod,'String','');
for i=1:length(handles.calibration.man_mp_group),
    contents = get(handles.popupmenuMod,'String');
    set(handles.popupmenuMod,'String',strvcat(contents,sprintf(' Manual MPPCA Model %d',i)));
end
        
for i=1:length(handles.calibration.mp_group),
    contents = get(handles.popupmenuMod,'String');
    set(handles.popupmenuMod,'String',strvcat(contents,sprintf(' MPPCA Model %d',i)));
end
     
for i=1:length(handles.calibration.mp_group2),
    contents = get(handles.popupmenuMod,'String');
    set(handles.popupmenuMod,'String',strvcat(contents,sprintf(' Merged MPPCA Model %d',i)));
end

handles.calibration.alph95=0.05;
handles.calibration.alph=0.01;
handles.calibration.alpr95=0.05;
handles.calibration.alpr=0.01;

% Matrix containing the test batch to be projected onto the latent
% structure.
handles.calibration.test_batch = [];
% Matrix containig the test batch to be check in the fault diagnosis step. 
handles.calibration.test_batchFD = 1;
% Flag indicating what mode must be used for fault diagnosis (by default
% off-line, i.e. a value equal to 0.
handles.modeMonitoring = 0;
% Flag indicating what multivariate statistic should be investigated in
% fault diagnosis. By default SPE is selected, i.e. a value equal to 0.
handles.statFaultDiagnosis = 0;
% Initialize the batch sampling point to take into account in the
% contribution estimation to the first one.
handles.batchTime = 0;
% Inititialize the counter of the selected data set. 
handles.selectedDataSet = 1;

set(handles.popupmenuVar,'Enable','off');
set(handles.textVar,'Enable','off');
set(handles.popupmenuBat,'Enable','off');
set(handles.textBat,'Enable','off');
%set(handles.pushbuttonRef,'Enable','off');
set(handles.pushbuttonPlo,'Enable','off');
set(handles.ImportMenuItem,'Enable','off');

popupmenuMod_Callback(handles.popupmenuMod, eventdata, handles);


% UIWAIT makes Monitoring wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Monitoring_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  STRUCTURE
%% 
%%  1.- MODEL PANEL
%%  2.- MONITORING SYSTEM PANEL
%%  3.- FAULT DIAGNOSIS PANEL
%%  4.- MENU ITEMS
%%  5.- AUXILIAR FUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                   1.- MODEL PANEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on selection change in popupmenuMod.
function popupmenuMod_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenuMod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuMod

model = get(hObject,'Value');
if model <= length(handles.calibration.man_mp_group),             
    model = handles.calibration.man_mp_group{model};
elseif model  <= length(handles.calibration.man_mp_group)+length(handles.calibration.mp_group),
     model = handles.calibration.mp_group{model-length(handles.calibration.man_mp_group)};
else
    model = handles.calibration.mp_group2{model-length(handles.calibration.man_mp_group)-length(handles.calibration.mp_group)};
end

handles.calibration.model = model;

set(handles.radiobuttonCV,'Enable','off');
set(handles.radiobuttonCV,'Value',0);

guidata(hObject,handles);
%pushbuttonRef_Callback(handles.pushbuttonRef, eventdata, handles);


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


% --- Executes on button press in pushbuttonCV.
function pushbuttonCV_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(gcf,'pointer','watch'); pause(.001);

try
    [handles.calibration.alph, handles.calibration.alpr, handles.calibration.alph95, handles.calibration.alpr95, handles.calibration.alpoh, handles.calibration.alpor, handles.calibration.alpoh95, handles.calibration.alpor95] = plot_distcv2(handles.calibration.x, handles.calibration.model.phases, 2, 1, handles.axes1, handles.axes2, handles.axesPostBatchT2, handles.axesPostBatchSPE);
catch err
   errordlg(err.message); 
   set(gcf,'pointer','arrow');
end

cla(handles.axes3);
cla(handles.axesPostBatchWI);
set(handles.axes3,'Visible','off');
set(handles.axesPostBatchWI,'Visible','off');

try
    if strcmp(handles.alignment.synchronization{handles.alignment.stages,1}.methodsyn,'dtw') || strcmp(handles.alignment.synchronization{handles.alignment.stages,1}.methodsyn,'multisynchro')
        plot_onwarp(handles.alignment.synchronization{handles.alignment.stages,1}.warp, handles.alignment.synchronization{handles.alignment.stages,1}.band,[],[],[],true,handles.axes3);
         set(handles.axes3,'Visible','on');
    end
catch err
   errordlg(err.message); 
   set(gcf,'pointer','arrow');
end
set(gcf,'pointer','arrow')


% Enable the monitoring option for the calibration data set
handles.calibration.test{1,1} = handles.alignment.synchronization{1,1}.nor_batches';

set(handles.popupmenuVar,'String',' '); 
set(handles.popupmenuVar,'String','calibration'); 
set(handles.textData,'String','calibration'); 

handles.calibration.test_batch = handles.calibration.test{handles.selectedDataSet}{1};
handles.calibration.test_batchFD = 1;

set(handles.popupmenuBat,'String','');
set(handles.popupmenuBatFD,'String','');

% Fill the popmenu with the number of batches
set(handles.popupmenuBat,'String','');
set(handles.popupmenuBatFD,'String','');
for i=1:length(handles.calibration.test{handles.selectedDataSet})
    contents = get(handles.popupmenuBat,'String');
    set(handles.popupmenuBat,'String',strvcat(contents,[' ',num2str(i)]));
    set(handles.popupmenuBatFD,'String',strvcat(contents,[' ',num2str(i)]));
end

set(handles.popupmenuBat,'Value',1);
set(handles.popupmenuBatFD,'Value',1);

% Fill the popmenu with the number of sampling points
set(handles.popupmenuTimePoint,'String','');
for i=1:size(handles.calibration.x,1)
    contents = get(handles.popupmenuTimePoint,'String');
    set(handles.popupmenuTimePoint,'String',strvcat(contents,[' ',num2str(i)]));
end

set(handles.pbContribution,'Enable','on');
set(handles.popupmenuVar,'Enable','on');   
set(handles.textVar,'Enable','on');
set(handles.popupmenuBat,'Enable','on');
set(handles.textBat,'Enable','on');
set(handles.pushbuttonPlo,'Enable','on');
set(handles.ImportMenuItem,'Enable','off');
set(handles.radiobuttonCV,'Enable','on');
set(handles.ImportMenuItem,'Enable','on');

handles.ParentFigure.track(5) = 1;

% Store the values of the CV 
handles.ParentFigure.s_monitoring = handles.calibration;

% Update handles structure
guidata(handles.ParentsWindow, handles.ParentFigure);
guidata(hObject,handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                  2.- MONITORING SYSTEM PANEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in popupmenuVar.
function popupmenuVar_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenuVar contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuVar

vars = get(handles.popupmenuVar,'String');
handles.selectedDataSet = get(handles.popupmenuVar,'Value');
set(handles.textData,'String',vars(handles.selectedDataSet,:));
set(handles.popupmenuBat,'String','');
set(handles.popupmenuBatFD,'String','');

for i=1:length(handles.calibration.test{handles.selectedDataSet})
    contents = get(handles.popupmenuBat,'String');
    set(handles.popupmenuBat,'String',strvcat(contents,[' ',num2str(i)]));
    set(handles.popupmenuBatFD,'String',strvcat(contents,[' ',num2str(i)]));
end

set(handles.popupmenuBat,'Value',1);
handles.calibration.test_batch = handles.calibration.test{handles.selectedDataSet}{1};

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenuVar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenuBat.
function popupmenuBat_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuBat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenuBat contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuBat

indx=get(handles.popupmenuBat,'Value');
handles.calibration.test_batch = handles.calibration.test{handles.selectedDataSet}{indx};

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenuBat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuBat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonPlo.
function pushbuttonPlo_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPlo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[synTestBatch,nsamplesToPlot,warptest] = onlineSynchronization(handles, handles.calibration.test_batch);

vars = cellstr(get(handles.popupmenuVar,'String'));
posTest = length(handles.alignment.x) + 1;
if strcmp(vars{get(handles.popupmenuVar,'Value')},'calibration')
    posTest = get(handles.popupmenuBat,'Value');
end

try 
    if get(handles.radiobuttonCV,'Value')
        plot_onstat(handles.calibration.x, synTestBatch, handles.calibration.model.phases, handles.calibration.model.arg.prep, 1, handles.calibration.alph, handles.calibration.alpr, handles.calibration.alph95, handles.calibration.alpr95,nsamplesToPlot, handles.axes1, handles.axes2, handles.calibration.alpoh, handles.calibration.alpor, handles.calibration.alpoh95, handles.calibration.alpor95, handles.axesPostBatchT2, handles.axesPostBatchSPE, posTest);
    else
        plot_onstat(handles.calibration.x, synTestBatch, handles.calibration.model.phases, handles.calibration.model.arg.prep, 1, 0.01, 0.01, 0.05, 0.05, nsamplesToPlot, handles.axes1, handles.axes2, 0.01, 0.01, 0.05, 0.05, handles.axesPostBatchT2, handles.axesPostBatchSPE, posTest);
    end
catch err
   errordlg(err.message); 
end

if strcmp(handles.alignment.synchronization{handles.alignment.stages}.methodsyn,'dtw') || strcmp(handles.alignment.synchronization{handles.alignment.stages}.methodsyn,'multisynchro')
    try
        plot_onwarp(handles.alignment.synchronization{handles.alignment.stages}.warp, handles.alignment.synchronization{handles.alignment.stages}.band,warptest,[],[],false,handles.axes3);
    catch err
       errordlg(err.message); 
    end
end

set(handles.pbContribution,'Enable','on');

function radiobuttonCV_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonCV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonCV


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                  3.- FAULT DIAGNOSIS PANEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in popupmenuBatFD.
function popupmenuBatFD_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuBatFD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuBatFD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuBatFD

handles.calibration.test_batchFD = get(handles.popupmenuBatFD,'Value');

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenuBatFD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuBatFD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenuTimePoint.
function popupmenuTimePoint_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuTimePoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuTimePoint contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuTimePoint

handles.batchTime = get(handles.popupmenuTimePoint,'Value');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenuTimePoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuTimePoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pbContribution.
function pbContribution_Callback(hObject, eventdata, handles)
% hObject    handle to pbContribution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

test_batch = handles.calibration.test{handles.selectedDataSet}{handles.calibration.test_batchFD};

[test_batchSyn] = onlineSynchronization(handles,test_batch);
    
if ~handles.modeMonitoring && handles.calibration.model.phases(1,3) ~= size(handles.alignment.alg_batches,1)-1, 
    warndlg('Overall contributions cannot be displayed because the model is not a batch-wise PCA. Only the option of contributions at a specific sampling time is possible.');
    return;
end

try 
    plot_contributions(handles.calibration.x, test_batchSyn, handles.calibration.model.phases, handles.modeMonitoring, handles.statFaultDiagnosis, 2, handles.ParentFigure.varNames(2:end,:) ,handles.batchTime);
catch err
   errordlg(err.message); 
end

% --------------------------------------------------------------------
function uipanelMonMode_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanelMonMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes when selected object is changed in uipanelMonMode.
function uipanelMonMode_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelMonMode 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobuttonOffline'
        % Code for when radiobutton5 is selected.
        handles.modeMonitoring = 0; % OFF-LINE
        set(handles.popupmenuTimePoint,'Enable','off');
    case 'radiobuttonOnline'
        % Code for when radiobutton6 is selected.
        handles.modeMonitoring = 1; % ON-LINE
        set(handles.popupmenuTimePoint,'Enable','on');
    % Continue with more cases as necessary.
    otherwise
        errordlg('An error has ocurred when the process monitoring mode has changed');
        
end
guidata(hObject,handles);

% --- Executes when selected object is changed in uipanelStatistics.
function uipanelStatistics_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelStatistics 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobuttonHotelling'
        % Code for when radiobuttonHotelling is selected.
        handles.statFaultDiagnosis = 0;
    case 'radiobuttonSPE'
        % Code for when radiobuttonSPE is selected.
        handles.statFaultDiagnosis = 1;
    otherwise
        % Code for when there is no match.
        errordlg('An error has ocurred when the multivariate statistic has changed');
end
guidata(hObject,handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                    4.- MENU ITEMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------------------------------------------------------
function SaveMenu_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function SaveMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, pathname] = uiputfile('*.mat','Save to File');
if ~isequal(file, 0)
    dataMon=handles.calibration;
    eval(['save ' fullfile(pathname, file) ' dataMon']);
end

% --------------------------------------------------------------------
function SaveWMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to SaveWMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

assignin('base','dataMon',handles.calibration);

% --------------------------------------------------------------------
function ImportMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to ImportMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = [];
eq = 0;
[filename, pathname] = uigetfile('*.mat', 'Select a MATLAB file');

if ~isequal(filename,0)
    try 
       aux = load(strcat(pathname,filename), 'test');
       xtest = aux.test;
    catch 
        % Give more information for mismatch.
        errordlg('An expected problem trying to load the test set from the selected file has occurred.');
    end 
else
    return;
end
   
if isvector(xtest) && size(xtest,1)>0
    x = cell(size(xtest,1),1);
    for j=1:size(xtest,1)
        if size(xtest(j).data,1) ~=  size(handles.ParentFigure.s_screening.batch_data(1,1).data,1)
            errordlg('The batches of the test set do not have the same frequency sampling as the calibration data set');
            return;
        end
     
        for z=1:size(xtest(j).data,2)
            if size(xtest(j).data{z},2)~=size(handles.ParentFigure.s_screening.batch_data(1,1).data{z},2)
                errordlg('The test data set does not contain the same number of variables with common sampling frequency as the calibration data set.');
                return;
            end
        end
        
        if size(xtest(j).data,1) == 1
            varsIn = find(handles.ParentFigure.dataset.VariablesIn==1); varsIn= varsIn(2:end);
            x{j} = xtest(j).data{1,1}(:,varsIn+ones(size(varsIn)).*1);
        else
            eq=1;
        end
    end
else
   errordlg('Recall that the variable trajectories must be contained in a Matlab column-wise cell array.'); 
end

if eq
    x = arrange2D(xtest,handles.ParentFigure.s_alignment.equalization.inter,handles.ParentFigure.s_alignment.equalization.units,handles.ParentFigure.s_alignment.equalization.method_interp);
    varsIn = find(handles.ParentFigure.dataset.VariablesIn==1);
    for i=1:length(x),
        x{i} =  x(:,varsIn+ones(size(varsIn)).*1);
    end
end
   
handles.calibration.test = [handles.calibration.test; {x}];

handles.calibration.test_batch = x{1};

handles.calibration.test_batchFD = 1;

contents = strvcat(get(handles.popupmenuVar,'String'),filename); 
set(handles.popupmenuVar,'String',contents); 
set(handles.popupmenuVar,'Value',size(get(handles.popupmenuVar,'String'),1));
set(handles.textData,'String',contents(size(contents,1),:));

handles.selectedDataSet = size(contents,1);

set(handles.popupmenuBat,'String','');
for i=1:length(x)
    contents = get(handles.popupmenuBat,'String');
    set(handles.popupmenuBat,'String',strvcat(contents,[' ',num2str(i)]));
    set(handles.popupmenuBatFD,'String',strvcat(contents,[' ',num2str(i)]));
end
set(handles.popupmenuBat,'Value',1);
set(handles.popupmenuBatFD,'Value',1);

  
set(handles.popupmenuVar,'Enable','on'); 
set(handles.textVar,'Enable','on');
set(handles.popupmenuBatFD,'Enable','on');
set(handles.textData,'Enable','on');
set(handles.popupmenuBat,'Enable','on');
set(handles.textBat,'Enable','on');
set(handles.pushbuttonPlo,'Enable','on');

guidata(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                    5.- AUXILIAR FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [synTestBatch,nsamplesToPlot,warptest] = onlineSynchronization(handles, test_batch)

warptest = [];

switch handles.alignment.synchronization{handles.alignment.stages}.methodsyn
     case 'iv'
        uBn =  handles.calibration.test_batch;
        if handles.alignment.synchronization{handles.alignment.stages,1}.cut
            indm = find(handles.alignment.synchronization{handles.alignment.stages,1}.max_ep>=handles.calibration.test_batch(:,handles.alignment.synchronization{handles.alignment.stages,1}.var),1);
            indM = find(handles.alignment.synchronization{handles.alignment.stages,1}.min_ep<=handles.calibration.test_batch(:,handles.alignment.synchronization{handles.alignment.stages,1}.var),1);
            if isempty(indm) || isempty(indM)
                errordlg('The batch selected cannot be synchronized because the IV has not the same initial and/or end point than the NOC IV.','Error Dialog','modal');
                return;
            end
            uBn =  test_batch(min(indm,indM):max(indm,indM),2:end);
        end

             if length(find(isnan(uBn(:,handles.alignment.synchronization{handles.alignment.stages,1}.var)))) <= 0.25*length(uBn(:,handles.alignment.synchronization{handles.alignment.stages,1}.var)) % Aling every batch where the iv was measured
                synTestBatch = align_IV(uBn,handles.alignment.synchronization{handles.alignment.stages,1}.var,handles.alignment.synchronization{handles.alignment.stages,1}.steps,handles.alignment.synchronization{handles.alignment.stages,1}.method);
                nsamplesToPlot = size(synTestBatch,1);
             else
                 errordlg('Too much missing data in the indicator variable.','Error Dialog','modal');
             end
         
     case 'dtw'
        Bref = scale_(handles.alignment.synchronization{handles.alignment.stages,1}.Xref,handles.alignment.synchronization{handles.alignment.stages,1}.rng);
        test = scale_(test_batch(:,2:end),handles.alignment.synchronization{handles.alignment.stages,1}.rng);
        [synTestBatch, warptest] = DTW(test,Bref,diag(handles.alignment.synchronization{handles.alignment.stages}.W));      
        
        for j=1:size(handles.alignment.synchronization{handles.alignment.stages}.nor_batches{1},2)
            synTestBatch(:,j)=synTestBatch(:,j).*handles.alignment.synchronization{handles.alignment.stages}.rng(j);
        end

        %[synTestBatch warptest] = onSyn(test_batch,Bref, handles.alignment.synchronization{handles.alignment.stages,1}.dtw.band,diag(handles.alignment.synchronization{handles.alignment.stages,1}.dtw.W), handles.alignment.dtw.zeta, handles.alignment.dtw.Xrng);

        st = size(synTestBatch);
        sr = size(Bref);
        nsamplesToPlot = size(synTestBatch,1);


case 'multisynchro'
        Xref = scale_(handles.alignment.synchronization{handles.alignment.stages,1}.Xref,handles.alignment.synchronization{handles.alignment.stages,1}.rng);
        test{1,1} = scale_(test_batch(:,1:end),handles.alignment.synchronization{handles.alignment.stages,1}.rng);

        [~,asynDetection] = high_multisynchro(test,Xref,handles.alignment.synchronization{handles.alignment.stages}.W,handles.alignment.synchronization{handles.alignment.stages}.Wconstr,handles.alignment.synchronization{handles.alignment.stages}.param.k,handles.alignment.synchronization{handles.alignment.stages}.param.psih,handles.alignment.synchronization{handles.alignment.stages}.param.psiv,1);
        
        test{1,1} = test_batch(:,1:end);
        [synTestBatch,warptest] = low_multisychro(test,handles.alignment.synchronization{handles.alignment.stages,1}.Xref,asynDetection,handles.alignment.synchronization{handles.alignment.stages}.Wconstr,handles.alignment.synchronization{handles.alignment.stages}.param.pcsMon,[],[],handles.alignment.synchronization{handles.alignment.stages}.specSynchronization);
        
        nsamplesToPlot = size(synTestBatch,1);

    otherwise
        errordlg('An error has ocurred in the selection of the synchronization method.','Synchronization Error');
end


