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

% Last Modified by GUIDE v2.5 16-Jun-2017 11:02:45

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
handles.batchTime = 1;

if length(varargin)<1, error('Error in the number of arguments.'); end;


% Obtain the parent windoe handles
handles.ParentsWindow=varargin{1};
handles.ParentFigure = guidata(handles.ParentsWindow);

handles.calibration = varargin{2};
handles.alignment = varargin{3};
 
% Populate the model popmenu
handles.models = []; 
handles.monitoringFlag = 1; % flag in case the user imports data sets with missing data to be treated internally
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

% Matrix containing the test batch to be projected onto the latent
% structure.
handles.calibration.test_batch = [];
% Flag indicating what mode must be used for fault diagnosis (by default
% off-line, i.e. a value equal to 0.
handles.modeMonitoring = 0;
% Flag indicating what multivariate statistics should be investigated in
% fault diagnosis. By default SPE is selected, i.e. a value equal to 0.
handles.statFaultDiagnosis = 0;
% Inititialize the counter of the selected data set. 
handles.selectedDataSet = 1;

% Enable the monitoring option for the calibration data set
handles.calibration.test{1,1} = handles.alignment.synchronization{1,1}.nor_batches';

set(handles.popupmenuVar,'String',' '); 
set(handles.popupmenuVar,'String','calibration'); 
set(handles.textData,'String','calibration'); 

handles.calibration.test_batch = handles.calibration.test{handles.selectedDataSet}{1};

set(handles.popupmenuBat,'String','');

% Fill the popmenu with the number of batches
set(handles.popupmenuBat,'String','');
for i=1:length(handles.calibration.test{handles.selectedDataSet})
    contents = get(handles.popupmenuBat,'String');
    set(handles.popupmenuBat,'String',strvcat(contents,[' ',num2str(i)]));
end

set(handles.popupmenuBat,'Value',1);

% Fill the popmenu with the number of sampling points
set(handles.popupmenuTimePoint,'String','');
for i=1:size(handles.calibration.x,1)
    contents = get(handles.popupmenuTimePoint,'String');
    set(handles.popupmenuTimePoint,'String',strvcat(contents,[' ',num2str(i)]));
end

% Update the Monitoring options
set(handles.pbContribution,'Enable','off');
set(handles.popupmenuVar,'Enable','on');   
set(handles.textVar,'Enable','on');
set(handles.popupmenuBat,'Enable','off');
set(handles.textBat,'Enable','on');
set(handles.pushbuttonPlo,'Enable','off');
set(handles.ImportMenuItem,'Enable','off');
set(handles.radiobuttonCV,'Enable','on');
set(handles.ImportMenuItem,'Enable','on');


handles.ParentFigure.track(:) = 1;

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

% Update handles structure
guidata(hObject, handles);

popupmenuMod_Callback(handles.popupmenuMod, eventdata, handles);

% UIWAIT makes Monitoring wait for user response (see UIRESUME)
% uiwait(handles.figureMonitoring);

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

nmodel = get(hObject,'Value');
if nmodel <= length(handles.calibration.man_mp_group),             
    model = handles.calibration.man_mp_group{nmodel};
elseif nmodel  <= length(handles.calibration.man_mp_group)+length(handles.calibration.mp_group),
     model = handles.calibration.mp_group{nmodel-length(handles.calibration.man_mp_group)};
else
    model = handles.calibration.mp_group2{nmodel-length(handles.calibration.man_mp_group)-length(handles.calibration.mp_group)};
end

handles.calibration.model = model;

% Check whether the monitoring systems of the selected model has been
% cross-validated
if isempty(handles.calibration.LVmodels(nmodel).latent_structure.cvD) % If no cross-validation has been performed
    set(handles.radiobuttonCV,'Enable','off');
    set(handles.radiobuttonCV,'Value',0);
else
    set(handles.radiobuttonCV,'Enable','on');
    set(handles.radiobuttonCV,'Value',1);
end

% Check whether the current test set has been projected on the latent model
if isempty(handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).evolD)
   set(handles.popupmenuBat,'Enable','off');
   set(handles.pushbuttonPlo,'Enable','off');
   set(handles.pbContribution,'Enable','off');
else
   set(handles.popupmenuBat,'Enable','on');
   set(handles.pushbuttonPlo,'Enable','on');
   set(handles.pbContribution,'Enable','on');
end

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

% Retrieve the model ID
nmodel = get(handles.popupmenuMod,'Value');


if isinf(handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_alpha.alpd95cv)
    % cross-validate the monitoring systems, yielding the cross-validated control limits     
     [handles.calibration.LVmodels(nmodel).latent_structure.cvevolD,...  
     handles.calibration.LVmodels(nmodel).latent_structure.cvevolQ,...
     handles.calibration.LVmodels(nmodel).latent_structure.cvD,...
     handles.calibration.LVmodels(nmodel).latent_structure.cvQ,...
     handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limd95cv,...
     handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limd99cv,...
     handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limq95cv,...
     handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limq99cv,...
     handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_alpha.alpd95cv,...
     handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_alpha.alpd99cv,...
     handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_alpha.alpq95cv,...
     handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_alpha.alpq99cv,...
     handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limd95cv,...
     handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limd99cv,...
     handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limq95cv,...
     handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limq99cv,...
     handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_alpha.alpd95cv,...
     handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_alpha.alpd99cv,...
     handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_alpha.alpq95cv,...
     handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_alpha.alpq99cv] = crossvalMVstatistics(handles.calibration.x,handles.calibration.LVmodels(nmodel).latent_structure.phases);
end
 
%% Display the control charts for the cross-validated multivariate statistics
visualizeStatEvolving(handles.calibration.LVmodels(nmodel).latent_structure.cvevolD, ...
handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limd95,...
handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limd99,...
handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limd95cv,...
handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limd99cv,...
[],...
'D-statistic',...
handles.axes1);

visualizeStatEvolving(handles.calibration.LVmodels(nmodel).latent_structure.cvevolQ, ...
handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limq95,...
handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limq99,...
handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limq95cv,...
handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limq99cv,...
[],...
'Q-statistic',...
handles.axes2);

cla(handles.axes3);
set(handles.axes3,'Visible','off');

try
    if strcmp(handles.alignment.synchronization{handles.alignment.stages,1}.methodsyn,'dtw') || strcmp(handles.alignment.synchronization{handles.alignment.stages,1}.methodsyn,'multisynchro')
        visualizeWarping(handles.alignment.synchronization{handles.alignment.stages,1}.warp, handles.alignment.synchronization{handles.alignment.stages,1}.band,[],[],[],true,handles.axes3);
         set(handles.axes3,'Visible','on');
    end
catch err
   errordlg(err.message); 
end

% If the model is batch-wise, visualize the overall D and Q statistic values
% for all batches
if ~isnan(handles.calibration.LVmodels(nmodel).latent_structure.cvD)
    
    visualizeStatGlobal(handles.calibration.LVmodels(nmodel).latent_structure.cvD,...
                        handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limd95,...
                        handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limd99,...
                        handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limd95cv,...
                        handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limd99cv,...
                        [],...
                        [],...
                        'D-statistic',handles.axesPostBatchT2);

    visualizeStatGlobal(handles.calibration.LVmodels(nmodel).latent_structure.cvQ,...
        handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limq95,...
        handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limq99,...
        handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limq95cv,...
        handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limq99cv,...
        [],...
        [],...
        'Q-statistic',handles.axesPostBatchSPE)
else
    cla(handles.axesPostBatchT2);
    set(handles.axesPostBatchT2,'Visible','off');
    cla(handles.axesPostBatchSPE);
    set(handles.axesPostBatchSPE,'Visible','off');
end

% Enable the option of using the control limits adjusted by
% cross-validation
set(handles.radiobuttonCV,'Enable','on');

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


% Populate the popmenu Batches for the selected data set
for i=1:length(handles.calibration.test{handles.selectedDataSet})
    contents = get(handles.popupmenuBat,'String');
    set(handles.popupmenuBat,'String',strvcat(contents,[' ',num2str(i)]));
end

% Visualize the first batch ID of the data set
set(handles.popupmenuBat,'Value',1);
handles.calibration.test_batch = handles.calibration.test{handles.selectedDataSet}{1};

% Check whether the current test set has been projected on the latent model
nmodel = get(handles.popupmenuMod,'Value');
if isempty(handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).evolD)
   set(handles.popupmenuBat,'Enable','off');
   set(handles.pushbuttonPlo,'Enable','off');
   set(handles.pbContribution,'Enable','off');
else
   set(handles.popupmenuBat,'Enable','on');
   set(handles.pushbuttonPlo,'Enable','on');
   set(handles.pbContribution,'Enable','on');
end

% Update handles structure
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

% Update handles structure
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

% Obtain the batch ID to monitor
nbatch = get(handles.popupmenuBat,'Value');

% Retrieve the index of the current LV model
nmodel = get(handles.popupmenuMod,'Value');

% Synchronize the corresponding batch
[synTestBatch,nsamplesToPlot,warptest] = onlineSynchronization(handles, handles.calibration.test_batch);

vars = cellstr(get(handles.popupmenuVar,'String'));
posTest = length(handles.alignment.x) + 1;
if strcmp(vars{get(handles.popupmenuVar,'Value')},'calibration')
    posTest = get(handles.popupmenuBat,'Value');
end

% Clear up the axes for the post-batch monitoring
cla(handles.axesPostBatchT2);
cla(handles.axesPostBatchSPE);
set(handles.axesPostBatchT2,'Visible','off');
set(handles.axesPostBatchSPE,'Visible','off');


%% Visualize the multivariate statistics

if get(handles.radiobuttonCV,'Value')
    % cross-validated online limits
    olimd95 = handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limd95cv;
    olimd99 = handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limd99cv;
    olimq95 = handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limq95cv;
    olimq99 = handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limq99cv;
    % cross-validated offline limits
    offlimd95 = handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limd95cv;
    offlimd99 = handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limd99cv;
    offlimq95 = handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limq95cv;
    offlimq99 = handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limq99cv;
else
    % online limits
    olimd95 = handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limd95;
    olimd99 = handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limd99; 
    olimq95 = handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limq95;
    olimq99 = handles.calibration.LVmodels(nmodel).latent_structure.control_limits.online_cl.limq99;
    % offline limits
    offlimd95 = handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limd95;
    offlimd99 = handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limd99;
    offlimq95 = handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limq95;
    offlimq99 = handles.calibration.LVmodels(nmodel).latent_structure.control_limits.offline_cl.limq99;
end

% Online D statistics
visualizeStatEvolving(handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).evolD(:,nbatch), ...
olimd95,...
olimd99,...
[],...
[],...
1,...
'D-statistic',...
handles.axes1);

% Online Q statistic
visualizeStatEvolving(handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).evolQ(:,nbatch), ...
olimq95,...
olimq99,...
[],...
[],...
1,...
'Q-statistic',...
handles.axes2);

if strcmp(handles.alignment.synchronization{handles.alignment.stages}.methodsyn,'dtw') || strcmp(handles.alignment.synchronization{handles.alignment.stages}.methodsyn,'multisynchro')
    try
        visualizeWarping(handles.alignment.synchronization{handles.alignment.stages}.warp, handles.alignment.synchronization{handles.alignment.stages}.band,warptest,[],[],false,handles.axes3);
    catch err
       errordlg(err.message); 
    end
end

% If the model is batch-wise, visualize the overall D and Q statistic values
% for all batches
if ~isempty(offlimd95) && numel(find(isnan(offlimd95))) ~= length(offlimd95)
    % Offline D statistics
    visualizeStatGlobal(handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).D,...
                        offlimd95,...
                        offlimd99,...
                        [],...
                        [],...
                        1,...
                        nbatch,...
                        'D-statistic',handles.axesPostBatchT2);

    % Offline Q statistics
    visualizeStatGlobal(handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).Q,...
                        offlimq95,...
                        offlimq99,...
                        [],...
                        [],...
                        1,...
                        nbatch,...
                        'Q-statistic',handles.axesPostBatchSPE)
end


set(handles.pbContribution,'Enable','on');

function radiobuttonCV_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonCV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbuttonMonitorBatches.
function pushbuttonMonitorBatches_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonMonitorBatches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Retrieve the model ID
nmodel = get(handles.popupmenuMod,'Value');

if ~isempty(handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).evolD)
    warndlg('Test data set already projected on the latent variable model','Message');
    return;
end

% Check whether the test batches have already been projected
% Total number of test batches
ntestbatches = length(handles.calibration.test{handles.selectedDataSet});

% Check whether the selected test set has already been projected. 
if isempty(handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).evolD)        
    h = waitbar(0/ntestbatches,sprintf('Projecting Batch #%d onto the latent structure',0),'Name','Model exploitation');
    try
        for b=1:ntestbatches
            [synTestBatch,nsamplesToPlot,warptest] = onlineSynchronization(handles, handles.calibration.test{handles.selectedDataSet}{b});

            %  Load the dataset for theright  [handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet). evolD...,
             [evolD,evolQ,D,Q,cont_evolD,cont_evolQ] = multiphaseProjection(handles.calibration.x,synTestBatch,...
                                                                              handles.calibration.LVmodels(nmodel).latent_structure.phases,...
                                                                              handles.calibration.LVmodels(nmodel).latent_structure.P,...
                                                                              handles.calibration.LVmodels(nmodel).latent_structure.T,...
                                                                              handles.calibration.LVmodels(nmodel).latent_structure.mn,...
                                                                              handles.calibration.LVmodels(nmodel).latent_structure.stnd,1);

             handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).evolD = [handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).evolD,evolD];
             handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).evolQ = [handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).evolQ,evolQ];
             handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).D     = [handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).D;D];
             handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).Q     = [handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).Q;Q];
             handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).cont_evolD = [handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).cont_evolD,cont_evolD];
             handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).cont_evolQ = [handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).cont_evolQ,cont_evolQ];
             waitbar(b/ntestbatches,h,sprintf('Projecting Batch #%d onto the latent structure',b),'Name','Model exploitation');
        end
    catch err
       errordlg(sprintf('Error when projecting test batches on the latent structure. %s',err.message)); 
       close(h);
       return;
    end
    close(h);
end
% Update handles structure
guidata(hObject,handles);

% Enable fault diagnosis
set(handles.pbContribution,'Enable','on');
% Enable the rest of functions
set(handles.popupmenuBat,'Enable','on');
set(handles.pushbuttonPlo,'Enable','on');

% Update handles structure
guidata(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                  3.- FAULT DIAGNOSIS PANEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in popupmenuTimePoint.
function popupmenuTimePoint_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuTimePoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuTimePoint contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuTimePoint

handles.batchTime = get(handles.popupmenuTimePoint,'Value');
% Update handles structure
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

ibatch = get(handles.popupmenuBat,'Value');
test_batch = handles.calibration.test{handles.selectedDataSet}{ibatch};

test_batchSyn = onlineSynchronization(handles,test_batch);
    
if ~handles.modeMonitoring && handles.calibration.model.phases(1,3) ~= size(handles.alignment.alg_batches,1)-1, 
    warndlg('Overall contributions cannot be displayed because the model is not a batch-wise PCA. Only contributions at a specific sampling time can be depicted. Please, enable the on-line fault diagnosis.');
    return;
end

try 
    % Retrieve the index of the current LV model
    nmodel = get(handles.popupmenuMod,'Value');

    % Retrieve the position of the process variables in the two-way array for the selected batch of the test set 
    indx = (handles.calibration.LVmodels(nmodel).latent_structure.dimensions.nvariables*(ibatch-1))+1:handles.calibration.LVmodels(nmodel).latent_structure.dimensions.nvariables*ibatch;

    % Pull out the contributions of the whole batch selected
    EvolvingDcontribution = handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).cont_evolD(:,indx);
    EvolvingQcontribution = handles.calibration.LVmodels(nmodel).latent_structure.monitoring_statistics(handles.selectedDataSet).cont_evolQ(:,indx);

    s=size(handles.calibration.x);
    
    % Retrieve the two-array containing the evolving contribution to the
    % selected statistic
    switch handles.statFaultDiagnosis
            case 0
            % D statistic
            contToStat = EvolvingDcontribution;
            label = 'Contribution to D statistic ';
            case 1
            % Q statistic
            contToStat = EvolvingQcontribution;
            label = 'Contribution to Q statistic ';
    end
   
    % Arrange the contributions based on the type (0: overall, 1:
    % intantaneous)
    
    [~,av] = preprocess3D(handles.calibration.x,handles.calibration.model.arg.prep);
    
    switch handles.modeMonitoring  
        case 0
        % Case 1: Evolving contribution to the multivariate statistics at time
        % point k
        ovContrib = nan(s(1)*s(2),1);
        for j=1:s(2)
            ovContrib((j-1)*s(1)+1:j*s(1)) =  contToStat(:,j);
        end 
        jContrib = nan(s(2),1);
        for j=1:s(2)
            jContrib(j,1) = nansum(contToStat(:,j));
        end
        kContrib = nan(s(1),1);
        for k=1:s(1)
            kContrib(k,1)= nansum(contToStat(k,:));
        end
        % Call function to display the overall contribution, and the
        % contribution per variable and batch time
        
        Diagnosis(handles.calibration.x,test_batchSyn,av,ovContrib,jContrib,kContrib,handles.calibration.model.phases(:,4),handles.ParentFigure.varNames,label);
        
        case 1
        % Case 2: Evolving contribution to the multivariate statistics for the
        % whole batch
        switch handles.statFaultDiagnosis
            case 0
            % D statistic
            contToStat = contToStat(handles.batchTime,:);
            
            case 1
            % Q statistic
            contToStat = contToStat(handles.batchTime,:);
            
        end
        createFigure(handles.calibration.x,test_batchSyn,av,handles.ParentFigure.varNames(:,1));
        bar(contToStat);
        xlabel('Variables','FontSize',12,'FontWeight','bold');
        ylabel(strcat(label,sprintf(' (k=%d)',handles.batchTime)),'FontSize',12,'FontWeight','bold');
    end
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
% Update handles structure
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
% Update handles structure
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
       aux = load(strcat(pathname,filename),'data');
       xtest = aux.data;
    catch 
        % Give more information for mismatch.
        errordlg('An expected problem trying to load the test set from the selected file has occurred.');
    end 
else
    return;
end
   
if isvector(xtest) && size(xtest,1)>0
    for j=1:size(xtest,1)
        if size(xtest{j},2)~=size(handles.ParentFigure.s_screening.data{1,1},2)
            errordlg('The test data set does not contain the same number of variables with common sampling frequency as the calibration data set.');
            return;
        end

    end
else
   errordlg('Recall that the variable trajectories must be contained in a Matlab column-wise cell array.'); 
end
   
% Identify those batches that contain missing data
md = zeros(length(xtest),1);

for i=1:length(xtest)
    if numel(find(isnan(xtest{i})))>0, md(i)=1; end
end

handles.BatcheslbIn = 1:length(xtest);
handles.filename = filename;
handles.xtest = xtest;

% Update handles structure
guidata(hObject,handles);

% In case that that there exist missing values in any of the batches, call
% the GUI for missing data imputation
if numel(find(md))>0
    MissingDataImputation(handles.output,md);
else
    handles.calibration.test = [handles.calibration.test; {xtest}];
    handles.calibration.test_batch = xtest{1};
    contents = strvcat(get(handles.popupmenuVar,'String'),filename);
    set(handles.popupmenuVar,'String',contents);
    set(handles.popupmenuVar,'Value',size(get(handles.popupmenuVar,'String'),1));
    set(handles.textData,'String',contents(size(contents,1),:));
    
    handles.selectedDataSet = size(contents,1);
    
    set(handles.popupmenuBat,'String','');
    for i=1:length(xtest)
        contents = get(handles.popupmenuBat,'String');
        set(handles.popupmenuBat,'String',strvcat(contents,[' ',num2str(i)]));
    end
    
    % Create a new instance of the Monitoring Parameters class object in each
    % model
    
    nmodels = length(handles.calibration.LVmodels);
    
    for n=1:nmodels
        ntestsets = length(handles.calibration.LVmodels(n).latent_structure.monitoring_statistics);
        handles.calibration.LVmodels(n).latent_structure.monitoring_statistics(ntestsets+1) = LatentStructure.MonitoringParameters(filename);
    end
    
    
    % Update the user interface
    set(handles.popupmenuBat,'Value',1);
    set(handles.popupmenuVar,'Enable','on');
    set(handles.textVar,'Enable','on');
    set(handles.textData,'Enable','on');
    set(handles.popupmenuBat,'Enable','off');
    set(handles.textBat,'Enable','on');
    set(handles.pushbuttonPlo,'Enable','off');
    
    % Update handles structure
    guidata(hObject,handles);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                    5.- AUXILIAR FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hf] = createFigure(xini,test,av,tagnames)

if nargin < 1, type = 0; end

[select_cdata,~,alpha] = imread('selectdata.png');
select_cdata = double(select_cdata) / 255;
select_cdata(~alpha) = NaN;

[tags_cdata,~,alpha] = imread('tags.png');
tags_cdata = double(tags_cdata) / 255;
tags_cdata(~alpha) = NaN;

[Trace_cdata,~,alpha] = imread('Trace.png');
Trace_cdata = double(Trace_cdata) / 255;
Trace_cdata(~alpha) = NaN;

% Create the toolbar
hf = figure;
tbh = uitoolbar(hf);


% Add a toggle tool to the toolbar
uipushtool(tbh,'CData',select_cdata,'Separator','on',...
'TooltipString','Select',...
'HandleVisibility','off','ClickedCallback',@hSelectCallback);  

uipushtool(tbh,'CData',Trace_cdata,'Separator','on',...
'TooltipString','Trace',...
'HandleVisibility','off','ClickedCallback',{@VisualizeVariablesCallback,xini,test,av,tagnames});

uipushtool(tbh,'CData',tags_cdata,...
'TooltipString','Tagname',...
'HandleVisibility','off','ClickedCallback',{@hTagsCallback,tagnames});
    
function hSelectCallback(hObject, eventdata)
% Callback function run when the Update button is pressed 
brush on

function hTagsCallback(hObject, eventdata, tagnames)
% Callback function run when the Update button is pressed

hLine = get(gca,'Children');
gname(tagnames,hLine(end));

function VisualizeVariablesCallback(hObject,eventdata,xini,test,av,tagnames)

% Due to changes in the new version of Matlab as of 2015, we need to check
% for special coding to indetify the bar

versionMatlab = version('-release');
versionMatlab = str2num(versionMatlab(1:4));

if versionMatlab < 2015
    hBrushLine = findall(gca, 'tag', 'Brushing');
    brushedData = get(hBrushLine, {'Xdata', 'Ydata'});
    if numel(find(~cellfun(@isempty,brushedData)))==0,
        warndlg('Please, press the ''Select'' icon, select a bar and then press the ''Trace'' icon.')
        return;
    end
    if numel(find(brushedData{1,2}~=0))==0
        errordlg('No variable has been selected, please select a variable and try again.','Error')
        return
    end
    Variable = find(brushedData{1,2});
    if numel(Variable)>1
        errordlg('You have selected multiple variables. Please, select a single variable.','Error')
        return
    end
else
    hLine = get(gca,'Children');
    selection = logical(hLine(end).BrushData);
    if numel(find(selection)) == 0
        errordlg('No variable has been selected, please select a variable and try again.','Error')
        return
    end
    if numel(find(selection)) > 1
        errordlg('Only one variable at a time can be selected.','Error')
        return
    end
    Variable = find(selection);
end
brush off;

figure;
cal = squeeze(xini(:,Variable,:));
h1 = plot(cal(:,1),'Color',[0.662745 0.662745 0.662745],'LineWidth',1.5);hold on;
plot(cal,'Color',[0.662745 0.662745 0.662745],'LineWidth',1.5);
h2 = plot(squeeze(test(:,Variable,:)),'r-','LineWidth',1.5);
h3 = plot(av(:,Variable),'g-','LineWidth',1.5); 
xlabel('Time','FontSize',14)
ylabel(tagnames{Variable,1},'FontSize',14)
legend([h1,h2,h3],'Historical batches','Test batch','Mean trajectory','Location','best');
axis tight;

function [synTestBatch,nsamplesToPlot,warptest] = onlineSynchronization(handles, test_batch)

warptest = [];

varIn = find(handles.ParentFigure.s_screening.VariablesIn);

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
            uBn =  test_batch(min(indm,indM):max(indm,indM),varIn);
        end

             if length(find(isnan(uBn(:,handles.alignment.synchronization{handles.alignment.stages,1}.var)))) <= 0.25*length(uBn(:,handles.alignment.synchronization{handles.alignment.stages,1}.var)) % Aling every batch where the iv was measured
                synTestBatch = align_IV(uBn,handles.alignment.synchronization{handles.alignment.stages,1}.var,handles.alignment.synchronization{handles.alignment.stages,1}.steps,handles.alignment.synchronization{handles.alignment.stages,1}.method);
                synTestBatch = synTestBatch(:,2:end); % remove time warping profile;
                nsamplesToPlot = size(synTestBatch,1);
             else
                 errordlg('Too much missing data in the indicator variable.','Error Dialog','modal');
             end
         
     case 'dtw'
        Bref = scale_(handles.alignment.synchronization{handles.alignment.stages,1}.Xref,handles.alignment.synchronization{handles.alignment.stages,1}.rng);
        test = scale_(test_batch(:,varIn),handles.alignment.synchronization{handles.alignment.stages,1}.rng);
        [synTestBatch,warptest] = DTW(test,Bref,diag(handles.alignment.synchronization{handles.alignment.stages}.W));      
        
        for j=1:size(handles.alignment.synchronization{handles.alignment.stages}.nor_batches{1},2)
            synTestBatch(:,j)=synTestBatch(:,j).*handles.alignment.synchronization{handles.alignment.stages}.rng(j);
        end

        %[synTestBatch warptest] = onSyn(test,Bref, handles.alignment.synchronization{handles.alignment.stages,1}.band,diag(handles.alignment.synchronization{handles.alignment.stages,1}.W), handles.alignment.synchronization{handles.alignment.stages,1}.zeta, handles.alignment.synchronization{handles.alignment.stages,1}.rng);
        nsamplesToPlot = size(synTestBatch,1);


case 'multisynchro'
        Xref = scale_(handles.alignment.synchronization{handles.alignment.stages,1}.Xref,handles.alignment.synchronization{handles.alignment.stages,1}.rng);
        test{1,1} = scale_(test_batch(:,varIn),handles.alignment.synchronization{handles.alignment.stages,1}.rng);

        [~,asynDetection] = high_multisynchro(test,Xref,handles.alignment.synchronization{handles.alignment.stages}.W,handles.alignment.synchronization{handles.alignment.stages}.Wconstr,handles.alignment.synchronization{handles.alignment.stages}.param.k,handles.alignment.synchronization{handles.alignment.stages}.param.psih,handles.alignment.synchronization{handles.alignment.stages}.param.psiv,1);
        
        {1,1} = test_batch(:,varIn);
        [synTestBatch,warptest] = low_multisychro(test,handles.alignment.synchronization{handles.alignment.stages,1}.Xref,asynDetection,handles.alignment.synchronization{handles.alignment.stages}.Wconstr,handles.alignment.synchronization{handles.alignment.stages}.param.pcsMon,handles.alignment.synchronization{handles.alignment.stages,1}.maxIter,[],[],handles.alignment.synchronization{handles.alignment.stages}.specSynchronization);
        
        nsamplesToPlot = size(synTestBatch,1);

    otherwise
        errordlg('An error has ocurred in the selection of the synchronization method.','Synchronization Error');
end



