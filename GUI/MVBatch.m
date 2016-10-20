function varargout = MVBatch(varargin)
% MVBATCH MATLAB code for MVBatch.fig
%      MVBATCH, by itself, creates a new MVBATCH or raises the existing
%      singleton*.
%
%      H = MVBATCH returns the handle to a new MVBATCH or the handle to
%      the existing singleton*.
%
%      MVBATCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MVBATCH.M with the given input arguments.
%
%      MVBATCH('Property','Value',...) creates a new MVBATCH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MVBatch_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MVBatch_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MVBatch

% Last Modified by GUIDE v2.5 07-Oct-2016 13:00:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MVBatch_OpeningFcn, ...
                   'gui_OutputFcn',  @MVBatch_OutputFcn, ...
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

% --- Executes just before MVBatch is made visible.
function MVBatch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MVBatch (see VARARGIN)

warning('OFF', 'ALL')
% Choose default command line output for MVBatch
handles.output = hObject;

% Load images of the front-page of the main window
load imagesBatchTools.mat
handles.images = batchToolsImages;
axes(handles.main_window);
image(handles.images{1});
axis off;
axis image;

% Initialize variables for the data matrix, data and models.
handles.x=[];
handles.data=[];
handles.s_screening = [];
handles.s_alignment = [];
handles.s_calibration = [];
handles.s_monitoring = [];
% Initialize system variables for keeping track of the process variables
% and batches selected in the data cleaning step.
handles.dataset.BatchesIn = [];
% System valriables keeping information about the process variables
handles.dataset.VariablesIn = [];
% Setting a binary vector to track the modelling procedure.
handles.track = zeros(5,1);

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = MVBatch_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pbCleaning.
function pbCleaning_Callback(hObject, eventdata, handles)
% hObject    handle to pbCleaning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.MVBatch, 'pointer', 'watch');
drawnow;
try
    Screening(handles.output);
catch err
    set(handles.MVBatch, 'pointer', 'arrow');
    errordlg(err.message);
end
set(handles.MVBatch, 'pointer', 'arrow');

% --- Executes on button press in pbAlignment.
function pbAlignment_Callback(hObject, eventdata, handles)
% hObject    handle to pbAlignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if numel(find(handles.track)) > 2
    set(handles.pbCalibration,'Enable','off');
    set(handles.pbMonitoring,'Enable','off');
    image(handles.images{3});
    axis off;
    axis image;
end
        
set(handles.MVBatch, 'pointer', 'watch');
drawnow;
try
    Alignment(handles.output);
catch err
    set(handles.MVBatch, 'pointer', 'arrow');
    errordlg(err.message);
end
    
%build3way(handles.output);
set(handles.MVBatch, 'pointer', 'arrow');

% --- Executes on button press in pbCalibration.
function pbCalibration_Callback(hObject, eventdata, handles)
% hObject    handle to pbCalibration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.MVBatch, 'pointer', 'watch');
try
    Modeling(handles.output);
catch err
   errordlg(err.message); 
   set(handles.MVBatch, 'pointer', 'arrow');
end

set(handles.MVBatch, 'pointer', 'arrow');

% --- Executes on button press in pbMonitoring.
function pbMonitoring_Callback(hObject, eventdata, handles)
% hObject    handle to pbMonitoring (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.main_window)
image(handles.images{6});
axis off;
axis image;
set(handles.MVBatch, 'pointer', 'watch');
try
    Monitoring(handles.output,handles.s_calibration, handles.s_alignment);
catch err
    errordlg(err.message);
    set(handles.MVBatch, 'pointer', 'arrow');
end
set(handles.MVBatch, 'pointer', 'arrow');
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pbLoad.
function pbLoad_Callback(hObject, eventdata, handles)
% hObject    handle to pbLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile('*.mat', 'Select a MATLAB file');

if ~isequal(filename,0)
    if  not(isempty(find(handles.track==1)))
        % Construct a questdlg with three options
        choice = questdlg('A data set is already loaded. The current batch data will be removed. Do you still want to proceed with the bilinear modeling?', ...
        'Yes', 'No');
        % Handle response
        switch choice
          case 'Yes'
              ;
          case 'No'
              return;
        end
    end
    
    try 
       S = load(strcat(pathname,filename), 'calibration');
       handles.s_screening = S.calibration;
    catch err
        % Give more information for mismatch.
        errordlg('An expected problem trying to load the structure ''calibration'' from the selected file has occurred');
        return;
    end 
    axes(handles.main_window);
    image(handles.images{1});
    axis off;
    axis image;
    % Reset the track vector
    handles.track(1) = 1;
    handles.track(2:end) = 0;
else
    return;
end

axes(handles.main_window);
image(handles.images{2});
axis off;
axis image;
set(handles.pbCleaning,'Enable','on');
set(handles.pbAlignment,'Enable','off');
set(handles.pbCalibration,'Enable','off');
set(handles.pbMonitoring,'Enable','off');

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menuFile_Callback(hObject, eventdata, handles)
% hObject    handle to menuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menuOpenProject_Callback(hObject, eventdata, handles)
% hObject    handle to menuOpenProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,pathname] = uigetfile('*.mat','Open MVBatch file');

if ~isequal(file, 0)
    
    if  not(isempty(find(handles.track==1)))
        % Construct a questdlg with three options
        choice = questdlg('Are you sure that you want to close out the current project?', ...
        'Open Project',...
        'Yes', 'No','Yes');
        % Handle response
        switch choice
          case 'Yes'
              ;
          case 'No'
              return;
        end
    end
    try
        eval(['load ''' fullfile(pathname, file) ''' dataStruct track dataset']);
    catch err
        errordlg(err.message); 
    end
    switch numel(find(track==1))
        case 1
            handles.data = dataStruct.data;
            set(handles.pbLoad,'Enable','on');
            set(handles.pbCleaning,'Enable','on');
            set(handles.pbAlignment,'Enable','off');
            set(handles.pbCalibration,'Enable','off');
            set(handles.pbMonitoring,'Enable','off');
            axes(handles.main_window)
            image(handles.images{2});
            axis off;
            axis image;
        case 2
            handles.data = dataStruct.data;
            handles.s_screening = dataStruct.s_screening;
            set(handles.pbLoad,'Enable','on');
            set(handles.pbCleaning,'Enable','on');
            set(handles.pbAlignment,'Enable','on');
            set(handles.pbCalibration,'Enable','off');
            set(handles.pbMonitoring,'Enable','off');
            axes(handles.main_window)
            image(handles.images{3});
            axis off;
            axis image;
        case 3
            handles.data = dataStruct.data;
            handles.s_screening = dataStruct.s_screening;
            handles.s_alignment = dataStruct.s_alignment;
            set(handles.pbLoad,'Enable','on');
            set(handles.pbCleaning,'Enable','on');
            set(handles.pbAlignment,'Enable','on');
            set(handles.pbCalibration,'Enable','on');
            set(handles.pbMonitoring,'Enable','off');
            axes(handles.main_window)
            image(handles.images{4});
            axis off;
            axis image;
        case 4
            handles.data = dataStruct.data;
            handles.s_screening = dataStruct.s_screening;
            handles.s_alignment = dataStruct.s_alignment;
            handles.s_calibration = dataStruct.s_calibration;
            set(handles.pbLoad,'Enable','on');
            set(handles.pbCleaning,'Enable','on');
            set(handles.pbAlignment,'Enable','on');
            set(handles.pbCalibration,'Enable','on');
            set(handles.pbMonitoring,'Enable','on');
            axes(handles.main_window)
            image(handles.images{5});
            axis off;
            axis image;
        case 5
            handles.data = dataStruct.data;
            handles.s_screening = dataStruct.s_screening;
            handles.s_alignment = dataStruct.s_alignment;
            handles.s_calibration = dataStruct.s_calibration;
            handles.s_monitoring = dataStruct.s_monitoring;
            set(handles.pbLoad,'Enable','on');
            set(handles.pbCleaning,'Enable','on');
            set(handles.pbAlignment,'Enable','on');
            set(handles.pbCalibration,'Enable','on');
            set(handles.pbMonitoring,'Enable','on');
            axes(handles.main_window)
            image(handles.images{6});
            axis off;
            axis image;
        otherwise
            errordlg('An error has been produced opening the project.');
    end
end

% --------------------------------------------------------------------
function menuSaveProject_Callback(hObject, eventdata, handles)
% hObject    handle to menuSaveProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(find(handles.track==1))
    warndlg('Empy project. No modeling step has been completed');
    return
end

[file, pathname] = uiputfile('*.mat','Save to File');

if ~isequal(file, 0)
    switch numel(find(handles.track==1))
        case 1
            dataStruct.data = handles.data;
        case 2
            dataStruct.data = handles.data;
            dataStruct.s_screening = handles.s_screening;
        case 3
            dataStruct.data = handles.data;
            dataStruct.s_screening = handles.s_screening;
            dataStruct.s_alignment = handles.s_alignment;
        case 4
            dataStruct.data = handles.data;
            dataStruct.s_screening = handles.s_screening;
            dataStruct.s_alignment = handles.s_alignment;
            dataStruct.s_calibration = handles.s_calibration;
        case 5
            dataStruct.data = handles.data;
            dataStruct.s_screening = handles.s_screening;
            dataStruct.s_alignment = handles.s_alignment;
            dataStruct.s_calibration = handles.s_calibration;
            dataStruct.s_monitoring = handles.s_monitoring;
        otherwise
            errordlg('An error has been produced saving the project.');
    end
        track = handles.track;
        dataset = handles.dataset;
        eval(['save ''' fullfile(pathname, file) ''' dataStruct track dataset']);
end
