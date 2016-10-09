function varargout = Screening(varargin)
% SCREENING M-file for Screening.fig
%      SCREENING, by itself, creates a new SCREENING or raises the existing
%      singleton*.
%
%      H = SCREENING returns the handle to a new SCREENING or the handle to
%      the existing singleton*.
%
%      SCREENING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCREENING.M with the given input arguments.
%
%      SCREENING('Property','Value',...) creates a new SCREENING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Screening_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Screening_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Screening

% Last Modified by GUIDE v2.5 26-Sep-2014 12:01:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Screening_OpeningFcn, ...
                   'gui_OutputFcn',  @Screening_OutputFcn, ...
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


% --- Executes just before Screening is made visible.
function Screening_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Screening (see VARARGIN)

% Choose default command line output for Screening
handles.output = hObject;

handles.ParentsWindow=varargin{1};
handles.ParentFigure = guidata(handles.ParentsWindow);

% Start collecting all the measurements related to 
if length(varargin) ~=1, error('Error in the number of arguments.'); end;
handles.s_screening = handles.ParentFigure.s_screening;

% Setting the variable for handles of subplots;
handles.handles_subplots = [];
% Set the first batch as the starting batch for visualization
handles.selectedBatch = 1;
% Initialize all the objects from the GUI
handles=Initialize_visualization(handles);

set(handles.figure1, 'pointer', 'watch')
drawnow;
handles.rng_variables = 1:min(9,length(handles.VariableslbIn));
vars = handles.VariableslbIn(handles.rng_variables);
obs = handles.BatcheslbIn;
if length(handles.s_screening.batch_data) > 0
    [handles.auxx, handles.test] = prepareData(handles.s_screening.batch_data,handles.selectedBatch,handles.BatcheslbIn,vars);
    handles.handles_subplots=plot3D_batchtools(handles.auxx,[],obs,vars,handles.s_screening.varNames,handles.test,handles.uipanelPlots);
else
    errodlg('The data set selected is empty','!!Error!!');
end
set(handles.figure1, 'pointer', 'arrow');



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Screening wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = Screening_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  STRUCTURE
%% 
%%  1.- BATCHES PANEL
%%  2.- VARIABLES PANEL
%%  3.- BATCH SELECTION
%%  4.- VARIABLE SELECTION
%%  5.- BUTTONS
%%  6.- MENU BAR
%%  7.- AUXILIAR FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                               1.- BATCHES PANEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on selection change in lb_BatchesOut.
function lb_BatchesOut_Callback(hObject, eventdata, handles)
% hObject    handle to lb_BatchesOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lb_BatchesOut contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lb_BatchesOut

handles.selectedBatchesOut = get(handles.lb_BatchesOut,'Value');
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function lb_BatchesOut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb_BatchesOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.selectedBatchesOut=1;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pb_BatchesOut.
function pb_BatchesOut_Callback(hObject, eventdata, handles)
% hObject    handle to pb_BatchesOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if isempty(handles.selectedBatchesIn), warndlg('No batch ID was selected. No action will be taken.','!! Warning !!'); return
end

indx = ones(size(handles.BatcheslbInprev,1),1);
indx(handles.selectedBatchesIn) = 0;

handles.BatcheslbInprev = handles.BatcheslbInprev(find(indx==1));
handles.BatchesIn = zeros(handles.nBatches,1);
handles.BatchesIn(handles.BatcheslbInprev)=1;

% Clean and new setting of the listbox IN
set(handles.lb_BatchesIn,'String',' ');
set(handles.lb_BatchesIn,'String',handles.s_screening.batchNames(find(handles.BatchesIn==1),1));

% Clean and new setting of the listbox OUT
set(handles.lb_BatchesOut,'String',' ');
set(handles.lb_BatchesOut,'String',handles.s_screening.batchNames(find(handles.BatchesIn==0),1));

set(handles.lb_BatchesIn,'Value',1);
set(handles.lb_BatchesOut,'Value',1);


handles.selectedBatchesOut=1;
if size(handles.BatcheslbInprev) == 1
    handles.selectedBatchesIn=[]; 
    set(handles.lb_BatchesIn,'Value',0);
else
    handles.selectedBatchesIn=1;
end


% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in lb_BatchesIn.
function lb_BatchesIn_Callback(hObject, eventdata, handles)
% hObject    handle to lb_BatchesIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lb_BatchesIn contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lb_BatchesIn

handles.selectedBatchesIn = get(handles.lb_BatchesIn,'Value');
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function lb_BatchesIn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb_BatchesIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.selectedBatchesOut=1;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pb_BatchesIn.
function pb_BatchesIn_Callback(hObject, eventdata, handles)
% hObject    handle to pb_BatchesIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.selectedBatchesOut), warndlg('No batch ID was selected. No action will be taken.','!! Warning !!'); return 
end

handles.BatcheslbOut = find(handles.BatchesIn==0);
indx = ones(size(handles.BatcheslbOut,1),1);
indx(handles.selectedBatchesOut) = 0;

handles.BatchesIn(handles.BatcheslbOut(find(indx==0)))= 1;
handles.BatcheslbInprev = find(handles.BatchesIn==1);

% Clean and new setting of the listbox IN
set(handles.lb_BatchesIn,'String',' ');
set(handles.lb_BatchesIn,'String',handles.s_screening.batchNames(find(handles.BatchesIn==1),1));

% Clean and new setting of the listbox OUT
set(handles.lb_BatchesOut,'String',' ');
set(handles.lb_BatchesOut,'String',handles.s_screening.batchNames(find(handles.BatchesIn==0),1));

set(handles.lb_BatchesIn,'Value',1);
set(handles.lb_BatchesOut,'Value',1);
handles.selectedBatchesIn=1;

if size(handles.BatcheslbOut,1) == 1,
    handles.selectedBatchesOut=[];
    set(handles.lb_BatchesOut,'Value',0);
else
    handles.selectedBatchesOut=1;
end

% Update handles structure
guidata(hObject, handles);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                               2.- VARIABLES PANEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on selection change in lb_VariablesIn.
function lb_VariablesIn_Callback(hObject, eventdata, handles)
% hObject    handle to lb_VariablesIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lb_VariablesIn contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lb_VariablesIn

handles.selectedVariablesIn = get(handles.lb_VariablesIn,'Value');

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function lb_VariablesIn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb_VariablesIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.selectedVariablesIn=1;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pb_VariablesIn.
function pb_VariablesIn_Callback(hObject, eventdata, handles)
% hObject    handle to pb_VariablesIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.selectedVariablesOut), warndlg('No process variable was selected. No action will be taken.','!! Warning !!'); return
end

handles.VariableslbOut = find(handles.VariablesIn==0);
indx = ones(size(handles.VariableslbOut,1),1);
indx(handles.selectedVariablesOut) = 0;

handles.VariablesIn(handles.VariableslbOut(find(indx==0)))= 1;
handles.VariableslbInprev = find(handles.VariablesIn==1);

% Clean and new setting of the listbox IN
set(handles.lb_VariablesIn,'String',' ');
set(handles.lb_VariablesIn,'String',handles.s_screening.varNames(find(handles.VariablesIn==1),1));

% Clean and new setting of the listbox OUT
set(handles.lb_VariablesOut,'String',' ');
set(handles.lb_VariablesOut,'String',handles.s_screening.varNames(find(handles.VariablesIn==0),1));

set(handles.lb_VariablesIn,'Value',1);
set(handles.lb_VariablesOut,'Value',1);
handles.selectedVariablesIn=1;

if size(handles.VariableslbOut) == 1
    handles.selectedVariablesOut=[]; 
    set(handles.lb_VariablesOut,'Value',0);
else
    handles.selectedVariablesOut=1;
end


% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in lb_VariablesOut.
function lb_VariablesOut_Callback(hObject, eventdata, handles)
% hObject    handle to lb_VariablesOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lb_VariablesOut contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lb_VariablesOut

handles.selectedVariablesOut = get(handles.lb_VariablesOut,'Value');
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function lb_VariablesOut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb_VariablesOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.selectedVariablesOut=1;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pb_VariablesOut.
function pb_VariablesOut_Callback(hObject, eventdata, handles)
% hObject    handle to pb_VariablesOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.selectedVariablesIn), warndlg('No process variable was selected. No action will be taken.','!! Warning !!'); return
end

indx = ones(size(handles.VariableslbInprev,1),1);
indx(handles.selectedVariablesIn) = 0;

handles.VariableslbInprev = handles.VariableslbInprev(find(indx==1));
handles.VariablesIn = zeros(handles.nVariables,1);
handles.VariablesIn(handles.VariableslbInprev)=1;

% Clean and new setting of the listbox IN
set(handles.lb_VariablesIn,'String',' ');
set(handles.lb_VariablesIn,'String',handles.s_screening.varNames(find(handles.VariablesIn==1),1));

% Clean and new setting of the listbox OUT
set(handles.lb_VariablesOut,'String',' ');
set(handles.lb_VariablesOut,'String',handles.s_screening.varNames(find(handles.VariablesIn==0),1));

set(handles.lb_VariablesIn,'Value',1);
set(handles.lb_VariablesOut,'Value',1);


handles.selectedVariablesOut=1;
if size(handles.VariableslbInprev) == 1
    handles.selectedVariablesIn=[]; 
    set(handles.lb_VariablesIn,'Value',0); 
else
    handles.selectedVariablesIn=1;
end

handles.selectedVariablesOut=1;

% Update handles structure
guidata(hObject, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                3.- BATCH SELECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in pb_backward_batch.
function pb_backward_batch_Callback(hObject, eventdata, handles)
% hObject    handle to pb_backward_batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pos = find(handles.BatcheslbIn == handles.selectedBatch);

if pos == 1
    errordlg('No batch backward is available','File Error'); return;
else
    handles.selectedBatch = handles.BatcheslbIn(pos-1);
    set(handles.e_batch,'String',handles.selectedBatch);
end

for z=1:size(handles.s_screening.batch_data(1).data,2)
    for j=2:size(handles.auxx{1}.data{z},2);
        plot(handles.handles_subplots(j-1),handles.auxx{pos}.data{z}(:,1),handles.auxx{pos}.data{z}(:,j),'-','Color',[0.466667 0.533333 0.68]);
        plot(handles.handles_subplots(j-1),handles.auxx{pos-1}.data{z}(:,1),handles.auxx{pos-1}.data{z}(:,j),'r-','LineWidth',1);
    end
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pb_forwad_batch.
function pb_forwad_batch_Callback(hObject, eventdata, handles)
% hObject    handle to pb_forwad_batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pos = find(handles.BatcheslbIn == handles.selectedBatch);

if pos == numel(handles.BatcheslbIn), errordlg('No batch forward is available','File Error')
else
    handles.selectedBatch = handles.BatcheslbIn(pos+1);
    set(handles.e_batch,'String',handles.selectedBatch);
end

for z=1:size(handles.s_screening.batch_data(1).data,2)
    for j=2:size(handles.auxx{1}.data{z},2);
        plot(handles.handles_subplots(j-1),handles.auxx{pos}.data{z}(:,1),handles.auxx{pos}.data{z}(:,j),'-','Color',[0.466667 0.533333 0.68]);
        plot(handles.handles_subplots(j-1),handles.auxx{pos+1}.data{z}(:,1),handles.auxx{pos+1}.data{z}(:,j),'r-','LineWidth',1);
    end
end



% Update handles structure
guidata(hObject, handles);

function e_batch_Callback(hObject, eventdata, handles)
% hObject    handle to e_batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function e_batch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                              4.- VARIABLE SELECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function e_variables_Callback(hObject, eventdata, handles)
% hObject    handle to e_variables (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function e_variables_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_variables (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pb_forward_variables.
function pb_forward_variables_Callback(hObject, eventdata, handles)
% hObject    handle to pb_forward_variables (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.current_plot_window == handles.n_plots_windows, warndlg('No variable forward is available','!!Warning!!'); return;
else
   handles.current_plot_window = handles.current_plot_window + 1;
end

set(handles.figure1, 'pointer', 'watch')
drawnow;

if handles.current_plot_window == handles.n_plots_windows,handles.rng_variables = (handles.current_plot_window-1)*9+1:length(handles.VariableslbIn);
else
    handles.rng_variables = (handles.current_plot_window-1)*9+1:handles.current_plot_window*min(9,numel(handles.VariableslbIn));
end
    
[handles.auxx handles.test] = prepareData(handles.s_screening.batch_data,handles.selectedBatch,handles.BatcheslbIn,handles.VariableslbIn(handles.rng_variables));
delete(handles.handles_subplots);
handles.handles_subplots=plot3D_batchtools(handles.auxx,[],handles.BatcheslbIn,handles.VariableslbIn(handles.rng_variables),handles.s_screening.varNames,handles.test,handles.uipanelPlots);
set(handles.e_variables,'String',[num2str(handles.rng_variables(1)) '-' num2str(handles.rng_variables(end)) '/' num2str(max(9,length(handles.VariableslbIn)))]);

set(handles.figure1, 'pointer', 'arrow')

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pb_backward_variables.
function pb_backward_variables_Callback(hObject, eventdata, handles)
% hObject    handle to pb_backward_variables (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.current_plot_window == 1, warndlg('No variable backward is available','!!Warning!!'); return;
else
   handles.current_plot_window = handles.current_plot_window - 1;
end

set(handles.figure1, 'pointer', 'watch')
drawnow;
if handles.current_plot_window == 1,handles.rng_variables = 1:min(9,length(handles.VariableslbIn));
else
    handles.rng_variables = (handles.current_plot_window-1)*9+1:handles.current_plot_window*9;
end
    
[handles.auxx handles.test] = prepareData(handles.s_screening.batch_data,handles.selectedBatch,handles.BatcheslbIn,handles.VariableslbIn(handles.rng_variables));
delete(handles.handles_subplots);
handles.handles_subplots=plot3D_batchtools(handles.auxx,[],handles.BatcheslbIn,handles.BatcheslbInhandles.VariableslbIn(handles.rng_variables),handles.s_screening.varNames,handles.test,handles.uipanelPlots);
set(handles.e_variables,'String',[num2str(handles.rng_variables(1)) '-' num2str(handles.rng_variables(end)) '/' num2str(max(9,length(handles.VariableslbIn)))]);

set(handles.figure1, 'pointer', 'arrow')

% Update handles structure
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                   5.- BUTTONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function pb_DataSet_Callback(hObject, eventdata, handles)
% hObject    handle to pb_DataSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.rng_variables = handles.VariableslbInprev;
handles.ParentFigure.varNames = handles.ParentFigure.s_screening.varNames(handles.rng_variables,:);
handles.ParentFigure.x = prepareData(handles.s_screening.batch_data,handles.selectedBatch,handles.BatcheslbIn,handles.rng_variables)';

axes(handles.ParentFigure.main_window)
image(handles.ParentFigure.images{3});
axis off;
axis image;
set(handles.ParentFigure.pbAlignment,'Enable','on');

% Save information related to the unit, process variables and batches
% selected for bilinear process modelling
% System variables keeping information about batches
handles.ParentFigure.dataset.BatchesIn = handles.BatchesIn;

% System valriables keeping information about the process variables
handles.ParentFigure.dataset.VariablesIn = handles.VariablesIn;

handles.ParentFigure.track(2) = 1;
handles.ParentFigure.track(3:end) = 0;
% Save the information of the parent handle    
guidata(handles.ParentsWindow,handles.ParentFigure)
delete(handles.figure1);

% --- Executes on button press in pb_Refresh.
function pb_Refresh_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.handles_subplots);
set(handles.figure1, 'pointer', 'watch')
drawnow;

handles.BatcheslbIn= handles.BatcheslbInprev;
handles.VariableslbIn = handles.VariableslbInprev;
handles.selectedBatch = handles.BatcheslbIn(1);

% Set the batch to plot
set(handles.e_batch,'String',handles.selectedBatch);

% Update the number of window plots
handles.n_plots_windows=ceil(numel(handles.VariableslbIn)/9);
handles.current_plot_window = 1;
set(handles.e_variables,'String',['1-' num2str(min(9,length(handles.VariableslbIn))) '/' num2str(numel(handles.VariableslbIn))]);

handles.rng_variables = 1:min(9,length(handles.VariableslbIn));
[handles.auxx handles.test] = prepareData(handles.s_screening.batch_data,handles.selectedBatch,handles.BatcheslbIn,handles.VariableslbIn(handles.rng_variables));
handles.handles_subplots=plot3D_batchtools(handles.auxx,[],handles.BatcheslbIn,handles.VariableslbIn(handles.rng_variables),handles.s_screening.varNames,handles.test,handles.uipanelPlots);

set(handles.figure1, 'pointer', 'arrow')

% Update handles structure
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                   6.- MENU BAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------------------------------------------------------
function icon_Open_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to icon_Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile('*.mat', 'Select a MATLAB file');

if ~isequal(filename,0)
    S = load(strcat(pathname,filename), 'calibration');
    handles.s_screening = S.calibration;
    delete(handles.handles_subplots); 
        
    handles=Initialize_visualization(handles);
    set(handles.pm_Units,'Value',1);
else
    return;
end

set(handles.figure1, 'pointer', 'watch')
drawnow;
handles.rng_variables = 1:min(9,length(handles.VariableslbIn));
if ~isempty(handles.s_screening.batch_data)
    [x handles.test] = prepareData(handles.s_screening.batch_data,handles.selectedBatch,handles.BatcheslbIn,handles.rng_variables);
    handles.handles_subplots=plot3D_batchtools(x,[],handles.VariableslbIn(handles.rng_variables),handles.s_screening.varNames,handles.test,handles.uipanelPlots);
else
    errodlg('The data set selected is empty','!!Error!!');
end
set(handles.figure1, 'pointer', 'arrow')
set(handles.uipushtool_save,'Enable','on');

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function uipushtool_save_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.rng_variables = handles.VariableslbIn(1:min(9,length(handles.VariableslbIn)));
x = prepareData(handles.s_screening.batch_data,handles.selectedBatch,handles.BatcheslbIn,handles.rng_variables);
[file,path] = uiputfile('*.mat','Save Workspace As');

if isequal(file,0) || isequal(path,0)
   return;
else
   save([path file],'x');
end

% --------------------------------------------------------------------
function uipushtool_open_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                 7.- AUXILIAR FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,test] = prepareData(xtest,ybatch,obs,vars)

s = size(xtest);
x = cell(size(obs,1),1); 
counterb = 0;
test{1} = struct('data',[]);
for i=1:s(1), % batches
    x{counterb+1} = struct('data',cell(size(xtest(i).data,2),1));
    for j=1:size(xtest(i).data,2) % sampling frequencies
        x{counterb+1}.data{j} = [];
        if ~isempty(find(obs==i)),counterb = counterb + 1; x{counterb}.data{j} = xtest(i).data{j}(:,1);
        else
        break;    
        end
        if ybatch == i, test{1}.data{j} = xtest(i).data{j}(:,1); end
        for z=1:size(xtest(i).data{j},2) % variables
            [auxv, auxpos] = find(vars' == z);
            if ~isempty(auxv)
                 if ~isempty(find(obs==i))
                    x{counterb}.data{j} = [x{counterb}.data{j} xtest(i).data{j}(:,1+vars(auxpos))];
                 end
                 if ybatch == i
                    test{1}.data{j} = [test{1}.data{j} xtest(i).data{j}(:,1+vars(auxpos))];
                 end
            end
        end
    end 
end


function handles = Initialize_visualization(handles)

%     if ~iscell(handles.s_screening.batch_data)
%         handles.nBatches = size(handles.s_screening.batch_data,3);
%     else
        handles.nBatches = size(handles.s_screening.batch_data,1);
%     end
    for i=1:handles.nBatches, handles.nBatchesID(i,1) = size(handles.s_screening.batch_data,1);end

if handles.ParentFigure.track(2) == 0 

    % Setting the variables for the number of batches
    handles.BatchesIn = ones(handles.nBatches,1);

    % Number of process variables
    handles.nVariables = size(handles.s_screening.varNames,1);
    handles.VariablesIn = ones(handles.nVariables,1);   
else 
    % When a dataset is already set in the system, the window will only
    % show the unit, process variables and batches selected previously.
    
    % Setting the variables for the number of batches
    handles.BatchesIn = handles.ParentFigure.dataset.BatchesIn;

    % Number of process variables
    handles.VariablesIn = handles.ParentFigure.dataset.VariablesIn;
end

    % Setting the variables for the number of batches
    handles.BatcheslbIn = find(handles.BatchesIn==1); handles.BatcheslbInprev = handles.BatcheslbIn;
    handles.BatcheslbOut = [];

    % Number of process variables
    handles.nVariables = size(handles.s_screening.varNames,1);
    handles.VariableslbIn = find(handles.VariablesIn==1); handles.VariableslbInprev = handles.VariableslbIn;
    handles.VariableslbOut = find(handles.VariablesIn==0);

    % Creating the varibales to track the selected and unselected process
    % variables and batches
    handles.selectedVariablesOut = [];
    if ~isempty(handles.VariableslbOut), handles.selectedVariablesOut=1;end
    handles.selectedVariablesIn = 1;
    handles.selectedBatchesOut = [];
    handles.selectedBatchesIn = 1;    

    % Cleaning the list menus of the GUI
    set(handles.lb_VariablesIn,'String',num2str(handles.VariableslbIn));
    set(handles.lb_VariablesOut,'String',handles.VariableslbOut);
    set(handles.lb_BatchesIn,'String',handles.BatcheslbIn);
    set(handles.lb_BatchesOut,'String',handles.BatcheslbOut);
    
    % Clean and new setting of the listbox IN
    set(handles.lb_VariablesIn,'String',' ');
    set(handles.lb_VariablesIn,'String',handles.s_screening.varNames(find(handles.VariablesIn==1),1));
    % Clean and new setting of the listbox OUT
    set(handles.lb_VariablesOut,'String',' ');
    set(handles.lb_VariablesOut,'String',handles.s_screening.varNames(find(handles.VariablesIn==0),1));

    % Clean and new setting of the listbox IN
    set(handles.lb_BatchesIn,'String',' ');
    set(handles.lb_BatchesIn,'String',handles.s_screening.batchNames(find(handles.BatchesIn==1),1));

    % Clean and new setting of the listbox OUT
    set(handles.lb_BatchesOut,'String',' ');
    set(handles.lb_BatchesOut,'String',handles.s_screening.batchNames(find(handles.BatchesIn==0),1));

    set(handles.lb_BatchesIn,'Value',1);
    set(handles.lb_BatchesOut,'Value',1);
    set(handles.lb_VariablesIn,'Value',1);
    set(handles.lb_VariablesOut,'Value',1);
    
    % Put the first batch of the list to plot
    set(handles.e_batch,'String',handles.selectedBatch);
    % Set the parameter of the number of plots windows in the GUI
    set(handles.e_variables,'String',['1-' num2str(min(9,length(handles.VariableslbIn))) '/' num2str(max(9,length(handles.VariableslbIn)))]);

    %Selected batch        
    handles.selectedBatch = handles.BatcheslbIn(1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setting the number of plots windows
    handles.n_plots_windows=ceil(handles.nVariables/9);
    handles.current_plot_window = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
