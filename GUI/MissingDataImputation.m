function varargout = MissingDataImputation(varargin)
% MISSINGDATAIMPUTATION MATLAB code for MissingDataImputation.fig
%      MISSINGDATAIMPUTATION, by itself, creates a new MISSINGDATAIMPUTATION or raises the existing
%      singleton*.
%
%      H = MISSINGDATAIMPUTATION returns the handle to a new MISSINGDATAIMPUTATION or the handle to
%      the existing singleton*.
%
%      MISSINGDATAIMPUTATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MISSINGDATAIMPUTATION.M with the given input arguments.
%
%      MISSINGDATAIMPUTATION('Property','Value',...) creates a new MISSINGDATAIMPUTATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MissingDataImputation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MissingDataImputation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MissingDataImputation

% Last Modified by GUIDE v2.5 16-Oct-2017 16:53:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MissingDataImputation_OpeningFcn, ...
                   'gui_OutputFcn',  @MissingDataImputation_OutputFcn, ...
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


% --- Executes just before MissingDataImputation is made visible.
function MissingDataImputation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MissingDataImputation (see VARARGIN)

% Choose default command line output for MissingDataImputation
handles.output = hObject;

% Retrieve the handles from the Screening GUI
handles.ParentsWindow=varargin{1};
handles.ParentFigure = guidata(handles.ParentsWindow);

% Populate the listbox of batches with missing values
handles.md = find(varargin{2});
handles.mdmaster = find(varargin{2});

handles = Initialize(handles);

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

% UIWAIT makes MissingDataImputation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MissingDataImputation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listboxImputedbatches.
function listboxImputedbatches_Callback(hObject, eventdata, handles)
% hObject    handle to listboxImputedbatches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxImputedbatches contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxImputedbatches

% --- Executes during object creation, after setting all properties.
function listboxImputedbatches_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxImputedbatches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listboxMDbatches.
function listboxMDbatches_Callback(hObject, eventdata, handles)
% hObject    handle to listboxMDbatches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxMDbatches contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxMDbatches

handles.selectedBatch = get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function listboxMDbatches_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxMDbatches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.selectedBatch = 1;
% Update handles structure
guidata(hObject, handles);


function editLags_Callback(hObject, eventdata, handles)
% hObject    handle to editLags (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLags as text
%        str2double(get(hObject,'String')) returns contents of editLags as a double

handles.lags = str2double(get(hObject,'String'));     
if handles.lags <0, errordlg('The number of lags cannot be lower than 0.'); end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editLags_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLags (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.lags = 0;
% Update handles structure
guidata(hObject, handles)


function editTolerance_Callback(hObject, eventdata, handles)
% hObject    handle to editTolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTolerance as text
%        str2double(get(hObject,'String')) returns contents of editTolerance as a double

handles.tolerance = str2double(get(hObject,'String'));     
if handles.tolerance <0, errordlg('Tolerance cannot be lower than 0.'); end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editTolerance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.tolerance = 1e-10;
set(hObject,'String',num2str(handles.tolerance));
% Update handles structure
guidata(hObject, handles)


function editIterations_Callback(hObject, eventdata, handles)
% hObject    handle to editIterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editIterations as text
%        str2double(get(hObject,'String')) returns contents of editIterations as a double

handles.iterations = str2double(get(hObject,'String'));     
if handles.tolerance <1, errordlg('The number of maximum iterations cannot be lower than 1.'); end
set(handles.editIterations,'String',num2str(handles.iterations));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editIterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editIterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.iterations = 5000;
set(hObject,'String',num2str(handles.iterations));
% Update handles structure
guidata(hObject, handles)


function editPCs_Callback(hObject, eventdata, handles)
% hObject    handle to editPCs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPCs as text
%        str2double(get(hObject,'String')) returns contents of editPCs as a double

handles.pcs = str2double(get(hObject,'String'));     
if handles.pcs <0, errordlg('The number of PCs cannot be lower than 1.'); end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editPCs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPCs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.pcs = 1;
set(hObject,'String',num2str(handles.pcs));
% Update handles structure
guidata(hObject, handles)


% --- Executes on button press in pushbuttonPCs.
function pushbuttonPCs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPCs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

batchID = handles.md(handles.selectedBatch);
cal = handles.dataset{batchID};

s = size(cal);
if handles.lags>s(1)-1, errordlg('The number of lags introduced for the selected batch is not feasible');end

% Disable the rest of actions and parameters
set(handles.pushbuttonImpute,'Enable','off');
set(handles.pushbuttonApply,'Enable','off');
% Disable option to save imputed data into data set
set(handles.pushbuttonSave,'Enable','off');
% Disable option to impute missing data
set(handles.pushbuttonImpute,'Enable','off');

try 
    % Lag the two-way array
    lagcal = lagmatrix(cal,handles.lags);
    % Estimate the pair-wise variance-covariance matrix
    var_cov=S_pairwise(lagcal);
    % Adjust a PCA model via singular value decomposition
    [U,S,V] = svd(var_cov);
    % Retrieve the eigenvalue
    eig_values = diag(S);
    % Estimate the optimum number of PCs
    cumpress = ckf(var_cov,U*S,V,0);

    nbatchid = handles.ParentFigure.BatcheslbIn(batchID);

    % Display the figures of merit
    axes(handles.axes3)
    plot(eig_values,'b.-','MarkerSize',16,'LineWidth',2);
    xlabel('Number of Components','FontSize',12);
    ylabel('Eigenvalues','FontSize',12);
    title(sprintf('Batch #%d',nbatchid),'FontSize',12);
    axes(handles.axes1)
    plot(cumpress,'b.-','MarkerSize',16,'LineWidth',2);
    xlabel('Number of Components','FontSize',12);
    ylabel('PRESS (ckf)','FontSize',12);
    
    % Suggest number of PCs
    [~,handles.pcs] = min(cumpress);
    set(handles.editPCs,'String',num2str(handles.pcs));
catch err
    errordlg(err.message);
    return
end

% Enable imputation
set(handles.pushbuttonImpute,'Enable','on');


% Update parameters
handles.lagcal = lagcal;
handles.cal = cal;

% Update handles structure
guidata(hObject, handles)


% --- Executes on button press in pushbuttonImpute.
function pushbuttonImpute_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonImpute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable option to save imputed data into data set
set(handles.pushbuttonSave,'Enable','off');

if ~get(handles.checkboxAllImputation,'Value')
    try 
        batchID = handles.md(handles.selectedBatch);
        cal = handles.dataset{batchID};
        handles.cal = cal;
        handles.lagcal = lagmatrix(cal,handles.lags);
        s = size(handles.lagcal);
        % Reconstruct the lagged two way array using the PCA
        % modeling building procedure to impute missing data
        if handles.pcs>s(2) || handles.pcs>s(2)>s(1), errordlg('The number of PCS exceeds the rank of the matrix'); return; end
        % Reconstruct the lagged two way array using the PCA
        % modeling building procedure to impute missing data
        lagrec=pcambtsr(handles.lagcal,handles.pcs,handles.iterations,handles.tolerance);  
        % Reconstruct the original two-way array with the imputed
        % values keeping the original dimensions.
        handles.rec = reclagmatrix(lagrec,handles.lags);
        % Plot the original and imputed two-way arrays
        plot3D(handles.rec,[],handles.cal)
        % Enable option to save imputed data into data set
        set(handles.pushbuttonSave,'Enable','on');
        % Update handles structure
        guidata(hObject, handles)
    catch err
        errordlg(err.message);
        return;
    end
else
    nbatches = length(handles.md);
    % Populate the listbox of batches with missing values
    if nbatches==0, warndlg('All the batches with missing values have already been computed. Please, close this user interface and continue with the bilinear modeling.');end
    prevSet = handles.dataset;
    for g=1:nbatches
        i = handles.md(1);
        cal = handles.dataset{i};
        % Lag the two-way array
        lagcal = lagmatrix(cal,handles.lags);
        s=size(lagcal);
        % Reconstruct the lagged two way array using the PCA
        % modeling building procedure to impute missing data
        if handles.pcs>s(2) || handles.pcs>s(2)>s(1), errordlg('The number of PCS exceeds the rank of the matrix'); return; end
        lagrec=pcambtsr(lagcal,handles.pcs,handles.iterations,handles.tolerance);  
        % Reconstruct the original two-way array with the imputed
        % values keeping the original dimensions.
        rec = reclagmatrix(lagrec,handles.lags);
        % Save imputed data set
        handles.dataset{i} = rec;
        
        % Retrieve the batchID from the listbox
        contentsImputed = cellstr(get(handles.listboxImputedbatches,'String'));
        % Retrieve the list of available batch IDs
        current_batchid = handles.ParentFigure.BatcheslbIn(handles.md);
        if isempty(contentsImputed{1,1})
            set(handles.listboxImputedbatches,'String',['' num2str(current_batchid(1))]);
        else
            set(handles.listboxImputedbatches,'String',[contentsImputed; num2str(current_batchid(1))]);
        end
        current_batchid(1) = [];
        handles.md(1) = [];
        set(handles.listboxMDbatches,'String',arrayfun(@num2str,current_batchid,'unif',0));
        % Update handles structure
        guidata(hObject, handles)
    end
    % Enable option to save imputed data into data set
    set(handles.pushbuttonApply,'Enable','on');
    % Enable option to reset everything
    set(handles.pushbuttonReset,'Enable','on');
    % Plot all the imputed batches
    nbatches = length(handles.mdmaster);
    plot3D(handles.dataset,[],prevSet);
end

% Update handles structure
guidata(hObject, handles)
    

% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Replace the batch trajectories with the imputed ones
batchID = handles.md(handles.selectedBatch);
handles.dataset{batchID} = handles.rec;

handles.rec = [];

% clean axes
axes(handles.axes3)
cla
set(handles.axes3,'Visible','off');
axes(handles.axes1)
cla
set(handles.axes1,'Visible','off');


% Retrieve the batchID from the listbox
contentsImputed = cellstr(get(handles.listboxImputedbatches,'String'));
% Retrieve the list of available batch IDs
current_batchid = handles.ParentFigure.BatcheslbIn(handles.md);

if isempty(contentsImputed{1,1})
    set(handles.listboxImputedbatches,'String',['' num2str(current_batchid(handles.selectedBatch))]);
else
    set(handles.listboxImputedbatches,'String',[contentsImputed; num2str(current_batchid(handles.selectedBatch))]);
end
current_batchid(handles.selectedBatch) = [];
handles.md(handles.selectedBatch) = [];
set(handles.listboxMDbatches,'String',arrayfun(@num2str,current_batchid,'unif',0));

% Update GUI
handles.selectedBatch = 1;
set(handles.listboxMDbatches,'Value',1);
set(handles.pushbuttonImpute,'Enable','on');
set(handles.pushbuttonSave,'Enable','off');
set(handles.pushbuttonPCs,'Enable','on');

if isempty(handles.md)
    % Enable option to save imputed data into data set
    set(handles.pushbuttonApply,'Enable','on');
    % Enable option to reset everything
    set(handles.pushbuttonReset,'Enable','on');
end

% Update handles structure
guidata(hObject, handles)

% --- Executes on button press in pushbuttonApply.
function pushbuttonApply_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonApply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% If there are not more elements in the MD list, enable the botton Close
% and transfer data to the main user interface

if isfield(handles.ParentFigure,'monitoringFlag')
    handles.ParentFigure.xtest = handles.dataset;
    handles.ParentFigure.calibration.test = [handles.ParentFigure.calibration.test; {handles.dataset}];
    handles.ParentFigure.calibration.test_batch = handles.dataset{1};
    handles.ParentFigure.OK = [];
    
    contents = strvcat(get(handles.ParentFigure.popupmenuVar,'String'),handles.ParentFigure.filename);
    set(handles.ParentFigure.popupmenuVar,'String',contents);
    set(handles.ParentFigure.popupmenuVar,'Value',size(get(handles.ParentFigure.popupmenuVar,'String'),1));
    set(handles.ParentFigure.textData,'String',contents(size(contents,1),:));
    
    handles.ParentFigure.selectedDataSet = size(contents,1);
    
    set(handles.ParentFigure.popupmenuBat,'String','');
    for i=1:length(handles.dataset)
        contents = get(handles.ParentFigure.popupmenuBat,'String');
        set(handles.ParentFigure.popupmenuBat,'String',strvcat(contents,[' ',num2str(i)]));
    end
    
    % Create a new instance of the Monitoring Parameters class object in each
    % model
    
    nmodels = length(handles.ParentFigure.calibration.LVmodels);
    
    for n=1:nmodels
        ntestsets = length(handles.ParentFigure.calibration.LVmodels(n).latent_structure.monitoring_statistics);
        handles.ParentFigure.calibration.LVmodels(n).latent_structure.monitoring_statistics(ntestsets+1) = LatentStructure.MonitoringParameters(handles.ParentFigure.filename);
    end
    
    
    % Update the user interface
    set(handles.ParentFigure.popupmenuBat,'Value',1);
    set(handles.ParentFigure.popupmenuVar,'Enable','on');
    set(handles.ParentFigure.textVar,'Enable','on');
    set(handles.ParentFigure.textData,'Enable','on');
    set(handles.ParentFigure.popupmenuBat,'Enable','off');
    set(handles.ParentFigure.textBat,'Enable','on');
    set(handles.ParentFigure.pushbuttonPlo,'Enable','off');
    
    % Save the information of the parent handle    
    guidata(handles.ParentsWindow,handles.ParentFigure)
    
else
    
    handles.ParentFigure.ParentFigure.x = handles.dataset;

    handles.ParentFigure.ParentFigure.track(2) = 1;
    handles.ParentFigure.ParentFigure.track(3:end) = 0;

    % Enable the following bilinear modeling step
    axes(handles.ParentFigure.ParentFigure.main_window)
    image(handles.ParentFigure.ParentFigure.images{3});
    axis off;
    axis image;
    set(handles.ParentFigure.ParentFigure.pbAlignment,'Enable','on');
    
    % Save the information of the parent handle    
    guidata(handles.ParentFigure.ParentsWindow,handles.ParentFigure.ParentFigure)
    
end



delete(handles.figure1);



% --- Executes on button press in pushbuttonPlotMD.
function pushbuttonPlotMD_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPlotMD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonPlotImputed.
function pushbuttonPlotImputed_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPlotImputed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkboxAllImputation.
function checkboxAllImputation_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAllImputation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxAllImputation

if get(hObject,'Value')
    % Disable the rest of actions and parameters
    set(handles.pushbuttonImpute,'Enable','on');
    set(handles.pushbuttonPCs,'Enable','off');
    set(handles.pushbuttonSave,'Enable','off');
    set(handles.pushbuttonApply,'Enable','off');
else
    set(handles.pushbuttonImpute,'Enable','on');
    set(handles.pushbuttonPCs,'Enable','on');
end

function handles = Initialize(handles)

batchid = handles.ParentFigure.BatcheslbIn(handles.mdmaster);
set(handles.listboxMDbatches,'String',arrayfun(@num2str,batchid,'unif',0));
set(handles.listboxImputedbatches,'String','');

% Default parameters in GUI
set(handles.editLags,'String','0');
handles.lags = 0;
set(handles.editIterations,'String','5000');
handles.iterations = 5000;
set(handles.editPCs,'String','1');
handles.pcs = 1;
set(handles.editTolerance,'String','1e-10');
handles.tolerance = 1e-10;

% Extract data set to impute
if isfield(handles.ParentFigure,'monitoringFlag')
    handles.dataset = handles.ParentFigure.xtest;
else
    handles.dataset = handles.ParentFigure.ParentFigure.x;
end

handles.lagcal = [];
handles.rec = [];
handles.selectedBatch = 1;
handles.md = handles.mdmaster;

% --- Executes on button press in pushbuttonReset.
function pushbuttonReset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = Initialize(handles);

set(handles.pushbuttonImpute,'Enable','on');
set(handles.pushbuttonPCs,'Enable','off');
set(handles.pushbuttonSave,'Enable','off');
set(handles.pushbuttonReset,'Enable','off');
set(handles.pushbuttonApply,'Enable','off');


% Update handles structure
guidata(hObject, handles)

