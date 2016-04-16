function varargout = Diagnosis(varargin)
% DIAGNOSIS MATLAB code for Diagnosis.fig
%      DIAGNOSIS, by itself, creates a new DIAGNOSIS or raises the existing
%      singleton*.
%
%      H = DIAGNOSIS returns the handle to a new DIAGNOSIS or the handle to
%      the existing singleton*.
%
%      DIAGNOSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIAGNOSIS.M with the given input arguments.
%
%      DIAGNOSIS('Property','Value',...) creates a new DIAGNOSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Diagnosis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Diagnosis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Diagnosis

% Last Modified by GUIDE v2.5 23-Sep-2014 17:52:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Diagnosis_OpeningFcn, ...
                   'gui_OutputFcn',  @Diagnosis_OutputFcn, ...
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


% --- Executes just before Diagnosis is made visible.
function Diagnosis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Diagnosis (see VARARGIN)

% Choose default command line output for contribution
handles.output = hObject;

if length(varargin) ~= 9, error('The number of input parameters are not expected.'); end

handles.x = varargin{1};
handles.test = varargin{2};
handles.goldenI = varargin{3};
handles.contrb = varargin{4};
handles.contribJ = varargin{5};
handles.contribK = varargin{6};
handles.phases = varargin{7};
handles.varNames = varargin{8};
handles.label = varargin{9};

% Plot the contributions of the corresponding statistic
plotcontrbGUI(handles)

handles.selectedVariable = 1;

Variable_popmenu_Callback(hObject, eventdata, handles)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes contribution wait for user response (see UIRESUME)
% uiwait(handles.contribution);

% --- Outputs from this function are returned to the command line.
function varargout = Diagnosis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Variable_popmenu.
function Variable_popmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Variable_popmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Variable_popmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Variable_popmenu

handles.selectedVariable  = get(handles.Variable_popmenu,'Value');
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Variable_popmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Variable_popmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Display_pushbutton.
function Display_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Display_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axesTrajectories)
cla
for i = 1:size(handles.x,3)
    p1=plot(handles.x(:,handles.selectedVariable,i),'Color',[0.392157 0.584314 0.929412],'LineWidth',1); hold on;
end
   p2=plot(handles.test(:,handles.selectedVariable),'LineStyle','-','Color',[1 0 0],'LineWidth',2);

p3=plot(handles.goldenI(:,handles.selectedVariable),'Color','y','LineWidth',2);  
xlabel('Batch Time','FontSize',12,'FontWeight','bold');
ylabel(handles.varNames{handles.selectedVariable},'FontSize',12,'FontWeight','bold');
legend([p1 p2 p3],'calibration batches','test batch','mean trajectory');
v=axis;
if not(isempty(handles.phases))
   for p=1:numel(handles.phases)
       line([handles.phases(p) handles.phases(p)],[v(3) v(4)],'LineStyle','--','Color','w');
   end
end
axis tight;
hold off;


% --- Executes on button press in Cancel_pushbutton.
function Cancel_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1);


function plotcontrbGUI(handles)

s = size(handles.x);

axes(handles.axesOverallContrib)


bar(handles.contrb,'FaceColor',[0 0 1]); hold on;

for j=1:s(2)
    line([j*s(1) j*s(1)],[min(handles.contrb) max(handles.contrb)],'LineStyle','--','Color','k');
end

xlabel('Variables','FontSize',10,'Color','k','fontweight','b');
ylabel(handles.label,'FontSize',10,'Color','k','fontweight','b');
axis tight


axes(handles.axesVarContrib);
bar(handles.contribJ,'FaceColor',[0 0 1]); axis tight;
xlabel('Original Process Variables','FontSize',10,'Color','k','fontweight','b');
ylabel([handles.label ' per variable'],'FontSize',10,'Color','k','fontweight','b');

axes(handles.axesTimeContrib);
bar(handles.contribK,'FaceColor',[0 0 1]);axis tight;
if not(isempty(handles.phases))
   for p=1:numel(handles.phases)
       line([handles.phases(p) handles.phases(p)],[min(handles.contribK) max(handles.contribK)],'LineStyle','--','Color','k');
   end
end
xlabel('Batch Time','FontSize',10,'Color','k','fontweight','b');
ylabel([handles.label ' per time'],'FontSize',10,'Color','k','fontweight','b');


% Fill the popmenu with tagnames
for i=1:length(handles.varNames)
    contents = get(handles.Variable_popmenu,'String');
    set(handles.Variable_popmenu,'String',strvcat(contents,strcat(strcat(num2str(i),'.- '),handles.varNames{i,1})));
end
