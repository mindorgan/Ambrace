function varargout = Ambrace(varargin)
% AMBRACE MATLAB code to find the best unnatural amino acid incorporation
% position.

%   Ambrace calculates incorporation scores of all possible DNA sequence
%   for a given amino acid sequence. The result will be sorted and saved 
%   in output directory.
%
%   Scores are calculated based on
%   1. Nine neucleotides including the amber codon: NNNTAGNNN
%       Ref:    Mao et al. - 2004, Miller - 1983, Bossi - 1983, Pedersen, Curran - 1991 
%   2. The last amino acid before the amber codon
%       Ref:    Zhang, Rydén-Aulin, Isaksson - 1996, Mottagui-Tabar, Isaksson - 1997
%   3. The second last amino acid before the amber codon
%       Ref:    Mottagui-tabar, Björnsson, Isaksson - 1994 , Björnsson,
%       Mottagui-Tabar, Isaksson - 1996, Mottagui-Tabar, Isaksson - 1997
%
% See also: AA2_score, AA1_score, Mao_score

% Edit the above text to modify the response to help Ambrace

% Last Modified by GUIDE v2.5 28-Apr-2017 20:58:48
% janghyun.yoo@nih.gov

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Ambrace_OpeningFcn, ...
                   'gui_OutputFcn',  @Ambrace_OutputFcn, ...
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


% --- Executes just before Ambrace is made visible.
function Ambrace_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Ambrace (see VARARGIN)

% Choose default command line output for Ambrace
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Ambrace wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Ambrace_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in push_calc.
function push_calc_Callback(hObject, eventdata, handles)
% hObject    handle to push_calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

if ~isfield(handles,'method')
    handles.method=0;
    guidata(hObject,handles);
end

% Haven't you set the codon table?
if strcmpi(get(handles.codon_text, 'String'),'Codon table')
    codon_select_Callback(hObject, eventdata, handles);
    handles=guidata(hObject);
    %h=errordlg('Select a condon table first');
    %return;
end

% Check the options
skip=get(handles.lastAAcheckbox,'value')+2*get(handles.slastAAcheckbox,'value');
% this is our protein sequence
AA_input = upper(get(handles.edit1,'String'));
if numel(AA_input)<4-handles.method
    h=errordlg(['Input at least ' num2str(4-handles.method) ' AAs']);
    return;
end

% Start scoring
try
    curPath=pwd;
    [codePath,name,ext] = fileparts(mfilename('fullpath'));
    addpath(fullfile(codePath,'Scores')); % add scoring fuctions path
    
    sequences=cell(0,1);
    scores=zeros(0,1);
    for ii=3:numel(AA_input)-1+handles.method
        AA2=AA_input(ii-2); % -2 position
        AA1=AA_input(ii-1); % -1 position
        AA_3=AA_input(ii+1-handles.method); % +1 position
        codons1=handles.codon_table.(AA1);
        codons_3=handles.codon_table.(AA_3);

        for j=1:numel(codons1)
            for k=1:numel(codons_3)
                nnSequence=[codons1{j} 'TAG' codons_3{k}];
                curScore=AA2_score([AA2 AA1],nnSequence,skip); % AA2_score ~ termination probability
                sequences{end+1}=[AA_input(1:ii-2) '...' codons1{j} '(' AA1 ')' ' TAG ' codons_3{k} '(' AA_3 ')' '...' AA_input(ii+1-handles.method+1:end)];
                scores(end+1)=curScore;
            end
        end
    end
    
    % Now sort the sequences by their scores
    [sortedScores, sortIndex]=sort(scores,'descend');
    
    outDir=fullfile(codePath,'Output');
    if exist(outDir,'dir')==0 %dir not exist
        mkdir(outDir);
    end
    cd(outDir);
    
    today_str=strcat(datestr(clock,'yyyy-mm-dd-HHMM'),'m',datestr(clock,'ss'),'s');
    result_file=['Ambrace_' today_str '.txt'];
    fileID=fopen(result_file,'w');
    
    fprintf(fileID,'%s\n',today_str);
    fprintf(fileID,'Ambrace result for input AA sequence:\n%s\n',AA_input);
    if handles.method==0
        methodString='Substiution';
    elseif handles.method==1
        methodString='Insertion';
    else
        methodString='Unknown';
    end
    
    if fix(skip/2)==0
        slastAAString='No';
    elseif fix(skip/2)==1
        slastAAString='Yes';
    else
        slastAAString='Unknown';
    end
    
    if rem(skip,2)==0
        lastAAString='No';
    else
        lastAAString='Yes';
    end
    fprintf(fileID,'Method: %s\nInclude the last AA: %s\nInclude the second last AA: %s\n\n',methodString,lastAAString,slastAAString);
    fprintf(fileID,'Rank\tScore\tSequence\n');
    for ii=1:numel(sortedScores)
        fprintf(fileID,'%d\t%f\t%s\n',ii,sortedScores(ii),sequences{sortIndex(ii)});
    end
    fclose(fileID);
    winopen(result_file);
    cd(curPath);
catch ME
    if strcmp(ME.identifier,'MATLAB:nonExistentField')
        h=errordlg(['Bad AA code; ' AA1 ' or ' AA_3 '?']);
        msg ='Unknown amino acid in the input string.';
        causeException = MException('MATLAB:nonExistentField',msg);
        ME = addCause(ME,causeException);
    end
    rethrow(ME)
end
guidata(hObject,handles);

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in button_subs.
function button_subs_Callback(hObject, eventdata, handles)
% hObject    handle to button_subs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
handles.method=0;
set(handles.uipanel1, 'Title', ['Amino acid sequence (single character notation, ' num2str(4-handles.method) ' AAs for minimum)']);
guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of button_subs


% --- Executes on button press in button_ins.
function button_ins_Callback(hObject, eventdata, handles)
% hObject    handle to button_ins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
handles.method=1;
set(handles.uipanel1, 'Title', ['Amino acid sequence (single character notation, ' num2str(4-handles.method) ' AAs for minimum)']);
guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of button_ins


% --- Executes on button press in codon_select.
function codon_select_Callback(hObject, eventdata, handles)
% hObject    handle to codon_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
curPath=pwd;
[codePath,name,ext] = fileparts(mfilename('fullpath'));
cd(fullfile(codePath,'Codon'));

[codon_file, foo, opened] = uigetfile('*.txt');

if opened==0
    h=errordlg('You must select a condon table');
    cd(curPath);
    return;
end

set(handles.codon_text, 'String', codon_file);
fileID=fopen(codon_file);

tline = fgetl(fileID);

assign_begin=0; % you can assign codons to amino acids

handles.codon_table={};
while ischar(tline)
    tline=strtrim(tline); % remove white spaces
    
    if numel(tline)<1
        tline = fgetl(fileID);
        continue;
    end
    
    if strcmp(tline(1),'/')
        tline = fgetl(fileID);
        continue;        
    end
    
    %Assign codons
    if isstrprop(tline(1),'digit')
        for ii=2:numel(tline)
            AA=tline(ii);
            if isstrprop(AA,'alpha')
                assign_begin=1;
                handles.codon_table.(AA)={};
                break;
            end
        end
    elseif isstrprop(tline(1),'alpha') && numel(tline)==3
        if assign_begin==1 % you know where to add these codons
            handles.codon_table.(AA){end+1}=tline;
        end
    end
    
    tline = fgetl(fileID);
end

fclose(fileID);
cd(curPath);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function codon_select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to push_calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function push_calc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to push_calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in lastAAcheckbox.
function lastAAcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to lastAAcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lastAAcheckbox


% --- Executes on button press in slastAAcheckbox.
function slastAAcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to slastAAcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of slastAAcheckbox


% --- Executes on button press in monobutton.
function monobutton_Callback(hObject, eventdata, handles)
% hObject    handle to monobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of monobutton


% --- Executes on button press in dibutton.
function dibutton_Callback(hObject, eventdata, handles)
% hObject    handle to dibutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dibutton
