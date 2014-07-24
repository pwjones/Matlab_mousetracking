function varargout = reviewTracking(varargin)
% REVIEWTRACKING MATLAB code for reviewTracking.fig
%      REVIEWTRACKING, by itself, creates a new REVIEWTRACKING or raises the existing
%      singleton*.
%
%      H = REVIEWTRACKING returns the handle to a new REVIEWTRACKING or the handle to
%      the existing singleton*.
%
%      REVIEWTRACKING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REVIEWTRACKING.M with the given input arguments.
%
%      REVIEWTRACKING('Property','Value',...) creates a new REVIEWTRACKING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before reviewTracking_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to reviewTracking_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help reviewTracking

% Last Modified by GUIDE v2.5 23-Jul-2014 16:58:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @reviewTracking_OpeningFcn, ...
                   'gui_OutputFcn',  @reviewTracking_OutputFcn, ...
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


% --- Executes just before reviewTracking is made visible.
function reviewTracking_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to reviewTracking (see VARARGIN)

% Choose default command line output for reviewTracking
handles.output = hObject;

% load the MouseTracker object to edit
if ~isempty(varargin)
    handles.tracker = varargin{1}; %the MouseTracker object 
    mt = handles.tracker;
else
    return;
end
% get the frames to zero in on
frames = 1:mt.nFrames;
handles.distThresh = 10;
handles.jumpFrames = handles.tracker.findNoseJumps(handles.distThresh, frames);
if (isempty(handles.jumpFrames))
    handles.jumpFrames = 1;
end
handles.currFrame = handles.jumpFrames(1);
handles.jumpFrameI = 1;
% Now tend to the displayed UI
handles.dispTypeList = {'diff', 'bin', 'orig'};
handles.dispType = 'diff';
handles.anatomyFlag = 1;
handles = changeFrame(hObject, handles, handles.currFrame);
set(hObject, 'Toolbar', 'figure');
set(hObject, 'KeyPressFcn', @(hObject, evt)keyResponder(hObject, evt));
% The propogation section
handles.corrThresh = .6;
if isempty(mt.blobID) % This is in case there are no blobIDs.  Necessary for easy propogation of changes
    mt.blob_num = 1;
    for ii=1:mt.nFrames
        mt.assignBlobIDs(ii);
    end
end
     
guidata(hObject, handles); % Update handles structure
% UIWAIT makes reviewTracking wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function keyResponder(hObject, evt)
% Makes sense of the key presses happening for the GUI
handles = guidata(hObject);
cf = handles.currFrame;
switch evt.Key
    case 'rightarrow'
        handles = changeFrame(hObject, handles, cf + 1);
    case 'leftarrow'
        handles = changeFrame(hObject, handles, cf - 1);
    case 'downarrow'
        nf = find(handles.jumpFrames > cf, 1, 'first');
        if isempty(nf)
            nf = handles.tracker.nFrames;
        end
        handles = changeFrame(hObject, handles, nf);
    case 'uparrow'
        nf = find(handles.jumpFrames < cf, 1, 'last');
        if isempty(nf)
            nf = 1;
        end
        handles = changeFrame(hObject, handles, nf);
    case 'slash'
        propogate_btn_Callback(hObject, evt, handles);
        handles = changeFrame(hObject, handles, cf);
end
guidata(hObject, handles);


% ---- changes the active Frame -------------------
function handles = changeFrame(hObject, handles, frame)
%
% Just what needs to happen to load and display a new frame
if frame < 1 % a little error checking so that calling functions don't have to worry
    frame = 1;
elseif frame > handles.tracker.nFrames
    frame = handles.tracker.nFrames;
end
handles.currFrame = frame;
axes(handles.im_ax);
cla(handles.im_ax); %clears the axis
handles.tracker.showFrame(handles.currFrame, handles.dispType, []);
line('Parent', handles.im_ax, 'Xdata', handles.tracker.nosePos(:,1), 'Ydata', handles.tracker.nosePos(:,2), 'Marker', 'none', 'Color', 'c');
% In order to register clicks, must set the ButtonDownFcn of the objects plotted
child_h = get(gca, 'Children');
for ii = 1:length(child_h)
    set(child_h, 'ButtonDownFcn', @changeAnatomyID);
end
set(handles.currFrame_edit, 'String', num2str(handles.currFrame));

%guidata(hObject, handles);
    

% ---Changes what is considered the nose or tail on click -------
function changeAnatomyID(hObject, eventdata, handles)
% Used as buttonDownFcn in order to set new nose or tail 
% positions on clicks in the main axis
%disp('In changeAnatomyID');
handles = guidata(hObject);
cp = get(gca, 'CurrentPoint');
clickP = [cp(1,1) cp(1,2)]

switch (handles.anatomyFlag)
    case(1) % the nose
        handles.tracker.setNosePosition(handles.currFrame, clickP);
    case(2)
        handles.tracker.setTailPosition(handles.currFrame, clickP);
end
hObject = handles.figure1;
handles = changeFrame(hObject, handles, handles.currFrame);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = reviewTracking_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.tracker;



function currFrame_edit_Callback(hObject, eventdata, handles)
% hObject    handle to currFrame_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currFrame_edit as text
%        str2double(get(hObject,'String')) returns contents of currFrame_edit as a double
newFrame = str2double(get(hObject,'String'));
handles = changeFrame(hObject, handles, newFrame);
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function currFrame_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currFrame_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%set(hObject, 'String', num2str(handles.currFrame));


% --- Executes on button press in prev_frame_btn.
function prev_frame_btn_Callback(hObject, eventdata, handles)
% hObject    handle to prev_frame_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newFrame = max(handles.currFrame-1, 1);
handles = changeFrame(hObject, handles, newFrame);
guidata(hObject, handles);

% --- Executes on button press in next_frame_btn.
function next_frame_btn_Callback(hObject, eventdata, handles)
% hObject    handle to next_frame_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newFrame = min(handles.currFrame+1, handles.tracker.nFrames);
handles = changeFrame(hObject, handles, newFrame);
guidata(hObject, handles);

% --- Executes on button press in adv10_btn.
function adv10_btn_Callback(hObject, eventdata, handles)
% hObject    handle to adv10_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newFrame = min(handles.currFrame+10, handles.tracker.nFrames);
handles = changeFrame(hObject, handles, newFrame);
guidata(hObject, handles);


% --- Executes on button press in back10_btn.
function back10_btn_Callback(hObject, eventdata, handles)
% hObject    handle to back10_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newFrame = max(handles.currFrame-10, 1);
handles = changeFrame(hObject, handles, newFrame);
guidata(hObject, handles);


% --- Executes on button press in next_jump_btn.
function next_jump_btn_Callback(hObject, eventdata, handles)
% hObject    handle to next_jump_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.jumpFrameI == length(handles.jumpFrames)
    nexti = handles.jumpFrameI;
else
    nexti = find(handles.jumpFrames > handles.currFrame, 1, 'first');
    if isempty(nexti)
        nexti = handles.jumpFrameI + 1;
    end
end
handles.jumpFrameI = nexti;
handles = changeFrame(hObject, handles, handles.jumpFrames(nexti));
guidata(hObject, handles);


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
contents = cellstr(get(hObject,'String'));
handles.dispType = contents{get(hObject,'Value')};
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.dispTypeList = {'diff', 'bin', 'orig'};
set(hObject, 'String', handles.dispTypeList);
set(hObject, 'Value', 1);

% --- Executes on selection change in anatomy_menu.
function anatomy_menu_Callback(hObject, eventdata, handles)
% hObject    handle to anatomy_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns anatomy_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from anatomy_menu
contents = cellstr(get(hObject,'String'));
handles.anatomyFlag = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function anatomy_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to anatomy_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', {'nose', 'tail'}, 'Value', 1);


% --- Executes on button press in save_btn.
function save_btn_Callback(hObject, eventdata, handles)
% hObject    handle to save_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.tracker.save();


% --- Executes on button press in find_jumps_btn.
function find_jumps_btn_Callback(hObject, eventdata, handles)
% hObject    handle to find_jumps_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in propogate_btn.
function propogate_btn_Callback(hObject, eventdata, handles)
% hObject    handle to propogate_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currFrame = handles.currFrame;
handles.tracker.propogateNosePosition(currFrame);
% noseID = mt.blobID(currFrame, mt.noseblobs(currFrame));
% % propogate forward
% nf = currFrame + 1; prop = 1;
% while (nf <= mt.nFrames) && prop
%     ids = mt.blobID(nf,:);
%     match = find(noseID == ids);
%     if ~isempty(match)
%         mt.noseblob(nf) = match;
%         mt.nosePos(nf,:) = mt.areas(nf,match).Centroid;
%         prop = 1;
%     else
%         prop = 0;
%     end
% end



function prop_thresh_edit_Callback(hObject, eventdata, handles)
% hObject    handle to prop_thresh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prop_thresh_edit as text
%        str2double(get(hObject,'String')) returns contents of prop_thresh_edit as a double


% --- Executes during object creation, after setting all properties.
function prop_thresh_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prop_thresh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in followingGaps_btn.
function followingGaps_btn_Callback(hObject, eventdata, handles)
% hObject    handle to followingGaps_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

paths = cat(1, handles.tracker.paths(1).PixelList, handles.tracker.paths(2).PixelList);
bodyCOM = handles.tracker.bodyCOM;
dm = ipdm(paths,bodyCOM);
min_dist = min(dm);
close = min_dist < 30;
missing = isnan(handles.tracker.nosePos(:,1));
handles.missing_nose_inds = find(close(:) & missing(:));
handles.missing_i = 1;
if ~isempty(handles.missing_nose_inds)
    handles = changeFrame(hObject, handles, handles.missing_nose_inds(handles.missing_i));
    guidata(hObject, handles);
end

% --- Executes on button press in advance_gap_btn.
function advance_gap_btn_Callback(hObject, eventdata, handles)
% hObject    handle to advance_gap_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'missing_nose_inds') && ~isempty(handles.missing_nose_inds)
    i = mod(handles.missing_i, length(handles.missing_nose_inds)-1) + 1;
    newFrame = handles.missing_nose_inds(i);
    handles.missing_i = i;
    handles = changeFrame(hObject, handles, newFrame);
    guidata(hObject, handles);
end