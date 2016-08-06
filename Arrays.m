function varargout = Arrays(varargin)
% ARRAYS MATLAB code for Arrays.fig
% Last Modified by GUIDE v2.5 11-Apr-2016 22:14:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Arrays_OpeningFcn, ...
                   'gui_OutputFcn',  @Arrays_OutputFcn, ...
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


% --- Executes just before Arrays is made visible.
function Arrays_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
%hide warnings of existing folder
warning ('off','all');

%Initialize to disabled.
set(handles.edit_elements, 'Enable', 'Off');
set(handles.edit_distance, 'Enable', 'Off');
set(handles.edit_waveLength, 'Enable', 'Off');
set(handles.edit_thetaD,'Enable','Off');


elements = str2double(get(handles.edit_elements,'string'));
handles.elements = elements;
distance = str2double(get(handles.edit_distance,'string'));
handles.distance = distance;
waveLength = str2double(get(handles.edit_waveLength,'string'));
handles.waveLength = waveLength;
thetaD = str2double(get(handles.edit_thetaD,'string'));
handles.thetaD = thetaD;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Arrays wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = Arrays_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%---------------------------------------------------------------
%---------------------------------------------------------------
%-----------------  >> MY FUNCTIONS <<    ----------------------
%---------------------------------------------------------------
%---------------------------------------------------------------
function enable(hObject, eventdata, handles)
handles.output = hObject;
set(handles.edit_elements, 'Enable', 'On');
set(handles.edit_distance, 'Enable', 'On');
set(handles.edit_waveLength, 'Enable', 'On');
guidata(hObject,handles);

function[ArrayFactor] =  Array_Factor_Cartesian(elements, psi)
numerator = sin((elements.*psi)/2);
denominator = sin(1/2.*psi);
ArrayFactor = (1/elements).*(numerator./denominator);

function printToFile(file,rows,info0,info1,info2,info3,info4,info5,info6)
fprintf(file,rows);
fprintf(file,'\n');
fprintf(file,info0);
fprintf(file,info1);
fprintf(file,info2);
fprintf(file,info3);
fprintf(file,info4);
fprintf(file,info5);
fprintf(file,info6);
fclose(file);

function printParToFile(file, elements, distance, wave)
fprintf(file,'\n');
fprintf(file,'Elements : ');
fprintf(file, elements);
fprintf(file,'\n');
fprintf(file,'Distance : ');
fprintf(file, distance);
fprintf(file,'\n');
fprintf(file,'Wave Length : ');
fprintf(file, wave);


function createFile(hObject, eventdata, handles,DateNow,rows,info0,info1,info2,info3,info4,info5,info6)
handles.output = hObject;
handles.elements
Radiation_Pattern
if strcmp(handles.currentData, 'Broadside Array')
          type = 'BroadsideArray _ ';
          t = datetime('now');
          time = datestr(t,'HH;MM;SS PM');
          fileName = [type  time '.txt'];
          file = fopen(fileName,'wt');
          fprintf(file,DateNow);
          fprintf(file,'\n');
          fprintf(file,handles.currentData);
          fprintf(file,'\n');
          fprintf(file,rows);
          printParToFile(file,num2str(handles.elements), num2str(handles.distance),...
                    num2str(handles.waveLength));
          fprintf(file,'\n');

          printToFile(file,rows,info0,info1,info2,info3,info4,info5,info6)
          cd ..

elseif strcmp(handles.currentData, 'End-Fire Array')
          type = 'End-Fire Array _ ';
          t = datetime('now');
          time = datestr(t,'HH;MM;SS PM');
          fileName = [type  time '.txt'];
          file = fopen(fileName,'wt');
          fprintf(file,DateNow);
          fprintf(file,'\n');
          fprintf(file,handles.currentData);
          fprintf(file,'\n');
          fprintf(file,rows);
          printParToFile(file,num2str(handles.elements), num2str(handles.distance),...
                    num2str(handles.waveLength));
          fprintf(file,'\n');

          printToFile(file,rows,info0,info1,info2,info3,info4,info5,info6)
          cd ..

elseif strcmp(handles.currentData, 'Hansen Woodyard')
          type = 'Hansen Woodyard _ ';
          t = datetime('now');
          time = datestr(t,'HH;MM;SS PM');
          fileName = [type  time '.txt'];
          file = fopen(fileName,'wt');
          fprintf(file,DateNow);
          fprintf(file,'\n');
          fprintf(file,handles.currentData);
          fprintf(file,'\n');
          fprintf(file,rows);
          printParToFile(file,num2str(handles.elements), num2str(handles.distance),...
                    num2str(handles.waveLength));
          fprintf(file,'\n');

          printToFile(file,rows,info0,info1,info2,info3,info4,info5,info6)
          cd ..

elseif strcmp(handles.currentData, 'Phased Array')
          type = 'Phased Array _ ';
          t = datetime('now');
          time = datestr(t,'HH;MM;SS PM');
          fileName = [type  time '.txt'];

end
guidata(hObject,handles);

function ArrayInfo(hObject, eventdata, handles)
global displayDate
global rows
global info0
global info1
global info2
global info3
global info4
global info5
global info6
handles.output = hObject;
displayDate = cellstr(datetime('now'));
displayDate = num2str(cell2mat(displayDate));
newLine = sprintf('\n');
rows = '*---------------------*';
if strcmp(handles.currentData, 'Broadside Array')
          Length  = (handles.elements-1)*handles.distance;
          Directivity = 2*handles.elements*(handles.distance/handles.waveLength);
          DirectivityDB = 10*log10(Directivity);
          FNBW = 2*(pi/2-acos(handles.waveLength/(handles.elements*handles.distance)));
          FSLBW = 2*(pi/2-acos((3*handles.waveLength)/(2*handles.elements*handles.distance)));
          HPBW = 2*(pi/2-acos((1.392*handles.waveLength)/(pi*handles.elements*handles.distance)));
          SLL = -13.46;
          
          info0 = sprintf(['Array Length : ', num2str(Length),'\n']);
          info1 = sprintf(['Directivity : ', num2str(Directivity),'\n']);
          info2 = sprintf(['Directivity in dB : ', num2str(DirectivityDB) ' dB','\n']);
          info3 = sprintf([ 'FNBW : ', num2str(FNBW), 'rad = ', num2str(rad2deg(FNBW)), 'degrees','\n']);
          info4 = sprintf(['FSLBW : ', num2str(FSLBW), 'rad = ', num2str(rad2deg(FSLBW)), 'degrees','\n']);
          info5 = sprintf(['HPBW : ', num2str(HPBW), 'rad = ', num2str(rad2deg(HPBW)), 'degrees','\n']);
          info6 = sprintf(['SLL : ', num2str(SLL), 'dB','\n']);
          info = [displayDate newLine rows newLine handles.currentData newLine rows newLine...
                    info0 info1 info2 info3...
                    info4 info5 info6];
          set(handles.directivity,'String', info)
          
elseif strcmp(handles.currentData, 'End-Fire Array')
          Length  = (handles.elements-1)*handles.distance;
          Directivity = (4*handles.elements*(handles.distance/handles.waveLength));
          DirectivityDB = 10*log10(Directivity);
          FNBW = 2*acos(1-(handles.waveLength/(handles.elements*handles.distance)));
          FSLBW = 2*acos(1-(3*handles.waveLength)/(2*handles.elements*handles.distance));
          HPBW = 2*acos(1-(1.391*handles.waveLength)/(pi*handles.elements*handles.distance));
          SLL = -13.46;
          
          info0 = sprintf(['Array Length : ', num2str(Length),'\n']);
          info1 = sprintf(['Directivity : ', num2str(Directivity),'\n']);
          info2 = sprintf(['Directivity in dB : ', num2str(DirectivityDB) ' dB','\n']);
          info3 = sprintf([ 'FNBW : ', num2str(FNBW), 'rad = ', num2str(rad2deg(FNBW)), 'degrees','\n']);
          info4 = sprintf(['FSLBW : ', num2str(FSLBW), 'rad = ', num2str(rad2deg(FSLBW)), 'degrees','\n']);
          info5 = sprintf(['HPBW : ', num2str(HPBW), 'rad = ', num2str(rad2deg(HPBW)), 'degrees','\n']);
          info6 = sprintf(['SLL : ', num2str(SLL), 'dB','\n']);
          info = [displayDate newLine rows newLine handles.currentData newLine rows newLine...
                    info0 info1 info2 info3...
                    info4 info5 info6];          
          set(handles.directivity,'String', info)
          
elseif strcmp(handles.currentData, 'Hansen Woodyard')
          Length  = (handles.elements-1)*handles.distance;
          Directivity = 1.805*(4*handles.elements*(handles.distance/handles.waveLength));
          DirectivityDB = 10*log10(Directivity);
          FNBW = 2*acos(1-(handles.waveLength/(2*handles.elements*handles.distance)));
          FSLBW = 2*acos(1-(handles.waveLength)/(handles.elements*handles.distance));
          HPBW = 2*acos(1-(0.1398*handles.waveLength)/(handles.elements*handles.distance));
          SLL = -13.46;

          Length  = (handles.elements-1)*handles.distance;
          Directivity = (4*handles.elements*(handles.distance/handles.waveLength));
          DirectivityDB = 10*log10(Directivity);
          FNBW = 2*acos(1-(handles.waveLength/(handles.elements*handles.distance)));
          FSLBW = 2*acos(1-(3*handles.waveLength)/(2*handles.elements*handles.distance));
          HPBW = 2*acos(1-(1.391*handles.waveLength)/(pi*handles.elements*handles.distance));
          SLL = -13.46;
          
          info0 = sprintf(['Array Length : ', num2str(Length),'\n']);
          info1 = sprintf(['Directivity : ', num2str(Directivity),'\n']);
          info2 = sprintf(['Directivity in dB : ', num2str(DirectivityDB) ' dB','\n']);
          info3 = sprintf([ 'FNBW : ', num2str(FNBW), 'rad = ', num2str(rad2deg(FNBW)), 'degrees','\n']);
          info4 = sprintf(['FLSBW : ', num2str(FSLBW), 'rad = ', num2str(rad2deg(FSLBW)), 'degrees','\n']);
          info5 = sprintf(['HPBW : ', num2str(HPBW), 'rad = ', num2str(rad2deg(HPBW)), 'degrees','\n']);
          info6 = sprintf(['SLL : ', num2str(SLL), 'dB','\n']);
          info = [displayDate newLine rows newLine handles.currentData newLine rows newLine...
                    info0 info1 info2 info3...
                    info4 info5 info6];           
          set(handles.directivity,'String', info)
          
elseif strcmp(handles.currentData, 'Phased Array')
         
          set(handles.directivity,'String', 'NOT READY YET...')
end
guidata(hObject,handles);

function PolarPlot(theta,RadiationPattern,data)
if strcmp(data, 'Broadside Array')
          hold off
          polar(theta,RadiationPattern);
elseif strcmp(data, 'End-Fire Array')
          hold off
          polar(theta,RadiationPattern);
elseif strcmp(data, 'Hansen Woodyard')
          hold off
          polar(theta,RadiationPattern);
elseif strcmp(data, 'Phased Array')
          hold off
          polar(theta,RadiationPattern);
elseif strcmp(data, 'Customized Array')
          hold off
          polar(theta,RadiationPattern);
end
          
function CartetianPlot(theta,RadiationPattern)
          plot(theta,abs(RadiationPattern));
          axis([0 350 0 1])

          grid on
          grid minor

%---------------------------------------------------------------
%---------------------------------------------------------------
%----------------------  << END >>    --------------------------
%---------------------------------------------------------------
%---------------------------------------------------------------

% --- Executes on selection change in SelectElement.
function SelectElement_Callback(hObject, eventdata, handles)
element = get(hObject, 'Value');
switch element
    %select Array Element...
    case 1
        warndlg('Select array please...','Warning')
        disp('No selection...')
    case 2
        disp('User selects -> Isotropic Array') 
        handles.arrayElement = 1;
        enable(hObject, eventdata, handles)
    case 3
        disp('User selects -> ë/2')
        handles.arrayElement = handles.waveLength/2;
        enable(hObject, eventdata, handles)
    case 4
        disp('User selects -> ë/4')
        handles.arrayElement = handles.waveLength/4;
        enable(hObject, eventdata, handles)
    case 5
        disp('User selects -> ë/6')
        handles.arrayElement = handles.waveLength/6;
        enable(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function SelectElement_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in selectArray.
function selectArray_Callback(hObject, eventdata, handles)
contents = get(hObject, 'Value');

switch contents
          %Select Array...
          case 1
                    warndlg('Select array please...','Warning')
                    disp('No selection...')
          %Broadside Array   
          case 2
                    cla(handles.Cartetian);
                    cla(handles.Polar);
                    cla(handles.Three);
                    
                    disp('User selects -> Broadside Array ') 
                    handles.currentData = 'Broadside Array';
                    enable(hObject, eventdata, handles)
                    set(handles.edit_thetaD,'Enable','Off');
                    
                    %set default parameters
                    set(handles.edit_elements,'String','10');
                    set(handles.edit_distance,'String','0.25');
                    set(handles.edit_waveLength,'String','1');
                    
                    handles.output = hObject;
                    elements = str2double(get(handles.edit_elements,'String'));
                    handles.elements = elements;
                    distance = str2double(get(handles.edit_distance,'String'));
                    handles.distance = distance;
                    waveLength = str2double(get(handles.edit_waveLength,'String'));
                    handles.waveLength = waveLength;
                    guidata(hObject,handles);
                    

          %End-Fire Array          
          case 3
                    cla(handles.Cartetian);
                    cla(handles.Polar);
                    cla(handles.Three);
                    
                    disp('User selects -> End-Fire Array ') 
                    handles.currentData = 'End-Fire Array';
                    enable(hObject, eventdata, handles)
                    set(handles.edit_thetaD,'Enable','Off');
                    
                    %set default parameters
                    set(handles.edit_elements,'String','10');
                    set(handles.edit_distance,'String','0.25');
                    set(handles.edit_waveLength,'String','1');
                    
                    elements = str2double(get(handles.edit_elements,'String'));
                    handles.elements = elements;
                    distance = str2double(get(handles.edit_distance,'String'));
                    handles.distance = distance;
                    waveLength = str2double(get(handles.edit_waveLength,'String'));
                    handles.waveLength = waveLength;


          %Hansen Woodyard
          case 4
                    cla(handles.Cartetian);
                    cla(handles.Polar);
                    cla(handles.Three);
                    
                    disp('User selects -> Hansen Woodyard ')
                    handles.currentData = 'Hansen Woodyard';
                    enable(hObject, eventdata, handles)
                    set(handles.edit_thetaD,'Enable','Off');
                    
                    %set default parameters
                    set(handles.edit_elements,'String','15');
                    set(handles.edit_distance,'String','0.25');
                    set(handles.edit_waveLength,'String','1');
                    
                    elements = str2double(get(handles.edit_elements,'String'));
                    handles.elements = elements;
                    distance = str2double(get(handles.edit_distance,'String'));
                    handles.distance = distance;
                    waveLength = str2double(get(handles.edit_waveLength,'String'));
                    handles.waveLength = waveLength;

          %Phased Array          
          case 5
                    cla(handles.Cartetian);
                    cla(handles.Polar);
                    cla(handles.Three);
                    
                    disp('User selects -> Phased Array ')
                    handles.currentData = 'Phased Array';
                    enable(hObject, eventdata, handles)                    
                    set(handles.edit_thetaD,'Enable','On');
                    
                    %set default parameters
                    set(handles.edit_elements,'String','10');
                    set(handles.edit_distance,'String','0.25');
                    set(handles.edit_waveLength,'String','1');
                    set(handles.edit_thetaD,'String','60');
                    
                    elements = str2double(get(handles.edit_elements,'String'));
                    handles.elements = elements;
                    distance = str2double(get(handles.edit_distance,'String'));
                    handles.distance = distance;
                    waveLength = str2double(get(handles.edit_waveLength,'String'));
                    handles.waveLength = waveLength;
                    thetaD = str2double(get(handles.edit_thetaD,'String'));
                    handles.thetaD = thetaD;
          %Customized          
          case 6
                    cla(handles.Cartetian);
                    cla(handles.Polar);
                    cla(handles.Three);
                    
                    disp('User selects -> Customized Array')
                    handles.currentData = 'Customized Array';
                    enable(hObject, eventdata, handles)
                    set(handles.edit_thetaD,'Enable','On');
                    
                    elements = str2double(get(handles.edit_elements,'String'));
                    handles.elements = elements;
                    distance = str2double(get(handles.edit_distance,'String'));
                    handles.distance = distance;

                    waveLength = str2double(get(handles.edit_waveLength,'String'));
                    handles.waveLength = waveLength;
                    thetaD = str2double(get(handles.edit_thetaD,'String'));
                    handles.thetaD = thetaD;
          
end
% Save the handles structure.
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function selectArray_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_elements_Callback(hObject, eventdata, handles)
elements = str2double(get(handles.edit_elements,'string'));
handles.elements = elements;
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_elements_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_distance_Callback(hObject, eventdata, handles)
distance = str2double(get(handles.edit_distance,'string'));
handles.distance = distance;
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_distance_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function edit_waveLength_Callback(hObject, eventdata, handles)
waveLength = str2double(get(handles.edit_waveLength,'string'));
handles.waveLength = waveLength;
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_waveLength_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_thetaD_Callback(hObject, eventdata, handles)
thetaD = str2double(get(handles.edit_thetaD,'string'));
handles.thetaD = thetaD;
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_thetaD_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%*****************************************
%************* > CARTESIAN < *************
%*****************************************
% --- Executes on button press in CartesianDiagram.
function CartesianDiagram_Callback(hObject, eventdata, handles)
handles.output = hObject;
k = 2*pi/handles.waveLength;
theta = 1:360;
if strcmp(handles.currentData, 'Broadside Array')
          psi = k*handles.distance.*cosd(theta)  ;
          handles.arrayFactor = Array_Factor_Cartesian(handles.elements, psi);
          axes(handles.Cartetian)
          AE= handles.arrayElement;
          hold off
          RadiationPattern = AE.*abs(handles.arrayFactor);
          CartetianPlot(theta,RadiationPattern(:));
          hold on
          plot(theta,ones(size(theta)) * AE)
          legend('AF','AE')

elseif strcmp(handles.currentData, 'End-Fire Array')
          psi = k*handles.distance.*cosd(theta) + deg2rad(90);
          handles.arrayFactor = Array_Factor_Cartesian(handles.elements, psi);
          axes(handles.Cartetian)
          AE = handles.arrayElement;
          hold off
          RadiationPattern = AE.*abs(handles.arrayFactor);
          CartetianPlot(theta,RadiationPattern(:));
          hold on
          plot(theta,ones(size(theta)) * AE)
          legend('AF','AE')          
          
elseif strcmp(handles.currentData, 'Hansen Woodyard')
          psi = k*handles.distance.*cosd(theta) + deg2rad(108);
          handles.arrayFactor = Array_Factor_Cartesian(handles.elements, psi);
          axes(handles.Cartetian)
          AE = handles.arrayElement;
          hold off
          RadiationPattern = AE.*abs(handles.arrayFactor);
          RadiationPattern = RadiationPattern/max(RadiationPattern(:));
          plot(theta/pi*180,RadiationPattern);
          temp = max(theta)/pi*180;
          axis([0 temp 0 1])
          grid on
          grid minor
          hold on
          plot(theta/pi*180,ones(size(theta/pi*180)) * AE)
          legend('AF','AE')
          
elseif strcmp(handles.currentData, 'Phased Array')
          psi = k*handles.distance.*(cosd(handles.thetaD)-cosd(theta)) ;
          handles.arrayFactor = Array_Factor_Cartesian(handles.elements, psi);
          axes(handles.Cartetian)
          AE= handles.arrayElement;
          hold off
          RadiationPattern = AE.*abs(handles.arrayFactor);
          CartetianPlot(theta,RadiationPattern(:));
          hold on
          plot(theta,ones(size(theta)) * AE)
          legend('AF','AE')

else
          disp('none..')
end

ArrayInfo(hObject, eventdata, handles)

guidata(hObject, handles);


%*****************************************
%************** > POLLAR < ***************
%*****************************************
% --- Executes on button press in PolarDiagram.
function PolarDiagram_Callback(hObject, eventdata, handles)
handles.output = hObject;
k = 2*pi/handles.waveLength;
theta = 1:360;
if strcmp(handles.currentData, 'Broadside Array')
          psi = k*handles.distance.*cosd(theta)  ;
          handles.arrayFactor = Array_Factor_Cartesian(handles.elements, psi);
          AE= handles.arrayElement;
          theta = linspace(0,2*pi,360);
          axes(handles.Polar)
          RadiationPattern = AE.*abs(handles.arrayFactor);
          PolarPlot(theta,RadiationPattern,handles.currentData);

elseif strcmp(handles.currentData, 'End-Fire Array')
          psi = k*handles.distance.*cosd(theta) + deg2rad(180);
          handles.arrayFactor = Array_Factor_Cartesian(handles.elements, psi);
          AE= handles.arrayElement;
          theta = linspace(0,2*pi,360);
          axes(handles.Polar)
          RadiationPattern = AE.*abs(handles.arrayFactor);
          PolarPlot(theta,RadiationPattern,handles.currentData);
          
elseif strcmp(handles.currentData, 'Hansen Woodyard')
          psi = k*handles.distance.*cosd(theta) + deg2rad(108);
          handles.arrayFactor = Array_Factor_Cartesian(handles.elements, psi);
          AE = handles.arrayElement;
          theta = linspace(0,2*pi,360);
          axes(handles.Polar)
          RadiationPattern = AE.*abs(handles.arrayFactor)/max(abs(handles.arrayFactor));
          PolarPlot(theta,RadiationPattern,handles.currentData)
          
elseif strcmp(handles.currentData, 'Phased Array')
          psi = k*handles.distance.*abs(cosd(handles.thetaD)-cosd(theta));
          handles.arrayFactor = Array_Factor_Cartesian(handles.elements, psi);
          AE= handles.arrayElement;
          theta = linspace(-0,2*pi,360);
          axes(handles.Polar)
          RadiationPattern = AE.*abs(handles.arrayFactor);
          PolarPlot(theta,RadiationPattern,handles.currentData)
else
          disp('none..')
end
ArrayInfo(hObject, eventdata, handles)

guidata(hObject, handles);


%*****************************************
%**************** > 3D < *****************
%*****************************************
% --- Executes on button press in ThreeD.
function ThreeD_Callback(hObject, eventdata, handles)
handles.output = hObject;
k = 2*pi/handles.waveLength;

if strcmp(handles.currentData, 'Broadside Array')   
          [theta,rho] = meshgrid(linspace(0,2*pi,360));
          num = ((1./handles.elements).*sin(handles.elements.*(k.*handles.distance./2).*cos(theta)));
          den = sin((k*handles.distance/2).*cos(theta));
          AF  = num./den; 
          AE= handles.arrayElement;
          RadiationPattern = AE.*abs(AF);
          axes(handles.Three)
          
          [x,y,z] = sph2cart(rho,theta,RadiationPattern);
          surf(x,y,z,'FaceLighting','phong',...
                     'LineStyle','none');
          colormap winter
          axis image

          light('Style','local',...
                'Position',[-10.162701816704 -0.924193626363743 14.9951905283833]);

elseif strcmp(handles.currentData, 'End-Fire Array')
          [theta,rho] = meshgrid(linspace(0,2*pi,360));
          num = ((1./handles.elements).*(sin(handles.elements.*(k.*handles.distance./2).*(cos(theta)+1)))) ;
          den = sin((k*handles.distance/2).*(cos(theta)+1));
          AF  = num./den; 
          AE= handles.arrayElement;
          RadiationPattern = AE.*abs(AF);
          axes(handles.Three)
          
          [x,y,z] = sph2cart(rho,theta,RadiationPattern);
          surf(x,y,z,'FaceLighting','phong',...
                     'LineStyle','none');
          colormap winter
          axis image

          light('Style','local',...
                'Position',[-10.162701816704 -0.924193626363743 14.9951905283833]);
          
elseif strcmp(handles.currentData, 'Hansen Woodyard')
          [theta,rho] = meshgrid(linspace(0,2*pi,360));
          num = ((1./handles.elements).*sin(handles.elements.*(k.*handles.distance./2).*(cos(theta)+deg2rad(108))));
          den = sin((k*handles.distance/2).*(cos(theta)+deg2rad(108)));
          AF  = num./den; 
          AE= handles.arrayElement;
          RadiationPattern = AE.*abs(AF);
          RadiationPattern = RadiationPattern/max(RadiationPattern(:));
          axes(handles.Three)
          
          [x,y,z] = sph2cart(rho,theta,RadiationPattern);
          surf(x,y,z,'FaceLighting','phong',...
                     'LineStyle','none');
          colormap winter
          axis image

          light('Style','local',...
                'Position',[-10.162701816704 -0.924193626363743 14.9951905283833]);
          
elseif strcmp(handles.currentData, 'Phased Array')
          [theta,rho] = meshgrid(linspace(0,2*pi,360));
          num = ((1./handles.elements).*sin(handles.elements.*(k.*handles.distance/2).*(abs(cosd(handles.thetaD)-cosd(theta)))));
          den = sin((k*handles.distance/2).*abs(cosd(handles.thetaD)-cosd(theta)));
          AF  = num./den; 
          AE= handles.arrayElement;
          RadiationPattern = AE.*abs(AF);
          axes(handles.Three)
          
          [x,y,z] = sph2cart(rho,theta,RadiationPattern);
          surf(x,y,z,'FaceLighting','phong',...
                     'LineStyle','none');
          colormap winter
          axis image

          light('Style','local',...
                'Position',[-10.162701816704 -0.924193626363743 14.9951905283833]);
else
          disp('none..')
end
ArrayInfo(hObject, eventdata, handles)

guidata(hObject, handles);



function Three_CreateFcn(hObject, eventdata, handles)

function directivity_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function directivity_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SavePlots.
function SavePlots_Callback(hObject, eventdata, handles)

% --- Executes on button press in SaveResults.
function SaveResults_Callback(hObject, eventdata, handles)
global displayDate
global rows
global info0
global info1
global info2
global info3
global info4
global info5
global info6
createFile(hObject, eventdata, handles,displayDate,rows,info0,info1,info2,info3,info4,info5,info6)
% --- Executes on button press in SaveAll.
function SaveAll_Callback(hObject, eventdata, handles)
