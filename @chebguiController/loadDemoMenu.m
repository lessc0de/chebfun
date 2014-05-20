function loadDemoMenu(handles)

% TODO:  Documentation.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Begin by checking whether we have already loaded the demos
if ( ~isempty(get(handles.menu_demos, 'UserData')) )
    return
end

% Set up ODEs, PDEs and EIGs demos separately

% Find the folders which demos are stored in. The chebguiDemos folder lives in
% the trunk folder, find the path of the Chebfun trunk.
trunkPath = fileparts(which('chebguiwindow'));

% Append directory information
bvppath = fullfile(trunkPath, 'chebguiDemos', 'bvpdemos');
pdepath = fullfile(trunkPath, 'chebguiDemos', 'pdedemos');
eigpath = fullfile(trunkPath, 'chebguiDemos', 'eigdemos');

% Setup ODEs
D = dir(bvppath);
for demoCounter = 1:length(D)
    demoPath = fullfile(bvppath, D(demoCounter, :).name);
    if ( isempty(strfind(demoPath, '.guifile')) )
        % Only want to load files ending in .guifile
        continue
    end

    % Need to obtain the name and type of the demo as well
    fid = fopen(demoPath);

    % Throw away ' at the ends of the string
    demoName = fgetl(fid);
    demoName = demoName(2:end-1);
    demoType = fgetl(fid);
    demoType = demoType(2:end-1);

    fclose(fid);

    demoFun = @(hObject, eventdata) ...
        hOpenMenuitemCallback(hObject, eventdata, handles, demoPath);
    switch demoType
        case 'bvp'
            hDemoitem  =  uimenu('Parent', handles.menu_bvps, ...
                'Label', demoName, ...
                'Separator', 'off', ...
                'HandleVisibility', 'callback',  ...
                'Callback', demoFun);
        case 'ivp'
            hDemoitem  =  uimenu('Parent', handles.menu_ivps, ...
                'Label', demoName, ...
                'Separator','off', ...
                'HandleVisibility','callback', ...
                'Callback', demoFun);
        case 'system'
            hDemoitem  =  uimenu('Parent', handles.menu_systems, ...
                'Label', demoName, ...
                'Separator', 'off', ...
                'HandleVisibility','callback', ...
                'Callback', demoFun);
    end
end

% Setup PDEs
D = dir(pdepath);
for demoCounter = 1:length(D) % First two entries are . and ..
    demoPath = fullfile(pdepath, D(demoCounter,:).name);

    % Only want to load files ending in .guifile
    if ( isempty(strfind(demoPath, '.guifile')) )
        continue
    end

    % Need to obtain the name and type of the demo as well
    fid = fopen(demoPath);

    % Throw away ' at the ends of the string
    demoName = fgetl(fid);
    demoName = demoName(2:end-1);
    demoType = fgetl(fid);
    demoType = demoType(2:end-1);

    fclose(fid);

    demoFun = @(hObject, eventdata) ...
        hOpenMenuitemCallback(hObject, eventdata, handles, demoPath);
    switch demoType
        case 'scalar'
            hDemoitem  =  uimenu('Parent', handles.menu_pdesingle, ...
                'Label', demoName, ...
                'Separator', 'off', ...
                'HandleVisibility', 'callback', ...
                'Callback', demoFun);
        case 'system'
            hDemoitem  =  uimenu('Parent', handles.menu_pdesystems, ...
                'Label', demoName, ...
                'Separator', 'off', ...
                'HandleVisibility', 'callback', ...
                'Callback', demoFun);
    end
end

% Setup EIGs
D = dir(eigpath);
for demoCounter = 1:length(D) % First two entries are . and ..
    demoPath = fullfile(eigpath,D(demoCounter,:).name);
    if isempty(strfind(demoPath,'.guifile')) % Only want to load files ending in .guifile
        continue
    end
    % Need to obtain the name and type of the demo as well
    fid = fopen(demoPath);

    % Throw away ' at the ends of the string
    demoName = fgetl(fid);
    demoName = demoName(2:end-1);
    demoType = fgetl(fid);
    demoType = demoType(2:end-1);

    fclose(fid);

    demoFun = @(hObject, eventdata) ...
        hOpenMenuitemCallback(hObject, eventdata, handles, demoPath);
    switch demoType
        case 'scalar'
            hDemoitem  =  uimenu('Parent', handles.menu_eigsscalar, ...
                'Label', demoName, ...
                'Separator', 'off', ...
                'HandleVisibility', 'callback', ...
                'Callback', demoFun);
        case 'system'
            hDemoitem  =  uimenu('Parent', handles.menu_eigssystem, ...
                'Label', demoName, ...
                'Separator', 'off', ...
                'HandleVisibility', 'callback', ...
                'Callback', demoFun);
    end
end

% Notify that we have loaded demos to prevent reloading
set(handles.menu_demos, 'UserData', 1);

end

function hOpenMenuitemCallback(hObject, eventdata, handles, demoPath)

% TODO:  Documentation.

% Callback function run when the Open menu item is selected
handles.guifile = loaddemos(handles.guifile, demoPath);
initSuccess = chebguiController.populate(handles, handles.guifile);

if ( initSuccess )
    switchModeCM = 'demo';
else
    switchModeCM = 'notdemo';
end

% Switch the mode of the GUI according to the type of the problem.
chebguiController.switchMode(handles, handles.guifile.type, switchModeCM);

% We no longer have a solution.
handles.hasSolution = 0;
set(handles.button_exportsoln, 'Enable', 'off');

% Update handle structure
guidata(hObject, handles);

end
