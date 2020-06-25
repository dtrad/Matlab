function edit(varargin)
echo on
%EDIT Edit M-file.
%   MYEDIT FUN opens the file FUN.M in a text editor.  FUN must be the
%   name of an m-file or a MATLABPATH relative partial pathname (see
%   PARTIALPATH).
%
%   EDIT FILE.EXT opens the specified text file.
%
%   EDIT FUN1 IN FUN2 opens the file FUN1 in the context of m-file
%   FUN2.  
%
%   EDIT FUN(A,B,C) opens the file FUN which matches its the 
%   given input arguments.  For example, edit feval(g), when
%   g=inline('sin(x)'), opens inline/feval.m.
%
%   EDIT, by itself, opens up a new editor window.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.43 $  $Date: 1997/11/21 23:31:45 $

% Reject anything but string arguments
if ~iscellstr(varargin)
   disp('The input must be a string.');
   return;
end

% Branch for Mac editing tools
if strncmp(computer, 'MAC', 3),
   if nargin == 0,
      mactools('meditor')
   elseif nargin == 1,
      name = varargin{:};
      if evalin('caller', ['exist(''' name ''') == 1']),
         disp(sprintf('Can''t edit the variable ''%s''.', name));
         return;   
      end
      mactools('meditor', name);
   else
      disp('Too many input arguments.');
   end
   
   return;
end


%
% Decide which file, if any, to edit.
%
% Pass the arguments through to WHICH.
% If there are no matching files in the WHICH's scope, assume
% that the file isn't on the path, and try opening it directly.
%
fcnname = '';
if nargin > 0
   fcnname = varargin{1};
   fullpath = [];
   whichargs = [ '''' fcnname ''''];
   
   % If there is more than one arg, we're using a special form of WHICH
   if nargin > 1
      whichargs = [whichargs sprintf(',''%s''', varargin{2:end})];
   end
   
   % Pass the incoming arguments through to WHICH
   fullpath = evalin('caller', ['which(' whichargs ')']);
   
   % Workaround for g34404: x=which(f), where f is an MDL file, has extra text in front.
   mdltag = 'loaded model ';
   if ~isempty(strmatch(mdltag, fullpath)) & length(fullpath) > length(mdltag)
      fullpath = fullpath(length(mdltag)+1:end);
   end
   
   %
   % Eliminate cases where WHICH returns a value that isn't a legitimate pathname.
   % Can't catch this using EXIST on fcnname below, because for some forms of EDIT
   % (e.g., edit subfun in fun), fcnname isn't enough information for EXIST.
   %
   if strcmp(fullpath, 'variable')
      disp(sprintf('Can''t edit the variable ''%s''.', fcnname));
      return;
   elseif strcmp(fullpath, 'built-in')
      disp(sprintf('Can''t edit the built-in function ''%s''.', fcnname));
      return;
   end
   
   % If we don't find a match, assume we're trying to edit a file
   % that isn't on the path.
   if isempty(fullpath) & nargin == 1
      fullpath = varargin{1};
   end
   
   % If we still don't have a match, give up.
   if isempty(fullpath)
      disp(sprintf('File ''%s'' not found.', fcnname));
      return;
   end
      
   % Make sure the file we're opening is of an openable type
   exists = evalin('caller', ['exist(''' fullpath ''')']);
   
   switch exists
   case 0
      disp(sprintf('File ''%s'' not found.', fcnname));
      return;
   case 1   % Variable: this is caught above.
      disp('Internal error in edit.m.');
      return;
   case 3
      disp(sprintf('Can''t edit the MEX-file ''%s''.', fcnname));
      return;
   case 4
      % Editing an MDL-file is allowed only when the user specifies the file
      % extension .mdl explicitly, making their intent clear.  If you edit
      % an MDL-file accidentally you can corrupt it and trash the model.  Only
      % MathWorks developers and sophisticated users edit MDL-files.
      % (Note: it's OK to do case-insensitive comparison against '.mdl' on UNIX
      % here even though the UNIX file system is case-sensitive, because exist
      % has already told us that this must be an MDL-file.)
      if (length(fcnname) < 4) | (strcmp(lower(fcnname(end-3:end)), '.mdl') ~= 1),
         disp(sprintf('Can''t edit the MDL-file ''%s'' unless you include the ''.mdl'' file extension.', fcnname));
         return;
      end
   case 5   % Built-in: this is caught above.
      disp('Internal error in edit.m');
      return;
   case 7
      disp(sprintf('Can''t edit the directory ''%s''.', fcnname));
      return;
   end
   
   
   % Map a P-file to the corresponding M-file.  Note: we can't do this by
   % checking exist(name) == 6, because name might point to a P-file that
   % is not on the path.
   if length(fullpath) >= 2,
      ext = fullpath(end-1:end);
      if ~isunix, ext = lower(ext); end  % file system on PC is case-tolerant
      if strcmp(ext, '.p') == 1,
         if exist([fullpath(1:end-2) '.m']) ~= 2,  % will return 2 if the M-file exists
            disp(sprintf('The P-file ''%s'' has no matching M-file ''%s'' to edit.', fullpath, [fullpath(1:end-2) '.m']));
            return;   
         end
         fullpath(end-1:end) = '.m';
      end
   end
end

% Now determine which editor to run
% Assume builtin editor to begin with
builtinEd = 1;

% Get the MATLAB editor preference.
% Note: On UNIX, the user's editor preference is read by miedit
% from the .Xdefaults file.
if strncmp(computer, 'PC', 2),
   if getprofl('Matlab Settings', 'Built-in Editor', '1', 'matlab.ini') == '0'
      editor = getprofl('Matlab Settings', 'Editor Pref', 'medit.exe', 'matlab.ini');
      builtinEd = 0;
   end
end
 %builtinEd = 0;

if builtinEd == 1,
   % Use the built-in editor
   if nargin == 0,
      miedit('editor');
   else
      miedit('editor', fullpath);
   end
else
   % User-specified editor
   if nargin == 0,
      eval(['!"' editor '" &'])
   else
      eval(['!"' editor '" "' fullpath '" &'])
   end
end
