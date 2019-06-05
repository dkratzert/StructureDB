; Script generated by the Inno Setup Script Wizard.
; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!

#define MyAppName "StructureFinder"
#define MyAppVersion "43"
#define MyAppPublisher "Daniel Kratzert"
#define pzipfile "Python3.6.1-32.7z"


[Setup]
; NOTE: The value of AppId uniquely identifies this application.
; Do not use the same AppId value in installers for other applications.
; (To generate a new GUID, click Tools | Generate GUID inside the IDE.)
AppId={{FD3791DD-E642-47A6-8434-FBD976271019}}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
DefaultDirName={pf}\StructureFinder
OutputBaseFilename=StructureFinder-setup-x32-v{#MyAppVersion}
Compression=lzma2/fast
SolidCompression=yes
SetupLogging=True
CloseApplications=False
RestartApplications=False
ShowLanguageDialog=no
ChangesAssociations=True
RestartIfNeededByRun=False
ChangesEnvironment=True
DisableFinishedPage=True
DisableReadyPage=True
DisableReadyMemo=True
DisableWelcomePage=True
AlwaysShowDirOnReadyPage=True
InternalCompressLevel=fast
EnableDirDoesntExistWarning=True
DirExistsWarning=no
UninstallLogMode=new
VersionInfoVersion={#MyAppVersion}
MinVersion=0,6.1
DefaultGroupName=StructureFinder
DisableProgramGroupPage=yes
AppendDefaultGroupName=True
AppContact=dkratzert@gmx.de
AppCopyright=Daniel Kratzert
AppSupportPhone=+49 761 203 6156
VersionInfoProductName=StructureFinder
AlwaysShowComponentsList=False
ShowComponentSizes=False
SetupIconFile="..\icons\strf.ico"


[UninstallRun]


[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

; adds a new page to the setup where you can choose if the path should be added
;Excludes: "*.pyc"

[Files]
Source: "..\apex\*"; DestDir: "{app}\apex"; Flags: ignoreversion createallsubdirs recursesubdirs; Excludes: "*.pyc"
Source: "..\displaymol\*"; DestDir: "{app}\displaymol"; Flags: ignoreversion createallsubdirs recursesubdirs; Excludes: "*.pyc"
Source: "..\icons\*"; DestDir: "{app}\icons"; Flags: ignoreversion createallsubdirs recursesubdirs; Excludes: "*.pyc"
Source: "..\lattice\*"; DestDir: "{app}\lattice"; Flags: ignoreversion createallsubdirs recursesubdirs; Excludes: "*.pyc"
Source: "..\misc\*"; DestDir: "{app}\misc"; Flags: ignoreversion createallsubdirs recursesubdirs; Excludes: "*.pyc"
Source: "..\pg8000\*"; DestDir: "{app}\pg8000"; Flags: ignoreversion createallsubdirs recursesubdirs; Excludes: "*.pyc"
Source: "..\pymatgen\*"; DestDir: "{app}\pymatgen"; Flags: ignoreversion createallsubdirs recursesubdirs; Excludes: "*.pyc"
Source: "..\searcher\*"; DestDir: "{app}\searcher"; Flags: ignoreversion createallsubdirs recursesubdirs; Excludes: "*.pyc"
Source: "..\ccdc\*"; DestDir: "{app}\ccdc"; Flags: ignoreversion createallsubdirs recursesubdirs; Excludes: "*.pyc"
Source: "..\shelxfile\*"; DestDir: "{app}\shelxfile"; Flags: ignoreversion createallsubdirs recursesubdirs; Excludes: "*.pyc"
Source: "..\p4pfile\*"; DestDir: "{app}\p4pfile"; Flags: ignoreversion createallsubdirs recursesubdirs; Excludes: "*.pyc"
Source: "..\gui\*"; DestDir: "{app}\gui"; Flags: ignoreversion createallsubdirs recursesubdirs; Excludes: "*.pyc"
Source: "..\strf.py"; DestDir: "{app}"; 
Source: "..\six.py"; DestDir: "{app}";
Source: "..\strf_cmd.py"; DestDir: "{app}"; 
Source: "win\strf_win_32.bat"; DestDir: "{app}"; DestName: "strf.bat"
; Caution change accordingly:
;Source: "F:\Programmieren\StructureFinder_distrib\{#PzipFile}"; DestDir: "{app}"; Flags: deleteafterinstall;
Source: "C:\tools\{#pzipfile}"; DestDir: "{app}"; Flags: deleteafterinstall;

[Run]
Filename: "{app}\misc\7z.exe"; Parameters: "x ""{app}\{#pzipfile}"" -o""{app}"" * -r -aoa"; Flags: runascurrentuser postinstall; 

[Icons]
Name: "{group}\{cm:UninstallProgram,{#MyAppName}}"; Filename: "{uninstallexe}"
Name: "{group}\StructureFinder"; Filename: "{app}\strf.bat"; WorkingDir: "{app}"; IconFilename: "{app}\icons\strf.ico"

[UninstallDelete]
Type: files; Name: "{app}\*.pyc"
Type: files; Name: "{app}\*.*"
Type: filesandordirs; Name: "{app}\*"

[Dirs]
Name: "{app}\displaymol"; Permissions: authusers-full
Name: "{app}\."; Permissions: authusers-full

[InstallDelete]

[Tasks]

[Code]


