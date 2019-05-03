;--------------------------------
;Include Modern UI

  !include "MUI2.nsh"

;--------------------------------
;General
  !define APPNAME "ssNake v1.0"
  !define DESCRIPTION "NMR processing utility"
  !define INSTALLSIZE 236500
  ;Name and file
  Name "${APPNAME}"
  OutFile "../ssNake v1.0 Installer.exe"

  ;Default installation folder
  InstallDir "$LOCALAPPDATA\ssNake"
  SetCompressor /SOLID lzma
  
  ;Request application privileges for Windows Vista
  RequestExecutionLevel user

;--------------------------------
;Interface Settings

  !define MUI_ABORTWARNING
  !define MUI_ICON "ssNake\Icons\logo.ico"
  !define MUI_UNICON "ssNake\Icons\logo.ico"

;--------------------------------
;Pages

  !insertmacro MUI_PAGE_LICENSE "gpl.rtf"
  ;!insertmacro MUI_PAGE_COMPONENTS
  !insertmacro MUI_PAGE_DIRECTORY
   Var StartMenuFolder
  !define MUI_STARTMENUPAGE_REGISTRY_ROOT "HKCU" 
  !define MUI_STARTMENUPAGE_REGISTRY_KEY "Software\ssNake v1.0" 
  !define MUI_STARTMENUPAGE_REGISTRY_VALUENAME "Start Menu Folder"
  
  !insertmacro MUI_PAGE_STARTMENU Application $StartMenuFolder
  !insertmacro MUI_PAGE_INSTFILES

    !define MUI_FINISHPAGE_NOAUTOCLOSE
    !define MUI_FINISHPAGE_RUN "$INSTDIR\ssNake\ssNake.exe"
    !define MUI_FINISHPAGE_RUN_TEXT "Start ssNake"
  !insertmacro MUI_PAGE_FINISH
  
  !insertmacro MUI_UNPAGE_CONFIRM
  !insertmacro MUI_UNPAGE_INSTFILES
  
;--------------------------------
;Languages
 
  !insertmacro MUI_LANGUAGE "English"

;--------------------------------
;Installer Sections

Section "Main Files"

  SetOutPath "$INSTDIR"
  file /r "ssNake"
  file /r "Tutorial"
  file "ReferenceManual.pdf"
  
  ;Store installation folder
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APPNAME}" "DisplayName" "${APPNAME} - ${DESCRIPTION}"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APPNAME}" "UninstallString" "$\"$INSTDIR\uninstall.exe$\""
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APPNAME}" "QuietUninstallString" "$\"$INSTDIR\uninstall.exe$\" /S"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APPNAME}" "InstallLocation" "$\"$INSTDIR$\""
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APPNAME}" "DisplayIcon" "$\"$INSTDIR\ssNake\Icons\logo.ico$\""
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APPNAME}" "Publisher" "ssNake"
  # There is no option for modifying or repairing the install
  WriteRegDWORD HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APPNAME}" "NoModify" 1
  WriteRegDWORD HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APPNAME}" "NoRepair" 1
  # Set the INSTALLSIZE constant (!defined at the top of this script) so Add/Remove Programs can accurately report the size
  WriteRegDWORD HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APPNAME}" "EstimatedSize" ${INSTALLSIZE}
  
  ;Create uninstaller
  WriteUninstaller "$INSTDIR\Uninstall.exe"
  SetOutPath "$INSTDIR\ssNake"
  
    !insertmacro MUI_STARTMENU_WRITE_BEGIN Application
    
    ;Create shortcuts
    CreateDirectory "$SMPROGRAMS\$StartMenuFolder"
	CreateShortcut "$SMPROGRAMS\$StartMenuFolder\${APPNAME}.lnk" "$INSTDIR\ssNake\ssNake.exe"
    CreateShortcut "$SMPROGRAMS\$StartMenuFolder\Uninstall.lnk" "$INSTDIR\Uninstall.exe"
  
    !insertmacro MUI_STARTMENU_WRITE_END

SectionEnd


;--------------------------------
;Uninstaller Section

Section "Uninstall"

  ;FILES

  Delete "$INSTDIR\Uninstall.exe"
  RMDir /r "$INSTDIR\ssNake"
  RMDir /r "$INSTDIR\Tutorial"
  Delete "$INSTDIR\ReferenceManual.pdf"
  RMDir "$INSTDIR"
  
  ;Registry
  DeleteRegKey HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APPNAME}"
  
  ; Remove Start Menu launcher
  !insertmacro MUI_STARTMENU_GETFOLDER Application $StartMenuFolder
  Delete "$SMPROGRAMS\$StartMenuFolder\Uninstall.lnk"
  Delete "$SMPROGRAMS\$StartMenuFolder\${APPNAME}.lnk"
  RMDir "$SMPROGRAMS\$StartMenuFolder"
  
  DeleteRegKey /ifempty HKCU "Software\ssNake v1.0"

SectionEnd
