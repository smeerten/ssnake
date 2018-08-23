Set oWS = WScript.CreateObject("WScript.Shell") 
set fso = CreateObject("Scripting.FileSystemObject")
strDesktop= oWS.SpecialFolders("Desktop") 
ssNakeDir = fso.GetAbsolutePathName(".")
sLinkFile = strDesktop + "\ssNake.lnk"  
Set oLink = oWS.CreateShortcut(sLinkFile) 
Target = ssNakeDir + "\WindowsRun.bat"
oLink.TargetPath = """"& Target &"""" 
oLink.IconLocation = ssNakeDir + "\Icons\logo.ico"
oLink.WorkingDirectory = ssNakeDir
oLink.Save 

' StartMenu
strStartMenu= oWS.SpecialFolders("Programs")
StartLocation = strStartMenu+"\ssNake.lnk"
 
Set oLink = oWS.CreateShortcut(StartLocation) 
oLink.TargetPath = """"& Target &"""" 
oLink.IconLocation = ssNakeDir + "\Icons\logo.ico"
oLink.WorkingDirectory = ssNakeDir
oLink.Save 
