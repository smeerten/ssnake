Set oWS = WScript.CreateObject("WScript.Shell") 
set fso = CreateObject("Scripting.FileSystemObject")
strDesktop= oWS.SpecialFolders("Desktop") 
ssNakeDir = fso.GetAbsolutePathName(".")
sLinkFile = strDesktop + "\ssNake.lnk"  
Set oLink = oWS.CreateShortcut(sLinkFile) 
oLink.TargetPath = "python" 
Arguments = ssNakeDir + "\ssNake.py"
oLink.Arguments = """"& Arguments &""""
oLink.IconLocation = ssNakeDir + "\logo.ico"
oLink.Save 

' StartMenu
strStartMenu= oWS.SpecialFolders("Programs")
StartLocation = strStartMenu+"\ssNake.lnk"
 
Set oLink = oWS.CreateShortcut(StartLocation) 
oLink.TargetPath = "python" 
Arguments = ssNakeDir + "\ssNake.py"
oLink.Arguments = """"& Arguments &""""
oLink.IconLocation = ssNakeDir + "\logo.ico"
oLink.Save 
