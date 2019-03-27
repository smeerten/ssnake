copy ssnake\ReferenceManual.pdf dist\ReferenceManual.pdf
xcopy /E ssnake\Tutorial dist\Tutorial\
rmdir /S/Q .\dist\Tutorial\src\
rmdir /S/Q .\dist\ssNake\Icons\Sources\
copy gpl.rtf dist
copy modern_ssNake.nsi dist
"C:\Program Files (x86)\NSIS\makensis.exe" dist\modern_ssNake.nsi
