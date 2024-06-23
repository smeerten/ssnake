# -*- mode: python -*-

block_cipher = None

basePath = 'c:\\Users\\Emulator\\Desktop\\ssNakeBuild\\ssnake\\src\\'

#Needed:
# scipy.sparse (for optimize it seems)
#scipy.spatial

#not needed from mpl-data:
#fonts: only need all in 'ttf'
#images: not needed
#sample data: not needed

#qt:
#pyqt5 --> Qt --> Translations (all except 'en' can be removed)
#plugins --> imageformats:  qjpeg, qtiff, qwebp, qicns,  not needed (only need gif ico svg)
# plugins --> pltforms: only need qwindows

a = Analysis([basePath + 'ssNake.py'],
             pathex=['C:\\Users\\Emulator\\AppData\\Local\\Programs\\Python\\Python38\\Lib\\site-packages\\scipy\\.libs','C:\\Users\\Emulator\\Desktop'],
             binaries=[],
             datas=[(basePath + 'IsotopeProperties', '.'),(basePath + 'licenseHtml.txt','.'),
					(basePath + '\\Icons','Icons'),
					(basePath + 'References.txt','.'),
					(basePath + 'Tooltips','.')],
             hiddenimports=['spectrum','fitting', 'updateWindow','saveFigure','loadIsotopes','scipy.optimize', 'scipy._lib.messagestream',
			 'matplotlib.backends.backend_svg','matplotlib.backends.backend_ps','matplotlib.backends.backend_pdf'],
             hookspath=[],
             runtime_hooks=[],
             excludes=['tkinter','PyQt5.QtWebEngineWidgets','PyQt5.QtWebEngineCore','PyQt5.QtWebEngine','PyQt5.QtDesigner',
			 'PyQt5.QtXmlPatterns','PyQt5.QtXml','PyQt5.QtWebSockets','PyQt5.QtWebChannel','PyQt5.QtTest',
			 'PyQt5.QtSql','PyQt5.QtSerialPort','PyQt5.QtSensors','PyQt5.QtQuickWidgets','PyQt5.QtQuick','PyQt5.QtQml',
			 'PyQt5.QtWinExtras','PyQt5.QtSvg','PyQt5.QtPrintSupport','PyQt5.QtPositioning','PyQt5.QtOpenGL',
			 'PyQt5.QtNfc','PyQt5.QtNetworkAuth','PyQt5.QtNetwork','PyQt5.QtMultimediaWidgets','PyQt5.QtMultimedia',
			 'PyQt5.QtLocation','PyQt5.QtHelp','PyQt5.QtDBus','PyQt5.QtBluetooth',
			 'PyQt4','h5py.tests','pytz',
			 'matplotlib.backends._tkagg','lib2to3','numba'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='ssNake',
          debug=False,
          strip=False,
          upx=True,
          console=False, icon = basePath + '\Icons\logo.ico' )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='ssNake')
