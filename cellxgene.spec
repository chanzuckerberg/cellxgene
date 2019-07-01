# -*- mode: python -*-

block_cipher = None


a = Analysis(['server/gui/main.py'],
             pathex=['/Users/charlotteweaver/Documents/Git/cellxgene'],
             binaries=[('/System/Library/Frameworks/Tk.framework/Tk', 'tk'), ('/System/Library/Frameworks/Tcl.framework/Tcl', 'tcl')],
             datas=[('server/app/web/templates/', 'server/app/web/templates/'), ('server/app/web/static/', 'server/app/web/static/')],
             hiddenimports=[],
             hookspath=['server/gui/'],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='cellxgene',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=False , icon='server/gui/images/icon.ico')
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='cellxgene')
app = BUNDLE(coll,
             name='cellxgene.app',
             icon='server/gui/images/icon.ico',
             bundle_identifier=None)
