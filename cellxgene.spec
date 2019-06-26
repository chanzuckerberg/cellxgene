# -*- mode: python -*-

block_cipher = None


a = Analysis(['server/gui/main.py'],
             pathex=['/Users/charlotteweaver/Documents/Git/cellxgene'],
             binaries=[],
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
          a.binaries,
          a.zipfiles,
          a.datas,
          [],
          name='cellxgene',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          runtime_tmpdir=None,
          console=False , icon='server/gui/images/icon.ico')
app = BUNDLE(exe,
             name='cellxgene.app',
             icon='server/gui/images/icon.ico',
             bundle_identifier=None)
