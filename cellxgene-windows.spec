# -*- mode: python -*-

block_cipher = None


a = Analysis(['server\\gui\\main.py'],
             pathex=['C:\\Users\\Charlotte\\Documents\\git\\cellxgene'],
             binaries=[],
             datas=[('server/app/web/templates/', 'server/app/web/templates'), ('server/app/web/static/', 'server/app/web/static')],
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
          console=True , icon='server\gui\images\icon.ico')
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='cellxgene')
