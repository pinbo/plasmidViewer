// Minimal Electron shell: loads the single-file app, adds native Save / Open-With support.
const { app, BrowserWindow, dialog, ipcMain, Menu, shell } = require('electron');
const fs = require('fs');
const path = require('path');

let win;
let pendingFiles = process.argv.slice(app.isPackaged ? 1 : 2).filter(f => /\.(gb|gbk|genbank|gbff|fa|fasta|fna|dna|seq|txt)$/i.test(f) && fs.existsSync(f));

function sendFile(file) {
  try {
    const data = fs.readFileSync(file).toString('base64');
    win.webContents.send('open-file', { name: path.basename(file), data, path: /\.dna$/i.test(file) ? null : file });
  } catch (e) { dialog.showErrorBox('Could not open file', String(e.message || e)); }
}

function createWindow() {
  win = new BrowserWindow({
    width: 1400, height: 900, minWidth: 900, minHeight: 600, title: 'Plasmid Viewer',
    webPreferences: { preload: path.join(__dirname, 'preload.js'), contextIsolation: true, nodeIntegration: false },
  });
  win.loadFile(path.join(__dirname, 'index.html'));
  win.webContents.setWindowOpenHandler(({ url }) => { shell.openExternal(url); return { action: 'deny' }; });
  win.webContents.on('did-finish-load', () => { for (const f of pendingFiles) sendFile(f); pendingFiles = []; });
}

ipcMain.handle('save-file', async (_e, { name, text, existingPath }) => {
  let target = existingPath;
  if (!target) {
    const r = await dialog.showSaveDialog(win, { defaultPath: name, filters: [{ name: 'GenBank', extensions: ['gb', 'gbk'] }] });
    if (r.canceled) return null;
    target = r.filePath;
  }
  fs.writeFileSync(target, text, 'utf8');
  return target;
});

// macOS "Open With" / double-click on an associated file
app.on('open-file', (e, file) => { e.preventDefault(); if (win && win.webContents && !win.webContents.isLoading()) sendFile(file); else pendingFiles.push(file); });

app.whenReady().then(() => {
  Menu.setApplicationMenu(Menu.buildFromTemplate([
    ...(process.platform === 'darwin' ? [{ role: 'appMenu' }] : []),
    { role: 'editMenu' }, { role: 'viewMenu' }, { role: 'windowMenu' },
  ]));
  createWindow();
  app.on('activate', () => { if (!BrowserWindow.getAllWindows().length) createWindow(); });
});
app.on('window-all-closed', () => { if (process.platform !== 'darwin') app.quit(); });
