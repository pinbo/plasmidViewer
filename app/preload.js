const { contextBridge, ipcRenderer, webUtils } = require('electron');
contextBridge.exposeInMainWorld('plasmidNative', {
  saveFile: (name, text, existingPath) => ipcRenderer.invoke('save-file', { name, text, existingPath }),
  // full path of a File chosen in the page (Open button / drag & drop), so Save can write back in place
  pathForFile: file => { try { return webUtils.getPathForFile(file); } catch (e) { return ''; } },
  onOpenFile: cb => ipcRenderer.on('open-file', (_e, payload) => cb(payload)),
});
