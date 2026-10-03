const { contextBridge, ipcRenderer } = require('electron');
contextBridge.exposeInMainWorld('plasmidNative', {
  saveFile: (name, text, existingPath) => ipcRenderer.invoke('save-file', { name, text, existingPath }),
  onOpenFile: cb => ipcRenderer.on('open-file', (_e, payload) => cb(payload)),
});
