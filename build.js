#!/usr/bin/env node
// Bundles src/ into ONE self-contained file: dist/PlasmidViewer.html (also copied to app/index.html for Electron).
const fs = require('fs');
const path = require('path');
const src = path.join(__dirname, 'src');
let html = fs.readFileSync(path.join(src, 'index.html'), 'utf8');
html = html.replace(/<link rel="stylesheet" href="([^"]+)">/g, (_, f) => `<style>\n${fs.readFileSync(path.join(src, f), 'utf8')}\n</style>`);
html = html.replace(/<script src="([^"]+)"><\/script>\s*/g, (_, f) => `<script>\n${fs.readFileSync(path.join(src, f), 'utf8').replace(/<\/script/gi, '<\\/script')}\n</script>\n`);
fs.mkdirSync(path.join(__dirname, 'dist'), { recursive: true });
fs.writeFileSync(path.join(__dirname, 'dist', 'PlasmidViewer.html'), html);
if (fs.existsSync(path.join(__dirname, 'app'))) fs.writeFileSync(path.join(__dirname, 'app', 'index.html'), html);
console.log(`Built dist/PlasmidViewer.html (${(html.length / 1024).toFixed(0)} KB)`);
