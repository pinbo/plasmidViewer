// Runs after `npm install`. Electron's own installer relies on `extract-zip`, which on some very new Node versions
// (e.g. Node 24+) finishes WITHOUT unpacking and without an error, leaving "Electron failed to install correctly".
// If that happened, unpack the (already downloaded) zip with the operating system's own tools instead.
const fs = require('fs'), path = require('path'), cp = require('child_process');
const dir = path.dirname(require.resolve('electron/package.json'));
if (fs.existsSync(path.join(dir, 'path.txt')) && fs.existsSync(path.join(dir, 'dist'))) process.exit(0);

const { downloadArtifact } = require(require.resolve('@electron/get', { paths: [dir, __dirname] }));
const { version } = require(path.join(dir, 'package.json'));
(async () => {
  console.log(`Electron ${version} was not unpacked by its installer – unpacking it manually…`);
  const zip = await downloadArtifact({ version, artifactName: 'electron', platform: process.platform, arch: process.arch, checksums: require(path.join(dir, 'checksums.json')) });
  const dist = path.join(dir, 'dist');
  fs.rmSync(dist, { recursive: true, force: true }); fs.mkdirSync(dist, { recursive: true });
  if (process.platform === 'darwin') cp.execFileSync('ditto', ['-x', '-k', zip, dist]);
  else if (process.platform === 'win32') cp.execFileSync('powershell', ['-NoProfile', '-Command', `Expand-Archive -Force -LiteralPath '${zip}' -DestinationPath '${dist}'`]);
  else cp.execFileSync('unzip', ['-q', '-o', zip, '-d', dist]);
  const exe = { darwin: 'Electron.app/Contents/MacOS/Electron', win32: 'electron.exe', linux: 'electron' }[process.platform];
  fs.writeFileSync(path.join(dir, 'path.txt'), exe);
  console.log('Electron unpacked.');
})().catch(e => { console.error(e); process.exit(1); });
