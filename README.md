# Plasmid Viewer

A SnapGene-Viewer-style plasmid map viewer **and editor** that runs entirely offline.
No server, no dependencies, no install.

## Run it (any OS)

Open **`dist/PlasmidViewer.html`** in Chrome, Edge, Firefox or Safari. That single file *is* the app —
copy it to a USB stick, email it, put it on another machine. Double-click and go.
(Chrome/Edge additionally let "Save" overwrite the original file in place; other browsers download a `.gb`.)

## Optional: native desktop app (Mac / Windows / Linux)

Requires Node.js 18+ once, on the machine you build from:

```bash
cd app
npm install
npm start            # try it
npm run dist:mac     # → app/dist/*.dmg      (build on a Mac)
npm run dist:win     # → app/dist/*.exe      (build on Windows, or Linux with wine)
npm run dist:linux   # → app/dist/*.AppImage / *.deb
```
The installers give you double-click-to-open for `.gb`, `.fasta` and `.dna`, and native save dialogs.

## Features

* **Open**: GenBank, FASTA, SnapGene `.dna` (read), drag & drop, multiple tabs. **Save**: GenBank (+ FASTA export).
* **Maps**: circular and linear, feature arrows, labels, tick marks, click/drag to select, zoomable linear map.
* **Sequence view**: both strands, feature bars, amino-acid translation in CDS bars, ORF finder, restriction sites with cut marks.
* **Auto-annotation**: 36 common features built in (promoters, tags, resistance genes, ori, polyA, FPs …), detected on both strands, tolerant to a few mismatches; protein tags (6xHis, FLAG, HA, Myc, V5 …) are found in any reading frame.
* **Feature library**: select a region → *Add feature* → tick “add to library”. Or *Library → Add a new feature* (DNA or protein motif). Import/export as JSON to share with colleagues.
* **Restriction enzymes**: 233 commercial enzymes incl. Type IIS (cross-checked against NEB with `node tools/check-enzymes.js`); hover any enzyme name for a pop-up with its recognition site, overhang type, **number of sites in the current plasmid**, and NEB data (supplied buffer, % activity in each NEBuffer, incubation/heat-inactivation temperature, Dam/Dcm/CpG sensitivity, star-activity and ligation notes, HF versions); unique / dual / all cutters, a **Golden Gate** set (BsaI, BsmBI/Esp3I, BbsI/BpiI, SapI/BspQI, AarI, PaqCI …) and a **My selected enzymes** set (tick enzymes in the Enzymes tab); add your own enzymes too.
* **Editing**: type bases, insert (with reverse-complement option), delete, replace, undo/redo,
  **copy / cut / paste with features** (fully-contained features travel with the sequence; paste reverse-complement flips them), reverse-complement a selection or the whole plasmid, set origin, find (both strands, IUPAC), go-to.

## Developing

Sources are in `src/` (plain JS, no framework, no build tooling beyond a 20-line bundler).

```bash
node build.js        # regenerates dist/PlasmidViewer.html (and app/index.html)
```

## Notes / limits

* Built-in library sequences were written from memory of canonical sequences; matching tolerates ~4 % mismatches for long elements, but
  please spot-check annotations on your own constructs and correct/add entries in the Library as needed.
* SnapGene `.dna` is read-only; save as GenBank (SnapGene opens it).
* Selection is a single contiguous range (not across the origin); features may span the origin.
