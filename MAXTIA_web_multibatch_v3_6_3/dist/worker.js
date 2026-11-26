"use strict";
// One worker owns one run. Terminating it releases the runtime and its virtual filesystem.
self.onmessage = async ({data}) => {
  const stage = message => self.postMessage({stage: message});
  try {
    stage("Starting the browser analysis engine…");
    const assets = new URL("./runtime/", self.location.href).href;
    importScripts(new URL("pyodide.js", assets).href);
    const [runtime, ...sources] = await Promise.all([
      loadPyodide({indexURL: assets}),
      ...["maxitia_mapping.py", "structure_explorer.py", "bridge.py"].map(async name => {
        const url = new URL(name, self.location.href);
        url.search = self.location.search;
        const response = await fetch(url, {cache: "no-cache"});
        if (!response.ok) throw Error("Could not load the analysis code. Reload the page and try again.");
        return [name, await response.text()];
      })
    ]);
    const home = "/home/pyodide";
    for (const [name, source] of sources) runtime.FS.writeFile(home + "/" + name, source);
    if (data.demo || data.files.cif) {
      importScripts("./explorer-assets.js" + self.location.search);
      await stageExplorerAssets(runtime, home, self.location.href);
    }
    const paths = {};
    for (const [role, file] of Object.entries(data.files)) {
      if (!file) continue;
      const dir = home + "/" + role;
      runtime.FS.mkdir(dir);
      const name = file.name.split(/[\\/]/).pop();
      if (!name || name === "." || name === "..") throw Error("Check the input filename and select the file again.");
      paths[role] = dir + "/" + name;
      runtime.FS.writeFile(paths[role], new Uint8Array(await file.arrayBuffer()));
    }
    const release = new URL(self.location.href).searchParams.get("v");
    if (release && runtime.runPython("import maxitia_mapping; maxitia_mapping.MODULE_VERSION") !== release) {
      throw Error("The analysis engine and page have different versions. Open the address printed by the launcher in a new tab.");
    }
    stage(data.files.cif ? "Running QC and mapping residues to the structure…" : "Matching plate wells and running QC…");
    runtime.globals.set("request_json", JSON.stringify({paths, demo:data.demo, ...data.settings}));
    const result = JSON.parse(await runtime.runPythonAsync("from bridge import analyze, Path, json\nanalyze(json.loads(request_json), Path('/home/pyodide'))"));
    const zip = runtime.FS.readFile(home + "/results.zip");
    self.postMessage({result, zip}, [zip.buffer]);
  } catch (err) {
    const lines = String(err.message || err).trim().split("\n");
    let message = lines.at(-1).replace(/^\w+(?:Error|Exception):\s*/, "");
    if (message.includes("Set ngs_plate_block")) message = "The NGS file contains multiple plate blocks. Select the primer barcode block for this experiment and run again.";
    self.postMessage({error:message || "Check the input files and run again."});
  }
};
