"use strict";

self.onmessage = async ({data}) => {
  const stage = text => self.postMessage({stage: text});
  try {
    stage("Starting the local analysis engine…");
    const assets = new URL("./runtime/", self.location.href).href;
    importScripts(new URL("pyodide.js", assets).href);
    const runtime = await loadPyodide({indexURL: assets});
    const home = "/home/pyodide";
    for (const name of ["maxitia_mapping.py", "structure_explorer.py", "bridge.py", "multibatch.py"]) {
      const url = new URL(name, self.location.href);
      url.search = self.location.search;
      const response = await fetch(url, {cache: "no-cache"});
      if (!response.ok) throw Error("Could not load the analysis engine. Reload this page and try again.");
      runtime.FS.writeFile(home + "/" + name, await response.text());
    }
    if (data.demo || data.shared?.cif) {
      importScripts("./explorer-assets.js" + self.location.search);
      await stageExplorerAssets(runtime, home, self.location.href);
    }
    const release = new URL(self.location.href).searchParams.get("v");
    if (release && (runtime.runPython("import maxitia_mapping; maxitia_mapping.MODULE_VERSION") !== release ||
                    runtime.runPython("import multibatch; multibatch.VERSION") !== release)) {
      throw Error("The analysis engine and page have different versions. Open the address printed by the launcher in a new tab.");
    }
    async function writeFiles(files, dir) {
      runtime.FS.mkdir(dir);
      const paths = {};
      for (const [role, file] of Object.entries(files)) {
        if (!file) continue;
        const basename = file.name.split(/[\\/]/).pop();
        if (!basename || basename === "." || basename === "..") throw Error("Select a file with a valid name.");
        runtime.FS.mkdir(dir + "/" + role);
        paths[role] = dir + "/" + role + "/" + basename;
        runtime.FS.writeFile(paths[role], new Uint8Array(await file.arrayBuffer()));
      }
      return paths;
    }
    const request = {demo:data.demo, partner:data.partner, wt:data.wt, chain:data.chain, batches:[],
      shared:await writeFiles(data.shared || {}, home + "/shared")};
    for (const [index, batch] of (data.batches || []).entries()) {
      const {files, ...settings} = batch;
      request.batches.push({...settings, paths:await writeFiles(files, home + "/input" + index)});
    }
    stage("Checking each batch and comparing WT-normalized effects…");
    runtime.globals.set("request_json", JSON.stringify(request));
    const result = JSON.parse(await runtime.runPythonAsync(
      "from multibatch import analyze_study, json, Path\nanalyze_study(json.loads(request_json), Path('/home/pyodide'))"));
    const zip = runtime.FS.readFile(home + "/results.zip");
    self.postMessage({result, zip}, [zip.buffer]);
  } catch (error) {
    self.postMessage({error:String(error.message || error).trim().split("\n").pop().replace(/^\w*Error:\s*/, "")});
  }
};
