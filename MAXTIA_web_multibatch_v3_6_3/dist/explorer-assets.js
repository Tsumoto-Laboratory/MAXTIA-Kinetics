"use strict";

// Include the bundled renderer in a portable HTML result only for mapped runs.
async function stageExplorerAssets(runtime, home, workerURL) {
  const names = ["style.css", "molecular-view.js", "explorer.js", "vendor/3Dmol-min.js", "vendor/3Dmol-LICENSE.txt",
    "vendor/3Dmol-min.js.LICENSE.txt"];
  runtime.FS.mkdirTree(home+"/vendor");
  await Promise.all(names.map(async name=>{
    const url=new URL(name,workerURL);
    url.search=new URL(workerURL).search;
    const response=await fetch(url);
    if (!response.ok) throw Error("Could not load the offline structure viewer. Reload the page and try again.");
    runtime.FS.writeFile(home+"/"+name,new Uint8Array(await response.arrayBuffer()));
  }));
}
