"use strict";
const assetQuery = new URL(document.currentScript.src).search;
const $ = id => document.getElementById(id);
const form = $("form"), fileIds = ["lsa","ngs","cif","full","summary","policy"];
// Keep the on-page picker text in English while using a native file input.
for (const id of fileIds) {
  const input = $(id);
  const picker = document.createElement("span");
  const button = document.createElement("span");
  const filename = document.createElement("span");
  picker.className = "file-picker";
  button.className = "file-button";
  button.textContent = "Choose file";
  filename.className = "filename";
  button.setAttribute("aria-hidden", "true");
  filename.setAttribute("aria-hidden", "true");
  input.before(picker);
  picker.append(button, filename, input);
  function showSelection() {
    const name = input.files[0]?.name;
    filename.textContent = name || "No file selected";
    input.title = name || "Choose file";
    const requiredMessage = id === "lsa" ? "Select the Carterra LSA workbook." : "Select the NGS library CSV.";
    input.setCustomValidity(input.required && !name ? requiredMessage : "");
  }
  input.addEventListener("change", showSelection);
  showSelection();
}
let worker = null;
const urls = [];
function clearResults() {
  $("results").hidden = true;
  structureExplorer.clear();
  $("download").removeAttribute("href");
  while (urls.length) URL.revokeObjectURL(urls.pop());
}
function blobURL(data, type) {
  const url = URL.createObjectURL(new Blob([data], {type})); urls.push(url); return url;
}
function status(message, error=false) {
  $("state").textContent = message;
  $("state").classList.toggle("error", error);
}
function step(index) {
  [1,2,3].forEach(i => $("step"+i).classList.toggle("on",i<=index));
}
function unlock() {
  if (worker) worker.terminate();
  worker = null;
  [...form.elements].forEach(el => {el.disabled=false;});
  $("cancel").hidden = true;
}
function mapping() {
  $("plateMap").replaceChildren();
  const start = Number($("block").value.slice(1,3));
  if (!start) return;
  for (let i=0;i<4;i++) {
    const chip = document.createElement("span");
    chip.textContent = "P"+String(start+i).padStart(2,"0")+" → P0"+(i+1);
    $("plateMap").append(chip);
  }
}
function begin(demo) {
  if (worker) return;
  clearResults();
  const files = demo ? {} : Object.fromEntries(fileIds.map(key=>[key,$(key).files[0]]));
  if (!demo && (files.full || files.summary) && !files.cif) {
    status("Select the structure CIF file that matches the AlphaFold 3 JSON files.",true); return;
  }
  const settings = demo ? {} : {block:$("block").value, chain:$("chain").value.trim(), wt:$("wt").value.trim()};
  status("Loading the analysis engine…"); step(1);
  try {
    worker = new Worker("./worker.js" + assetQuery);
    [...form.elements].forEach(el => {el.disabled=el.id!=="cancel";});
    $("cancel").hidden = false;
    worker.onerror = () => {unlock(); status("The analysis engine could not start. Reload the page and try again.",true);};
    worker.onmessage = ({data}) => {
      if (data.stage) {status(data.stage); step(2); return;}
      unlock();
      if (data.error) {status(data.error,true); step(1); return;}
      const r = data.result;
      $("resultTitle").textContent = demo ? "Example results" : "Analysis results";
      $("metrics").replaceChildren();
      const items = [["Matched LSA rows",r.joined+" / "+r.rows],["Rows eligible for quantitation",r.quantitative],["Inactive rows retained",r.inactive]];
      if (r.mapped!=null) items.push(["Residues mapped to structure",r.mapped]);
      for (const [label,value] of items) {
        const div=document.createElement("div"), b=document.createElement("b"), s=document.createElement("span");
        b.textContent=value; s.textContent=label; div.append(b,s); $("metrics").append(div);
      }
      $("note").textContent=(r.block ? "NGS plate block "+r.block : "Using original NGS plate numbers")+
        " · "+r.warnings.length+(r.warnings.length===1 ? " parser warning" : " parser warnings")+(r.warnings.length ? " · "+r.warnings.slice(0,2).join(" / ") : "");
      $("report").textContent = r.report;
      if (typeof showAnalysisReport === "function") showAnalysisReport($("report"), r.report);
      $("download").href=blobURL(data.zip,"application/zip");
      $("download").download=demo ? "MAXTIA-demo-results.zip" : "MAXTIA-results.zip";
      $("results").hidden=false; step(3); status("Analysis complete");
      structureExplorer.show(r.viewer);
      $("results").scrollIntoView({block:"start"});
    };
    worker.postMessage({files,settings,demo});
  } catch (err) {unlock();status(String(err.message||err),true);}
}
form.addEventListener("submit",e=>{e.preventDefault();begin(false);});
form.addEventListener("input",()=>{clearResults();mapping();step($("lsa").files.length && $("ngs").files.length ? 1 : 0);status("Awaiting input");});
$("demo").onclick=()=>begin(true);
$("cancel").onclick=()=>{unlock();step(0);status("Run canceled.");};
