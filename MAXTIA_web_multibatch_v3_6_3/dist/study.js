"use strict";
const assetQuery = new URL(document.currentScript.src).search;

const $ = id => document.getElementById(id);
const form = $("form");
let worker = null;
let result = null;
let downloadURL = "";
let nextBatch = 1;


function clearResults() {
  result = null;
  structureExplorer.clear();
  $("results").hidden = true;
  $("download").removeAttribute("href");
  if (downloadURL) URL.revokeObjectURL(downloadURL);
  downloadURL = "";
  $("effects").replaceChildren();
  $("report").textContent = "";
}

function setStatus(text, error=false) {
  $("state").textContent = text;
  $("state").classList.toggle("error", error);
}

function updateRemove() {
  const buttons = $("batches").querySelectorAll(".remove");
  buttons.forEach(button => {button.hidden = buttons.length <= 2;});
}

function setupPickers(scope) {
  for (const input of scope.querySelectorAll('input[type="file"]')) {
    const wrapper = document.createElement("span");
    wrapper.className = "file-picker";
    const button = document.createElement("span");
    button.className = "file-button";
    button.textContent = "Choose file";
    const name = document.createElement("span");
    name.className = "filename";
    button.setAttribute("aria-hidden", "true");
    name.setAttribute("aria-hidden", "true");
    input.before(wrapper);
    wrapper.append(button, name, input);
    function selected() {
      name.textContent = input.files[0]?.name || "No file selected";
      input.title = input.files[0]?.name || "Choose file";
      input.setCustomValidity(input.required && !input.files.length ? "Select this batch input file." : "");
    }
    input.addEventListener("change", selected);
    selected();
  }
}

function addBatch() {
  const card = $("batchTemplate").content.firstElementChild.cloneNode(true);
  const name = "Batch " + nextBatch++;
  card.querySelector("h2").textContent = name;
  const id = card.querySelector('[data-key="id"]');
  id.value = name;
  id.addEventListener("input", () => {card.querySelector("h2").textContent=id.value || name;});
  card.querySelector(".remove").onclick = () => {
    card.remove(); updateRemove(); clearResults(); setStatus("Awaiting input");
  };
  setupPickers(card);
  $("batches").append(card);
  updateRemove(); clearResults(); setStatus("Awaiting input");
}

function unlock() {
  if (worker) worker.terminate();
  worker = null;
  for (const el of form.elements) el.disabled = false;
  $("cancel").hidden = true;
}

function display(value) {
  if (value === null || value === undefined) return "—";
  if (typeof value === "number") return Number(value.toPrecision(4)).toString();
  return String(value);
}

function drawTable() {
  if (!result) return;
  const view = $("view").value;
  const columns = [["Variant","variant"]];
  if (view === "batch_rows") columns.push(["Batch","batch_id"]);
  if (view !== "overview") columns.push(["Group","condition_id"],["Condition label","condition"]);
  if (view === "batch_rows") {
    columns.push(["Observed ROIs","observed_rois"],["ROIs used for medians","median_rois"],
      ["Median koff (s⁻¹)","median_koff_s_inv"],["Median kon (M⁻¹ s⁻¹)","median_kon_M_inv_s"],["Median Rmax (RU)","median_Rmax_RU"],
      ["PASS active","n_pass_active"],["Inactive","n_inactive"],["WT QC eligible","usable_wt"],
      ["KD fold","KD_fold"],["kon fold","kon_fold"],["koff fold","koff_fold"],["Variant QC","variant_qc_status"],["Phenotype","phenotype"]);
  } else {
    columns.push(["Observed batches","batches_observed"],["Eligible batches","batches_with_fold"],
      ["KD fold","KD_fold_center"],["KD min","KD_fold_min"],["KD max","KD_fold_max"],
      ["kon fold","kon_fold_center"],["koff fold","koff_fold_center"],["Interpretation","interpretation"]);
  }
  const table = $("effects");
  table.replaceChildren();
  const head = table.createTHead().insertRow();
  for (const [title] of columns) {
    const th = document.createElement("th"); th.scope="col"; th.textContent=title; head.append(th);
  }
  const body = table.createTBody();
  const search = $("find").value.trim().toLowerCase();
  for (const record of result[view]) {
    if (!record.variant.toLowerCase().includes(search)) continue;
    const row = body.insertRow();
    for (const [,key] of columns) row.insertCell().textContent=display(record[key]);
  }
  if (!body.rows.length) {
    const cell=body.insertRow().insertCell(); cell.colSpan=columns.length;
    cell.textContent="No matching variants.";
  }
}

function start(demo) {
  if (worker) return;
  clearResults();
  const shared = {};
  for (const role of ["cif","full","summary"]) shared[role]=$(role).files[0];
  if (!demo && (shared.full || shared.summary) && !shared.cif) {
    setStatus("Select the structure CIF from the same AlphaFold model.",true); return;
  }
  const batches = [...$("batches").children].map(card => {
    const batch = {files:{}};
    for (const input of card.querySelectorAll("[data-key]")) {
      if (input.type === "file") batch.files[input.dataset.key]=input.files[0];
      else batch[input.dataset.key]=input.value.trim();
    }
    return batch;
  });
  try {
    worker = new Worker("./study-worker.js" + assetQuery);
    for (const el of form.elements) el.disabled = el.id !== "cancel";
    $("cancel").hidden=false;
    setStatus("Loading the analysis engine…");
    worker.onerror=() => {unlock(); setStatus("Could not start analysis. Reload the page and try again.",true);};
    worker.onmessage=({data}) => {
      if (data.stage) {setStatus(data.stage); return;}
      unlock();
      if (data.error) {setStatus(data.error,true); return;}
      result=data.result;
      $("resultTitle").textContent = demo ? "Example study results" : "Study results";
      $("metrics").replaceChildren();
      for (const [key,label] of [["batches","Batches"],["wells","LSA rows"],["variants","Distinct variants"],["repeated","Quantified in multiple batches"]]) {
        const div=document.createElement("div"), n=document.createElement("b"), s=document.createElement("span");
        n.textContent=result[key]; s.textContent=label; div.append(n,s); $("metrics").append(div);
      }
      downloadURL=URL.createObjectURL(new Blob([data.zip],{type:"application/zip"}));
      $("download").href=downloadURL;
      $("report").textContent = result.report;
      if (typeof showAnalysisReport === "function") showAnalysisReport($("report"), result.report);
      $("find").value=""; $("view").value="overview"; drawTable();
      $("results").hidden=false; setStatus("Analysis complete");
      structureExplorer.show(result.viewer);
      $("results").scrollIntoView({block:"start"});
    };
    worker.postMessage(demo ? {demo:true} : {demo:false, shared, batches, partner:$("partner").value.trim(),
      chain:$("chain").value.trim(), wt:$("wt").value.trim()});
  } catch (e) {unlock(); setStatus(e.message || "Could not start analysis.",true);}
}

form.addEventListener("submit", e=>{e.preventDefault(); start(false);});
form.addEventListener("input", ()=>{clearResults(); setStatus("Awaiting input");});
$("add").onclick=addBatch;
$("demo").onclick=()=>start(true);
$("cancel").onclick=()=>{unlock(); setStatus("Run canceled.");};
$("view").onchange=drawTable;
$("find").oninput=drawTable;
setupPickers($("form"));
addBatch(); addBatch();
