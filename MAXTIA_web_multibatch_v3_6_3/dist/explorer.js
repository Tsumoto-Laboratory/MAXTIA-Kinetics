"use strict";

// The legend, residue inspector and renderer share the same polymer-to-WT mapping.
window.structureExplorer = (() => {
  const panel = document.getElementById("modelPanel");
  panel.innerHTML = `
    <div class="map-heading"><h3>Structure explorer</h3><div class="map-options"><label>Representation<select class="representation"><option>Ribbon</option><option>Atoms and bonds</option></select></label><label>Map<select class="scope-picker"></select></label></div></div>
    <div class="map-layout"><div class="map-left">
      <div class="molecule" aria-label="Interactive three-dimensional protein structure"></div>
      <div class="map-tools"><span>Drag to rotate · Scroll or pinch to zoom</span><button type="button" class="secondary reset-camera">Reset view</button></div>
      <p class="render-note" role="status" hidden></p>
      <div class="map-legend"><span class="gradient"></span><span>Normalized ΔΔG &nbsp; −1 blue · 0 white · +1 red</span></div>
      <div class="map-keys"></div>
    </div><div class="map-detail">
      <label>WT residue<select class="residue-picker"></select></label>
      <div class="inspection"></div><button type="button" class="secondary unpin" hidden>Resume hover</button>
    </div></div>`;
  const get = selector => panel.querySelector(selector);
  const viewport=get(".molecule"), pick=get(".residue-picker"), scopePick=get(".scope-picker");
  const details=get(".inspection"), unpin=get(".unpin"), reset=get(".reset-camera");
  const names={NO_DATA:"No data", REVIEW:"Review", INACTIVE_ONLY:"Inactive only",
    NO_PASS_DATA:"No PASS data", QUANTITATIVE:"Quantitative"};
  for (const [key,color] of Object.entries({NO_DATA:"#808080",REVIEW:"#800080",INACTIVE_ONLY:"#ffff00",NO_PASS_DATA:"#ff00ff"})) {
    const item=document.createElement("span"), chip=document.createElement("i");
    chip.style.backgroundColor=color; item.append(chip,names[key]); get(".map-keys").append(item);
  }
  let view=null, payload=null, pinned=false, selected=null;
  const note=get(".render-note");
  let byModel=new Map(), byPosition=new Map();
  function add(parent, tag, value) {
    const node=document.createElement(tag); node.textContent=value; parent.append(node); return node;
  }
  function inspect(r, pin=false) {
    if (!r || (pinned && !pin)) return;
    if (pin) pinned=true;
    unpin.hidden=!pinned;
    if (selected===r && !pin) return;
    selected=r; pick.value=String(r.pos); details.replaceChildren();
    add(details,"h4",r.aa+r.pos);
    add(details,"p",names[r.status]+" · Model "+r.author);
    if (r.selected) add(details,"p",`Color from ${r.selected} · ΔΔG ${Number(r.energy.toPrecision(3))} kcal/mol`);
    if (!r.effects.length) add(details,"p","No sequence-verified substitution was observed here in this map.");
    for (const e of r.effects) {
      const entry=document.createElement("div"); entry.className="variant-card"; details.append(entry);
      add(entry,"strong",e.variant);
      if (e.multi) add(entry,"p","Multi-substitution clone · Individual residue effect unresolved.");
      const ratios=e.ratios.map(x=>x==null ? "—" : Number(x.toPrecision(3))+"×");
      add(entry,"p",`KD ${ratios[0]} · kon ${ratios[1]} · koff ${ratios[2]}`);
      add(entry,"small",[names[e.status], e.support, e.note?.replaceAll("_"," ")].filter(Boolean).join(" · "));
    }
    if (r.effects.length) add(details,"small","Ratios are relative to WT. A dash marks an unavailable estimate.");
    view?.mark(r);
  }
  function defaultDetail() {
    pinned=false; selected=null; unpin.hidden=true; view?.mark(null); pick.value="";
    details.replaceChildren(); add(details,"p","Hover over a residue to see its variants. Click a residue or choose one above to pin its measurements.");
  }
  function switchScope() {
    const previous=pinned ? selected?.pos : null;
    const scope=payload.scopes[Number(scopePick.value)];
    byModel=new Map(scope.residues.map(r=>[r.model,r]));
    byPosition=new Map(scope.residues.map(r=>[r.pos,r]));
    pick.replaceChildren(new Option("Choose a residue…",""));
    for (const r of scope.residues) pick.add(new Option(r.aa+r.pos+" · "+names[r.status],String(r.pos)));
    defaultDetail();
    const atomic=get(".representation").value==="Atoms and bonds";
    view?.style(byModel,atomic);
    if (previous && byPosition.has(previous)) inspect(byPosition.get(previous),true);
  }
  scopePick.onchange=switchScope;
  get(".representation").onchange=switchScope;
  pick.onchange=()=>{if (pick.value) inspect(byPosition.get(Number(pick.value)),true); else defaultDetail();};
  unpin.onclick=defaultDetail;
  reset.onclick=()=>view?.reset();
  function clear() {
    view?.destroy(); view=null; payload=null; pinned=false; selected=null;
    panel.hidden=true; byModel.clear(); byPosition.clear(); details.replaceChildren();
    note.hidden=true; note.textContent="";
  }
  function show(data) {
    clear();
    if (!data || !data.scopes.length) return;
    payload=data; panel.hidden=false;
    scopePick.replaceChildren(...data.scopes.map((s,i)=>new Option(s.label,String(i))));
    scopePick.parentElement.hidden=data.scopes.length===1;
    reset.disabled=false;
    try {
      view=window.createMolecularView(viewport,data.geometry,(id,pin)=>inspect(byModel.get(id),pin));
    } catch (err) {
      reset.disabled=true; note.hidden=false;
      note.textContent="The 3D image could not be drawn. Choose a WT residue to inspect its measurements.";
      console.error("Structure rendering failed",err);
    }
    switchScope();
  }
  return {show,clear};
})();
