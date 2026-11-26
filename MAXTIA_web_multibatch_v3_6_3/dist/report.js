"use strict";

// Render the headings, lists and pipe tables emitted by the analysis engines.
// All experiment-derived content is inserted as text.
function showAnalysisReport(container, markdown) {
  const lines = markdown.replace(/\r\n?/g, "\n").split("\n");
  const report = document.createDocumentFragment();
  let heading = "Analysis report";

  function appendText(element, value) {
    for (const part of value.split(/(`[^`]+`)/g)) {
      if (part.length > 2 && part.startsWith("`") && part.endsWith("`")) {
        const code = document.createElement("code");
        code.textContent = part.slice(1, -1);
        element.append(code);
      } else {
        element.append(part);
      }
    }
  }

  function tableRow(section, line, header) {
    const row = section.insertRow();
    const fields = line.trim().slice(1, -1).split(/(\\\||\|)/);
    let value = "";
    function addCell() {
      const cell = document.createElement(header ? "th" : "td");
      if (header) cell.scope = "col";
      appendText(cell, value.trim());
      row.append(cell);
      value = "";
    }
    for (const field of fields) {
      if (field === "|") addCell();
      else value += field === "\\|" ? "|" : field;
    }
    addCell();
  }

  let i = 0;
  while (i < lines.length) {
    const line = lines[i++].trim();
    if (!line) continue;
    const title = /^(#{1,6})\s+(.+)$/.exec(line);
    if (title) {
      heading = title[2];
      const node = document.createElement("h" + Math.min(6, title[1].length + 2));
      appendText(node, heading);
      report.append(node);
    } else if (line.startsWith("|") && /^\|(?:\s*:?-{3,}:?\s*\|)+$/.test((lines[i] || "").trim())) {
      const scroll = document.createElement("div");
      scroll.className = "scroll-table";
      scroll.tabIndex = 0;
      scroll.setAttribute("role", "region");
      scroll.setAttribute("aria-label", heading);
      const table = document.createElement("table");
      tableRow(table.createTHead(), line, true);
      i++;
      const body = table.createTBody();
      while (i < lines.length && lines[i].trim().startsWith("|")) tableRow(body, lines[i++], false);
      scroll.append(table);
      report.append(scroll);
    } else if (line.startsWith("- ")) {
      const list = document.createElement("ul");
      let item = line;
      while (true) {
        const li = document.createElement("li");
        appendText(li, item.slice(2));
        list.append(li);
        if (i >= lines.length || !lines[i].trim().startsWith("- ")) break;
        item = lines[i++].trim();
      }
      report.append(list);
    } else {
      const p = document.createElement("p");
      appendText(p, line);
      report.append(p);
    }
  }
  container.replaceChildren(report);
}
