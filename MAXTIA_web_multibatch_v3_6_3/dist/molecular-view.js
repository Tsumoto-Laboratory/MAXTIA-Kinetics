"use strict";

// Build the cartoon with 3Dmol's actual secondary-structure and mesh routines.
// A depth-buffered Canvas 2D image uses the same surface, camera and light on
// machines with or without WebGL. Each visible pixel also records its residue.
window.createMolecularView = function(host, geometry, onPick) {
  const model = new window.$3Dmol.GLModel(0, {cartoonQuality: 6});
  model.addMolData(geometry, "cif");
  const atoms = model.selectedAtoms({});
  if (!atoms.length) throw Error("No drawable coordinates");

  const ids = [null], slotOf = new Map(), anchors = new Map();
  const center = [0, 0, 0];
  for (const a of atoms) {
    if (!slotOf.has(a.resi)) {
      slotOf.set(a.resi, ids.length);
      ids.push(a.resi);
    }
    if (!anchors.has(a.resi) || a.atom === "CA") anchors.set(a.resi, a.index);
    center[0] += a.x;
    center[1] += a.y;
    center[2] += a.z;
  }
  for (let k = 0; k < 3; k++) center[k] /= atoms.length;
  const radius = Math.max(
      2,
      atoms.reduce((r, a) => Math.max(r, Math.hypot(a.x - center[0], a.y - center[1], a.z - center[2])), 0));
  const points = new Float32Array(3 * atoms.length);
  const slot = new Uint32Array(atoms.length), atomRadius = new Float32Array(atoms.length);
  const alphaCarbons = [], bonds = [], atomIndex = new Map(atoms.map((a, i) => [a.index, i]));
  atoms.forEach((a, i) => {
    const j = 3 * i;
    points[j] = a.x - center[0];
    points[j + 1] = a.y - center[1];
    points[j + 2] = a.z - center[2];
    slot[i] = slotOf.get(a.resi);
    atomRadius[i] = .23 * (window.$3Dmol.GLModel.vdwRadii?.[a.elem] || 1.5);
    if (a.atom === "CA") {
      alphaCarbons.push(i);
      atomRadius[i] = .55;
    }
    for (const other of a.bonds || []) {
      const b = atomIndex.get(other);
      if (b > i) bonds.push(i, b);
    }
  });

  // Encode the polymer picking key as an integer color while extracting the
  // mesh. The experimental palette can then change without rebuilding geometry.
  model.setStyle({}, {cartoon: {colorfunc: a => slotOf.get(a.resi)}});
  const meshes = [];
  function collect(node) {
    for (const g of node.geometry?.geometryGroups || []) {
      if (!g.faceidx) continue;
      const n = g.vertices, xyz = g.vertexArray.slice(0, n * 3), normals = g.normalArray.slice(0, n * 3);
      const slots = new Uint32Array(n);
      for (let i = 0; i < n; i++) {
        const j = i * 3, colors = g.colorArray;
        xyz[j] -= center[0];
        xyz[j + 1] -= center[1];
        xyz[j + 2] -= center[2];
        slots[i] = (Math.round(colors[j] * 255) << 16) | (Math.round(colors[j + 1] * 255) << 8) |
            Math.round(colors[j + 2] * 255);
        const norm = Math.hypot(normals[j], normals[j + 1], normals[j + 2]) || 1;
        normals[j] /= norm;
        normals[j + 1] /= norm;
        normals[j + 2] /= norm;
      }
      meshes.push({
        xyz,
        normals,
        slots,
        indices: g.faceArray.slice(0, g.faceidx),
        screen: new Float32Array(3 * n),
        light: new Float32Array(n)
      });
    }
    for (const child of node.children || []) collect(child);
  }
  model.globj({add: collect, remove() {}}, {supportsImposters: false, supportsAIA: false});

  const canvas = document.createElement("canvas"), label = document.createElement("span");
  canvas.setAttribute("role", "img");
  canvas.setAttribute("aria-label", "Interactive three-dimensional protein structure");
  canvas.className = "molecular-canvas";
  label.className = "molecular-label";
  label.hidden = true;
  const context = canvas.getContext("2d", {alpha: false});
  if (!context) throw Error("Canvas 2D unavailable");
  host.append(canvas, label);
  const identity = () => [1, 0, 0, 0, 1, 0, 0, 0, 1];
  let matrix = identity(), zoom = 1, atomic = false, selected = null, frame = 0, disposed = false;
  let w = 0, h = 0, scale = 1, image, pixels, depth, hits;
  const screenAtoms = new Float32Array(points.length);
  const palette = new Uint32Array(ids.length * 256);
  const pack = (r, g, b) => (255 << 24) | (b << 16) | (g << 8) | r;
  const background = pack(242, 246, 247);
  function setPalette(map) {
    ids.forEach((id, idx) => {
      const color = parseInt((map.get(id)?.color || "#808080").slice(1), 16);
      const r = color >> 16, g = (color >> 8) & 255, b = color & 255;
      for (let l = 0; l < 256; l++) {
        const light = l / 255, diffuse = .39 + .58 * light, shine = 21 * Math.pow(light, 18);
        const channel = base => Math.min(255, Math.round(base * diffuse + shine));
        palette[idx * 256 + l] = pack(channel(r), channel(g), channel(b));
      }
    });
  }
  // Right-handed model coordinates, positive Z toward the observer. Using the
  // same scale on all three axes also gives atoms and bonds round cross sections.
  function project(source, dest) {
    const m = matrix;
    for (let j = 0; j < source.length; j += 3) {
      const x = source[j], y = source[j + 1], z = source[j + 2];
      dest[j] = w * .5 + scale * (m[0] * x + m[1] * y + m[2] * z);
      dest[j + 1] = h * .5 - scale * (m[3] * x + m[4] * y + m[5] * z);
      dest[j + 2] = scale * (m[6] * x + m[7] * y + m[8] * z);
    }
  }
  function lightIndex(x, y, z) {
    return Math.max(0, Math.min(255, Math.round(255 * (.28 * x + .50 * y + .82 * z))));
  }
  function drawMesh(mesh) {
    const {xyz, normals, slots, indices, screen, light} = mesh, m = matrix;
    project(xyz, screen);
    for (let i = 0, j = 0; i < slots.length; i++, j += 3) {
      const nx = normals[j], ny = normals[j + 1], nz = normals[j + 2];
      const z = m[6] * nx + m[7] * ny + m[8] * nz, sign = z < 0 ? -1 : 1;
      light[i] = lightIndex(
          sign * (m[0] * nx + m[1] * ny + m[2] * nz), sign * (m[3] * nx + m[4] * ny + m[5] * nz),
          Math.abs(z));
    }
    for (let f = 0; f < indices.length; f += 3) {
      const a = indices[f], b = indices[f + 1], c = indices[f + 2], ia = 3 * a, ib = 3 * b, ic = 3 * c;
      const ax = screen[ia], ay = screen[ia + 1], bx = screen[ib], by = screen[ib + 1], cx = screen[ic],
            cy = screen[ic + 1];
      const dx1 = bx - ax, dy1 = by - ay, dx2 = cx - ax, dy2 = cy - ay, det = dx1 * dy2 - dy1 * dx2;
      if (Math.abs(det) < 1e-6) continue;
      const left = Math.max(0, Math.ceil(Math.min(ax, bx, cx) - .5)),
            right = Math.min(w - 1, Math.floor(Math.max(ax, bx, cx) - .5));
      const top = Math.max(0, Math.ceil(Math.min(ay, by, cy) - .5)),
            bottom = Math.min(h - 1, Math.floor(Math.max(ay, by, cy) - .5));
      const du = dy2 / det, dv = -dy1 / det, az = screen[ia + 2], dz1 = screen[ib + 2] - az,
            dz2 = screen[ic + 2] - az;
      for (let y = top; y <= bottom; y++) {
        const px = left + .5 - ax, py = y + .5 - ay;
        let u = (px * dy2 - py * dx2) / det, v = (py * dx1 - px * dy1) / det;
        for (let x = left, i = y * w + left; x <= right; x++, i++, u += du, v += dv) {
          const wa = 1 - u - v;
          if (u < 0 || v < 0 || wa < 0) continue;
          const z = az + u * dz1 + v * dz2;
          if (z <= depth[i]) continue;
          const id = slots[wa >= u && wa >= v ? a : u >= v ? b : c] || slots[a] || slots[b] || slots[c];
          const l = Math.round(wa * light[a] + u * light[b] + v * light[c]);
          depth[i] = z;
          hits[i] = id;
          pixels[i] = palette[id * 256 + l];
        }
      }
    }
  }
  function sphere(i, radiusInAngstroms) {
    const j = 3 * i, cx = screenAtoms[j], cy = screenAtoms[j + 1], cz = screenAtoms[j + 2],
          r = Math.max(.65, radiusInAngstroms * scale);
    const left = Math.max(0, Math.ceil(cx - r - .5)), right = Math.min(w - 1, Math.floor(cx + r - .5));
    const top = Math.max(0, Math.ceil(cy - r - .5)), bottom = Math.min(h - 1, Math.floor(cy + r - .5));
    for (let y = top; y <= bottom; y++)
      for (let x = left, p = y * w + left; x <= right; x++, p++) {
        const nx = (x + .5 - cx) / r, ny = (cy - y - .5) / r, q = 1 - nx * nx - ny * ny;
        if (q < 0) continue;
        const nz = Math.sqrt(q), z = cz + r * nz;
        if (z <= depth[p]) continue;
        depth[p] = z;
        hits[p] = slot[i];
        pixels[p] = palette[slot[i] * 256 + lightIndex(nx, ny, nz)];
      }
  }
  function bond(a, b) {
    const i = 3 * a, j = 3 * b, ax = screenAtoms[i], ay = screenAtoms[i + 1], az = screenAtoms[i + 2];
    const bx = screenAtoms[j], by = screenAtoms[j + 1], dx = bx - ax, dy = by - ay,
          dz = screenAtoms[j + 2] - az, len2 = dx * dx + dy * dy;
    if (len2 < .01) return;
    const r = Math.max(.65, .17 * scale), left = Math.max(0, Math.floor(Math.min(ax, bx) - r)),
          right = Math.min(w - 1, Math.ceil(Math.max(ax, bx) + r));
    const top = Math.max(0, Math.floor(Math.min(ay, by) - r)),
          bottom = Math.min(h - 1, Math.ceil(Math.max(ay, by) + r));
    for (let y = top; y <= bottom; y++)
      for (let x = left, p = y * w + left; x <= right; x++, p++) {
        const ex = x + .5 - ax, ey = y + .5 - ay, t = (ex * dx + ey * dy) / len2;
        if (t < 0 || t > 1) continue;
        const nx = (ex - t * dx) / r, ny = (t * dy - ey) / r, q = 1 - nx * nx - ny * ny;
        if (q < 0) continue;
        const nz = Math.sqrt(q), z = az + t * dz + r * nz;
        if (z <= depth[p]) continue;
        const id = t < .5 ? slot[a] : slot[b];
        depth[p] = z;
        hits[p] = id;
        pixels[p] = palette[id * 256 + lightIndex(nx, ny, nz)];
      }
  }
  function draw() {
    frame = 0;
    const cw = host.clientWidth, ch = host.clientHeight;
    if (!cw || !ch || disposed) return;
    // Fixed sampling density for consistent edge quality on Retina and standard displays.
    const density = Math.min(1.6, Math.sqrt(1100000 / (cw * ch)));
    const width = Math.max(1, Math.round(cw * density)), height = Math.max(1, Math.round(ch * density));
    if (width !== w || height !== h) {
      canvas.width = w = width;
      canvas.height = h = height;
      image = context.createImageData(w, h);
      pixels = new Uint32Array(image.data.buffer);
      depth = new Float32Array(w * h);
      hits = new Uint32Array(w * h);
    }
    pixels.fill(background);
    hits.fill(0);
    depth.fill(-Infinity);
    scale = Math.min(w, h) * .46 * zoom / radius;
    project(points, screenAtoms);
    if (atomic) {
      for (let i = 0; i < bonds.length; i += 2) bond(bonds[i], bonds[i + 1]);
      for (let i = 0; i < atoms.length; i++) sphere(i, atomRadius[i]);
    } else {
      meshes.forEach(drawMesh);
      if (alphaCarbons.length)
        for (const i of alphaCarbons) sphere(i, .55);
      else if (!meshes.length)
        for (let i = 0; i < atoms.length; i++) sphere(i, atomRadius[i]);
    }
    context.putImageData(image, 0, 0);
    placeLabel();
  }
  function schedule() {
    if (!frame && !disposed) frame = requestAnimationFrame(draw);
  }
  function placeLabel() {
    const i = selected ? atomIndex.get(anchors.get(selected.model)) : undefined;
    label.hidden = i === undefined || !w || !h;
    if (label.hidden) return;
    label.textContent = selected.aa + selected.pos;
    label.style.left = (100 * screenAtoms[3 * i] / w) + "%";
    label.style.top = (100 * screenAtoms[3 * i + 1] / h) + "%";
  }
  function hit(e, pin) {
    if (!hits) return;
    const r = canvas.getBoundingClientRect();
    const x = Math.floor((e.clientX - r.left) * w / r.width),
          y = Math.floor((e.clientY - r.top) * h / r.height);
    if (x >= 0 && y >= 0 && x < w && y < h && hits[y * w + x]) onPick(ids[hits[y * w + x]], pin);
  }
  function rotate(dx, dy) {
    // Screen Y grows downward, so positive pitch moves the front surface down.
    const ax = dy * .012, ay = dx * .012, cx = Math.cos(ax), sx = Math.sin(ax), cy = Math.cos(ay),
          sy = Math.sin(ay);
    const old = matrix.slice();
    for (let i = 0; i < 3; i++) {
      const x = old[i], y = old[i + 3], z = old[i + 6];
      matrix[i] = cy * x + sy * z;
      matrix[i + 3] = cx * y - sx * (-sy * x + cy * z);
      matrix[i + 6] = sx * y + cx * (-sy * x + cy * z);
    }
    schedule();
  }
  const pointers = new Map();
  let moved = false, start = null;
  function span() {
    const [a, b] = pointers.values();
    return b ? Math.hypot(b.x - a.x, b.y - a.y) : 0;
  }
  function magnify(factor) {
    zoom = Math.max(.22, Math.min(12, zoom * factor));
    schedule();
  }
  canvas.onpointerdown = e => {
    if (e.button !== 0) return;
    const p = {x: e.clientX, y: e.clientY};
    if (!pointers.size) {
      start = p;
      moved = false;
    } else
      moved = true;
    pointers.set(e.pointerId, p);
    canvas.setPointerCapture(e.pointerId);
  };
  canvas.onpointermove = e => {
    const old = pointers.get(e.pointerId);
    if (!old) {
      if (e.pointerType !== "touch") hit(e, false);
      return;
    }
    const before = span(), p = {x: e.clientX, y: e.clientY};
    pointers.set(e.pointerId, p);
    if (before) {
      magnify(span() / before);
      return;
    }
    if (Math.hypot(p.x - start.x, p.y - start.y) > 3) moved = true;
    if (moved) rotate(p.x - old.x, p.y - old.y);
  };
  canvas.onpointerup = e => {
    if (!pointers.has(e.pointerId)) return;
    if (!moved && pointers.size === 1) hit(e, true);
    pointers.delete(e.pointerId);
  };
  canvas.onpointercancel = e => {
    pointers.delete(e.pointerId);
    moved = true;
  };
  canvas.addEventListener("wheel", e => {
    e.preventDefault();
    magnify(Math.exp(-Math.sign(e.deltaY) * .1));
  }, {passive: false});
  const resize = new ResizeObserver(schedule);
  resize.observe(host);
  return {
    style(map, showAtoms) {
      setPalette(map);
      atomic = showAtoms;
      schedule();
    },
    mark(r) {
      selected = r;
      placeLabel();
    },
    reset() {
      matrix = identity();
      zoom = 1;
      schedule();
    },
    destroy() {
      disposed = true;
      cancelAnimationFrame(frame);
      resize.disconnect();
      canvas.remove();
      label.remove();
    }
  };
};
