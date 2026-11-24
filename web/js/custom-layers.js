// Bilimsel veri katmanları: web/data/layers/ altındaki "*.layer.json" tanım dosyalarından
// haritaya otomatik katman, panel düğmesi, lejant ve açıklama penceresi üretir.
// Biçim için web/README.md › "Veri katmanı ekleme" bölümüne bakın.

const INDEX_URL = 'data/layers/index.json';
const BASE = 'data/layers/';
const DEFAULT_COLOR = '#8aa0bb';
const MISSING_COLOR = '#3a4659';

const esc = (s) => String(s ?? '').replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

const resolve = (p) => (/^(https?:)?\//.test(p) || p.startsWith('data/') ? p : BASE + p);

// MapLibre'nin vektör karo kodlaması çok büyük sayılarda (≥1e15, ör. N·m cinsinden moment) katmanı
// hiç çizmiyor. Bu değerleri metne çeviriyoruz; ifadeler 'to-number' ile yine sayı olarak okuyabiliyor.
function sanitize(data) {
  for (const f of data.features || []) {
    for (const [k, v] of Object.entries(f.properties || {})) {
      if (typeof v === 'number' && Math.abs(v) >= 1e15) f.properties[k] = v.toExponential(4);
      else if (v !== null && typeof v === 'object') f.properties[k] = JSON.stringify(v);
    }
  }
  return data;
}

// ---------------------------------------------------------------- İfadeler
function valueExpr(prop) {
  const props = Array.isArray(prop) ? prop : [prop];
  return props.length === 1 ? ['get', props[0]] : ['coalesce', ...props.map((p) => ['get', p])];
}

function colorExpr(color = {}) {
  const type = color.type || (color.categories ? 'categorical' : color.stops ? 'interpolate' : 'fixed');
  if (type === 'fixed' || !color.property) return color.value || DEFAULT_COLOR;
  const raw = valueExpr(color.property);
  const missing = color.missing || MISSING_COLOR;
  if (type === 'categorical') {
    const m = ['match', ['to-string', raw]];
    for (const [k, c] of Object.entries(color.categories || {})) m.push(k, c);
    m.push(color.default || DEFAULT_COLOR);
    return m;
  }
  const num = ['to-number', raw, 0];
  const stops = color.stops || [[0, '#3f7be0'], [1, '#ff3b30']];
  let e;
  if (type === 'step') {
    e = ['step', num, stops[0][1]];
    stops.slice(1).forEach(([v, c]) => e.push(v, c));
  } else {
    e = ['interpolate', ['linear'], num, ...stops.flat()];
  }
  return ['case', ['==', raw, null], missing, e];
}

function sizeExpr(size, fallback) {
  if (!size?.property) return size?.value ?? fallback;
  return ['interpolate', ['linear'], ['to-number', valueExpr(size.property), 0], ...(size.stops || [[0, 3], [1, 12]]).flat()];
}

const GEOM = {
  fill: ['in', ['geometry-type'], ['literal', ['Polygon', 'MultiPolygon']]],
  line: ['in', ['geometry-type'], ['literal', ['LineString', 'MultiLineString']]],
  circle: ['in', ['geometry-type'], ['literal', ['Point', 'MultiPoint']]],
};

function geometryKinds(data) {
  const kinds = new Set();
  for (const f of data.features || []) {
    const t = f.geometry?.type || '';
    if (t.includes('Polygon')) kinds.add('fill');
    else if (t.includes('LineString')) kinds.add('line');
    else if (t.includes('Point')) kinds.add('circle');
  }
  return [...kinds];
}

// ---------------------------------------------------------------- Katman kurulumu
function addMapLayers(map, def, data) {
  const id = 'x-' + def.id;
  const st = def.style || {};
  const opacity = st.opacity ?? 0.85;
  const color = colorExpr(def.color);
  const below = (kind) => (kind === 'fill' || kind === 'heatmap' ? 'faults-glow' : 'quakes-glow');
  const ids = [];
  map.addSource(id, { type: 'geojson', data });

  const add = (layer, kind) => {
    map.addLayer({ ...layer, source: id, layout: { ...(layer.layout || {}), visibility: 'none' } }, map.getLayer(below(kind)) ? below(kind) : undefined);
    ids.push(layer.id);
  };

  if (st.type === 'heatmap') {
    const stops = def.color?.stops || [[0, '#3f7be0'], [1, '#ff3b30']];
    const lo = stops[0][0], hi = stops[stops.length - 1][0];
    const weight = def.color?.property
      ? ['interpolate', ['linear'], ['to-number', valueExpr(def.color.property), lo], lo, 0, hi, 1]
      : 1;
    add({
      id: id + '-heat', type: 'heatmap', filter: GEOM.circle,
      paint: {
        'heatmap-weight': weight,
        'heatmap-radius': ['interpolate', ['linear'], ['zoom'], 4, st.radius ?? 18, 10, (st.radius ?? 18) * 3],
        'heatmap-opacity': opacity,
        'heatmap-intensity': st.intensity ?? 2,
        'heatmap-color': ['interpolate', ['linear'], ['heatmap-density'], 0, 'rgba(0,0,0,0)',
          ...stops.flatMap(([, c], i) => [0.05 + (0.95 * i) / Math.max(1, stops.length - 1), c])],
      },
    }, 'heatmap');
    return ids;
  }

  const kinds = st.type && st.type !== 'auto' ? [st.type] : geometryKinds(data);
  for (const kind of kinds) {
    if (kind === 'fill') {
      add({ id: id + '-fill', type: 'fill', filter: GEOM.fill, paint: { 'fill-color': color, 'fill-opacity': st.opacity ?? 0.45 } }, 'fill');
      add({ id: id + '-outline', type: 'line', filter: GEOM.fill, paint: { 'line-color': st.outline || color, 'line-width': st.outlineWidth ?? 0.6, 'line-opacity': 0.7 } }, 'fill');
    } else if (kind === 'line') {
      const w = sizeExpr(def.size, st.width ?? 2.4);
      if (st.glow !== false) {
        add({ id: id + '-glow', type: 'line', filter: GEOM.line, layout: { 'line-cap': 'round' }, paint: { 'line-color': color, 'line-width': ['*', 3.5, w], 'line-blur': 6, 'line-opacity': 0.3 } }, 'line');
      }
      add({ id: id + '-line', type: 'line', filter: GEOM.line, layout: { 'line-cap': 'round', 'line-join': 'round' }, paint: { 'line-color': color, 'line-width': w, 'line-opacity': opacity } }, 'line');
    } else if (kind === 'circle') {
      add({
        id: id + '-circle', type: 'circle', filter: GEOM.circle,
        paint: {
          'circle-color': color,
          'circle-radius': sizeExpr(def.size, st.radius ?? 5),
          'circle-opacity': opacity,
          'circle-stroke-color': st.outline || 'rgba(255,255,255,0.7)',
          'circle-stroke-width': st.outlineWidth ?? 1,
        },
      }, 'circle');
    }
  }
  return ids;
}

// ---------------------------------------------------------------- Lejant
function fmt(v, digits = 2) {
  if (v === null || v === undefined || v === '' || v === 'null') return '–';
  const n = typeof v === 'number' ? v : (/^-?[\d.]+(e[+-]?\d+)?$/i.test(v) ? +v : NaN);
  if (Number.isNaN(n)) return String(v);
  if (n !== 0 && (Math.abs(n) >= 1e6 || Math.abs(n) < 1e-3)) {
    const [m, e] = n.toExponential(2).split('e');
    return `${m.replace('.', ',')}×10^${+e}`;
  }
  return n.toLocaleString('tr-TR', { maximumFractionDigits: digits });
}

function legendHtml(def) {
  const c = def.color || {};
  const lg = def.legend || {};
  const unit = lg.unit ? ` ${esc(lg.unit)}` : '';
  const title = esc(lg.title || def.label);
  let body = '';
  const type = c.type || (c.categories ? 'categorical' : c.stops ? 'interpolate' : 'fixed');
  if (type === 'categorical') {
    const labels = lg.categoryLabels || {};
    body = `<div class="keys">${Object.entries(c.categories || {}).map(([k, col]) =>
      `<span><i class="dot" style="background:${esc(col)}"></i>${esc(labels[k] || k)}</span>`).join('')}</div>`;
  } else if (type === 'fixed' || !c.property) {
    body = `<div class="keys"><span><i class="dot" style="background:${esc(c.value || DEFAULT_COLOR)}"></i>${esc(lg.itemLabel || def.label)}</span></div>`;
  } else {
    const stops = c.stops || [];
    const grad = type === 'step'
      ? stops.flatMap(([, col], i) => [`${col} ${(100 * i) / stops.length}%`, `${col} ${(100 * (i + 1)) / stops.length}%`]).join(',')
      : stops.map(([, col]) => col).join(',');
    const labels = lg.labels || [stops[0]?.[0], stops[Math.floor((stops.length - 1) / 2)]?.[0], stops[stops.length - 1]?.[0]].map((v) => fmt(v));
    body = `<div class="ramp" style="background:linear-gradient(90deg,${esc(grad)})"></div>
      <div class="ramp-lbl">${labels.map((l, i) => `<span>${esc(l)}${i === labels.length - 1 ? unit : ''}</span>`).join('')}</div>`;
  }
  const note = def.note ? `<p class="lg-note">${esc(def.note)}</p>` : '';
  return `<div class="lg"><h3>${title}</h3>${body}${note}</div>`;
}

// ---------------------------------------------------------------- Açıklama penceresi
function tooltipHtml(def, p) {
  const tt = def.tooltip || {};
  const titleProp = tt.title;
  const title = (titleProp && p[titleProp]) || p.name || p.id || def.label;
  const fields = tt.fields || Object.keys(p).filter((k) => k !== titleProp).slice(0, 8).map((k) => ({ property: k }));
  const rows = fields.map((f) => {
    const v = p[f.property];
    if (v === undefined) return '';
    return `<dt>${esc(f.label || f.property)}</dt><dd>${esc(fmt(v, f.digits ?? 2))}${f.unit ? ' ' + esc(f.unit) : ''}</dd>`;
  }).join('');
  return `<div class="pop-place" style="margin:0 0 6px">${esc(title)}</div>
    <div class="pop-tag">${esc(def.label)}</div>
    <dl class="pop-grid">${rows}</dl>${def.note ? `<p class="lg-note">${esc(def.note)}</p>` : ''}`;
}

// ---------------------------------------------------------------- Yükleyici
/**
 * Dizin dosyasındaki her katman tanımını yükler, haritaya ekler ve panel için OVERLAYS biçiminde girişler döndürür.
 * Dizin yoksa boş liste döner.
 */
export async function loadCustomLayers(map) {
  let index;
  try {
    const r = await fetch(INDEX_URL, { cache: 'no-cache' });
    if (!r.ok) return [];
    index = await r.json();
  } catch (_) {
    return [];
  }
  const files = Array.isArray(index) ? index : index.layers || [];
  const defs = (await Promise.all(files.map(async (file) => {
    try {
      const def = await (await fetch(resolve(file), { cache: 'no-cache' })).json();
      def.id = String(def.id || file.replace(/^.*\//, '').replace(/\.layer\.json$/, '')).replace(/[^\w-]/g, '_');
      return def;
    } catch (e) {
      console.warn('Katman tanımı okunamadı:', file, e);
      return null;
    }
  }))).filter(Boolean).sort((a, b) => (a.order ?? 100) - (b.order ?? 100) || String(a.label).localeCompare(String(b.label), 'tr'));

  const popup = new maplibregl.Popup({ closeButton: false, closeOnClick: false, offset: 8, maxWidth: '320px' });
  const out = [];
  for (const def of defs) {
    const entry = {
      id: 'x-' + def.id,
      label: def.label || def.id,
      desc: def.description || '',
      group: def.group || 'Model katmanları',
      visible: !!def.visible,
      soon: true,
      mapLayers: [],
      legend: () => legendHtml(def),
    };
    out.push(entry);
    if (!def.data) continue;
    try {
      const data = sanitize(await (await fetch(resolve(def.data), { cache: 'no-cache' })).json());
      def.note = def.note || data.metadata?.note || '';
      entry.mapLayers = addMapLayers(map, def, data);
      entry.soon = false;
      const n = data.features?.length ?? 0;
      entry.desc = [def.description, `${n.toLocaleString('tr-TR')} öğe`, def.updated && `güncelleme ${def.updated}`].filter(Boolean).join(' · ');
      if (def.tooltip !== false) {
        for (const lid of entry.mapLayers) {
          if (lid.endsWith('-glow') || lid.endsWith('-heat')) continue;
          map.on('mousemove', lid, (ev) => {
            map.getCanvas().style.cursor = 'pointer';
            popup.setLngLat(ev.lngLat).setHTML(tooltipHtml(def, ev.features[0].properties)).addTo(map);
          });
          map.on('mouseleave', lid, () => { map.getCanvas().style.cursor = ''; popup.remove(); });
        }
      }
    } catch (e) {
      console.warn('Katman verisi okunamadı:', def.data, e);
      entry.desc = 'Veri dosyası okunamadı';
    }
  }
  return out;
}
