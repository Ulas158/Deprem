import { fetchSince, connectLive, inRegion } from './feeds.js';
import {
  baseStyle, addDataLayers, OVERLAYS, AGE_STOPS, DEPTH_STOPS, FAULT_TYPES, SLIP_STOPS, MOMENT_STOPS,
  quakeColorExpr, faultColorExpr, ageColor, depthColor, faultTypeKey,
} from './layers.js';

const HOUR = 3600e3, DAY = 24 * HOUR;
const POLL_MS = 60e3;            // WebSocket'e ek olarak düzenli yoklama
const PULSE_MAX_AGE_H = 3;       // bu kadar yeni depremler sürekli "nabız" atar
const PULSE_MAX = 25;

const state = {
  events: new Map(),   // id -> olay
  loadedSince: Infinity,
  windowDays: 7,
  minMag: 0,
  colorMode: 'age',
  faultMode: 'type',
  live: 'connecting',
  source: 'EMSC',
  visible: Object.fromEntries(OVERLAYS.map((o) => [o.id, !!o.visible])),
  ready: false,
};

const $ = (s) => document.querySelector(s);
const fmtTime = new Intl.DateTimeFormat('tr-TR', { timeZone: 'Europe/Istanbul', day: '2-digit', month: 'short', hour: '2-digit', minute: '2-digit', second: '2-digit' });
const fmtDay = new Intl.DateTimeFormat('tr-TR', { timeZone: 'Europe/Istanbul', day: 'numeric', month: 'short' });
const fmtHour = new Intl.DateTimeFormat('tr-TR', { timeZone: 'Europe/Istanbul', hour: '2-digit', minute: '2-digit' });
const esc = (s) => String(s ?? '').replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

function ago(ms) {
  const s = Math.max(0, (Date.now() - ms) / 1000);
  if (s < 60) return `${Math.floor(s)} sn`;
  if (s < 3600) return `${Math.floor(s / 60)} dk`;
  if (s < 86400) return `${Math.floor(s / 3600)} sa`;
  return `${Math.floor(s / 86400)} g`;
}
const magStr = (m) => (m == null ? '?' : m.toFixed(1));
const colorOf = (e) => (state.colorMode === 'depth' ? depthColor(e.depth) : ageColor((Date.now() - e.time) / HOUR));

// ---------------------------------------------------------------- Harita
const map = new maplibregl.Map({
  container: 'map',
  style: baseStyle(),
  center: [35.2, 39.0],
  zoom: window.innerWidth < 820 ? 4.3 : 5.4,
  minZoom: 3.5,
  maxZoom: 12,
  maxBounds: [[5, 22], [66, 55]],
  attributionControl: false,
  dragRotate: false,
  pitchWithRotate: false,
});
map.touchZoomRotate.disableRotation();
map.addControl(new maplibregl.NavigationControl({ showCompass: false }), 'bottom-right');

function visibleEvents() {
  const since = Date.now() - state.windowDays * DAY;
  return [...state.events.values()]
    .filter((e) => e.time >= since && e.mag >= state.minMag)
    .sort((a, b) => b.time - a.time);
}

function toGeoJSON(list) {
  const now = Date.now();
  return {
    type: 'FeatureCollection',
    features: list.map((e) => ({
      type: 'Feature',
      geometry: { type: 'Point', coordinates: [e.lon, e.lat] },
      properties: { id: e.id, mag: e.mag, depth: e.depth, ageH: (now - e.time) / HOUR },
    })),
  };
}

// ---------------------------------------------------------------- Çizim
let renderQueued = false;
function render() {
  if (renderQueued) return;
  renderQueued = true;
  requestAnimationFrame(() => {
    renderQueued = false;
    const list = visibleEvents();
    if (state.ready) map.getSource('quakes').setData(toGeoJSON(list));
    renderStats(list);
    renderList(list);
    renderTimeline(list);
    renderPulses(list);
  });
}

function renderStats(list) {
  const label = { 1: 'son 24 saat', 7: 'son 7 gün', 30: 'son 30 gün' }[state.windowDays];
  $('#st-count').textContent = list.length.toLocaleString('tr-TR');
  $('#st-count-label').textContent = label;
  const max = list.reduce((m, e) => (e.mag > (m?.mag ?? -1) ? e : m), null);
  $('#st-max').textContent = max ? 'M' + magStr(max.mag) : '–';
  $('#st-max').parentElement.title = max ? `${max.place} · ${fmtTime.format(max.time)}` : '';
  $('#st-last').textContent = list[0] ? ago(list[0].time) : '–';
  $('#ev-count').textContent = list.length > 200 ? '· ilk 200' : '';
}

let lastListKey = '';
const freshIds = new Set();
function renderList(list) {
  const top = list.slice(0, 200);
  const key = state.colorMode + top.map((e) => e.id + e.mag).join();
  const ol = $('#event-list');
  if (key === lastListKey) {
    // yalnızca "önce" metinlerini güncelle
    ol.querySelectorAll('.ev-ago').forEach((n) => { n.textContent = ago(+n.dataset.t); });
    return;
  }
  lastListKey = key;
  ol.innerHTML = top.map((e) => `
    <li data-id="${esc(e.id)}" class="${freshIds.has(e.id) ? 'fresh' : ''}">
      <span class="mag-chip" style="--c:${colorOf(e)}">${magStr(e.mag)}</span>
      <div style="min-width:0">
        <div class="ev-place">${esc(e.place || '—')}</div>
        <div class="ev-meta">${fmtTime.format(e.time)} · ${e.depth.toFixed(0)} km · ${esc(e.auth)}</div>
      </div>
      <span class="ev-ago" data-t="${e.time}">${ago(e.time)}</span>
    </li>`).join('');
  freshIds.clear();
}

$('#event-list').addEventListener('click', (ev) => {
  const li = ev.target.closest('li');
  if (li) focusEvent(li.dataset.id);
});

function renderTimeline(list) {
  const svg = $('#timeline');
  const W = svg.clientWidth || 600, H = svg.clientHeight || 70;
  const t1 = Date.now(), t0 = t1 - state.windowDays * DAY;
  const maxM = Math.max(5, ...list.map((e) => e.mag));
  const x = (t) => ((t - t0) / (t1 - t0)) * (W - 8) + 4;
  const y = (m) => H - 12 - (Math.max(0, m) / maxM) * (H - 18);
  const parts = [];
  // yatay kılavuzlar
  for (let m = 1; m <= maxM; m += maxM > 6 ? 2 : 1) {
    parts.push(`<line class="ax" x1="0" x2="${W}" y1="${y(m)}" y2="${y(m)}"/>${H > 50 ? `<text x="2" y="${y(m) - 2}">${m}</text>` : ""}`);
  }
  // zaman etiketleri
  const stepH = state.windowDays === 1 ? 3 : state.windowDays === 7 ? 24 : 24 * 5;
  const start = Math.ceil(t0 / (stepH * HOUR)) * stepH * HOUR;
  for (let t = start; t <= t1; t += stepH * HOUR) {
    const lbl = state.windowDays === 1 ? fmtHour.format(t) : fmtDay.format(t);
    parts.push(`<line class="ax" x1="${x(t)}" x2="${x(t)}" y1="0" y2="${H - 10}"/><text x="${x(t) + 3}" y="${H}">${lbl}</text>`);
  }
  // olaylar: küçükler önce, büyükler üstte
  const sorted = [...list].sort((a, b) => a.mag - b.mag);
  for (const e of sorted) {
    const r = Math.max(1.4, 0.55 * 1.55 ** e.mag);
    parts.push(`<circle data-id="${esc(e.id)}" cx="${x(e.time).toFixed(1)}" cy="${y(e.mag).toFixed(1)}" r="${r.toFixed(1)}" fill="${colorOf(e)}" fill-opacity="0.8"><title>M${magStr(e.mag)} · ${esc(e.place)} · ${fmtTime.format(e.time)}</title></circle>`);
  }
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.innerHTML = parts.join('');
  $('#tl-range').textContent = `${fmtDay.format(t0)} – şimdi · ${state.source}`;
}

$('#timeline').addEventListener('click', (ev) => {
  const c = ev.target.closest('circle');
  if (c) focusEvent(c.dataset.id);
});

// Yeni depremler için DOM tabanlı dalga animasyonu
const pulseMarkers = new Map();
function makePulse(e, once = false) {
  const el = document.createElement('div');
  el.className = 'pulse' + (once ? ' once' : '');
  const size = Math.max(26, 9 * 1.5 ** e.mag);
  el.style.setProperty('--s', `${size}px`);
  el.style.setProperty('--c', colorOf(e));
  el.style.setProperty('--d', `${2 + e.mag * 0.25}s`);
  el.innerHTML = '<i></i><i></i><i></i>';
  return new maplibregl.Marker({ element: el }).setLngLat([e.lon, e.lat]);
}

function renderPulses(list) {
  const show = state.visible.quakes && state.ready;
  const recent = show ? list.filter((e) => Date.now() - e.time < PULSE_MAX_AGE_H * HOUR)
    .sort((a, b) => b.mag - a.mag).slice(0, PULSE_MAX) : [];
  const want = new Set(recent.map((e) => e.id));
  for (const [id, m] of pulseMarkers) if (!want.has(id)) { m.remove(); pulseMarkers.delete(id); }
  for (const e of recent) {
    if (!pulseMarkers.has(e.id)) pulseMarkers.set(e.id, makePulse(e).addTo(map));
    else pulseMarkers.get(e.id).getElement().style.setProperty('--c', colorOf(e));
  }
}

// Haritada bir depreme odaklan ve ayrıntılarını göster
let popup;
function focusEvent(id) {
  const e = state.events.get(id);
  if (!e) return;
  map.flyTo({ center: [e.lon, e.lat], zoom: Math.max(map.getZoom(), 7.2), speed: 1.4 });
  showPopup(e);
  if (window.innerWidth < 820) $('#panel').classList.remove('open');
}

function showPopup(e) {
  popup?.remove();
  popup = new maplibregl.Popup({ offset: 12, maxWidth: '300px' })
    .setLngLat([e.lon, e.lat])
    .setHTML(`
      <div class="pop-mag"><b style="color:${colorOf(e)}">M${magStr(e.mag)}</b><span>${esc(e.magType)} · ${ago(e.time)} önce</span></div>
      <div class="pop-place">${esc(e.place)}</div>
      <dl class="pop-grid">
        <dt>Zaman</dt><dd>${fmtTime.format(e.time)} TSİ</dd>
        <dt>Derinlik</dt><dd>${e.depth.toFixed(1)} km</dd>
        <dt>Konum</dt><dd>${e.lat.toFixed(3)}°K, ${e.lon.toFixed(3)}°D</dd>
        <dt>Kaynak</dt><dd>${esc(e.auth)} via ${esc(e.src)}</dd>
      </dl>
      ${e.url ? `<a class="pop-link" href="${esc(e.url)}" target="_blank" rel="noopener">Ayrıntılar ↗</a>` : ''}`)
    .addTo(map);
}

// ---------------------------------------------------------------- Bildirim
function toast(e) {
  const box = $('#toasts');
  while (box.children.length >= 4) box.firstElementChild.remove();
  const el = document.createElement('div');
  el.className = 'toast glass';
  el.innerHTML = `<span class="mag-chip" style="--c:${colorOf(e)}">${magStr(e.mag)}</span>
    <div><small>Yeni deprem</small><p>${esc(e.place)}</p><span class="ev-meta">${fmtHour.format(e.time)} · ${e.depth.toFixed(0)} km</span></div>`;
  el.onclick = () => focusEvent(e.id);
  box.appendChild(el);
  setTimeout(() => el.remove(), 10000);
  if (state.visible.quakes) {
    const m = makePulse(e, true).addTo(map);
    setTimeout(() => m.remove(), 12000);
  }
}

// ---------------------------------------------------------------- Veri
function ingest(events, { announce = false } = {}) {
  let changed = false;
  for (const e of events) {
    if (!inRegion(e)) continue;
    const prev = state.events.get(e.id);
    if (prev && prev.mag === e.mag && prev.lat === e.lat && prev.lon === e.lon && prev.depth === e.depth) continue;
    state.events.set(e.id, e);
    changed = true;
    if (!prev && announce && Date.now() - e.time < 2 * HOUR) {
      freshIds.add(e.id);
      toast(e);
    }
  }
  if (changed) render();
}

async function loadWindow() {
  const since = Date.now() - state.windowDays * DAY;
  if (since >= state.loadedSince) return render();
  try {
    const { source, events, snapshot } = await fetchSince(since);
    state.source = source;
    state.loadedSince = since;
    ingest(events);
    render();
    if (snapshot) setLive('snapshot');
    else if (state.live === 'connecting' || state.live === 'offline') setLive('poll');
  } catch (err) {
    console.error(err);
    setLive('offline');
  }
}

async function poll() {
  try {
    const { source, events, snapshot } = await fetchSince(Date.now() - 3 * HOUR);
    if (snapshot) return;
    state.source = source;
    ingest(events, { announce: true });
    if (state.live !== 'live') setLive('poll');
  } catch (err) {
    console.warn('Yoklama başarısız', err);
    if (state.live !== 'live') setLive('offline');
  }
}

function setLive(s) {
  state.live = s;
  const el = $('#live');
  el.dataset.state = s;
  $('#live-text').textContent = { live: 'Canlı', poll: 'Güncel', offline: 'Bağlantı yok', connecting: 'Bağlanıyor', snapshot: 'Önizleme' }[s];
  el.title = {
    live: 'EMSC canlı akışına bağlı; yeni depremler anında görünür',
    poll: `Canlı akış yok; veriler her ${POLL_MS / 1000} sn'de yenileniyor`,
    offline: 'Veri kaynaklarına ulaşılamıyor; yeniden deneniyor',
    connecting: 'Bağlanıyor…',
    snapshot: 'Canlı kaynaklara ulaşılamadı; kaydedilmiş anlık görüntü gösteriliyor',
  }[s];
}

// ---------------------------------------------------------------- Katmanlar
async function loadOverlays() {
  const res = await fetch('data/faults.geojson');
  const faults = await res.json();
  for (const f of faults.features) f.properties.kind = faultTypeKey(f.properties.slip_type);
  map.getSource('faults').setData(faults);
  OVERLAYS.find((o) => o.id === 'faults').desc = `GEM GAF-DB · ${faults.features.length} segment`;

  // İleride eklenecek bilimsel katmanlar: dosya varsa otomatik etkinleşir
  for (const o of OVERLAYS.filter((x) => x.probe)) {
    try {
      const r = await fetch(o.probe);
      if (!r.ok) continue;
      const data = await r.json();
      // MapLibre vektör karo kodlaması çok büyük tamsayı benzeri değerlerde (ör. 1e18 N·m) bozuluyor;
      // bu değerleri metne çevirerek koru.
      for (const f of data.features || []) {
        for (const [k, v] of Object.entries(f.properties || {})) {
          if (typeof v === 'number' && Math.abs(v) >= 1e15) f.properties[k] = v.toExponential(3);
        }
      }
      o.setup(map, data);
      o.soon = false;
      if (data.metadata?.note) o.desc = data.metadata.note;
      o.count = data.features?.length;
      state.visible[o.id] = false;
      applyVisibility(o);
    } catch (_) { /* henüz yok */ }
  }
  renderLayerList();
}

function applyVisibility(o) {
  const v = state.visible[o.id];
  for (const id of o.mapLayers || []) if (map.getLayer(id)) map.setLayoutProperty(id, 'visibility', v ? 'visible' : 'none');
  if (o.dom === 'cities') cityMarkers.forEach((m) => (m.getElement().style.display = v ? '' : 'none'));
  if (o.dom === 'pulses') render();
}

function renderLayerList() {
  $('#layer-list').innerHTML = OVERLAYS.map((o) => `
    <label class="layer ${o.soon ? 'soon' : ''}">
      <input type="checkbox" data-id="${o.id}" ${state.visible[o.id] ? 'checked' : ''} ${o.soon ? 'disabled' : ''}>
      <span class="sw"></span>
      <span class="txt"><b>${esc(o.label)}</b><span>${esc(o.desc)}</span></span>
      ${o.soon ? '<span class="tag">yakında</span>' : ''}
    </label>`).join('');
}

$('#layer-list').addEventListener('change', (ev) => {
  const o = OVERLAYS.find((x) => x.id === ev.target.dataset.id);
  if (!o) return;
  state.visible[o.id] = ev.target.checked;
  applyVisibility(o);
  renderLegend();
});

const cityMarkers = [];
async function loadCities() {
  const d = await (await fetch('data/cities.geojson')).json();
  for (const f of d.features) {
    const p = f.properties;
    const el = document.createElement('div');
    el.className = 'city' + (p.pop > 1e6 ? ' big' : '') + (p.tr ? '' : ' foreign');
    el.textContent = p.name;
    el.dataset.minzoom = p.pop > 1e6 || !p.tr ? 4 : p.pop > 4e5 ? 5.6 : 6.6;
    cityMarkers.push(new maplibregl.Marker({ element: el, anchor: 'left', offset: [-2, 0] }).setLngLat(f.geometry.coordinates).addTo(map));
  }
  updateCities();
}
function updateCities() {
  const z = map.getZoom();
  for (const m of cityMarkers) {
    const el = m.getElement();
    el.style.visibility = z >= +el.dataset.minzoom ? 'visible' : 'hidden';
  }
}
map.on('zoomend', updateCities);

// ---------------------------------------------------------------- Lejant
function renderLegend() {
  const parts = [];
  if (state.visible.quakes) {
    const sizes = [2, 3, 4, 5, 6].map((m) => {
      const r = Math.min(40, 2 * Math.max(2, 0.9 * 1.6 ** m));
      return `<div><i style="width:${r}px;height:${r}px"></i>${m}</div>`;
    }).join('');
    const ramp = state.colorMode === 'depth'
      ? `<div class="ramp" style="background:linear-gradient(90deg,${DEPTH_STOPS.map((s) => s[1]).join(',')})"></div>
         <div class="ramp-lbl"><span>0</span><span>30</span><span>150 km+</span></div>`
      : `<div class="ramp" style="background:linear-gradient(90deg,${AGE_STOPS.map((s) => s.c).join(',')})"></div>
         <div class="ramp-lbl"><span>şimdi</span><span>24 sa</span><span>30 g</span></div>`;
    parts.push(`<div class="lg"><h3>Büyüklük</h3><div class="sizes">${sizes}</div></div>`);
    parts.push(`<div class="lg"><h3>${state.colorMode === 'depth' ? 'Derinlik' : 'Geçen süre'}</h3>${ramp}</div>`);
  }
  const mo = OVERLAYS.find((o) => o.id === 'moment');
  if (state.visible.moment && !mo.soon) {
    parts.push(`<div class="lg"><h3>Moment açığı · eşdeğer Mw</h3>
      <div class="ramp" style="background:linear-gradient(90deg,${MOMENT_STOPS.map((s) => s[1]).join(',')})"></div>
      <div class="ramp-lbl"><span>5</span><span>6.5</span><span>7.5+</span></div>
      <p class="lg-note">Tahmin değildir. 1900 öncesi depremler dahil değildir.</p></div>`);
  }
  if (state.visible.faults) {
    const body = state.faultMode === 'rate'
      ? `<div class="ramp" style="background:linear-gradient(90deg,${SLIP_STOPS.map((s) => s[1]).join(',')})"></div>
         <div class="ramp-lbl"><span>0</span><span>4</span><span>22+ mm/yıl</span></div>`
      : `<div class="keys">${FAULT_TYPES.map((t) => `<span><i style="background:${t.c}"></i>${t.label}</span>`).join('')}</div>`;
    parts.push(`<div class="lg"><h3>Faylar · ${state.faultMode === 'rate' ? 'kayma hızı' : 'mekanizma'}</h3>${body}</div>`);
  }
  $('#legend').innerHTML = parts.join('');
  $('#legend').style.display = parts.length ? '' : 'none';
}

// ---------------------------------------------------------------- Kontroller
function seg(sel, fn) {
  $(sel).addEventListener('click', (ev) => {
    const b = ev.target.closest('button');
    if (!b) return;
    $(sel).querySelectorAll('button').forEach((x) => x.classList.toggle('on', x === b));
    fn(b.dataset.v);
  });
}
seg('#seg-window', (v) => { state.windowDays = +v; loadWindow(); });
seg('#seg-color', (v) => {
  state.colorMode = v;
  for (const id of ['quakes', 'quakes-glow']) map.setPaintProperty(id, 'circle-color', quakeColorExpr(v));
  lastListKey = '';
  render(); renderLegend();
});
seg('#seg-fault', (v) => {
  state.faultMode = v;
  for (const id of ['faults', 'faults-glow']) map.setPaintProperty(id, 'line-color', faultColorExpr(v));
  renderLegend();
});
$('#minmag').addEventListener('input', (ev) => {
  state.minMag = +ev.target.value;
  $('#minmag-out').textContent = `M ${state.minMag.toFixed(1)}+`;
  render();
});
$('#panel-toggle').addEventListener('click', () => $('#panel').classList.toggle('open'));

// Harita üzerinde etkileşim
let hoverId = null;
map.on('mousemove', 'quakes', (ev) => {
  map.getCanvas().style.cursor = 'pointer';
  const id = ev.features[0].properties.id;
  if (hoverId !== id) {
    if (hoverId) map.setFeatureState({ source: 'quakes', id: hoverId }, { hover: false });
    hoverId = id;
    map.setFeatureState({ source: 'quakes', id }, { hover: true });
  }
});
map.on('mouseleave', 'quakes', () => {
  map.getCanvas().style.cursor = '';
  if (hoverId) map.setFeatureState({ source: 'quakes', id: hoverId }, { hover: false });
  hoverId = null;
});
map.on('click', 'quakes', (ev) => {
  const e = state.events.get(ev.features[0].properties.id);
  if (e) showPopup(e);
});

const TYPE_TR = Object.fromEntries(FAULT_TYPES.map((t) => [t.key, t.label]));
const SLIP_TR = { Dextral: 'Sağ yanal doğrultu atımlı', Sinistral: 'Sol yanal doğrultu atımlı', Normal: 'Normal', Reverse: 'Ters', Subduction_Thrust: 'Dalma-batma bindirmesi' };
const faultTip = new maplibregl.Popup({ closeButton: false, closeOnClick: false, offset: 8, className: 'fault-tip' });
map.on('mousemove', 'faults', (ev) => {
  if (map.queryRenderedFeatures(ev.point, { layers: ['quakes'] }).length) return faultTip.remove();
  const p = ev.features[0].properties;
  const st = (p.slip_type || '').split('-').map((s) => SLIP_TR[s] || s).join(' + ') || TYPE_TR[p.kind];
  const rate = p.slip_rate != null && p.slip_rate !== 'null'
    ? `${(+p.slip_rate).toFixed(1)} mm/yıl${p.slip_rate_min !== 'null' && p.slip_rate_min != null ? ` <span style="color:var(--faint)">(${(+p.slip_rate_min).toFixed(1)}–${(+p.slip_rate_max).toFixed(1)})</span>` : ''}`
    : 'bilinmiyor';
  faultTip.setLngLat(ev.lngLat).setHTML(`
    <div class="pop-place" style="margin:0 0 6px">${esc(p.name || 'Adsız fay segmenti')}</div>
    <dl class="pop-grid"><dt>Tür</dt><dd>${esc(st)}</dd><dt>Kayma hızı</dt><dd>${rate}</dd><dt>Katalog</dt><dd>${esc(p.catalog)}</dd></dl>`).addTo(map);
});
map.on('mouseleave', 'faults', () => faultTip.remove());

// Moment açığı katmanı (veri varsa)
const fmtNum = (x, d = 1) => (x == null || x === 'null' ? '–' : (+x).toLocaleString('tr-TR', { maximumFractionDigits: d }));
map.on('mousemove', 'moment', (ev) => {
  const p = ev.features[0].properties;
  const mw = p.mw_equiv ?? p.mw_equiv_since_ref;
  faultTip.setLngLat(ev.lngLat).setHTML(`
    <div class="pop-place" style="margin:0 0 6px">${esc(p.name || p.fault_id || 'Fay segmenti')}</div>
    <dl class="pop-grid">
      <dt>Birikmiş moment</dt><dd>Mw ${fmtNum(mw, 2)} eşdeğeri</dd>
      ${p.m0_since_ref ? `<dt>M₀</dt><dd>${esc(String(p.m0_since_ref).replace('e+', '×10^'))} N·m</dd>` : ''}
      <dt>Kayma hızı</dt><dd>${fmtNum(p.slip_rate_mm_yr)} mm/yıl</dd>
      <dt>Uzunluk</dt><dd>${fmtNum(p.length_km, 0)} km</dd>
      <dt>Son büyük deprem</dt><dd>${p.last_big_event && p.last_big_event !== 'null' ? esc(p.last_big_event) : 'katalogda yok'}</dd>
    </dl>
    <p class="lg-note">1900'den bu yana katalog içi hesap; tahmin değildir.</p>`).addTo(map);
});
map.on('mouseleave', 'moment', () => faultTip.remove());

// ---------------------------------------------------------------- Başlat
map.on('load', async () => {
  addDataLayers(map);
  state.ready = true;
  renderLayerList();
  renderLegend();
  loadOverlays().catch((e) => console.error('Fay verisi yüklenemedi', e));
  loadCities().catch((e) => console.error('Şehirler yüklenemedi', e));
  await loadWindow();

  connectLive({
    onEvent: (e) => ingest([e], { announce: true }),
    onState: (s) => {
      if (state.live === 'snapshot' && s !== 'live') return;
      setLive(s === 'live' ? 'live' : (state.loadedSince < Infinity ? 'poll' : 'offline'));
    },
  });
  setInterval(poll, POLL_MS);
  setInterval(render, 30e3); // yaş renkleri ve "önce" süreleri tazelensin
  document.addEventListener('visibilitychange', () => { if (!document.hidden) poll(); });
  window.addEventListener('resize', render);
});

// Konsoldan hata ayıklama için
window.deprem = { state, map, ingest, focusEvent };
