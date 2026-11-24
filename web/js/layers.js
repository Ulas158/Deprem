// Harita katmanları ve katman kayıt defteri.
// Temel katmanlar (depremler, faylar, şehirler, ızgara). Bilimsel veri katmanları için
// web/data/layers/ altına bir *.layer.json tanımı eklemek yeterlidir (bkz. custom-layers.js, README).

// ---------- Renk ölçekleri ----------
export const AGE_STOPS = [
  { h: 1, c: '#ff3b30', label: '1 sa' },
  { h: 6, c: '#ff6a3d', label: '6 sa' },
  { h: 24, c: '#ff9f43', label: '24 sa' },
  { h: 72, c: '#ffd25a', label: '3 g' },
  { h: 168, c: '#a9b8cc', label: '7 g' },
  { h: Infinity, c: '#5f7290', label: '30 g' },
];

export const DEPTH_STOPS = [
  [0, '#fff07a'], [10, '#ffb347'], [30, '#ff6a3d'], [70, '#e0497a'], [150, '#8a5cf0'], [300, '#3f7be0'],
];

export const FAULT_TYPES = [
  { key: 'Dextral', label: 'Sağ yanal', c: '#ff8c42' },
  { key: 'Sinistral', label: 'Sol yanal', c: '#ffd25a' },
  { key: 'Normal', label: 'Normal', c: '#4fc3f7' },
  { key: 'Reverse', label: 'Ters / bindirme', c: '#ff4f8b' },
  { key: 'Oblique', label: 'Verev (karma)', c: '#b18cff' },
];

export const SLIP_STOPS = [[0, '#33415a'], [1, '#6f5bd0'], [4, '#e0497a'], [10, '#ff8c42'], [22, '#ffe680']];

export function ageColor(hours) {
  return AGE_STOPS.find((s) => hours < s.h).c;
}

export function depthColor(d) {
  let c = DEPTH_STOPS[0][1];
  for (const [v, col] of DEPTH_STOPS) if (d >= v) c = col;
  return c;
}

export const faultTypeKey = (t) => {
  if (!t) return 'Oblique';
  if (t.includes('-')) return 'Oblique';
  if (t === 'Subduction_Thrust' || t === 'Reverse') return 'Reverse';
  return t;
};

// ---------- Paint ifadeleri ----------
const magRadius = ['interpolate', ['exponential', 1.6], ['get', 'mag'], 0, 2, 2, 3.2, 3, 4.8, 4, 8, 5, 13, 6, 21, 7, 34];

export const quakeRadius = ['interpolate', ['linear'], ['zoom'], 4, magRadius, 9, ['*', 1.8, magRadius]];
const quakeGlowRadius = ['interpolate', ['linear'], ['zoom'], 4, ['*', 2.2, magRadius], 9, ['*', 4, magRadius]];

export function quakeColorExpr(mode) {
  if (mode === 'depth') {
    return ['interpolate', ['linear'], ['get', 'depth'], ...DEPTH_STOPS.flat()];
  }
  const e = ['step', ['get', 'ageH'], AGE_STOPS[0].c];
  AGE_STOPS.slice(0, -1).forEach((s, i) => e.push(s.h, AGE_STOPS[i + 1].c));
  return e;
}

export function faultColorExpr(mode) {
  if (mode === 'rate') {
    return ['interpolate', ['linear'], ['coalesce', ['get', 'slip_rate'], 0], ...SLIP_STOPS.flat()];
  }
  const m = ['match', ['get', 'kind']];
  FAULT_TYPES.forEach((t) => m.push(t.key, t.c));
  m.push('#8899aa');
  return m;
}

const faultWidth = ['interpolate', ['linear'], ['coalesce', ['get', 'slip_rate'], 0.3],
  0, 0.7, 1, 1.1, 5, 1.8, 15, 2.8, 30, 3.6];

// ---------- Taban harita stili (karo sunucusu gerektirmez) ----------
export function baseStyle() {
  return {
    version: 8,
    sources: {
      countries: { type: 'geojson', data: 'data/countries.geojson', tolerance: 0.3 },
    },
    layers: [
      { id: 'sea', type: 'background', paint: { 'background-color': '#070c14' } },
      {
        id: 'land', type: 'fill', source: 'countries',
        paint: { 'fill-color': ['case', ['==', ['get', 'iso'], 'TUR'], '#121c2a', '#0d141f'] },
      },
      {
        id: 'coast-glow', type: 'line', source: 'countries',
        paint: { 'line-color': '#2b5a8a', 'line-width': 6, 'line-blur': 6, 'line-opacity': 0.12 },
      },
      {
        id: 'borders', type: 'line', source: 'countries',
        paint: {
          'line-color': ['case', ['==', ['get', 'iso'], 'TUR'], '#36506f', '#223249'],
          'line-width': ['interpolate', ['linear'], ['zoom'], 4, 0.6, 9, 1.4],
        },
      },
    ],
  };
}

// Graticule (enlem/boylam ızgarası) – teknik bir görünüm katar.
export function graticule() {
  const f = [];
  for (let lon = 10; lon <= 60; lon += 2) f.push({ type: 'Feature', geometry: { type: 'LineString', coordinates: [[lon, 25], [lon, 50]] }, properties: {} });
  for (let lat = 26; lat <= 50; lat += 2) f.push({ type: 'Feature', geometry: { type: 'LineString', coordinates: [[10, lat], [60, lat]] }, properties: {} });
  return { type: 'FeatureCollection', features: f };
}

// ---------- Katman kayıt defteri ----------
// id: benzersiz ad · label/desc: panel metni · visible: başlangıç durumu
// mapLayers: aç/kapat ile görünürlüğü değişen MapLibre katman kimlikleri
// soon: henüz veri yoksa pasif gösterilir · group: paneldeki başlık
// Bilimsel veri katmanları burada değil, web/data/layers/*.layer.json tanımlarından gelir (custom-layers.js).
export const OVERLAYS = [
  {
    id: 'quakes', label: 'Depremler', desc: 'EMSC canlı akışı', visible: true, group: 'Temel katmanlar',
    mapLayers: ['quakes-glow', 'quakes'], dom: 'pulses',
  },
  {
    id: 'faults', label: 'Aktif faylar', desc: 'GEM küresel aktif fay veritabanı', visible: true,
    mapLayers: ['faults-glow', 'faults'],
  },
  {
    id: 'cities', label: 'Şehirler', desc: 'Nüfusu 150 binin üzerindeki iller', visible: true, dom: 'cities',
  },
  {
    id: 'grid', label: 'Koordinat ızgarası', desc: '2° aralıklı', visible: true, mapLayers: ['graticule'],
  },
];

export function addDataLayers(map) {
  map.addSource('graticule', { type: 'geojson', data: graticule() });
  map.addLayer({ id: 'graticule', type: 'line', source: 'graticule', paint: { 'line-color': '#1a2738', 'line-width': 0.6, 'line-dasharray': [2, 4] } }, 'land');

  map.addSource('faults', { type: 'geojson', data: { type: 'FeatureCollection', features: [] } });
  map.addLayer({
    id: 'faults-glow', type: 'line', source: 'faults', layout: { 'line-cap': 'round', 'line-join': 'round' },
    paint: { 'line-color': faultColorExpr('type'), 'line-width': ['*', 4, faultWidth], 'line-blur': 5, 'line-opacity': 0.22 },
  });
  map.addLayer({
    id: 'faults', type: 'line', source: 'faults', layout: { 'line-cap': 'round', 'line-join': 'round' },
    paint: {
      'line-color': faultColorExpr('type'),
      'line-width': ['interpolate', ['linear'], ['zoom'], 4, faultWidth, 10, ['*', 2.2, faultWidth]],
      'line-opacity': 0.85,
    },
  });

  map.addSource('quakes', { type: 'geojson', data: { type: 'FeatureCollection', features: [] }, promoteId: 'id' });
  map.addLayer({
    id: 'quakes-glow', type: 'circle', source: 'quakes',
    paint: {
      'circle-radius': quakeGlowRadius,
      'circle-color': quakeColorExpr('age'),
      'circle-blur': 1,
      'circle-opacity': ['interpolate', ['linear'], ['get', 'ageH'], 0, 0.55, 24, 0.3, 168, 0.12, 720, 0.06],
    },
  });
  map.addLayer({
    id: 'quakes', type: 'circle', source: 'quakes',
    layout: { 'circle-sort-key': ['get', 'mag'] },
    paint: {
      'circle-radius': quakeRadius,
      'circle-color': quakeColorExpr('age'),
      'circle-opacity': ['interpolate', ['linear'], ['get', 'ageH'], 0, 0.95, 24, 0.85, 168, 0.6, 720, 0.4],
      'circle-stroke-color': ['case', ['boolean', ['feature-state', 'hover'], false], '#ffffff', 'rgba(255,255,255,0.55)'],
      'circle-stroke-width': ['case', ['boolean', ['feature-state', 'hover'], false], 2, ['>=', ['get', 'mag'], 4], 1.2, 0.5],
    },
  });
}
