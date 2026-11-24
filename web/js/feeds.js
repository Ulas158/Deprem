// Deprem veri akışları: EMSC FDSN (birincil), USGS FDSN (yedek), EMSC WebSocket (canlı).
// Tüm kaynaklar ortak bir olay biçimine dönüştürülür:
//   { id, time (ms), lat, lon, depth (km), mag, magType, place, auth, src, url }

export const REGION = { minlat: 34, maxlat: 43.5, minlon: 24, maxlon: 46 };

const EMSC_FDSN = 'https://www.seismicportal.eu/fdsnws/event/1/query';
const USGS_FDSN = 'https://earthquake.usgs.gov/fdsnws/event/1/query';
const EMSC_WS = 'wss://www.seismicportal.eu/standing_order/websocket';

export function inRegion(e) {
  return e.lat >= REGION.minlat && e.lat <= REGION.maxlat && e.lon >= REGION.minlon && e.lon <= REGION.maxlon;
}

const titleCase = (s) => (s || '').toLocaleLowerCase('tr-TR').replace(/(^|[\s,\-(])(\p{L})/gu, (m, a, b) => a + b.toLocaleUpperCase('tr-TR'));

function fromEmsc(f) {
  const p = f.properties || {};
  const id = p.unid || f.id;
  return {
    id: 'emsc:' + id,
    time: Date.parse(p.time),
    lat: +p.lat, lon: +p.lon,
    depth: Math.abs(+p.depth || 0),
    mag: p.mag == null ? null : +p.mag,
    magType: (p.magtype || '').toUpperCase(),
    place: titleCase(p.flynn_region),
    auth: p.auth || 'EMSC',
    src: 'EMSC',
    url: `https://www.seismicportal.eu/eventdetails.html?unid=${encodeURIComponent(id)}`,
  };
}

function fromUsgs(f) {
  const p = f.properties || {};
  const [lon, lat, depth] = f.geometry.coordinates;
  return {
    id: 'usgs:' + f.id,
    time: p.time,
    lat, lon,
    depth: Math.abs(depth || 0),
    mag: p.mag == null ? null : +p.mag,
    magType: (p.magType || '').toUpperCase(),
    place: p.place || '',
    auth: (p.net || 'us').toUpperCase(),
    src: 'USGS',
    url: p.url,
  };
}

const iso = (ms) => new Date(ms).toISOString().slice(0, 19);

async function getJson(url, timeoutMs = 20000) {
  const ctl = new AbortController();
  const t = setTimeout(() => ctl.abort(), timeoutMs);
  try {
    const r = await fetch(url, { signal: ctl.signal });
    if (r.status === 204) return { features: [] }; // FDSN: sonuç yok
    if (!r.ok) throw new Error(`${r.status} ${url}`);
    return await r.json();
  } finally {
    clearTimeout(t);
  }
}

/** Belirtilen zamandan bu yana bölgedeki depremleri çeker. EMSC başarısız olursa USGS'e düşer. */
export async function fetchSince(sinceMs) {
  const q = new URLSearchParams({
    starttime: iso(sinceMs), orderby: 'time', limit: '10000',
    minlatitude: REGION.minlat, maxlatitude: REGION.maxlat,
    minlongitude: REGION.minlon, maxlongitude: REGION.maxlon,
  });
  try {
    const d = await getJson(`${EMSC_FDSN}?format=json&${q}`);
    return { source: 'EMSC', events: d.features.map(fromEmsc).filter(valid) };
  } catch (err) {
    console.warn('EMSC alınamadı, USGS deneniyor:', err);
    const d = await getJson(`${USGS_FDSN}?format=geojson&${q}`);
    return { source: 'USGS', events: d.features.map(fromUsgs).filter(valid) };
  }
}

const valid = (e) => Number.isFinite(e.time) && Number.isFinite(e.lat) && Number.isFinite(e.lon) && e.mag != null;

/**
 * EMSC canlı WebSocket akışına bağlanır, kopunca artan gecikmeyle yeniden bağlanır.
 * onEvent(event, action) — action: 'create' | 'update'
 * onState('live' | 'down')
 */
export function connectLive({ onEvent, onState }) {
  let ws, retry = 0, closed = false, timer;

  const open = () => {
    try {
      ws = new WebSocket(EMSC_WS);
    } catch (e) {
      schedule();
      return;
    }
    ws.onopen = () => { retry = 0; onState('live'); };
    ws.onmessage = (msg) => {
      try {
        const m = JSON.parse(msg.data);
        if (!m.data || !m.data.properties) return;
        const e = fromEmsc(m.data);
        if (valid(e) && inRegion(e)) onEvent(e, m.action || 'create');
      } catch (err) {
        console.warn('WebSocket mesajı çözülemedi', err);
      }
    };
    ws.onclose = () => { onState('down'); schedule(); };
    ws.onerror = () => { try { ws.close(); } catch (_) { /* yok say */ } };
  };

  const schedule = () => {
    if (closed) return;
    clearTimeout(timer);
    const delay = Math.min(60000, 2000 * 2 ** retry++);
    timer = setTimeout(open, delay);
  };

  open();
  return () => { closed = true; clearTimeout(timer); ws && ws.close(); };
}
