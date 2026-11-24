#!/usr/bin/env python3
"""web/data/layers/*.layer.json tanımlarını tarayıp web/data/layers/index.json dizinini yazar.

GitHub Pages yayın iş akışı bunu her yayında çalıştırır; yerelde denemek için:
  python3 web/tools/build_layer_index.py
Bozuk JSON içeren tanımlar atlanır ve uyarı verilir (yayını durdurmaz).
"""
import json
import sys
from pathlib import Path

LAYERS = Path(__file__).resolve().parent.parent / "data" / "layers"


def main() -> int:
    LAYERS.mkdir(parents=True, exist_ok=True)
    entries = []
    for path in sorted(LAYERS.glob("*.layer.json")):
        try:
            d = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as e:
            print(f"UYARI: {path.name} okunamadı, atlandı: {e}", file=sys.stderr)
            continue
        data = d.get("data")
        if data and not data.startswith(("http://", "https://", "/")):
            target = (LAYERS.parent.parent / data) if data.startswith("data/") else (LAYERS / data)
            if not target.exists():
                print(f"UYARI: {path.name} → veri dosyası yok: {data}", file=sys.stderr)
        entries.append(path.name)
    (LAYERS / "index.json").write_text(json.dumps({"layers": entries}, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(f"{len(entries)} katman tanımı → {LAYERS / 'index.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
