"""Bütün yerel pencerelerin ağını ve Coulomb önbelleğini (eğitim ve sınav çizelgesi) önceden kur."""
import sys
from multiprocessing import get_context
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))


def kur(ad):
    from deprem.mekanik import ozumseme as oz
    import mekanik_ikinci_tur as M
    oz.yerel_ag_yukle(ROOT, ad)
    for sinav in (False, True):
        _, ol, say = M._hazirla("kucuk", ad, sinav)
    print(ad, "ayar", sum(o.puan == "ayar" for o in ol), "sinav", sum(o.puan == "sinav" for o in ol),
          "coulomb", say.get("coulomb_eklenen"), flush=True)


if __name__ == "__main__":
    from deprem.mekanik import ozumseme as oz
    with get_context("fork").Pool(4) as p:
        list(p.imap_unordered(kur, list(oz.YEREL_AGLAR)))
