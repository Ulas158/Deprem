"""Türkiye ve çevresi için aletsel dönem öncesi (1900 öncesi) tarihsel deprem kataloğu.

Ham kaynak dosyaları ``tarihsel/data/raw`` altındadır. ``build`` modülü bunları
ortak şemaya çevirir, kataloglar arası yinelenen kayıtları gruplar ve her kayda
belirsizliğe dayalı eğitim ağırlığı yazar.
"""

from .config import Config

__all__ = ["Config"]
