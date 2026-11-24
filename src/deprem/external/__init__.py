"""Gaz, sıvı ve dış etkiler (gelgit, mevsim, uzak tetikleme, insan kaynaklı).

İçe aktarıldığında fay elemanı kaydına ve gizli parametre listesine yeni
alanları ekler; hücre × pencere özellik tablosuna yeni özellik grupları
kaydeder. Ayrıntı: ``docs/akiskan_ve_dis_etkiler.md``.

Buradaki hiçbir değişken doğrulanmış bir öncü değildir; hepsi geçmiş
veride referans modellere karşı sınanacak hipotezlerdir.
"""

from . import params  # noqa: F401  (kayıtları ekler)
from . import features  # noqa: F401
