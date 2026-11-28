"""Marmara (Kuzey Anadolu Fayı kuzey kolu) için fizik tabanlı fay kırılma modeli.

Modüller:
    geometri   fay izi, 1 km x 1 km hücreler, elastik etkileşim (Green) matrisi
    termal     ısı akısı -> sıcaklık -> sürtünme rejimi ve akışkanlık (viskozite)
    simulator  yükleme, yavaş kayma (krip), sünek akma, deprem kırılması
    tarih      tarihsel kırılmalar ve sırayla oynatma (veri asimilasyonu)
    dogrulama  her büyük depremden önce modelin gösterdiği yerin sınanması
"""
