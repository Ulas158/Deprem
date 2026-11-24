# Kaynaklar, lisanslar ve dosya özetleri

Ham dosyalar `data/raw/` altında, indirildikleri hâliyle (değiştirilmeden) saklanır.
İndirme tarihi: 2026-11-24.

| Dosya | Kaynak | Kapsam | Lisans | SHA-256 |
|---|---|---|---|---|
| `EPICA_v1.1.xlsx` | https://emidius.eu/epica/data/EPICA_v1.1.xlsx | Avrupa-Akdeniz, 1000-1899, Mw (makrosismik + katalog) | CC BY 4.0 (dosyanın `format` sayfasında) | `c51189997ff557eee8a1b50e2eb40ad4bee6255844f4f11538bfec46c11c0ba1` |
| `SHARE_CET.xls` | https://emidius.eu/SHEEC/docs/SHARE_CET.xls | Orta ve Doğu Türkiye (32-45°D), 1003-2006 | Dosyada ve sayfada lisans belirtilmemiş; **doğrulanmalı** | `08b1d799ddee9428240f5f305f39e7f9022529fd6109ba52a0c0958842a0f40b` |
| `GEM-GHEC-v1.txt` | https://emidius.eu/GEH/download/GEM-GHEC-v1.txt | Dünya, 1000-1903, M≥~7 | CC BY-SA 3.0 (paylaşımda aynı lisans şartı) | `22a6346cdc408c876d502c98776ae0029e75a7585f00d3b2dbd156f30101e67e` |
| `ncei_tr_pre1900_p1.json`, `_p2.json` | NOAA NCEI Hazel API, `earthquakes?minLatitude=34&maxLatitude=44&minLongitude=24&maxLongitude=46&maxYear=1899` | MÖ 2000 - 1899, önemli (hasarlı) depremler | ABD kamu malı | p1: `1e84d7c2cb9bd7b6c1623faf7fe61477ac82eec406e5a854458bb70bd913411d`<br>p2: `591306166cf2f4d1eff2a824ae77f7032308ce403cfd6b1675d008f93d9f270e` |

## Atıflar

- **EPICA**: Rovida A., Antonucci A. (2021). EPICA - European PreInstrumental earthquake CAtalogue, version 1.1. INGV. https://doi.org/10.13127/epica.1.1
  ve Rovida A., Antonucci A., Locati M. (2022). The European Preinstrumental Earthquake Catalogue EPICA, the 1000-1899 catalogue for the European Seismic Hazard Model 2020. *Earth Syst. Sci. Data* 14, 5213-5231. https://doi.org/10.5194/essd-14-5213-2022
- **SHEEC / SHARE**: Stucchi M. ve ark. (2013). The SHARE European Earthquake Catalogue (SHEEC) 1000-1899. *J. Seismology* 17, 523-544. https://doi.org/10.1007/s10950-012-9335-2
- **GEM GHEC**: Albini P., Musson R.M.W., Rovida A., Locati M., Gomez Capera A.A., Viganò D. (2014). The Global Earthquake History. *Earthquake Spectra* 30(2), 607-624. https://doi.org/10.1193/122013EQS297
  ve Albini P. ve ark. (2013). Global Historical Earthquake Archive and Catalogue (1000-1903). GEM Technical Report 2013-01.
- **NOAA NCEI**: National Geophysical Data Center / World Data Service: NCEI/WDS Global Significant Earthquake Database. NOAA NCEI. https://doi.org/10.7289/V5TD9V7K
- **Ms → Mw**: Scordilis E.M. (2006). Empirical global relations converting MS and mb to moment magnitude. *J. Seismology* 10, 225-236.
- **Tamlık**: Stepp J.C. (1972). Analysis of completeness of the earthquake sample in the Puget Sound area. NOAA Tech. Rep. ERL 267-ESL 30.

GHEC kayıtlarının her biri kendi birincil kaynağını (`reference` sütunu) taşır: Ambraseys & Jackson (1998),
Ambraseys (2009), Ambraseys & Melville (1982), Guidoboni & Comastri (2005), Shebalin & Tatevossian (1997),
Sbeinati ve ark. (2005) vb. EPICA kayıtlarında `reference` makrosismik veri ve bölgesel katalog kodlarıdır
(ör. `CPTI15`, `AMBR009`).

## Erişilemeyen / eklenmeyen kaynaklar

- **GEM GHEC v2015**: GEM sitesinde form arkasında; doğrudan indirme adresi bulunamadı. Bunun yerine v1.0 (2014'te NC maddesi kaldırılmış hâli) kullanıldı.
- **Ambraseys (2009)** *Earthquakes in the Mediterranean and Middle East* ve **Ambraseys & Finkel (1995)** *The Seismicity of Turkey and Adjacent Areas*: kitaplar; makine-okunur sürümleri yok. İçerikleri GHEC, SHARE-CET ve EPICA üzerinden dolaylı olarak katalogda.
- **Guidoboni & Comastri (2005)** ve **KOERI tarihsel kataloğu**: açık makine-okunur dosya bulunamadı.
