"""Sıcaklık ve kayaç davranışı.

Jeoterm (Turcotte & Schubert 2014, Bölüm 4): yüzey ısı akısı q, ısıl
iletkenlik k, derinlikle üstel azalan radyojenik ısı A(z) = A0 exp(-z/h):

    T(z) = T0 + (q - A0 h) z / k + A0 h^2 / k (1 - exp(-z/h))

Sürtünme rejimi, ıslak granit deneylerinden (Blanpied vd. 1995, JGR) basit
bir sıcaklık eşleştirmesiyle: 100 °C altı hıza göre güçlenen (yavaş kayan),
100-350 °C hıza göre zayıflayan (kilitlenip deprem üreten), 350-450 °C geçiş,
üstü sünek akma. Sünek kuşağın akma (Maxwell) zamanı Arrhenius yasasıyla
sıcaklığa bağlanır: t_M(T) = t_ref exp(Q/R (1/T - 1/T_ref)).

Isı akısı her 1 km için deprem derinliklerinden ters çözülür: kırılgan
tabanı (depremlerin %90'ının üstünde kaldığı derinlik, D90) 350 °C eş
sıcaklık eğrisine denk getirilir. Çıkan ısı akısı bağımsız ölçümlerle
(Marmara için literatürde ~50-100 mW/m², Pfister vd. 1998; Erkan 2015)
karşılaştırılır.
"""
from __future__ import annotations

import numpy as np

K_ILETKEN = 2.5      # W/m/K
A0 = 1.5e-6          # W/m^3
H_RAD = 10e3         # m
T0 = 15.0            # °C
T_KIRILGAN = 350.0   # °C, kırılgan-sünek geçiş (kuvars)
Q_AKTIF = 135e3      # J/mol, kuvars dislokasyon sünmesi (Hirth vd. 2001 mertebesi)
R_GAZ = 8.314


def sicaklik(z_km, q_wm2):
    z = np.asarray(z_km, float) * 1e3
    q = np.asarray(q_wm2, float)
    return T0 + (q - A0 * H_RAD) * z / K_ILETKEN + A0 * H_RAD ** 2 / K_ILETKEN * (1 - np.exp(-z / H_RAD))


def isi_akisi_derinlikten(d_km, T_hedef=T_KIRILGAN):
    """T(d) = T_hedef olacak yüzey ısı akısı (W/m²)."""
    z = np.asarray(d_km, float) * 1e3
    rad = A0 * H_RAD ** 2 / K_ILETKEN * (1 - np.exp(-z / H_RAD))
    return (T_hedef - T0 - rad) * K_ILETKEN / z + A0 * H_RAD


def a_eksi_b(T):
    """Hız-durum sürtünmesi (a-b), sıcaklığa göre parçalı doğrusal."""
    T = np.asarray(T, float)
    return np.interp(T, [0, 80, 120, 330, 370, 450, 2000],
                     [0.012, 0.008, -0.004, -0.004, 0.0, 0.02, 0.02])


def maxwell_zamani(T, t_ref_yil, T_ref=500.0):
    Tk, Tr = np.asarray(T, float) + 273.15, T_ref + 273.15
    return t_ref_yil * np.exp(Q_AKTIF / R_GAZ * (1 / Tk - 1 / Tr))
