# Wiener Linien

JSON File enthält 3 Keys: EA, E und SR.

E: Haltestellen (proper UTF8 pending) als Liste

EA: Distance matrix von Haltestellen ie. müssten ungefähr kilometer sein. Habe keine geo library genutzt (siehe ipynb). Index entspricht position in E.

SR: Set relations. Liste von Listen, Index entspricht position in E.

Beispiel: E[2] ist "Burggasse-Stadthalle". EA[2,:] hat km Distanzen von Burggasse-Stadthalle zu allen andern Haltestellen, zB EA[2,0] Burggasse-Stadthalle -- Alaudagasse. Kilometer, nicht Stationen mit Ubahn. SR[2] hat die Linien bei Burggasse-Stadthalle, also U6.
