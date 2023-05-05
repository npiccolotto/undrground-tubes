# Wiener Linien

JSON File enthält 5 Keys: EA, E, S, SA und SR.

E: Haltestellen als Liste

EA: Distance matrix von Haltestellen ie. müssten ungefähr kilometer sein. Habe keine geo library genutzt (siehe ipynb). Index entspricht position in E.

SR: Set relations, i.e., welches Element (Haltestelle) gehört in welche Sets (Linien). Liste von Listen, Index entspricht position in E.

S: Liste von Sets (Linien).

SA: Distance Matrix von Haltestellen bezogen auf gemeinsame Linien.

Beispiel: E[2] ist "Burggasse-Stadthalle". EA[2,:] hat km Distanzen von Burggasse-Stadthalle zu allen andern Haltestellen, zB EA[2,0] Burggasse-Stadthalle -- Alaudagasse. Kilometer, nicht Stationen mit Ubahn. SA[2,:] das gleiche aber für gemeinsame Linien (Jaccard Distance). SR[2] hat die Linien bei Burggasse-Stadthalle, also U6. S enthält alle Sets (Linien) die es gibt.
