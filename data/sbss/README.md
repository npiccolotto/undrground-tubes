# SBSS Dataset

Die Bilder (gross = hoher Wert, rot = negativ, blau = positiv) habe ich in ownCloud gegeben weil 294 davon [1]. Brauchst du aber nur für visuellen Check. Die für unsere Vis relevanten Daten sind in dem JSON File `ensemble_set.json`. Key "E" sind Namen für die Elemente und Dateinamen von den Bildern. Key "EA" ist die Distanzmatrix, Key "SR" enthält Set Informationen, speziell das "bin" Array ist ein Label fürs Set. Also z.B. wenn Index 3 (0-basiert) in "bin" "2" enthält, heisst es dass das vierte Element (E[3]) in Set "2" ist. Das kannst du dir also raussuchen was du brauchst, z.B. alle 294 Karten und Bin 1 vom Feature "moransi" und Bin 4 vom Feature "Kurtosis" als Sets.

[1] https://owncloud.tuwien.ac.at/index.php/s/QYG3WNrnGb2fOVq
