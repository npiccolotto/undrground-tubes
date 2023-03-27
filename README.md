# draw-ensemble-set-vis

TODO

- √ hex glyphs
- √ ordered bundles ("Edge routing with ordered bundles", Pupyrev et al., 2016)
- configurable spacing and sizing
- smooth curves instead of straight lines
- square glyphs
- kelpfusion drawing ("KelpFusion: A Hybrid Set Visualization Technique", Meulemans et al., 2013)
- edge bundling, e.g., ("Winding Roads: Routing edges into bundles", Lambert et al., 2010)


# Optimize MST choice

Observation: MSTs may not be unique, e.g., a 3-clique has 3 MSTs. The choice of MST
affects the drawing. We could want to minimize the number of edges used and the number of crossings.

Could this be an ILP? Thoughts based on hex diagram.

There are 78 edges in a glyph: Fully connected 6 side nodes, 6 corner nodes, 1 center node = 13*12/2 = 78 edges per glyph.
Each triangular face of the lattice contains 9 edges in the glyph margins: 3 neighbor edges, 3 to the anchor, 3 from the anchor to anchors of adjacent faces.

That seems like a lot, but couple more observations:

* Margin edges can only cross if they are in the same face and of different type (anchor vs neighbor)
* Vertex edges can only cross if they are in the same vertex
* Vertex edges connected to the glyph center cannot cross
* Vertex edges not connected to glyph center only cross if their (cw) endpoints are not contained in each other
  * Like, 2--5 crosses with 4--8, but not with 1--9, 3--4 or 2--5.
      123456789
    a  ----
    b    -----   <- crosses a
    c   --       <- crosses b
    d --------
* Margin and vertex edges cannot cross

If we can model those repeating structures, we can maybe save a lot of variables and constraints.

That's the graph. Then we have choices for all the MSTs. There are n sets, each using k connected components. We have m choices for each CC and l choices for the MST between components. `O(k*m * l)` choices per set, though value of variables differ per set. The choice for MST and each CC is independent. We just need to choose once per CC (k times in total) and once per connecting MST (1). Each choice entails using a set of edges = boolean vector (e1,e2,...,e). The bundle-edges objective is then to minimize non-zero entries (=sum) in the disjunction (OR) of all the choices of all the sets.

Crossings are even more difficult, naively need to compute them for all pairs of edges. However one could determine the set of potentially crossing edges quite quickly with a lookup for all edges (vertex # for vertex edges, negative face # for margin edges)?
Then we only need to look at those. But no idea how to pack this in an ILP.
