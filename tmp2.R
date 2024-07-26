# 1. Fixed CS3 processing: reference and expression samples used are now all from Vancouver
# 2. Fixed xsite normalization to use pools
# 3. Can I use xsite replicates in CS3_R (random1) as part of CS3_USC_X and CS3_AOC_X? No for now

library(DiagrammeR)
grViz(
  paste0(
    "digraph processing {

    # layout
    graph [layout = dot overlap = TRUE]

    # node definitions
    node [shape = box]
    a [label = 'Selected Cohorts from all CodeSets\nN= 3674']
    b [label = 'N=3441']
    c1 [label = 'CS1\nn=286']
    c2 [label = 'CS2\nn=882']
    c3 [label = 'CS3\nn=2273']
    d1x [label = 'CS1_X\nn=279']
    d1r [label = 'CS1_R\nn=7']
    d2x [label = 'CS2_X\nn=876']
    d2r [label = 'CS2_R\nn=6']
    d3x [label = 'CS3_X\nn=2264']
    d3r [label = 'CS3_R\nn=9']
    e3van [label = 'CS3_VAN_X\nn=2178']
    e3usc [label = 'CS3_USC_X\nn=45']
    e3aoc [label = 'CS3_AOC_X\nn=41']

    f1x [label = 'CS1_X\nn=263']
    f1r [label = 'CS1_R\nn=5']
    f2x [label = 'CS2_X\nn=827']
    f2r [label = 'CS2_R\nn=5']
    f3r [label = 'CS3_R\nn=5']

    f3van [label = 'CS3_VAN_X\nn=2094']
    f3usc [label = 'CS3_USC_X\nn=28']
    f3usc_pools [label = 'CS3_USC_R\nn=17']
    f3aoc [label = 'CS3_AOC_X\nn=22']
    f3aoc_pools [label = 'CS3_AOC_R\nn=19']

    g1x [label = 'CS1_T\nn=263']
    g2x [label = 'CS2_T\nn=827']

    g3usc [label = 'CS3_USC_T\nn=28']
    g3aoc [label = 'CS3_AOC_T\nn=22']

    h3 [label = 'CS3_T\nn=520']
    train_all [label = 'Normalized Expression from all CodeSets\nN=1610']
    train_hist [label = 'Major Histotypes\nN=1535']
    train_uniq [label = 'Final Training\nN=1243']

    # edge definitions
    a -> b [label = '  Quality Control']
    b -> {c1 c2}
    b -> c3 [label = 'Split by CodeSets']
    c1 -> {d1x d1r}
    c2 -> {d2x d2r}
    c3 -> d3x  [label = 'Split into Expression and Reference Samples']
    c3 -> d3r
    d1x -> f1x
    d1r -> f1r
    d2x -> f2x
    d2r -> f2r
    d3x -> {e3van e3usc}
    d3x -> e3aoc [label = '  Split by site']
    e3van -> f3van
    e3usc -> {f3usc f3usc_pools}
    e3aoc -> f3aoc
    e3aoc -> f3aoc_pools [label = 'Remove duplicates and separate pools']
    d3r -> f3r
    f1r -> g1x [xlabel = 'Normalization by Random1                 ']
    {f1x f3r} -> g1x
    {f2x f2r f3r} -> g2x
    {f3usc f3usc_pools} -> g3usc
    f3aoc -> g3aoc
    f3aoc_pools -> g3aoc [label = 'Normalization by Pools']
    {f3van g3usc} -> h3
    g3aoc -> h3 [label = 'Remove test sets and pools']
    {g1x g2x} -> train_all
    h3 -> train_all [label = 'Combine CodeSets']
    train_all -> train_hist [label = '  Remove other Histotypes']
    train_hist -> train_uniq [label = ' Remove cross CodeSet and cross site replicates']

    c1 -> c2 -> c3 [style = invis]

    # subgraph definitions
    subgraph {rank = same; c1; c2; c3}
    subgraph {rank = same; d1x; d1r; d2x; d2r; d3r; d3x}
    subgraph {rank = same; e3van; e3usc; e3aoc}
    subgraph {rank = same; f1x; f1r; f2x; f2r; f3r; f3usc; f3usc_pools; f3aoc; f3aoc_pools}
    subgraph {rank = same; f3van; g3usc; g3aoc}
    subgraph {rank = same; g1x; g2x; h3}
    }"
  )
)
