# Symmetry

Symmetries for each complex of the Benchmark5.0 (receptor and ligand) has been calculated using [AnAnaS v.1.1](https://team.inria.fr/nano-d/software/ananas/).

For each complex:

```bash
ananas 1K4C_r_u.pdb -y --json 1k4c_r.json
ananas 1K4C_l_u.pdb -y --json 1k4c_l.json
```

JSON output files have been parsed and summarized on [symmetries.csv](symmetries.csv)
