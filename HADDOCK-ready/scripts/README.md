# HADDOCK related scripts

## Setup scripts

* `setup-ana_scripts.csh` : when cloning this repo, run this script to correct the path in all analysis scripts (takes as argument the name of the complex directory)

## restraints-related scripts

* `create-true-interface-AIRs.csh` : generates an `ambig.tbl` restraint files from true interface contacts at 3.9A cutoff
* `generate-act-act.sh` : generates `ambig.tbl` from a list of active and passive residues
* `generate-restraint-bodies.csh` : detects and generate distance restraints to keep unconnected segments together
* `generate-restraint-ligands.csh` : detects and generate distance restraints to keep ligand in place
* `generate-restraint.csh` : combines bodies and ligand restraints into one file

## Analysis-related scripts

* `i-rmsd_unbounds_to_xray.csh` : calculate the iRMSD values of the unbound superimposed complex and its components
* `results-stats.csh` : extract stats from a HADDOCK run once analysis has been done
* `results-stats-topclust.csh` : extract stats from a HADDOCK run considering only the first acceptable cluster once analysis has been done
* `extract-stats-iRMSD.csh` : extract some stats from the stats file generated by the previous script

## HADDOCK running-related scripts

The following scripts will setup and launch HADDOCK runs. They will detected possible crashed runs and resubmit them (for a few scenarios).
Edit them to define the run name and the maximum number of runs submitted.

* `setup-run-grid.csh` : Setups and launches HADDOCK2.4 runs using the grid
* `setup-run-node.csh` : Setups and launches HADDOCK2.4 runs using a full node per run
* `setup-run-local.csh : Setups and launches HADDOCK2.4 runs starting HADDOCK locally and using a batch system 
* `setup-analysis.csh` : Setups and launches the analysis of a completed HADDOCK run

The following scripts will monitor and lauch launch HADDOCK runs. They call the setup scripts above.
Edit them to define the run name and the maximum number of runs submitted.

* `monitor-submit-grid-jobs.csh`
* `monitor-submit-node-jobs.csh`
* `monitor-submit-local-jobs.csh`
* `monitor-submit-analysis.csh`
