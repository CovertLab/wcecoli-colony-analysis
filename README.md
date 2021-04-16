# Analysis Code for Whole-Cell *E. coli* Colony Simulations

This repository contains code and data to accompany a paper
"Whole-Colony Modeling of *Escherichia coli*" that is currently in
preparation. This project is hosted on GitHub at
[https://github.com/CovertLab/wcecoli-colony-analysis](https://github.com/CovertLab/wcecoli-colony-analysis)
and licensed under the terms in [`LICENSE.txt`](LICENSE.txt).

We welcome comments, which you can leave by filing an issue on the
GitHub repository.

## Requirements

### Hardware Requirements

* 5 GB of free space because the un-compressed simulation data is about
  2.4 GB.
* At least 4 GB RAM because you will need to load at least one
  experiment's data into memory to generate the figures. This data can
  be over 1 GB. If you generate all the figures at once, the analysis
  code will need to load all 2.4 GB of simulation data into memory.
  Depending on what other processes consume RAM on your system, you may
  need more memory and/or swap space.

### Software Requirements

* If you are running on a headless system, you will need to install the
  X virtual framebuffer (XVFB) and prepend the figure generation
  commands with `xvfb-run -a`.
* To analyze phylogenetic data, you will need to have
  [R](https://www.r-project.org/) (version 4.0.4) with the `argparse`
  and `phytools` packages. Developers will also need the `lintr` package
  for linting.
* [Python](https://python.org) 3.8.3
* If you clone from GitHub, you'll also need these tools:
  * [Git](https://git-scm.com)
* [GPG](https://gnupg.org) to verify code integrity if you don't want to
  rely on GitHub's security checking.

## Setup

1. Clone the repository.

   ```console
   $ git clone https://github.com/CovertLab/wcecoli-colony-analysis.git
   ```

   Alternatively, you can simply extract an archive of the source code
   if you got the code in that format.

   **You should always validate the code's integrity and make sure it
   comes from a trusted source.** See the [Security section](#security)
   below for instructions.
2. (recommended but optional) Setup a Python virtual environment.
3. Install Python dependencies

   ```console
   $ pip install numpy
   $ pip install -r requirements.txt
   ```

   Note that you have to install numpy first because the `setup.py`
   script of one of our dependencies requires it.
4. Install R dependencies. In an R shell, run:

   ```r
   install.packages('argparse')
   install.packages('phytools')
   ```

   If you want to run the lint checks (only developers need to do this),
   also run:

   ```r
   install.packages('lintr')
   ```

## Reproduce Analyses

There are two ways to make the simulation data accessible to the
analysis code: reading from files or reading from a MongoDB database.
Here, we will only describe how to read from files since MongoDB is
over-kill for just reproducing our analyses.

1. The raw simulation data is not included in the repository due to its
   size. Instead, it is available under the DOI
   [10.5281/zenodo.4697519](https://doi.org/10.5281/zenodo.4697519).
   Follow the DOI link and download all the files that it contains to
   `data/`.

   **You should verify the integrity of the downloaded files**. You can
   do so like this:

   ```console
   $ cd data
   $ gpg --verify SHA512SUMS.txt.asc
   gpg: assuming signed data in 'SHA512SUMS.txt'
   gpg: Signature made Fri Apr 16 14:26:12 2021 EDT
   gpg:                using RSA key F76925D5D12B91104587678FC98CBB9C501917E0
   gpg: Good signature from "anonymous <cs.temporary@icloud.com>" [full]
   $ shasum -c SHA512SUMS.txt
   LICENSE.txt: OK
   README.md: OK
   archived_simulations.tar.gz: OK
   search.json: OK
   ```

2. The raw simulation data is now an archive at
   `data/archived_simulations.tar.gz`. You can extract it like this:

   ```console
   $ cd data
   $ tar -xf archived_simulations.tar.gz
   ```

   Now you should see the simulation data as a series of JSON files
   under `data/archived_simulations`.

   We also provide the `data/search.json` file, which specifies the
   theoretical boundary that appears in magenta in Figures 5G 5H, and
   5I. This is already un-compressed, so you don't need to do anything
   more with it except pass it to scripts as specified below.

3. Now you can generate figures using the `src/make_figures.py` script.
   Run `python -m src.make_figures -h` to see the available options:

   ```console
   $ python -m src.make_figures -h
   usage: make_figures.py [-h] [--atlas] [--port PORT]
                          [--host HOST]
                          [--database_name DATABASE_NAME]
                          [--data_path DATA_PATH] [--3A] [--3B]
                          [--3C] [--3D] [--3E] [--3F] [--3G]
                          [--5A] [--5B] [--5C] [--5D] [--5E]
                          [--5F] [--5G] [--5H] [--5I] [--X1]
                          [--all]
                          search_data

   Generate selected figures and associated stats from
   simulation data.

   positional arguments:
     search_data           Path to boundary search data.

   optional arguments:
     -h, --help            show this help message and exit
     --atlas, -a           Read data from an mongoDB Atlas
                           instead of a local mongoDB.
                           Credentials, cluster subdomain, and
                           database name should be specified in
                           secrets.json.
     --port PORT, -p PORT  Port at which to access local mongoDB
                           instance. Defaults to "27017".
     --host HOST, -o HOST  Host at which to access local mongoDB
                           instance. Defaults to "localhost".
     --database_name DATABASE_NAME, -d DATABASE_NAME
                           Name of database on local mongoDB
                           instance to read from. Defaults to
                           "simulations".
     --data_path DATA_PATH
                           Folder of JSON files to read data
                           from instead of Mongo
     --3A                  Generate figure & stats for fig 3A:
                           snapshots of growing colony consuming
                           glucose
     --3B                  Generate figure & stats for fig 3B:
                           environment cross-sections showing
                           glucose depletion
     --3C                  Generate figure & stats for fig 3C:
                           snapshots in the basal condition
     --3D                  Generate figure & stats for fig 3D:
                           snapshots in the anaerobic condition
     --3E                  Generate figure & stats for fig 3E:
                           colony mass on basal and anaerobic
                           media
     --3F                  Generate figure & stats for fig 3F:
                           snapshots showing expression
                           heterogeneity
     --3G                  Generate figure & stats for fig 3G:
                           distributions of protein
                           concentrations
     --5A                  Generate figure & stats for fig 5A:
                           parameter scan for tolerance
                           threshold
     --5B                  Generate figure & stats for fig 5B:
                           snapshot of final colony under
                           nitrocefin
     --5C                  Generate figure & stats for fig 5C:
                           box plot showing distances from
                           center
     --5D                  Generate figure & stats for fig 5D:
                           phylogenetic tree
     --5E                  Generate figure & stats for fig 5E:
                           dotplot of final [AmpC] colored by
                           survival
     --5F                  Generate figure & stats for fig 5F:
                           dotplot of final [AcrAB-TolC] colored
                           by survival
     --5G                  Generate figure & stats for fig 5G:
                           final [AmpC] and plotted against
                           [AcrAB-TolC]
     --5H                  Generate figure & stats for fig 5H:
                           paths of dead agents through
                           concentration space
     --5I                  Generate figure & stats for fig 5I:
                           paths of a lineage of agents through
                           concentration space
     --X1                  Generate figure & stats for fig X1:
                           death_snapshots_antibiotic
     --all                 Generate all figures and stats.
     ```

     For example, to generate Figure 3A from the paper:

     ```console
     $ python -m src.make_figures data/search.json --data_path data/archived_simulations --3A
     ```

     This will save the plots used in Figure 3A to `out/figs/`. Note
     that in this case there will be multiple plots generated, one for
     each simulation, even though only one of these plots was actually
     used for Figure 3A in the paper. Along with the figures, a
     `stats.json` file will also be created.

     You can also generate figures and associated statistics for all the
     figures, even those that don't appear in the paper (e.g. `X1`). You
     can do this like so:

     ```console
     $ python -m src.make_figures data/search.json --data_path data/archived_simulations --all
     ```

4. Calculate summary statistics using `src/analyze_stats.py`. You can
   view its help text by running:

   ```console
   $ python -m src.analyze_stats -h
   usage: analyze_stats.py [-h] [-o OUT] stats_json

   positional arguments:
     stats_json         Path to stats JSON file

   optional arguments:
     -h, --help         show this help message and exit
     -o OUT, --out OUT  Path to write summary stats to
   ```

   For example, you could generate summary statistics from the figure 3A
   statistics by running:

   ```console
   $ python -m src.analyze_stats out/figs/stats.json -o out/figs/summary_stats.json
   ```

   The out/figs/summary_stats.json file stores the summary statistics in
   a human-readable format.

5. Analyze phylogeny data using `src/analyze_phylogeny.r` like this:

   ```console
   $ Rscript src/analyze_phylogeny.r out/figs/phylogeny.nw out/figs/agent_survival.csv
   ```

   The analysis will be printed to the console.

## For Developers

### Create Simulation Data Archives

You can use the `src/archive_experiments.py` script to create an archive
of all the simulation data used by `src/make_figures.py`:

```console
$ python -m src.archive_experiments <connection args>
```

where `<connection args>` can include `-o IP` for MongoDB IP address
`IP` and `-p PORT` for MongoDB port `PORT`.

### Run Tests

We use `mypy` for type checking, `pylint` for linting, and `pytest` for
unit tests. You can run all these by executing the `test.sh` script.
There should be no errors.

## Security

**This code is not hardened against malicious inputs.** Therefore, you
should only run it on data you trust is not malformed. This trusted data
might include the archived simulation data we provide or the outputs of
simulations you ran.

All releases and commits are signed with the OpenPGP key
`0x12C1B5FF558317E5`, which has the fingerprint:

```text
1642C5C9C092F5AC1FE9222012C1B5FF558317E5
```

### Verifying Release Signatures

You can validate releases by verifying tag signatures:

```console
$ git verify-tag <tag name>
```

If you retrieved an archive of the code with an accompanying signature,
you can also verify that signature:

```console
$ gpg --verify <signature file>
```

### Verifying Commit Signatures

You can use `git log --show-signature` to verify these signatures,
though note that the signatures will show a fingerprint of the signing
sub-key:

```text
F76925D5D12B91104587678FC98CBB9C501917E0
```

We have provided a script to automatically do this check in
`check_signatures.py`, though of course you should make sure you
understand the script before you trust it.

### Verifying Using GitHub

If you trust GitHub, you may also rely on the `Verified` labels GitHub
adds to signed commits and releases. Click on the label to check which
key was used.
