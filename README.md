# Paperboy

A command line program to check for the most recently published articles from specified scientific journals.

## Features

What Paperboy does:
* Pulls all articles, news, editorials, etc. published since the last check (or within the last 4 days if used the first time)
* For the set of scientific journals specified from all journals curated in the NCBI Pubmed/Medline database.

Note:
**No publisher subscription** is required, but only the identifying information is printed, not the abstract or full text.

## Get started

Download/clone the code. Install required Python packages (see Requirements below).

#### Find your journal
Let say you want to add the journal *PLOS Computational Biology*. Find its NlmID. (The following command greps for the journal's Medline abbreviation, but you can grep for anything you want.)
```
$ python paperboy-cli.py journal-list-all|grep -e 'PLoS Comput Biol'
INFO: Found 33599 journals.
Journal: PLoS computational biology (PLoS Comput Biol) (NlmId: 101238922) (ISSN: 1553-734X, 1553-7358)
```

#### Mark your journal active
Add the NlmID to the active journal list.
```
$ python paperboy-cli.py journal-add 101238922
INFO: PLoS Comput Biol (NlmId: 101238922) was marked active.
```

#### (Optionally) Verify your journal was set active
```
$ python paperboy-cli.py journal-list-active
INFO: Found 1 active journals.
Journal: PLoS computational biology (PLoS Comput Biol) (NlmId: 101238922) (ISSN: 1553-734X, 1553-7358)
```

#### Check for new articles
Any new articles since the last check, or since 4 days ago if this is the first time, will be downloaded and shown.
```
$ python paperboy-cli.py 
INFO: Loading new articles since 2021/01/03.
INFO: Found 27 new articles from PLoS Comput Biol.
INFO: Displaying 26 articles, excluding 1 Letter, Editorial, News, Interview.
Probabilistic transmission models incorporating sequencing data for healthcare-associated Clostridioides difficile outperform.. Eyre, (4 more). PLoS Comput Biol (https://doi.org/10.1371/journal.pcbi.1008417) (PMID 33444378) (Journal Article)
Investigating the mitochondrial genomic landscape of Arabidopsis thaliana by long-read sequencing. Masutani, Arimura, Morishita. PLoS Comput Biol (https://doi.org/10.1371/journal.pcbi.1008597) (PMID 33434206) (Journal Article)
...
```
Once your set of journals is active, just repeat the last step.

## Requirements

The following python packages are required: `biopython`, `lxml`, `colorama`, `pyparsing`, `appdirs`.

One way to set up an environment is by using conda.
```
conda create -n paperboy python=3.8 biopython lxml colorama pyparsing appdirs
conda activate paperboy
```

## Usage

The usage help is:
```
paperboy-cli.py [-h] [-d] {update,update-show,show,show-last-update,set-last-update,journal-list-all,journal-list-active,journal-add,journal-remove}
```

The commands do the following:
* `update`: Check for new published articles and download their IDs since the last update.
* `show`: Show details about all new articles, including title, authors, and DOI link.
* `update-show`: Run `update` followed by `show`. This is the **default** command.
* `show-last-update`: Print the date of the last update.
* `set-last-update`: Set the date of the last update if you want to check older articles again.
* `journal-list-all`: List all journals available in the Pubmed/Medline database. This is a fairly large number (~30k at the time of writing) and can take a bit. This list is downloaded the first time and stored locally, and only updated periodically after that.
* `journal-list-active`: List all journals for which new articles are pulled.
* `journal-add`: Add a journal to the list of active journals; requires the journal's NlmID (listed by `journal-list-all`).
* `journal-remove`: Remove a journal from the list of active journals; requires the journal's NlmID.

Output is formatted to terminal width, unless the terminal width is small, or the output is redirected.

## Files

* The code is contained in just one file, `paperboy-cli.py`.
* One config file is stored in an OS-specific user-dependent location, e.g. on Linux in `~/.config/paperboy/paperboy.conf`.

## Contribute

**Any help is welcome**, bug reports, feature requests, or pull requests for code or documentation!
