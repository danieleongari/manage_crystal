# manage_crystal
A tool to convert crystal files (atoms coordinates + unit cell) into common files and extract some useful info

## Installation:
```
pip install -e .[pre-commit,testing]
```
## Usage:

- to get default info about the crystal:

```
$ manage_crystal.py inputfilename.inputformat [options]`
```

- to convert to another format:

```
$ manage_crystal.py inputfilename.inputformat [options] -o outputfilename.outputformat
```

or

```
$ manage_crystal.py inputfilename.inputformat [options] -o outputformat
```

- to get help and explore the functionalities:

```
$ manage_crystal.py --help
```

Tips:

- you may want to use `-silent` to suppress default verbose output: several options "skip -silent" so that you can print just that information on the screen (e.g., `-printatoms -silent` prints on the screen just the atom types on one line). This make easy to use bash loops for lists of structures.

## Development

To enable automatic code formatting for every commit, do

```
pip install pre-commit yapf prospector pylint
pre-commit install
```
