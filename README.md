## Requirements
1. Must have Python installed (version >= 3.6): https://www.python.org/

## Installation
1. Download the latest release of the package: https://github.com/cfe-lab/HAPLOID/releases
2. Decompress the archive and open a terminal window within
3. Create a virtual environment from your terminal:

```bash
python -m venv venv
```

4. Source the virtual environment:

```bash
# Windows
source ./venv/Scripts/activate
# Linux
source ./venv/bin/activate
```

5. Install the package (ensure you are inside the extracted folder, if you `ls` you should see the `README.md` file):

```bash
python -m pip install .
```

6. Run the software (with -h for help, specify other flags as needed):

```bash
haploid -h
```