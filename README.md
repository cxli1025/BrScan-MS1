<div align="left">

## BrScan-MS1: an automated scanning algorithm for Brominated Internal Standards (BrIS) identification

</div>

<p align="left">
<a href="https://www.python.org/"><img alt="python" src="https://img.shields.io/badge/Python-3.12.11-yellow.svg"/></a>
<a href="https://pypi.org/project/numpy/"><img alt="numpy" src="https://img.shields.io/badge/Numpy-2.3.2-teal.svg"/></a>
<a href="https://pypi.org/project/pandas/"><img alt="pandas" src="https://img.shields.io/badge/Pandas-2.3.2-brightgreen.svg"/></a>
<a href="https://pypi.org/project/pyteomics/"><img alt="pyteomics" src="https://img.shields.io/badge/Pyteomics-4.7.5-purple.svg"/></a>
<a href="https://pypi.org/project/lxml/"><img alt="lxml" src="https://img.shields.io/badge/lxml-6.0.1-purple.svg"/></a>
<a href="https://github.com/michalsta/opentims"><img alt="opentims" src="https://img.shields.io/badge/Opentims-1.0.16-red.svg"/></a>
<a href="https://pypi.org/project/brain-isotopic-distribution/"><img alt="brainpy" src="https://img.shields.io/badge/Brainpy-1.5.19-blue.svg"/></a>
</p>

## Algorithm workflow
BrScan-MS1 is a Python-based tool designed for the automated identification and quantification of Brominated Internal Standards (BrIS) using MS1 spectra. 
![workflow](./workflow.tiff)



```shell
git clone git@github.com:KCLabMTU/LM-OGlcNAc-Site.git
```
### 1. Prerequisites

### Python Environment

Ensure you have Python 3.8+ installed. You can install the required Python libraries using pip:

```shell
pip install numpy pandas pyteomics opentims_bruker_bridge opentimspy
```

### External Tools (For .raw files)
To process Thermo .raw files directly, **BrScan-MS1** requires the ThermoRawFileParser. Download [ThermoRawFileParser](https://github.com/CompOmics/ThermoRawFileParser). Unzip the compressed file into the `BrScan-MS1` folder and rename the folder to `third_party`.

**Structure:** `BrScan-MS1/third_party/ThermoRawFileParser.exe`

### 2. Example Data
To run **BrScan-MS1**, you can download an example data named `wgl_dda_30min_20240705_10_Slot1-36_1_27757.d` from [PXD000000](https://github.com/CompOmics/ThermoRawFileParser), and extract the .d folder into the `BrScan-MS1` folder.

### 3. Usage
Run **BrScan-MS1** via the following command line. Use the `--profile` flag if your data is in profile mode (**BrScan-MS1** will automatically perform centroiding).

```shell
python BrScan_MS1.py -f ./wgl_dda_30min_20240705_10_Slot1-36_1_27757.d -l ./BrIS_library.csv -o ./results
```

**Output Structure:**

Results are saved in the `./results/[wgl_dda_30min_20240705_10_Slot1-36_1_27757]/` directory, containing 6 key files:

`BrIS_isotope[0-2]_result.csv`: Detailed identification results for three isotopes (m/z, intensity, charge, etc.).

`BrIS_isotope0_xic.csv`: Extracted Ion Chromatogram (XIC) data for the first isotope across all scans.

`BrIS_all_isotopes_xic.csv`: Comparative XIC data for all isotopic peaks.

`BrIS_isotopes_ratios.csv`: Experimental vs. theoretical intensity ratios for isotopes.

### 4. Help Message

Run the following command line to print parameter descriptions of **BrScan-MS1**:

```shell
python BrScan_MS1.py --help
```
