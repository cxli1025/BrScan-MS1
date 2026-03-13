## BrScan-MS1: an automated scanning algorithm for Brominated Internal Standards (BrIS) identification

BrScan-MS1 is a Python-based tool designed for the automated identification and quantification of Brominated Internal Standards (BrIS) using MS1 spectra. 

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
