# Fitting torsional angle

In this repository, we can find 2 files. The **torsional_fitting.py**, where the program fitting is stored. The file
**xtbscan.log** contains the results from a **xtb mol --opt tight --input scan.inp** calculations performed using the 
molecule *hexane* stored in **mol.xyz**. This article presents the development of the process:

https://pubs.acs.org/doi/10.1021/ct200908r


## Remarks

The script runs with the python library symfit, this needs to be installed by typing:

```
pip install symfit
```

## How to use the code

The code can be used just by typping:

```
python torsional_fitting.py
```

The code has 4 different capabilities:

1. **XTB Template:** This capability creates a file **scan.inp** that can be subsequently used by **xtb** to produce the torsional angle as by typping this command:

```
xtb mol.xyz --opt --input scan.inp
```

where **mol.xyz** is the molecule to be reparametrised. The user **must** provide the atom numbers and the initial dihedral angle. This should be done using 
**AVOGADRO** or other visualisation program.


2. **XTB data collect:** This module reads the file **xtbscan.log** file produced after a successful calculation of a torsional angle. The information is stored in the
file **scan_data.dat**. The user must provide the same values as in the previous case for the variables scan_points, scan_start, scan_end, like:

```
xtb_data_collector.save_xtb_scan_data("xtbscan.log", "scan_data.dat", scan_points, scan_start, scan_end)
```
3. **Coordinate processor:** This module creates a new folder with all coordinates/frames created by **xtb** or **ORCA**. The many files are stored in a folder called **sp_frames**.
They can be used for subsequent calculation using **Gromacs** in a recursive manner. This module just needs the **xtbscan.lo** file. 

4. **FourierFitTool:** This module perform a Fourier series procedure:

$`x=y`$


```
if __name__ == "__main__":           
  #Example usage XTB_Template
  #template_manager = XTB_Template()
  #template_string = template_manager.create_template(
  #    force_constant=0.15,
  #    dihedral=[2, 3, 4, 5, 62.9],
  #    scan_start=0.0,
  #    scan_end=360.0,
  #    scan_points=100,
  #)
  #template_manager.save_template_to_file("scan.inp", template_string)

  #Example usage XTB data collect
  #xtb_data_collector = XTB_data_collect()
  #xtb_data_collector.save_xtb_scan_data("xtbscan.log", "scan_data.dat", 100, 0.0, 360.0)

  #Postprocessing xtbscan.log for creating many frames for subsequent Gromacs calculation.
  #processor = CoordinateProcessor("xtbscan.log")
  #processor.format_and_write_data()

  # Example usage for fitting the data AFTER xtb - gromacs energy substraction.
  fit = FourierFitTool()
  minimizers = [BFGS, BasinHopping, BFGS, BasinHopping, BFGS, BFGS]
  # Using full Fourier series
  #fit.fit_fourier_series(data_file="scan_data.dat", n_order=5, minimizer=minimizers, fourier_type="regular", plot=True)  
  # Or using cosine Fourier series
  fit.fit_fourier_series(data_file="scan_data.dat", n_order=3, minimizer=minimizers, fourier_type='cosine', output_file="fit_results.json")
```


The code prints the results and outputs a plot showing the fitting and the raw data superimposed. The parameters are also stored in 
a .json file called in this case **fit_results.json**. This can be opened by any text editor or by typing: 

```
head -n 90 fit_results.json
```

The fitting procedure can also be monitored by looking at 4 different variables, these are the **standard deviation** for each 
coefficient (this needs to be as low as possible), the **chi square** parameter (also as low as possible), the **objective value** (as low
as possible), and the **r_squared** parameter that should be very close to 1. 
