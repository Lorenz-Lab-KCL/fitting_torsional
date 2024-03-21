# Fitting torsional angle

In this repository, we can find 2 files. The **angle_fitting.py**, where the program fitting is stored. The file
**d_10_5_2_1.dat**, where the angles vs relaxed energies obtained by ORCA has been stored. This is for the molecule
*poly_t_cispt*.

## Remarks

The script runs with the python library symfit, this needs to be installed by typing:

```
pip install symfit
```

Similarly, there are two functions for fitting the data. The first one is **fourier_series**, that contains sine terms and the second
one is **fourier_series_cosine** that only contains the cosine terms. This last one is the default. One can change it by changing the
following part of the code (line 71):

```
 model_dict = {y: fourier_series_cosine(x, f=w, n=n_order)}
````

The code prints the results and outputs a plot showing the fitting and the raw data superimposed. The parameters are also stored in 
a .json file called in this case **fit_results.json**. This can be opened by any text editor or by typing: 

```
head -n 90 fit_results.json
```

The fitting procedure can also be monitored by looking at 4 different variables, these are the **standard deviation** for each 
coefficient (this needs to be as low as possible), the **chi square** parameter (also as low as possible), the **objective value** (as low
as possible), and the **r_squared** parameter that should be very close to 1. 
