from symfit import parameters, variables, sin, cos, Fit
from symfit.core.minimizers import BFGS, BasinHopping
import numpy as np
import matplotlib.pyplot as plt
import json

def fourier_series(x, f, n=0):
    """
    Returns a symbolic fourier series of order `n`.

    :param n: Order of the fourier series.
    :param x: Independent variable
    :param f: Frequency of the fourier series
    """
    # Make the parameter objects for all the terms
    a0, *cos_a = parameters(','.join(['a{}'.format(i) for i in range(0, n + 1)]))
    sin_b = parameters(','.join(['b{}'.format(i) for i in range(1, n + 1)]))
    # Construct the series
    series = a0 + sum(ai * cos(i * f * x) + bi * sin(i * f * x)
                     for i, (ai, bi) in enumerate(zip(cos_a, sin_b), start=1))
    return series

def fourier_series_cosine(x, f, n=0):
  """
  Returns a symbolic cosine Fourier series of order `n`.

  :param n: Order of the cosine Fourier series.
  :param x: Independent variable
  :param f: Frequency of the cosine Fourier series
  """
  # Make the parameter objects for all cosine terms
  a_n = parameters(','.join(['a{}'.format(i) for i in range(n + 1)]))

  # Construct the cosine series
  series = a_n[0] + sum(a_n[i] * cos(i * f * x) for i in range(1, n + 1))  # Access a_n elements directly
  return series


def read_dat_file(filename, conversion_factor=1.0):
  """
  Reads a .dat file with two space-separated floating point values per line 
  and returns them as separate NumPy arrays.

  Args:
      filename: The path to the .dat file.

  Returns:
      A tuple containing two NumPy arrays: one for the first column and 
      another for the second column.
  """
  with open(filename, 'r') as f:
    # Use np.fromstring to efficiently read lines and convert to floats
    data = np.fromstring(f.read(), sep=' ', dtype=np.float64)
  # Reshape data into two separate arrays
  return data.reshape(-1, 2)[:, 0], data.reshape(-1, 2)[:, 1]*np.float64(conversion_factor)

def fit_fourier_series(data_file, n_order, minimizer, plot=True, output_file=None):
  """
  Fits a Fourier series (cosine by default) to data from a .dat file.

  Args:
      data_file: Path to the .dat file.
      n_order: Order of the Fourier series.
      minimizers: List of minimizer objects from symfit.core.minimizers (default: [BFGS]).
      plot: Boolean flag to enable plotting (default: True).
      output_file: Path to a JSON file to save fit results (default: None).
  """
  xdata, ydata = read_dat_file(data_file)
  x, y = variables('x, y')
  w, = parameters('w')
  model_dict = {y: fourier_series_cosine(x, f=w, n=n_order)}
  print (model_dict)

  # Define a Fit object and fit the model
  fit = Fit(model_dict, x=xdata, y=ydata, minimizer=list(minimizer))
  fit_result = fit.execute(BFGS={'tol': 1e-5})
  print (fit_result)

  # Save fit results to JSON (if specified)
  if output_file:
    with open(output_file, 'w') as f:
      json.dump(fit_result.params, f, indent=4)  # Save only parameters

  # Plotting (optional)
  if plot:
      plt.plot(xdata, ydata)
      plt.plot(xdata, fit.model(x=xdata, **fit_result.params).y, ls=':')
      plt.xlabel('x')
      plt.ylabel('y')
      plt.show()


# Example Usage
data_file = "d_10_5_1_2.dat"
n_order = 5
minimizers = [BFGS, BasinHopping, BFGS, BasinHopping]  # List of minimizers from scipy 
output_file = "fit_results.json"  # Optional JSON output file (uncomment to use)
fit_fourier_series(data_file, n_order, minimizer=minimizers, plot=True, output_file=output_file)
