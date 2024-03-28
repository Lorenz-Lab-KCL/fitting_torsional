from jinja2 import Environment, Template
from symfit import parameters, variables, sin, cos, Fit
from symfit.core.minimizers import BFGS, BasinHopping
from itertools import count, groupby
import matplotlib.pyplot as plt
import numpy as np
import re, json, os


class XTB_Template():
  """
  This class manages the creation and saving of templates using Jinja2.
  """

  def __init__(self):
    self.env = Environment()

  def create_template(self, force_constant, dihedral, scan_start, scan_end, scan_points):
    """
    This method creates a template string using Jinja2 based on the provided variables.

    Args:
      force_constant: The force constant value (float).
      dihedral: A list of 5 elements (4 integers and 1 float).
      scan_start: The starting value for the scan (float).
      scan_end: The ending value for the scan (float).
      scan_points: The number of points for the scan (integer).

    Returns:
      A string containing the formatted template.
    """
    
    template = self.env.from_string("""$constrain
  force constant={{ force_constant }}
  dihedral: {% for item in dihedral %}{{ item }}{% if not loop.last %},{% endif %}{% endfor %}
$scan
  1: {{ scan_start }},{{ scan_end }},{{ scan_points }}
$end
""")
    return template.render(force_constant=force_constant, dihedral=dihedral, scan_start=scan_start, scan_end=scan_end, scan_points=scan_points)

  def save_template_to_file(self, filename, template_string):
    """
    This method saves a template string to a file.

    Args:
      filename: The name of the file to save the template to (string).
      template_string: The template string to save (string).
    """
    try:
      with open(filename, 'w') as f:
        f.write(template_string)
      print(f"Template saved successfully to file: {filename}")
    except FileNotFoundError:
      print(f"Error: Could not create file '{filename}'.")
    except Exception as e:
      print(f"An error occurred while saving the file: {e}")

class XTB_data_collect():
  """
  This class provides functions for collecting and saving XTB scan data.
  """

  def save_xtb_scan_data(self, filename, output_filename, scan_points, scan_start, scan_end):
    """
    Generates simulated XTB scan data and saves it to a text file.

    Args:
      filename: Path to the XTB scan log file (for reference, not used in this function).
      output_filename: Name of the text file to save the data. (provided by user)
      scan_points: Number of data points in the scan.
      scan_start: Starting value of the scan.
      scan_end: Ending value of the scan.
    """

    energies = self.read_xtb_file(filename)  # Assuming read_xtb_file is within the class

    # Round energies and calculate corresponding scan positions
    xdata = np.array([np.round(((scan_end - scan_start) / (scan_points - 1)) * i + scan_start, 2) for i in range(int(scan_points))])
    ydata = np.array([np.round(energy, 6) for energy in energies])

    # Combine xdata and ydata into a single array for easier saving
    data = np.stack((xdata, ydata), axis=1)

    # Save data to text file with user-provided name
    with open(output_filename, 'w') as f:
      # Save each data point on a separate line, separated by spaces
      np.savetxt(f, data, fmt='%.2f %.6f', delimiter=' ')

    print(f"XTBScan data saved to: {output_filename}")

  def read_energy(self, data):
    """
    This function reads an energy value from a string using regular expressions.

    Args:
      data: The string containing the energy value.

    Returns:
      The energy value as a float, or None if not found.
    """
    match = re.search(r"energy: ([-+]?\d+\.\d+)", data)
    if match:
      return float(match.group(1))
    else:
      return None

  def read_xtb_file(self, filename):
    """
    This function reads an energy value from each line in a file and returns them as a list.

    Args:
      filename: The path to the file to read.

    Returns:
      A list of energy values (floats), or an empty list if the file is not found or an error occurs.
    """
    energies = []
    try:
      with open(filename, 'r') as f:
        for line in f:
          energy = self.read_energy(line)
          if energy is not None:
            energies.append(energy)
    except FileNotFoundError:
      print(f"Error: File '{filename}' not found.")
    except Exception as e:
      print(f"An error occurred while reading the file: {e}")
    return energies

class FourierFitTool():
    """
    This class provides tools for fitting Fourier series to data.
    """

    def __init__(self):
        pass

    def read_dat_file(self, filename, conversion_factor=1.0):
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
        return data.reshape(-1, 2)[:, 0], data.reshape(-1, 2)[:, 1] * np.float64(conversion_factor)

    def fourier_series_cosine(self, x, f, n=0):
        """
        Returns a symbolic cosine Fourier series of order `n`.

        :param n: Order of the cosine Fourier series.
        :param x: Independent variable
        :param f: Frequency of the cosine Fourier series
        """
        
        # Make the parameter objects for all cosine terms
        a_n = parameters(','.join(['a{}'.format(i) for i in range(n + 1)]))

        # Construct the cosine series
        series = a_n[0] + sum(a_n[i]*(cos(i * f * x)) for i in range(1, n + 1))  # Access a_n elements directly
  
        return series

    def fourier_series(self, x, f, n=0):
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

    def fit_fourier_series(self, data_file, n_order, minimizer, fourier_type="cosine", plot=True, output_file=None):
        """
        Fits a Fourier series (cosine or regular) to data from a .dat file.

        Args:
            data_file: path to the .dat file to be parsed
            n_order: Order of the Fourier series.
            minimizer: List of minimizer objects from symfit.core.minimizers (default: [BFGS]).
            fourier_type: The type of Fourier series to use ('cosine' or 'regular').
            plot: Boolean flag to enable plotting (default: True).
            output_file: Path to a JSON file to save fit results (default: None).
        """
        xdata, ydata = self.read_dat_file(data_file)
        x, y = variables('x, y')
        w, = parameters('w')

        # Choose the model based on the specified Fourier type
        model_dict = {
            y: self.fourier_series(x, f=w, n=n_order) if fourier_type == "regular" else self.fourier_series_cosine(x, f=w, n=n_order)
        }

        print(model_dict)

        # Define a Fit object and fit the model
        fit = Fit(model_dict, x=xdata, y=ydata, minimizer=list(minimizer))
        fit_result = fit.execute(BFGS={'tol': 1e-5})
        print(fit_result)

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

class CoordinateProcessor():
  """Processes coordinates from a log file and creates formatted output files."""

  def __init__(self, logfile_path="xtbopt.log"):
    """Initializes the class with the path to the log file.

    Args:
      logfile_path (str, optional): Path to the log file. Defaults to "xtbopt.log".
    """
    self.logfile_path = logfile_path

  def group_coordinates(self):
    """Groups coordinates from the log file.

    Returns:
      dict: A dictionary containing grouped coordinates.
    """
    with open(self.logfile_path, 'r') as infile:
      firstline = infile.readline()
      num_lines = int(firstline) + 2

    grouped_coords = {}
    with open(self.logfile_path) as f:
      for g, group in groupby(f, key=lambda _, c=count(): next(c) // num_lines):
        grouped_coords[g] = list(group)

    return grouped_coords

  def format_and_write_data(self):
    """Formats grouped coordinates and writes them to individual files."""
    # Create the sp_frames folder if it doesn't exist
    os.makedirs("sp_frames", exist_ok=True)  # Ensure the folder exists

    grouped_coords = self.group_coordinates()

    for i, key in enumerate(grouped_coords):
      formatted_data = self._format_data(grouped_coords[key])
      output_name = f"sp_frames/frame_{i}.xyz"
      self._write_formatted_data(formatted_data, output_name)

  def _format_data(self, data):
    """Formats data for a single group."""
    formatted_data = [data[0].rstrip()]  # Add the first line
    formatted_data.extend(line.strip() for line in data[1:])  # Remove extra spaces
    return formatted_data

  def _write_formatted_data(self, formatted_data, output_name):
    """Writes formatted data to a file."""
    with open(output_name, "w") as f:
      for line in formatted_data:
        f.write(line + "\n")

    print("Data written to {} successfully!".format(output_name))

if __name__ == "__main__":           
  #Example usage XTB_Template
  template_manager = XTB_Template()
  template_string = template_manager.create_template(
      force_constant=0.15,
      dihedral=[2, 3, 4, 5, 62.9],
      scan_start=0.0,
      scan_end=360.0,
      scan_points=100,
  )
  template_manager.save_template_to_file("scan.inp", template_string)

  #Example usage XTB data collect
  #xtb_data_collector = XTB_data_collect()
  #xtb_data_collector.save_xtb_scan_data("xtbscan.log", "scan_data.dat", 100, 0.0, 360.0)

  #Postprocessing xtbscan.log for creating many frames for subsequent Gromacs calculation.
  #processor = CoordinateProcessor("xtbscan.log")
  #processor.format_and_write_data()

  # Example usage for fitting the data AFTER xtb - gromacs energy substraction.
  #fit = FourierFitTool()
  #minimizers = [BFGS, BasinHopping, BFGS, BasinHopping, BFGS, BFGS]
  # Using full Fourier series
  #fit.fit_fourier_series(data_file="scan_data.dat", n_order=5, minimizer=minimizers, fourier_type="regular", plot=True)  
  # Or using cosine Fourier series
  #fit.fit_fourier_series(data_file="scan_data.dat", n_order=3, minimizer=minimizers, fourier_type='cosine', output_file="fit_results.json")  

