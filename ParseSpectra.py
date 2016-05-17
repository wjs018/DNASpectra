"""Parses spectra text files from the Rogers spectrometer and calculates the
absorption spectra for DNA (260 nm).

A lot of the automation in this program is done through using filename
conventions. Files should be named according to the following:

For Temperature series experiments:

<any text but time>_<number>_Eth_Gly_<number>_mM_NaCl_<number>C.txt

For time series experiments:

<any text containing time>_<number>_Eth_Gly_<number>_mM_NaCl_<number>C_<timestamp from spectrometer>.txt

For a more detailed description of these, see the readme file.

"""

import os
import re
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt


class Spectrum:
    """A class containing all the data from a single spectrum."""

    def __init__(self, eth_gly, mM_NaCl, temperature, time=None):
        """The ethylene glycol, mM NaCl, and temperature need to be set in the
        call initializing the spectrum. Time is an optional parameter to be 
        used in time series experiments."""

        self.eth_gly = eth_gly
        self.mM_NaCl = mM_NaCl
        self.temperature = temperature
        self.time = time

        if temperature == 'blank':
            self.blank = True
        else:
            self.blank = False

    def add_data(self, lambdas, intensities):
        """Add the wavelength and intensity data to the spectrum. Should be in
        the form of numpy arrays."""

        self.lambdas = lambdas
        self.intensities = intensities


class Experiment:
    """A class containing all the spectra from a single experiment."""

    def __init__(self, exp_type):
        """Define the experiment type when initializing an Experiment. Can
        either be 'time' for a time series measurement, or 'temp' for a 
        temperature series experiment."""

        self.spectra_list = []
        self.abs_list = []
        self.blank_spectrum = False
        self.exp_type = exp_type

    def add_spectrum(self, spectrum):
        """Add one more spectrum to the experiment."""

        if spectrum.blank == True:
            self.blank_spectrum = spectrum
        else:
            self.spectra_list.append(spectrum)

    def get_temps(self):
        """Return the list of temperatures for the spectra in the Experiment."""

        temp_list = []

        for spec in self.spectra_list:

            temp_list.append(spec.temperature)

        return temp_list

    def calc_abs(self):
        """Calculate the absorbance spectra."""

        if self.blank_spectrum:

            for spec in self.spectra_list:

                trans = np.divide(
                    spec.intensities, self.blank_spectrum.intensities)
                trans = trans.clip(min=1e-10)

                absorb = - np.log10(trans)
                abs_spec = Spectrum(
                    spec.eth_gly, spec.mM_NaCl, spec.temperature)
                abs_spec.add_data(spec.lambdas, absorb)

                self.abs_list.append(abs_spec)

        else:
            print "No blank spectrum set!"

    def get_abs_maxes(self):
        """Get the 260 nm absorbance maximum from each spectrum."""

        abs_max_list = []

        if len(self.abs_list) == 0:
            self.calc_abs()

        if len(self.abs_list) == len(self.spectra_list):

            lambdas = self.blank_spectrum.lambdas

            ix = np.where((258 < lambdas) & (lambdas < 262))

            for spec in self.abs_list:

                abs_max_list.append(np.max(spec.intensities[ix]))

        return abs_max_list

    def plot_melting(self):
        """Plots Absorbance vs. Temperature curve."""

        if len(self.abs_list) == 0:
            self.calc_abs()

        if len(self.abs_list) == len(self.spectra_list):

            temps = self.get_temps()
            maxes = self.get_abs_maxes()

            plt.plot(temps, maxes, 'o')
            plt.title(str(self.spectra_list[0].eth_gly) + "% Ethylene Glycol")
            plt.show()

    def get_times(self):
        """Returns a list of times the spectra in the experiment were taken."""

        if self.exp_type != 'time':
            print "Experiment is wrong type for this function."
            return None

        time_list = []

        for spec in self.spectra_list:

            time_list.append(spec.time)

        return time_list

    def plot_time(self):
        """Plots absorption as a function of time."""

        if len(self.abs_list) == 0:
            self.calc_abs()

        if len(self.abs_list) == len(self.spectra_list):

            times = self.get_times()
            maxes = self.get_abs_maxes()

            plt.plot(times, maxes, 'o')
            plt.title(str(self.spectra_list[0].eth_gly) + "% Ethylene Glycol")
            plt.show()

    def save(self, output_file):
        """Saves the results of the experiment. For time series, it will save a
        csv file with two columns, time and absorbance. For temperature series,
        it will save a csv with two columns, temperature and absorbance."""
        
        import csv
        
        outfile = open(output_file, 'wb')
        writer = csv.writer(outfile)
        
        if self.exp_type == 'time':
            
            col1 = self.get_times()
            col2 = self.get_abs_maxes()
        
        elif self.exp_type == 'temp':
            
            col1 = self.get_temps()
            col2 = self.get_abs_maxes()
        
        for i in range(len(col1)):
            
            writer.writerow([col1[i], col2[i]])
        
        outfile.close()

def parse_folder(dir_path):
    """Parse the DNA spectra in the given directory. Returns a dictionary of 
    the Experiment objects discovered from the files in the directory."""

    source_dir = dir_path

    for root, dirs, files in os.walk(source_dir):

        files = [os.path.join(root, f) for f in files if f.endswith('.txt')]

    # Experiments distinguished according to ethylene glycol content

    experiment_dict = {}

    for spectrum_file in files:

        # First, get the filename separate from the directory structure

        os_test = os.path.join('test', 'path')

        split_char = os_test[4]

        filename = spectrum_file.split(split_char)[-1]
        filename = filename[0:-4]

        # Next, get the ethylene glycol content

        filename_parts = filename.split('_')
        chunks = len(filename_parts)

        # Determine if this spectrum is a blank

        blank_flag = False

        if 'blank' in filename_parts:

            temperature = 'blank'
            blank_flag = True

        for i in range(1, chunks - 1):

            # Get the ethylene glycol content

            if filename_parts[i] == 'Eth' and filename_parts[i + 1] == 'Gly' and filename_parts[i - 1].isdigit():
                eth_gly = float(filename_parts[i - 1])

            # Get the salt content

            if filename_parts[i] == 'mM' and filename_parts[i + 1] == 'NaCl' and filename_parts[i - 1].isdigit():
                mM_NaCl = float(filename_parts[i - 1])

        # Extract the temperature if it is not a blank

        if not blank_flag:
            temperature_inds = re.search("[0-9]C", filename)
            temperature = float(
                filename[temperature_inds.start() - 1:temperature_inds.end() - 1])

        # Actually read in the data from the text file (16 rows of header)

        data = np.loadtxt(spectrum_file, delimiter="\t", skiprows=16)

        lambdas = data[:, 0]
        intensities = data[:, 1]

        # Save to a Spectrum object

        spectrum_obj = Spectrum(eth_gly, mM_NaCl, temperature)
        spectrum_obj.add_data(lambdas, intensities)

        # Check whether this is a temperature or time series experiment

        if not any('time' in s.lower() for s in filename_parts):

            # This is a temperature series experiment

            exp_time = None
            exp_type = 'temp'
            exp_key = str(eth_gly) + '_' + str(mM_NaCl) + '_' + exp_type

        elif any('time' in s.lower() for s in filename_parts):

            # This is a time series experiment, we need to extract the timestamp
            # unless it is a blank

            if not blank_flag:

                time_str = filename_parts[-1]
                time_parts = time_str.split('-')

                # We need to convert strings into ints for the time object

                for i in range(len(time_parts)):
                    time_parts[i] = int(time_parts[i])

                exp_short_time = dt.time(
                    time_parts[0], time_parts[1], time_parts[2], time_parts[3])
                today_date = dt.date.today()
                exp_time = dt.datetime.combine(today_date, exp_short_time)

            exp_type = 'time'
            exp_key = str(eth_gly) + '_' + str(mM_NaCl) + '_' + exp_type

        # Save to a Spectrum object

        spectrum_obj = Spectrum(eth_gly, mM_NaCl, temperature, time=exp_time)
        spectrum_obj.add_data(lambdas, intensities)

        # Add the spectrum to an existing Experiment or create a new one

        if exp_key in experiment_dict:

            experiment_dict[exp_key].add_spectrum(spectrum_obj)

        else:

            if exp_time:
                experiment_dict[exp_key] = Experiment('time')
            else:
                experiment_dict[exp_key] = Experiment('temp')

            experiment_dict[exp_key].add_spectrum(spectrum_obj)

    # Return the dictionary of experiments

    return experiment_dict


if __name__ == '__main__':

    #=========================================================================
    # Change the line below if you want to specify a different directory.
    #=========================================================================
    source_dir = '/media/sf_M_DRIVE/DNA Spectra/20160513'

    experiment_dict = parse_folder(source_dir)

    # Plot results depending on type of experiment

    for key in experiment_dict:

        exp = experiment_dict[key]

        if exp.exp_type == 'temp':

            exp.plot_melting()

        else:

            exp.plot_time()
        
        # Save results to a csv file
        
        exp.save(os.path.join(source_dir, key + '.csv'))