#=========================================================================
# Parses spectra text files from the Rogers spectrometer and calculates the
# absorption spectra for DNA (260 nm).
#=========================================================================

import os
import re
import numpy as np
import matplotlib.pyplot as plt


class Spectrum:
    """A class containing all the data from a single spectrum."""

    def __init__(self, eth_gly, mM_NaCl, temperature):
        """The ethylene glycol, mM NaCl, and temperature need to be set in the
        call initializing the spectrum."""

        self.eth_gly = eth_gly
        self.mM_NaCl = mM_NaCl
        self.temperature = temperature

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


if __name__ == '__main__':

    source_dir = '/media/sf_M_DRIVE/DNA Spectra/20160513'

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

        for i in range(1, chunks - 1):
            if filename_parts[i] == 'Eth' and filename_parts[i + 1] == 'Gly' and filename_parts[i - 1].isdigit():
                eth_gly = float(filename_parts[i - 1])

            # Get the salt content

            if filename_parts[i] == 'mM' and filename_parts[i + 1] == 'NaCl' and filename_parts[i - 1].isdigit():
                mM_NaCl = float(filename_parts[i - 1])

        # Check if it is a blank or a sample

        if 'blank' in filename_parts:

            # This is our blank spectrum
            temperature = 'blank'

        else:

            # This is a sample measurement, so we need to get the temperature

            temperature_inds = re.search("[0-9]C", filename)
            temperature = float(
                filename[temperature_inds.start() - 1:temperature_inds.end() - 1])

        # print filename + ': ' + str(eth_gly) + ' Ethylene Glycol, ' +
        # str(mM_NaCl) + ' mM NaCl, ' + str(temperature) + ' C'

        data = np.loadtxt(spectrum_file, delimiter="\t", skiprows=16)

        lambdas = data[:, 0]
        intensities = data[:, 1]
        
        # Save to a Spectrum object

        spectrum_obj = Spectrum(eth_gly, mM_NaCl, temperature)
        spectrum_obj.add_data(lambdas, intensities)
        
        # Add the spectrum to an existing Experiment or create a new one
        
        if eth_gly in experiment_dict:
            
            experiment_dict[eth_gly].add_spectrum(spectrum_obj)
        
        else:
            
            experiment_dict[eth_gly] = Experiment('temp')
            experiment_dict[eth_gly].add_spectrum(spectrum_obj)
    
    for key in experiment_dict:
        
        exp = experiment_dict[key]
        exp.plot_melting()
