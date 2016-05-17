# DNA Spectra

This is a piece of code written to quickly analyze absorbance spectra of DNA samples obtained in the lab in which I work. The data is written to file by the spectrometer (mfd. by Ocean Optics). This program is very dependent upon filename conventions, so make sure that as the experiments proceed that all the data files are named appropriately. See below section for more info.

## Filename Conventions

There are two main types of experiments this program is made to handle: time series and temperature series. A time series experiment consists of setting a temperature on the spectrometer and measuring absorbance of your sample as a function of time. This basically measures the time it takes to heat up or cool down your sample. A temperature series experiment consists of changing the temperature of your sample, letting it equilibrate at the new temperature and then recording a spectrum. This allows for equilibrium measurements of absorbance rather than the kinetic measurements of a time series experiment.

The two types of experiment have their own file naming conventions. For time series, the format looks like this:

`<text A>_Time_<text B>_<number 1>_Eth_Gly_<number 2>_mM_NaCl_<number 3>C_<timestamp>.txt`

where

* `<text A>` optional and can be any text you want so long as it doesn't conflict with any other markers
* `<text B>` optional and can be any text you want so long as it doesn't conflict with any other markers
* `<number 1>` is the percent by volume of Ethylene Glycol in the sample
* `<number 2>` is the salt (NaCl) concentration of the sample in mM
* `<number 3>` is the temperature of the sample in Celcius
* `<timestamp>` is the automatically generated timestamp from the spectrometer

For Temperature series experiments, the filename convention looks like this:

`<text A>_<number 1>_Eth_Gly_<number 2>_mM_NaCl_<number 3>C_<timestamp>.txt`

where

* `<text A>` optional and can be any text you want so long as it doesn't conflict with any other markers. Especially make sure it does not contain the word time in it, otherwise it will be allocated to the wrong type of experiment.
* `<number 1>` is the percent by volume of Ethylene Glycol in the sample
* `<number 2>` is the salt (NaCl) concentration of the sample in mM
* `<number 3>` is the temperature of the sample in Celcius

## Dependencies

* Python 2.x

## Contact

Walter Schwenger, wjs018@gmail.com
