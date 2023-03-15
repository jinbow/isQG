
================================================
interior+surface Quasi-Geostrophic (isQG)
================================================

The latest documentation can be found here: https://isqg.readthedocs.io/en/latest/.

This software applies QG framework to retrieve ocean interior state using Sea Surface Temperature (SST), Sea Surface Height (SSH) anomaly fields and ocean stratification as the input. It has assumption that the considered ocean physics follow QG dynamics, which is a fair assumption for spatial scales larger than :math:`\mathcal{O}(10km)`.

Installation
============

To install isQG, simply clone this repository and run the setup file using the following command:

.. code-block:: sh

    $ git clone https://github.com/username/isqg.git
    $ cd isqg
    $ python setup.py install

Dependencies
============

isQG requires the following dependencies to be installed:

Usage
=====

isQG can be used to reconstruct the three-dimensional subsurface state of the ocean based on the Sea Surface Temperature (SST), Sea Surface Height (SSH) anomaly fields, and ocean stratification climatology.

To use isQG, first prepare the input data (SST, SSH anomaly, and ocean stratification climatology) in NetCDF format. Then, use the `isqg_reconstruction.py` script to perform the reconstruction. An example command is provided below:

.. code-block:: sh

    $ python isqg_reconstruction.py -i input_data.nc -o output_data.nc

This will output the reconstructed subsurface state in NetCDF format.

Contributing
============

Contributions to isQG are welcome! Please feel free to fork this repository and submit pull requests with your improvements.

If you encounter any bugs or have suggestions for new features, please submit an issue on the GitHub issue tracker.

License
=======

isQG is licensed under the MIT License. See the LICENSE file for more information.

Citation
===============

This software is based on the interior+surface Quasi-geostrophic (isQG) framework proposed by Wang et al. in 2013. 
