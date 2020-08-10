# Realistic Microstructure Simulator (RMS): Monte Carlo simulations of diffusion in 3D cells (CUDA C++)

## Part 1: Diffusivity time-dependence along realistic white matter axons
The code implements 3d Monte Carlo simulations originally developed in [Lee, et al., Journal of Neuroscience Methods, 2020](), demonstrating the power spectrum of realistic axonal shapes along white matter axons (Fig. 1) and the diffusivity and kurtosis time-dependence along axons (Fig. 2) in [Lee, et al., Communications Biology 2020](https://doi.org/10.1038/s42003-020-1050-x). The white matter axons were segmented from the corpus callosum of a mouse brain, with details in [Lee, et al., Brain Structure and Function](https://doi.org/10.1007/s00429-019-01844-6) and our another Github toolbox [Random Walker (RaW) segmentation](https://github.com/NYU-DiffusionMRI/RaW-seg).

* **Demo 1, power spectrum:** Calculate the power spectrum of realistic axonal shapes along axons (Fig. 1, Fig. 6, and Supplementary Fig. 1).
* **Demo 2, artificial shape generation:** Generation of artificially designed microgeometry based on realistic axons (Fig. 2a).
* **Demo 3, simulations:** Perform Monte Carlos simulations of diffusion 3d cell geometries. The code is implemented in CUDA C++, and you need an Nvidia GPU to run the code.
* **Demo 4, analysis:** Calculate diffusivity and kurtosis time-dependence based on displacement cumulants (Fig. 2b-h).

The simulation of wide pulse sequence and permeability-based exchange will be released in another repository.

## Part 2: Why elastic collision is the most reliable particle-membrane interaction? (This part will be updated soon.)
The code implements 1d, 2d, and 3d Monte Carlo simulations for the educational purpose, with details in the Appendices A and B in [Lee, et al., Journal of Neuroscience Methods, 2020](), demonstrating the bias caused by the following two particle-membrane interactions: equal-step-length random leap (ERL) and rejection sampling. This part justifies the choice of elastic collision in the RMS. The codes in this part are implemented in Matlab.

* **Demo 1, equal-step-length random leap:** Perform simple Monte Carlo simulations of diffusion between 1d, 2d, and 3d impermeable parallel planes, and show the inhomogeneous particle density around membranes and bias in the diffusivity transverse and parallel to the planes caused by the ERL.
* **Demo 2, rejection sampling:** Perform simple Monte Carlo simulations of diffusion between 1d, 2d, and 3d (1) impermeable parallel planes to show the bias in the diffusivity parallel to membranes caused by rejection sampling, and (2) permeable parallel planes to demonstrate the compatibility of rejection sampling with the simulation of water exchange.

These are good exercises if you just start your own MC simulation codes.
Some results can suprise you, even if you are well experienced!!

## References
* **Monte Carlo simulation**
  - [Fieremans, et al., NMR Biomed, 2010](https://doi.org/10.1002/nbm.1577)
  - [Novikov, et al., Nature Physics, 2011](https://doi.org/10.1038/nphys1936)
  - [Fieremans and Lee, NeuroImage 2018](https://doi.org/10.1016/j.neuroimage.2018.06.046)
  - [Lee, et al., Communications Biology 2020](https://doi.org/10.1038/s42003-020-1050-x)
  - [Lee, et al., Journal of Neuroscience Methods 2020]()

* **Realistic axonal shape**
  - [Lee, et al., Brain Structure and Function](https://doi.org/10.1007/s00429-019-01844-6)
  - [Electron microscopy data](https://www.cai2r.net/resources/software/intra-axonal-space-segmented-3d-scanning-electron-microscopy-mouse-brain-genu)
  - [Random Walker (RaW) segmentation](https://github.com/NYU-DiffusionMRI/RaW-seg)
  
* **Watson Distribution sampling**
  - [SphericalDistributionsRand](https://www.mathworks.com/matlabcentral/fileexchange/52398-sphericaldistributionsrand)

## Authors
* [Hong-Hsi Lee](http://www.diffusion-mri.com/people/hong-hsi-lee)
* [Dmitry S Novikov](http://www.diffusion-mri.com/people/dmitry-novikov)
* [Els Fieremans](http://www.diffusion-mri.com/people/els-fieremans)

## Acknowledgement
We would like to thank [Sune N Jespersen](https://pure.au.dk/portal/en/persons/sune-jespersen(f4d1a00c-677b-4aca-b9b0-c7ad14f1fddc).html) for the fruitful discussion of simulation implementation and the theory of diffusion.

## License
This project is licensed under the [LICENSE](https://github.com/NYU-DiffusionMRI/monte-carlo-simulation-3D-RMS/blob/master/LICENSE).
