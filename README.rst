.. |ech2o| replace:: EcH\ :sub:`2`\ O

|ech2o|-iso, a process-based ecohydrologic model
=================================================

|ech2o|-iso builds on the **process-based, spatially-distributed ecohydrologic model EcH**\ :sub:`2`\ **O** developed in C++ in the Regional Hydrology Lab at the University of Montana (Maneta & Silverman, 2013) (`link <http://hs.umt.edu/RegionalHydrologyLab/software/default.php>`_).

The general structure of the |ech2o| model is given in Fig. 1.

.. figure:: ./EcH2O_Model.png
   :width: 400px
   :align: center
   :height: 300px
   :figclass: align-center

   **Figure 1.** Conceptual diagram of the structure of the |ech2o| model. The grid cells of the simulation domain (a) are laterally connected via overland runoff, streamflow and lateral flow in the saturated profile (DEM-derived drainage network), with for each grid cell a process-oriented description of (b) the energy balance, (c) hydrological transfers and (d) vegetation growth and dynamics. Adapted from Kuppel et al. (2018a) and Douinot et al. (2019).


The specific features of |ech2o|-iso is the implementation of stable water isotopes (:sup:`2`\ H and :sup:`18`\ O), chloride and age tracking.
It is mostly based on an immediate, complete implicit scheme for mixing water transiting between compartments within time steps (Fig. 2).
Fractionation of isotopes during evaporation, evapoconcentration of chloride during evaporation and transpiration, are also included.


.. figure:: ./EcH2O-iso_Model.png
   :width: 400px
   :align: center
   :height: 300px
   :figclass: align-center

   **Figure 2.** Water compartments (black rectangles) and fluxes (coloured arrows) as represented in |ech2o| and used for isotope and age tracking in |ech2o|-iso, with the numbers between brackets reflecting the sequence of calculation within a time step. Note that water routing (steps [8] to [13]) differs between cells where a stream is present (◦) or not (∗). 


**Other EcH2O(-iso) repositories**

The present repository provides the version(s) of |ech2o|-iso primarily maintained by S. Kuppel and colleagues. Depending on the research questions addressed and thus the model features sought by potential users, other relevant repositories are:

- the |ech2o| repository maintained by M. Maneta and colleagues at the University of Montana (`link <http://bitbucket.org/maneta/ech2o>`_).
- the |ech2o|-iso repository maintained at the Leibniz Institute for Freshwater Ecology and Inland Fisheries (IGB, `link <http://bitbucket.igb-berlin.de:7990/users/ech2o/repos/ech2o_iso/browse>`_).
  
   
**Latest Version**

The latest stable version can be found in the *master* branch of the `source repository <https://bitbucket.org/scicirc/ech2o-iso/src/master/>`_. 


**Documentation**

The documentation for installing an runnnig |ech2o|-iso, available as of the date of this release, can be found on its `ReadTheDocs webpage <http://ech2o-iso.readthedocs.io/en/latest/>`_.



**Third-party dependencies**

|ech2o|-iso depends on the following third-party libraries with the following licenses:
  
- armadillo (Mozilla Public License 2.0) and dependencies therein 
- libcsf (BSD License)
  
For convenience, precompiled versions of the libcsf librairies for Linux, Windows 64 bit, and Mac architectures are distributed with the source code.   


**Compilation of source code**

Please see the fille called INSTALL.rst

**Data Preprocessing**

|ech2o|-iso uses the PCRASTER map format (a cross-system format) for data pre- and post-processing, and for visulalization. 
PCRASTER can be downloaded free of charge from http://pcraster.geo.uu.nl/downloads

**Licensing**

Please see the file called LICENSE.txt.

**Bugs**

Should you encounter any bug, please file a ticket `here <https://bitbucket.org/scicirc/ech2o-iso/issues>`_.
Known issues can be found there, as well as on the `main EcH2O page <https://bitbucket.org/maneta/ech2o/issues>`_.

**How to Cite**

Please, acknowledge the use of |ech2o|-iso by citing:

- Kuppel, S, Tetzlaff, D, Maneta, MP & Soulsby, C (2018). |ech2o|-iso 1.0: Water isotopes and age tracking in a process-based, distributed ecohydrological model, Geosci. Model Dev., 11, 3045-3069, `<https://doi.org/10.5194/gmd-11-3045-2018>`_.
  
Further references documenting previous/specific developments and/or application of the code:

- Maneta, MP & Silverman, N (2013). A spatially-distributed model to simulate water, energy and vegetation dynamics using information from regional climate models. Earth Interactions, 17, 1-44, [`link <https://doi.org/10.1175/2012EI000472.1>`_].
- Lozano-Parra, J, Maneta, MP & Schnabel, S (2014). Climate and topographic controls on simulated pasture production in a semiarid Mediterranean watershed with scattered tree cover. Hydrology and Earth System Sciences, 18, 1439 [`link <https://doi.org/10.5194/hess-18-1439-2014>`_].
|  
- Kuppel, S et al. (2018). What can we learn from multi-data calibration of a process-based ecohydrological model?. Environmental Modelling & Software, 101, 301–316 [`link <https://doi.org/10.1016/j.envsoft.2018.01.001>`_].
- Maneta, MP et al. (2018). Conceptualizing catchment storage dynamics and nonlinearities. Hydrological Processes, 32, 3299–3303 [`link <https://doi.org/10.1002/hyp.13262>`_].
|
- Douinot A et al. (2019). Ecohydrological modelling with EcH2O-iso to quantify forest and grassland effects on water partitioning and flux ages. Hydrological Processes 33 (16): 2174–2191 [`link <https://doi.org/10.1002/hyp.13480>`_].
- Simeone, C et al. (2019). Coupled ecohydrology and plant hydraulics modeling predicts ponderosa pine seedling mortality and lower treeline in the US Northern Rocky Mountains. New Phytologist, 221(4), 1814-1830, [`link <https://doi.org/10.1111/nph.15499>`_].
- Smith A et al. (2019). Assessing the influence of soil freeze–thaw cycles on catchment water storage–flux–age interactions using a tracer-aided ecohydrological model. Hydrology and Earth System Sciences 23 (8): 3319–3334, [`link <https://doi.org/10.5194/hess-23-3319-2019>`_].
|
- Knighton, J et al. (2020). Using isotopes to incorporate tree water storage and mixing dynamics into a distributed ecohydrologic modelling framework. Ecohydrology, 13(3), e2201 [`link <https://doi.org/10.1002/eco.2201>`_].
- Kuppel, S et al. (2020). Critical zone storage controls on the water ages of ecohydrological outputs. Geophysical Research Letters, 47, e2020GL088897 [`link <https://doi.org/10.1029/2020GL088897>`_].
- Neill, AJ et al. (2020). An agent-based model that simulates the spatio-temporal dynamics of sources and transfer mechanisms contributing faecal indicator organisms to streams. Part 1: Background and model description. Journal of environmental management, 270, 110903 [`link <https://doi.org/10.1016/j.jenvman.2020.110903>`_].
- Neill, AJ et al. (2020). An agent-based model that simulates the spatio-temporal dynamics of sources and transfer mechanisms contributing faecal indicator organisms to streams. Part 2: Application to a small agricultural catchment. Journal of environmental management, 270, 110905 [`link <https://doi.org/10.1016/j.jenvman.2020.110905>`_].
- Smith, A et al. (2020). Isotope‐aided modelling of ecohydrologic fluxes and water ages under mixed land use in central Europe: the 2018 drought and its recovery. Hydrological Processes, 34(16), 3406-3425 [`link <https://doi.org/10.1002/hyp.13838>`_].
|
- Gillefalk, M et al. (2021). Quantifying the effects of urban green space on water partitioning and ages using an isotope-based ecohydrological model. Hydrology and Earth System Sciences, 25, 3635-3652 [`link <https://doi.org/10.5194/hess-25-3635-2021>`_].
- Kleine, L et al (2021). Modelling ecohydrological feedbacks in forest and grassland plots under a prolonged drought anomaly in Central Europe 2018–2020. Hydrological Processes, 35(8), e14325 [`link <https://doi.org/10.1002/hyp.14325>`_].
- Neill, AJ et al. (2021). Structural changes to forests during regeneration affect water flux partitioning, water ages and hydrological connectivity: Insights from tracer-aided ecohydrological modelling. Hydrology and Earth System Sciences, 25, 4861-4886 [`link <https://doi.org/10.5194/hess-25-4861-2021>`_].
- Smith, A., et al. (2021). Quantifying the effects of land use and model scale on water partitioning and water ages using tracer-aided ecohydrological models. Hydrology and Earth System Sciences, 25(4), 2239-2259 [`link <https://doi.org/10.5194/hess-25-2239-2021>`_]..
- Yang, X, et al. (2021). Catchment Functioning Under Prolonged Drought Stress: Tracer‐Aided Ecohydrological Modeling in an Intensively Managed Agricultural Catchment. Water Resources Research, 57(3), e2020WR029094 [`link <https://doi.org/10.1029/2020WR029094>`_].
|
- Gillefalk, M., et al. (2022). Estimates of water partitioning in complex urban landscapes with isotope‐aided ecohydrological modelling. Hydrological Processes, 36(3), e14532, [`link <https://doi.org/10.1002/hyp.14532>`_].
- Smith, A, et al. (2022). Critical Zone Response Times and Water Age Relationships Under Variable Catchment Wetness States: Insights Using a Tracer‐Aided Ecohydrological Model. Water Resources Research, 58(4), e2021WR030584 [`link <https://doi.org/10.1029/2021WR030584>`_]..
- Smith, A., et al. (2022). Modelling temporal variability of in situ soil water and vegetation isotopes reveals ecohydrological couplings in a riparian willow plot. Biogeosciences, 19(9), 2465-2485, [`link <https://doi.org/10.5194/bg-19-2465-2022>`_].
|
- Li, K., et al. (2023). Parameterizing Vegetation Traits with a Process‐Based Ecohydrological Model and Xylem Water Isotopic Observations. Journal of Advances in Modeling Earth Systems, e2022MS003263 [`link <https://doi.org/10.1029/2022MS003263>`_].
- Yang, X., et al. (2023). Upscaling tracer-aided ecohydrological modeling to larger catchments: Implications for process representation and heterogeneity in landscape organization. Water Resources Research, 59, e2022WR033033. [`link <https://doi.org/10.1029/2022WR033033>`_].
- Wu, S., et al. (2023). Integrating tracers and soft data into multi-criteria calibration: Implications from distributed modeling in a riparian wetland. Water Resources Research, 59, e2023WR035509. [`link <https://doi.org/10.1029/2023WR035509>`_].
- Ackerer, J., et al. (2023).  Exploring the critical zone heterogeneity and the hydrological diversity using an integrated ecohydrological model in three contrasted long-term observatories. Water Resources Research, 59, e2023WR035672. [`link <https://doi.org/10.1029/2023WR035672>`_].

**Contacts**

If you have any questions, please contact sylvain[dot]kuppel[at]ird[dot]fr.
