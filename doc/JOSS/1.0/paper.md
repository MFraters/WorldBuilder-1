---
title: 'The Geodynamic World Builder: An initial conditions generator for the geosciences'
tags:
  - C++
  - CPP
  - C
  - Fortran
  - Python
  - Geodynamics
  - Earth science
  - Tectonics
  - Seismology
authors:
  - name: Menno R. T. Fraters
    orcid: 0000-0003-0035-7723
    affiliation: "1" # (Multiple affiliations must be quoted)
affiliations:
 - name: Universiy of Florida, USA
   index: 1
date: 9 September 2023
bibliography: paper.bib
---

# Summary

Geodynamic models can be broadly divided into three classes. The first class are synthetic process-oriented models, which feature simple geometries and are easy to set up and parameterize. The second class are data driven models, which derive their initial conditions (semi-)automatically from datasets like topography, or geophysical imaging like satellite images, seismic or electromagnetic tomography. The third and very common class of models are complex synthetic models. This class features complex geometrical features based on the modelers expectation of the geologic situation that are nevertheless parameterized and not data-based, because the purpose of the modeling study requires easy changes to the initial condition, for example to perform a parameter study. The Geodynamic World Builder (GWB) has been designed to simplify the creation of this third class of models, but can also be used to design simple synthetic models and several types of datasets are supported to be used to setup models in a way which aligns with the GWB design philosophies. Besides setting up initial conditions for geodynamic models, the GWB can also be used to make detailed 3D visualizations of geologic and geodynamic settings.

# Statement of need

The increased availabability of computational resources and new numerical methods and infrastructure, has enabled the creation of 3D numerical models which closer resemble or mimic actual location on the Earth with a much finer resolution. This also means that setting up these more detailed regional or global models, becomes increasingly difficult, both in 2D, but especially in 3D. Furthermore, to properly investigate these models, often many smaller and bigger changes to the model need to be tested, which requires a way to change these models quickly and easily. Although successful attempts have been published, such model setups often have one or many of the following issues:

1. The configuration is not human-readable.
2. The software is not easily modifiable and extendable.
3. The model setup is not portable to other computing systems or reproducible in other software frameworks.
4. The model setup is not shareable with other users.

These issues lead to a number of problems with the reproducibility and reliability of geodynamic modeling studies, which threaten to undermine the predictive power and usefulness of modeling results. The GWB was designed to solve this situation, by creating human readable, parameterized, portable, reproducible and shareable geodynamic model setups. As a side benefit, since the GWB comes with its own programs to visualize the constructed model through programs like Paraview, and creating the models requires no programming knowledge, it can be easily used to visualize tectonic and geodynamic settings for publications, teaching or public outreach. 

# Methods

The Geodynamic World Builder is designed to solve the problems stated in the previous section through implementing specific code and user design principles.

## GWB Code Principles
GWB is build around the following principles: 1. A single, text-based, and human readable input file using a standard conforming JSON file. Following an established format like JSON allows for the use of pre-existing readers and writers in many common programming languages. 2. Software, language, and platform independent software that supports **Linux**, **OSX** and **Windows** and has official interfaces for coupling to **C++**, **C**, **FORTRAN** and **Python** software. 3. Support for parallel execution so that GWB can be easily used from laptops to modern supercomputers. 4. Readable and extensible software design through modern modular programming paradigms. 5. Strict version numbering of source and input files to ensure reproducible results.

These principles are implemented in an object-oriented C++ code with the mentioned interfaces. All variable parts of the software are implemented as plugin systems using interface classes that decouple individual modules and allow to easily extend the code with new features. The code is extensively tested with unit and integration tests and continuously build and tested.

## GWB Use Principles

GWB's user philosophy is build around the idea that a complex model region can be split into individual tectonic features. These tectonic features can be parameterized by defining their areal extent or contour and their properties in a map view. For example, a continental plate can be defined as an area on a map. A fault is a line feature, so the user defines that fault line as a feature. The user also provides information on the properties of the feature, such as a fault or plate thickness, which can be spatially variable, or dip angles. These parameters are then used by the GWB world builder plugin to create the volume extent of the feature. Next users can attach one or many models to those volumes to define quantities like temperature or composition variables. These can be very simple models, such as a uniform temperature distribution, or a more complex distribution, such as a half space cooling model, or a McKenzie (ref) or mass conserving (ref) slab temperature model.

All these parameterizations are bundled in a single JSON styled input file, which is human readable, writeable and editable. The main idea behind this design of the GWB so that users can easily create and modify complex parameterized initial conditions for their geodynamic or tectonic setting, which has been illustrated in models as diverse as ... TODO list a number of GWB application cases with their paper here. 

## Example 
Below we show an example input file for a cartesian model of a subducting plate beneath an overriding ocean that specifies position, volume, temperature, and composition for both plates. This model can be easily extended to a spherical model with additional capabilities, such as shown in the example figure below, which was build using an input file of only 66 lines.

```json
{
  "version": "1.0",
  "coordinate system":{"model":"cartesian"},
  "features":
  [
    {
       "model":"oceanic plate", "name":"Overriding Plate", "max depth":100e3, 
       "coordinates":[[0,0],[0,1000e3],[1500e3,1000e3],[1600e3,350e3],[1500e3,0]],
       "temperature models":[{"model":"linear", "max depth":100e3}],
       "composition models":[{"model":"uniform", "compositions":[0]}]
    },
    { 
      "model":"subducting plate", "name":"Slab", "dip point":[0,0],
      "coordinates":[[1500e3,1000e3],[1600e3,350e3],[1500e3,0]],
      "segments":[{"length":300e3, "thickness":[100e3], "angle":[0,60]}],
      "temperature models":[{"model":"plate model", "plate velocity":0.02}],
    }
  ]
}
```

![An schematic example of what can be build with 66 lines of GWB inptu file formatted in the same way as in the above input file example. \label{fig:example}](../../sphinx/_static/images/user_manual/basic_starter_tutorial/BST_17.png)

# Acknowledgements

CIG, NSF projects Cascadia and FRES
see todo

# References


# TODO
1. Add citations to other initial conditions generators
2. Decide on co-authors
3. write acknowledgments for other contributors who are not co-authors and for funding agencies
4. Update the example with the new contours slab and maybe add a fault as well

