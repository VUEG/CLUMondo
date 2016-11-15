# CLUMondo

CLUMondo is based on the land systems approach. Land systems are socio-ecological systems that reflect land use in a spatial unit in terms of land cover composition, spatial configuration, and the management activities employed. The precise definition of land systems depends on the scale of analysis, the purpose of modelling, and the case study region. In contrast to land cover classifications the role of land use intensity and livestock systems are explicitly addressed. Each land system can be characterized in terms of the fractional land covers.

This repository contains the source code for the CLUMondo model, which has a command-line interface (CLI). For downloading the graphical user intergafe (GUI, non-OSS), see [the Environmental Geography Group homepage](http://www.environmentalgeography.nl/site/data-models/data/clumondo-model/).

### Core publications an example applications of CLUMondo:

+ Van Asselen, S. and Verburg, P.H., 2012. A land system representation for global assessments and land-use modeling. Global Change Biology 18, 3125-3148. http://dx.doi.org/10.1111/j.1365-2486.2012.02759.x
+ Van Asselen, S., Verburg, P.H., 2013. Land cover change or land-use intensification: simulating land system change with a global-scale land change model. Global Change Biology 19, 3648–3667. http://dx.doi.org/10.1111/gcb.12331
+ Eitelberg, D.A., van Vliet, J., Verburg, P.H., 2015. A review of global potentially available cropland estimates and their consequences for model-based assessments. Global Change Biology 21, 1236–1248. http://dx.doi.org/10.1111/gcb.12733
+ Ornetsmüller, C., Verburg, P. & Heinimann, A. (2016). Scenarios of land system change in the Lao PDR: Transitions in response to alternative demands on goods and services provided by the land. Applied Geography, 75, 1-11. http://dx.doi.org/10.1016/j.apgeog.2016.07.010

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisities

**NOTE**: currently CLUMondo reliably rund only on 32-bit Windows.

CLUMondo is built using pure `std` libraries. For building CLUMondo from source, you will need the following packages:

#### Linux

+ [CMake](https://cmake.org/)
+ [Make](https://www.gnu.org/software/make/)
+ [gcc-c++](https://gcc.gnu.org/projects/cxx-status.html)


### Installing

CLUMondo uses [CMake](https://cmake.org/) to genrerate the needed build files on different platforms. To build, do the following in the source code folder:

```
cmake .
make
```

## Authors

* **Peter Verburg** <peter.verburg@vu.nl> 

See also the list of [contributors](https://github.com/VUEG/CLUMondo/contributors) who participated in this project.

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE.md](LICENSE.md) file for details
