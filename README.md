# MDM-TexPPAC system simulation
TODO

## Prerequisites
- Geant4 v11 (for v10, switch to ``G4v10`` branch)
- C++17 (C++14 for v10)
- CMake>=3.16

## How to use
TODO. Below are not useful at all.
### Compile
Make a build directory and enter it.
```
$ mkdir build && cd build
```

cmake & make 
```
$ cmake ../ && make
```

### Configuration
TODO
<!-- This program takes JSON file as configuration (it is mandatory!). The example config file is ``config/config.json``. Here below shows the content of this JSON file.
```
{
    "GUI": true,
    "RunMac": "run1.mac",
    "Threads":10,
    "Foil": 2,
    "BeamEnergy": 0.5
}
```
- GUI: false - batch mode, true - interactive mode.
- RunMac: The .mac file defines how many particles you want to fire. Only valid for batch mode.
- Threads: 1 - serial mode, N - multithread mode with N threads (N>0).
- Foil: 0 - no foil, 1 - 304 stainless steel, 2 - Be.
- BeamEnergy: any non negative value. -->

### Run
```
$ ./exampleB1 ../config/config.json
```

### Analysis
The output data file is ``Stage#.root``. A simple analysis ROOT macro is provided as ``analyze.C``. 
