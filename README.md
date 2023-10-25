StripMaker: Perception-driven Learned Vector Sketch Consolidation
=================================================================

<strong>Chenxi Liu<sup>1</sup>, Toshiki Aoki<sup>2</sup>, Mikhail Bessmeltsev<sup>3</sup>, Alla Sheffer<sup>1</sup></strong>

<small><sup>1</sup>University of British Columbia, <sup>2</sup>University of Tokyo, <sup>3</sup>Université de Montréal</small>

<p align="center">
	<img src="https://www.cs.ubc.ca/labs/imager/tr/2023/stripmaker/teaser.png" width="600"/>
</p>

This repository contains the source code for the paper
[StripMaker: Perception-driven Learned Vector Sketch Consolidation](https://www.cs.ubc.ca/labs/imager/tr/2023/stripmaker/).

 - `src` contains the source code for *libstripmaker*, a C++ library for consolidating raw vector sketches.
 - `stroke_strip_src` contains a modified version of [StrokeStrip](https://www.cs.ubc.ca/labs/imager/tr/2021/StrokeStrip/), a C++ library for parameterizing and fitting stroke strips. The modification increases the default sampling rate of input strokes and improves the parameterization robustness and fitting speed.
 - `python` contains Python scripts that call C++ binaries and run training and inference experiments.
 - `cmake` contains cmake files for fetching dependencies.
 - `external` contains extra dependencies.

The data can be downloaded from
[StripMakerData](https://www.cs.ubc.ca/labs/imager/tr/2023/stripmaker/StripMakerData.zip).
The `data` folder contains the pre-processed line drawings and annotations used to train the model, as well as the drawings used to generate the results in the paper.
Put `data` under the repository root for the following instructions.

Installation
--------

### Dependencies

All dependencies except for Gurobi are either included in `external` folder (Cornucopia, SketchConnectivity, glm) or automatically downloaded by cmake (cli11, eigen3, nlohmannjson, pybind11, span-lite, spdlog, tbb).

The Gurobi optimization library must be installed and licensed on your machine. If you are at a university, a [free academic license can be obtained](https://www.gurobi.com/downloads/end-user-license-agreement-academic/).
This project was built and tested with Gurobi 9.1 and 10.0; if you are using a newer version of Gurobi, update `FindGUROBI.cmake` to reference your installed version (e.g. append `gurobi101` to `NAMES gurobi gurobi91 gurobi100`).

### Building

To build `libstripmaker` and all its binary entrances with debug symbols and optimizations on, run

    mkdir build && cd build
    cmake ..
    cmake --build . --config RelWithDebInfo

This will produce a static library and multiple executable binaries.

### Python Dependencies

Python dependencies are installed inside a new conda environment by running

    conda env create -f python/environment.yml

To activate this environment, run

    conda activate stripmaker

Running Our Method
------------------

After building, make sure to add the `python` directory under the project root and the `python/core/` directory under the build directory to your `PYTHONPATH` environment variable, or else Python will not be able to import scripts of this project.

The results seen in our [supplementary materials](https://www.cs.ubc.ca/labs/imager/tr/2023/stripmaker/supplementary_materials.zip) can be generated with

    python3 ./python/launching/run_prediction.py ./data/inputs_test.yml --exe [build_folder]/clustering-solve,[build_folder]/clustering-endend --output ./snapshot/stripmaker_test.pdf --spiral

`--spiral` option turns on special fitting handling of strips containing spiral strokes.
(Consult the output of `python3 python/launching/run_prediction.py --help` for details on optional flags)

This will create a file `./snapshot/stripmaker_test.pdf` summarizing the results, and a `snapshot` directory containing intermediate and final outputs within individual result folders.

To run on a single input, replace the yml file with a single input scap file.

### Multithreading

The code is sequential by default.
It supports multithreading, which can be enabled by setting number of threads in [this line](https://github.com/squidrice21/sketch_clustering/blob/cab35b16f4c701e7bd2e784902d89d2c42c79d12/stroke_strip_src/Context.cpp#L7) and recompiling.

Training the Classifiers
------------------------

We provide pre-trained models in two forms: [scikit-learn models provided as a Python pickle](python/sketching/resources/classify_junction_models-current.pickle) in `models` folder in [StripMakerData](https://www.cs.ubc.ca/labs/imager/tr/2023/stripmaker/StripMakerData.zip), and an auto-generated C++ version of [the local classifier](src/classifier/forest.cpp) and [the global classifier](src/classifier/forest_sec.cpp).
The rest of this section is for people who wish to retrain the classifiers themselves.

Build the project and make sure to add the `python` directory under the project root and the `python/core/` directory under the build directory to your `PYTHONPATH` environment variable.
Run

    python3 ./python/launching/run_batch.py ./data/inputs_train+validate.yml -e [build_folder]/clustering-secondary-sample,[build_folder]/clustering-secondary-feature -s ./snapshots/local_classifier/ --save

to train the local classifier.
The `--save` option tells the script to convert the resulting `./snapshots/local_classifier/model.sav` and `./snapshots/local_classifier/sc_model.sav` to `src/classifier/forest.cpp`; it requests confirmation before overwriting.

Run

    python3 ./python/launching/run_batch.py ./data/inputs_train+validate.yml -e [build_folder]/clustering-secondary-sample,[build_folder]/clustering-secondary-feature -s ./snapshots/global_classifier/ --secondary --save

to train the global classifier.
The `--save` option tells the script to convert the resulting `./snapshots/global_classifier/model_sec.sav` and `./snapshots/global_classifier/sc_model_sec.sav` to `src/classifier/forest_sec.cpp`; it requests confirmation before overwriting.

The training script also automatically runs leave-one-out cross-validation.
The statistics are saved in the `./snapshots/local_classifier/cv/` or `./snapshots/global_classifier/cv/` folder.

Data
----

The `data` directory contains pre-processed inputs (line drawings) in the form of SCAP files, which encode strokes as polylines.
Strokes are contained in pairs of braces `{ ... }`.
Each stroke has a unique stroke ID and a strip ID shared by all strokes that collectively make up one intended curve.
```
#[width]	[height]
@[thickness]
{
	#[stroke_id]	[strip_id]
	[x1]	[y1]	[time1]
	[x2]	[y2]	[time2]
	[x3]	[y3]	[time3]
	[...etc]
}
[...etc]
```
The strip IDs are ignored in the input files and overwritten by StripMaker in the output files.
The timestamps per polyline samples are omitted; only the drawing order indicated by the stroke IDs is considered by StripMaker.

The SCAP files can be visualized by converting them to SVG files using `python/visualization/scap_to_svg.py`.
New SCAP files can be generated from SVG files using the `python/visualization/svg_to_scap_width_exact.py`.

License
-------

The source code (everything under `src` and `python`) is licensed under [Version 2.0 of the Apache License](LICENSE).
The source code under `stroke_strip_src` is a modification of [StrokeStrip](https://www.cs.ubc.ca/labs/imager/tr/2021/StrokeStrip/), licensed under [MIT](stroke_strip_src/LICENSE.txt).
The drawings (which need to be downloaded separately from [StripMakerData](https://www.cs.ubc.ca/labs/imager/tr/2023/stripmaker/StripMakerData.zip)) are licensed under separate licenses.
Please refer to `data/inputs_train+validate.yml` and `data/inputs_test.yml` for license information for each drawing.

BibTeX
------

```
@article{stripmaker,
  title = {StripMaker: Perception-Driven Learned Vector Sketch Consolidation},
  author = {Liu, Chenxi and Aoki, Toshiki and Bessmeltsev, Mikhail and Sheffer, Alla},
  year = {2023},
  issue_date = {August 2023},
  publisher = {Association for Computing Machinery},
  address = {New York, NY, USA},
  volume = {42},
  number = {4},
  journal = {ACM Trans. Graph.},
  month = {jul},
  articleno = {55},
  numpages = {15},
  issn = {0730-0301},
  url = {https://doi.org/10.1145/3592130},
  doi = {10.1145/3592130}
}
```
