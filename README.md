# [Siggraph Asia 2024] Optimized shock-protecting microstructures

<h1 align="center">
<a href="https://cims.nyu.edu/~zh1476/research/shock.html"><img alt="shock" src="https://cims.nyu.edu/~zh1476/assets/img/research/shock/website-teaser.png" width="80%"></a>
</h1><br>

A repository of the data and code used in our work, ["Optimized shock-protecting microstructures" [Huang et al. 2024]](https://cims.nyu.edu/~zh1476/research/shock.html).


Compilation
-----------

The code is tested on Linux and Mac (including Apple Silicon). To run the code:

    git clone --recurse-submodules git@github.com:Huangzizhou/ShockProtection.git

    # Compile the inflator code
    cd ShockProtection/inflator
    mkdir build
    cd build
    cmake ..
    make -j16

    # Compile the polyfem code
    cd ../..
    mkdir build
    cd build
    cmake ..
    make -j16

    # Unit test
    ./tests/unit_tests "isosurface-inflator-periodic"

The inflator code is used to map a vector of shape parameters to a periodic mesh of the microstructure cell, and compute the shape velocity of the shape parameters needed in the optimizations. The inflator code is called using terminal commands during the optimizations, since the cmake setup of CGAL (which is used by the inflator) may break the cmake of PolyFEM.

Usage
-------------

To run a shape optimization and obtain the microstructure cell shape that corresponds to a specific target stress:

    cd scripts/
    python optimize.py 5000 \
    ../inflator/data/patterns/2D/topologies/0105.obj \
    --no_tile --strain 0.2 --n_samples 3

where `5000` is the desired stress, the optimization will sample 3 strains from `10%` to `20%` and optimize the stress on those samples to match the desired stress.

All 105 topologies used in the paper are in folder `inflator/data/patterns/2D/topologies/0105.obj`, one can provide a new `.obj` edge mesh to run the optimization on a custom cell topology. The edge mesh has to fit into the unit cube and be periodic in both directions.

The above command creates a folder `result/0105_0.2_5000.0` and set up the JSON files needed to run PolyFEM. Once the optimization finishes successfully, one can further optimize the cell shape for wider flat region by

    python optimize.py 5000 \
    ../inflator/data/patterns/2D/topologies/0105.obj \
    --no_tile --strain 0.25 --n_samples 4 \
    --params ../result/0105_0.2_5000.0/sol.txt 

which adds one more sample strain, and uses the previous optimized shape parameters as a starting point.

Note that in the paper, we run simulations on a 2x2 periodic tile of the cell, which takes significantly more time to simulate but produces more accurate results. To simulate with 2x2 tiles in the optimization, remove `--no_tile` from above commands.

Documentation
-------------

The full documentation of PolyFEM can be found at [https://polyfem.github.io/](https://polyfem.github.io/)

License
-------

The code of PolyFEM itself is licensed under [MIT License](LICENSE). However, please be mindful of third-party libraries which are used by PolyFEM and may be available under a different license.

Citation
--------

If you use this code in your project, please consider citing our work:

```bibtex
@misc{huang2023shock,
      title={Optimized shock-protecting microstructures}, 
      author={Zizhou Huang and Daniele Panozzo and Denis Zorin},
      year={2023},
      eprint={2310.08609},
      archivePrefix={arXiv},
      primaryClass={math.OC},
      url={https://arxiv.org/abs/2310.08609}, 
}
```

```bibtex
@article{huang2024diffipc,
  author = {Huang, Zizhou and Tozoni, Davi Colli and Gjoka, Arvi and Ferguson, Zachary and Schneider, Teseo and Panozzo, Daniele and Zorin, Denis},
  title = {Differentiable solver for time-dependent deformation problems with contact},
  year = {2024},
  issue_date = {June 2024},
  publisher = {Association for Computing Machinery},
  address = {New York, NY, USA},
  volume = {43},
  number = {3},
  issn = {0730-0301},
  url = {https://doi.org/10.1145/3657648},
  doi = {10.1145/3657648},
  journal = {ACM Trans. Graph.},
  month = {may},
  articleno = {31},
  numpages = {30},
  keywords = {Differentiable simulation, finite element method, elastodynamics, frictional contact}
}
```
