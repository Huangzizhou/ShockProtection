# [Siggraph Asia 2024] Optimized shock-protecting microstructures

<h1 align="center">
<a href="https://cims.nyu.edu/~zh1476/research/shock.html"><img alt="shock" src="https://cims.nyu.edu/~zh1476/assets/img/research/shock/website-teaser.png" width="80%"></a>
</h1><br>

A repository of the data and code used in our work, ["Optimized shock-protecting microstructures" [Huang et al. 2024]](https://cims.nyu.edu/~zh1476/research/shock.html).


Compilation
-----------

The code is only tested on Linux. To run the code:

    git clone --recurse-submodules git@github.com:Huangzizhou/ShockProtection.git

    # Compile the inflator code
    cd inflator
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
