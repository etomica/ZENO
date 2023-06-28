.. ZENO documentation master file, created by
   sphinx-quickstart on Wed May  2 17:31:44 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

ZENO
====

The ZENO software tool computes material, solution, and suspension properties for a specifed particle
shape or molecular structure using path-integral and Monte Carlo methods. These properties include:
capacitance, electric polarizability tensor, intrinsic conductivity, volume, gyration tensor, hydrodynamic
radius, intrinsic viscosity, friction coeffcient, diffusion coeffcient, sedimentation coeffcient, and related
quantities.

Run the code in a web browser
-----------------------------

The ZENO code can be run directly at https://zeno.nist.gov/zenoweb


How to get the code
-------------------

The source code can be found at https://github.com/usnistgov/ZENO

For users who prefer a graphical user interface, the code is also available at https://nanohub.org/resources/zeno

How to cite the code
--------------------

Derek Juba, Debra J. Audus, Michael Mascagni, Jack F. Douglas, and Walid Keyrouz. 
ZENO: Software for Calculating Hydrodynamic, Electrical, and Shape Properties of Polymer and Particle Suspensions
*Journal of Research of NIST*, 122:20, 2017. URL: https://doi.org/10.6028/jres.122.020
:download:`BibTeX <_static/ZENO.bib>` :download:`RIS <_static/ZENO.ris>`

Additionally, for virial-coefficient calculations:

Arpit Bansal, Andrew J. Schultz, Jack F. Douglas, and David A. Kofke. 
Probabilistic computations of virial coefficients of polymeric structures described by rigid configurations of spherical particles: A fundamental extension of the ZENO program
*Journal of Chemical Physics*, 157:224801, 2022. URL: https://doi.org/10.1063/5.0127465
:download:`BibTeX <_static/virial.bib>` :download:`RIS <_static/virial.ris>`


.. toctree::
    :maxdepth: 2
    :caption: Algorithms

    Calculations

.. toctree::
    :maxdepth: 2
    :caption: Compile the code

    Compilation

.. toctree::
    :maxdepth: 2
    :caption: Run the code

    RunCode
    Input
    Output
    Units

.. toctree::
    :maxdepth: 2
    :caption: Miscellaneous

    Contact
    Notes
    Validation
    License
    ZReferences
