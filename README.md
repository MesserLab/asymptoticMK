# asymptoticMK
A tool for computing alpha, the fraction of nucleotide substitutions in a given genomic region that were driven to fixation by positive selection.

**Citation:** B.C. Haller, P.W. Messer. (2017). asymptoticMK: A web-based tool for the asymptotic McDonald–Kreitman test. *G3: Genes, Genomes, Genetics 7*(5), 1569–1575.

---------------------
#### A web version of this software is available at [http://benhaller.com/messerlab/asymptoticMK.html](http://benhaller.com/messerlab/asymptoticMK.html).  This source repository is intended for users who wish to run asymptoticMK on their local machine, or wish to modify asymptoticMK for their own use.
---------------------
 
Files
-------
* asymptoticMK_local.R : an R script for running asymptoticMK locally
* asymptoticMK_run.html.R : the web-based version, using FastRWeb
* asymptoticMK.html : the HTML splash page for the web-based version
* LICENSE : the full GNU GPL v.3 license for asymptoticMK
* README.md : this readme file
* sample_polymorphism_levels.txt : a sample Drosophila dataset (use with d=59570, d0=159058)
* simulate_alpha.slim : a [SLIM](https://messerlab.org/slim/) model to generate simulated data for testing asymptoticMK

License
----------

Copyright (c) 2017 Philipp Messer.  All rights reserved.

asymptoticMK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

asymptoticMK is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with asymptoticMK.  If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).


Development & Feedback
-----------------------------------
If you have feedback or feature requests, or if you are interested in contributing to asymptoticMK, please contact Philipp Messer at [messer@cornell.edu](mailto:messer@cornell.edu). Please note that Philipp is also looking for graduate students and postdocs.