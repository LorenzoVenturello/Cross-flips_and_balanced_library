# Cross-flips_and_balanced_library

A SageMath implementation of cross-flips and database of balanced triangulations of manifolds.

Try it on [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/LorenzoVenturello/Cross-flips_and_balanced_library/master?filepath=demo.ipynb).

The goal of this code is to take as input a balanced triangulation of a combinatorial manifold with a high number of vertices, and return a balanced triangulation of the same manifold with fewer vertices. This is obtained by performing local moves, called "cross-flips" or flips for short, which preserve both the coloring and the PL-homeomorphism type. 

Due to the complexity of the flips, the software can currently deal with 2- and 3- dimensional combinatorial manifolds. The program makes two choices: which combinatorial type of move to pick and where (in the triangulation) to apply it. Already in dimension 3 it is not, in general, possible to find a sequence of flips which decreases the number of vertices at each step. When no such move is applicable, we simply produce a random sequence of flips increasing the number of vertices and then start again with the reduction.  

An extended abstract containing more details on the implementation and the results obtained appeared in the Séminaire Lotharingien de Combinatoire (82B), in [Proceedings of the FPSAC 2019](http://fpsac2019.fmf.uni-lj.si/resources/Proceedings/201.pdf). See also the [slides](http://fpsac2019.fmf.uni-lj.si/resources/Slides/203slides.pdf) of the software presentation delivered by Alexander Wang at [FPSAC 2019](http://fpsac2019.fmf.uni-lj.si/).


