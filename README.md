# Metaheuristics for Clojure

Metaheuristics for Clojure is a Clojure library for common (and some
uncommon) metaheuristic algorithms. It uses Clojure Agents to parallel
compute the fitness functions, so there will be a distinct speedup on
multi-core machines. It was tested on a *dual* 4-core Intel Xeon
(W5580), where it ran 6.3 times faster than on a single core machine.

These algorithms are implemented at the moment:
* Particle Swarm Optimization (PSO), including some extensions (velocity clamping, inertia weights, multi-start resetting)
* Genetic Algorithm (GA)
* Evolution Strategies (ES)
* Invasive Weed Optimization (IWO)

For more information about these algorithms and techniques, I
recommend these sources:

* [Introduction to Evolutionary Computing - Eiben/Smith]
  (http://www.cs.vu.nl/~gusz/ecbook/ecbook.html)
* [Fundamentals of Computational Swarm Intelligence - Engelbrecht](http://si.cs.up.ac.za/)
* [A novel numerical optimization algorithm inspired from weed
colonization - Mehrabiana et al. (paper)](http://dx.doi.org/10.1016/j.ecoinf.2006.07.003)

## Usage

The algorithms live in the namespace 'metaheuristics.*', where *
stands for one of the abbreviations from above. There is only one
public function to call for each variant (same abbrv.). They all
return a struct with the global best solution. The fitness value and
the position ca be accessed by (:position gbest) and/or (:fitness
gbest). See the code for more information. It's not that complicated.

## Known issues and limitations

For now, the functions can only handle real value solutions. Binary,
decimal and permutations will be added later.

Also, I am by no means a Clojure expert. Bugs, bad design...it's
possible that you will find all of that in my code. If not, phew!! :)

## Installation

* If you want to work with the source: It's a
  [Leiningen](http://github.com/technomancy/leiningen) project. From
  here on I think you know what to do.

* If you just want to use the library in your own project: Use the
  Leiningen (see above) to compile a .jar and put it to the path of
  your own project.

## License
Creative Commons - BY-NC-SA

Feel free to contact me for any questions.


