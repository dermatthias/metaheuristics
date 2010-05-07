# Metaheuristics for Clojure

Metaheuristics for Clojure is a Clojure library for common metaheuristic algorithms. 

## Usage

At the moment, there are only three (PSO, IWO and ES) algorithms implemented. They live in the namespace 'metaheuristics.pso', '.iwo' and '.es' respectively. There is only one public function to call for each variant, (iwo & params), (pso & params) and (es & params). They return a struct with the global best solution. The fitness value and the position ca be accessed by (:position gbest) and/or (:fitness gbest). 

## Known issues and limitations
For now, the functions can only handle real value solutions. Binary, decimal and permutations will be added later. 

## Installation
* If you want to work with the source:
  It's a [Leiningen](http://github.com/technomancy/leiningen) project. From here on I think you know what to do. 

* If you just want to use the library in your own project:
  Use the Leiningen (see above) to compile a .jar and put it to the path of your own project.

## License
Creative Commons - BY-NC-SA

Feel free to contact me. 
