(ns metaheuristics.ga
  (:use [clojure.set :only (intersection difference)])
  (:use [clojure.contrib.math])
  (:use metaheuristics.testfunctions))

; helper functions
; ----------------
(defn- scramble
  [l]
  (let [items (java.util.ArrayList. l)]
	scramble (do (java.util.Collections/shuffle items) (seq items))))

(defn- normal [mu sigma]
  (let [r (new java.util.Random)]
    (+ mu (* sigma (.nextGaussian r)))))

(defn- garand-int [max number]
  (map (fn [_] (int (* max (rand)))) (range number)))

(defn- euclidean
  [v1 v2]
  (Math/sqrt (areduce v1 i ret 0
	   (+ ret (Math/pow (- (aget v1 i) (aget v2 i)) 2)))))

; genetic algorithm
; -----------------
(defstruct individual :tag :chromosome :steps :fitness)
(defstruct population :poplist)

(defn- generate-chromosome
  [n bits]
  (let [max (int (- (Math/pow 2 bits) 1))]
    (garand-int max n)))
;(generate-chromosome 22 8)

(defn- init-individual
  [n bits]
  (struct individual 0 (generate-chromosome n bits)
	  (map #(* 10 %) (take 22 (repeatedly rand))) 0))

(defn- init-population
  [n dim bits]
  (let [poplist (map (fn [_] (init-individual dim bits)) (range n))]
    (struct population poplist)))

(defn chromo-to-phenotype
  [chromosome]
  (let [sum (reduce + chromosome)]
    (double-array (for [gene chromosome] (/ gene sum)))))

(defn- share
  [dist sigma alpha]
  (if (<= dist sigma)
    (- 1 (Math/pow (/ dist sigma) alpha))
    0))

(defn- fitness-sharing 
  [ind popu]
  (let [sigma 100 ;; 100 seems to be fine
	alpha 1
	dist-sum (reduce + (for [other (:poplist popu)]
		   (share (euclidean
			   (double-array (:chromosome ind)) 
			   (double-array (:chromosome other)))
			  sigma 
			  alpha)))]
    (assoc ind :fitness (/ (:fitness ind) dist-sum))))
;; (def foo (init-population 10 22 8))
;; (fitness-sharing (init-individual 22 8) foo)

(defn- evaluate-individual
  [ind fitness]
  (let [w-pheno (chromo-to-phenotype (:chromosome ind))	
	fitval (fitness w-pheno)]
    (assoc ind :fitness fitval)))

(defn- evaluate-all
  [popu fitness fs?]
  (let [agentlist (for [ind (:poplist popu) :when (= (:tag ind) 1)] (agent ind))]

    (dorun (map #(send %1 evaluate-individual fitness) agentlist))
    (apply await agentlist)

    ;; fitness sharing
    (if fs? (do 
	      (dorun (map #(send %1 fitness-sharing popu) agentlist))
	      (apply await agentlist)))

    (let [children (for [agent agentlist] @agent)
	  parents (for [p (:poplist popu) :when (= (:tag p) 0)] p)]
      (assoc popu :poplist (concat parents children)))))

(defn- evaluate-all-firstrun
  [popu fitness]
  (let [agentlist (for [ind (:poplist popu)] (agent ind))]
    (dorun (map #(send %1 evaluate-individual fitness) agentlist))
    (apply await agentlist)
    (assoc popu :poplist (for [agent agentlist] @agent))))

(defn- gene-crossover
  [gene1 gene2 len]
  (let [mask (int (- (Math/pow 2 len) 1))
	last-g1 (bit-and gene1 mask)
	last-g2 (bit-and gene2 mask)
	new-g1 (bit-or (bit-shift-left (bit-shift-right gene1 len) len) last-g2)
	new-g2 (bit-or (bit-shift-left (bit-shift-right gene2 len) len) last-g1)]
    (list new-g1 new-g2)))
;(gene-crossover 8 3 6)

(defn- crossover
  [chromo1 chromo2 bits]
  (let [split-pos (nth (range 1 bits) (rand-int (- bits 1)))
	split-length (- bits split-pos)
	result (map #(gene-crossover %1 %2 split-length) chromo1 chromo2)
	c1-new (map #(first %1) result)
	c2-new (map #(second %1) result)]
    (list c1-new c2-new)))
;(crossover '(255 12 18 238 210 199 88) '(4 8 15 16 23 42 108) 8)

(defn- mutation
  [chromosome bits percentage]
  (for [gene chromosome]
    (loop [g gene b 0]
      (if (= b bits)
	g
	(let [new-g (if (< (rand) 0.25) (bit-flip g b) g)]
	  (recur new-g (inc b)))))))
;;(mutation '(34 21 57 56 11 10) 8)

(defn- do-offspring
  [acc parents-pair]
  (let [percentage 0.3
	p1c (:chromosome (first parents-pair))
	p2c (:chromosome (second parents-pair))
	children (crossover p1c p2c 8)
	child1 {:tag 1 :chromosome (mutation (first children) 8 percentage) :steps (list) :fitness 0}
	child2 {:tag 1 :chromosome (mutation (second children) 8 percentage) :steps (list)  :fitness 0}
	old-pop-list (:poplist acc)]
    (assoc acc :poplist (conj old-pop-list child1 child2))))

(defn- adapted-crossover
  [p1 p2]
  (let [dim (count (:chromosome p1))
	dim2 (/ dim 2)
	p1ch (:chromosome p1)
	p1st (:steps p1)
	p2ch (:chromosome p2)
	p2st (:steps p2)
	childc (map (fn [c1 c2] (if (= (int (* 2 (rand))) 0) c1 c2)) p1ch p2ch)
	childs (map (fn [s1 s2] (if (= (int (* 2 (rand))) 0) s1 s2)) p1st p2st)]
    (list childc childs)))

(defn- adapted-mutation
  [chromosome steps]
  (let [dim 22
	bound 0.5
	tau (/ 1 (sqrt (* 2 (sqrt dim))))
	tauprime-rand (* (/ 1 (sqrt (* 2 dim))) (normal 0 1))
	new-steps-raw (for [s steps]
		    (* s (expt (Math/E) (+ tauprime-rand (* tau (normal 0 1))))))
	new-steps (map (fn [s] (if (< s bound) bound s)) new-steps-raw)
	new-pos (map (fn [c s] (+ c (* s (normal 0 1))))
		     chromosome new-steps)]
    (list new-pos new-steps)))

(defn- do-offspring-adapted
  [acc parents-pair]
  (let [p1 (first parents-pair)
	p2 (second parents-pair)
	crossed (adapted-crossover p1 p2)
	mutated (adapted-mutation (first crossed) (second crossed))
	child {:tag 1 :chromosome (first mutated)  :steps (second mutated) :fitness 0}
	old-pop-list (:poplist acc)]
    (assoc acc :poplist (conj old-pop-list child))))

(defn- generate-offspring
  [acc parents adapted?]
  (if adapted?
    (reduce do-offspring-adapted acc parents)
    (reduce do-offspring acc parents)))

(defn- survivor-selection
  [popu ftype percent popsize]
  (let [sorted-all (sort-by :fitness ftype (:poplist popu))
	size-all (count sorted-all)
	size-top (* percent size-all)
	size-rest (- popsize size-top)
	splitted (split-at size-top sorted-all)
	top-percent-members (nth splitted 0)
	rest-offspring (take size-rest (sort-by :fitness ftype
				  (for [i (nth splitted 1) :when (if (= (:tag i) 1) true)] i)))]
    (assoc popu :poplist (concat top-percent-members rest-offspring))))

(defn- parent-selection
  [popu ftype percent popsize adapted?]
  (let [sorted (sort-by :fitness ftype (:poplist popu))
	size-parents (int (* percent popsize))
	parents (take size-parents sorted)

	ext-parents (if adapted?
		     (scramble (take (* 2 popsize) (cycle parents)))
		     (scramble (take popsize (cycle parents))))
	splitted (if adapted? 
		   (split-at popsize ext-parents)
		   (split-at (/ popsize 2) ext-parents))]
    (map #(list %1 %2) (nth splitted 0) (nth splitted 1))))

(defn ga
  "Starts the GA algorithm.
  Parameter: 
  fitness - defines the fitness function used. The only argument is the position of the particle (a double-vector).
  ftype - defines the type of optimization problem (minimize (<) or maximize (>)).
  dim - number of dimensions in solution.
  popsize - size of the population.
  par-perc - percentage of parents chosen for parent selection.
  surv-perc - percentage of population chosen for next genertion.
  adpated? - use mutation step size adaption (boolean).
  max-iterations - maximum number of iterations.
  "  
  [fitness ftype dim popsize par-perc surv-perc adapted? max-iterations]
  (let [popu (init-population popsize 22 8)
	popu-evaluated (evaluate-all-firstrun popu fitness)]

    (loop [runs max-iterations popu popu-evaluated]
      (if DEBUG? (println "best:" (:fitness (first (sort-by :fitness ftype 
							    (:poplist popu))))))
      (if (zero? runs)
	(first (sort-by :fitness ftype (:poplist popu)))
	(let [parents (parent-selection popu ftype par-perc popsize adapted?)
	      popu-with-children (generate-offspring popu parents adapted?)
	      popu-eva (evaluate-all popu-with-children fitness adapted?)
	      popu-surv (survivor-selection popu-eva ftype surv-perc popsize)
	      new-pop (assoc popu-surv :poplist
			     (for [i (:poplist popu-surv)] (assoc i :tag 0)))]
	  (if DEBUG? (do (println "no. parent pairs:" (count parents))
			 (println "no. parents + children:" (count (:poplist popu-with-children)))
			 (println "no. survivors: " (count (:poplist popu-surv)))))
	  (recur (dec runs) new-pop))))))

(def DEBUG? false)

;; usage
;; (seq (chromo-to-phenotype (:chromosome (ga griewank < 4 30 0.4 0.4 false 50))))