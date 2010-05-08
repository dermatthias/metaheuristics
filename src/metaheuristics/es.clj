(ns metaheuristics.es
  (:use [clojure.set :only (intersection difference)])
  (:use [clojure.contrib.math])
  (:use metaheuristics.testfunctions))

; helper functions
; ----------------
(defn- scramble
  [l]
  (let [items (java.util.ArrayList. l)]
	scramble (do (java.util.Collections/shuffle items) (seq items))))

;; rand from normal distribution
(defn- normal [mu sigma]
  (let [r (new java.util.Random)]
    (+ mu (* sigma (.nextGaussian r)))))
(defn- esrand-int [max number]
  (map (fn [_] (* max (rand))) (range number)))

; genetic algorithm
; -----------------
(defstruct individual :tag :chromosome :steps :fitness)
(defstruct population :poplist)

(defn- generate-chromosome
  [n bits]
  (let [max (int (- (Math/pow 2 bits) 1))]
    (double-array (esrand-int max n))))
;(generate-chromosome 22 8)

(defn- init-individual
  [n bits]
  (struct individual 0 (generate-chromosome n bits) 
	  (double-array (map #(* 10 %) (take 22 (repeatedly rand)))) 
	  0))
;(init-individual 22 8)

(defn- init-population
  [n dim bits]
  (let [poplist (map (fn [_] (init-individual dim bits)) (range n))]
    (struct population poplist)))

(defn chromo-to-phenotype
  [chromosome]
  (let [sum (areduce chromosome i ret 0
		     (+ ret (aget chromosome i)))]
    (amap chromosome i ret
	  (/ (aget chromosome i) sum))))
;(chromo-to-phenotype (generate-chromosome 22 8) 0.2)

(defn- evaluate-individual
  [ind fitness]
  (let [w-pheno (chromo-to-phenotype (:chromosome ind))
	fitval (fitness w-pheno)]
    (assoc ind :fitness fitval)))

(defn- evaluate-all
  [popu fitness]
  (if DEBUG? (println "Evaluating population..."))
  (let [agentlist (for [ind (:poplist popu) :when (= (:tag ind) 1)] (agent ind))]
    
    (dorun (map #(send %1 evaluate-individual fitness) agentlist))
    (apply await agentlist)

    (let [children (for [agent agentlist] @agent)
	  parents (for [p (:poplist popu) :when (= (:tag p) 0)] p)]
      (assoc popu :poplist (concat parents children)))))

(defn- evaluate-all-firstrun
  [popu fitness]
  (let [agentlist (for [ind (:poplist popu)] (agent ind))]
    (dorun (map #(send %1 evaluate-individual fitness) agentlist))
    (apply await agentlist)
    (assoc popu :poplist (for [agent agentlist] @agent))))

(defn- crossover-inter
  "Intermediate recombination"
  [p1 p2]
  (let [dim (count (:chromosome p1))
	p1ch (:chromosome p1)
	p1st (:steps p1)
	p2ch (:chromosome p2)
	p2st (:steps p2)
	childc (amap p1ch i ret
		    (/ (+ (aget p1ch i) (aget p2ch i)) 2.0))
	childs (amap p1st i ret
		    (/ (+ (aget p1st i) (aget p2st i)) 2.0))]
    (list childc childs)))
;; (crossover (init-individual 22 8) (init-individual 22 8))

(defn- crossover-discrete
  "Discrete recombination"
  [p1 p2]
  (let [dim (count (:chromosome p1))
	p1ch (:chromosome p1)
	p1st (:steps p1)
	p2ch (:chromosome p2)
	p2st (:steps p2)
	childc (amap p1ch i ret
		     (if (= (int (* 2 (rand))) 0) (aget p1ch i) (aget p2ch i)))
	childs (amap p1st i ret
		     (if (= (int (* 2 (rand))) 0) (aget p1st i) (aget p2st i)))]
    (list childc childs)))
;; (crossover (init-individual 22 8) (init-individual 22 8))


(defn- adapted-mutation
  "Uncorrelated Mutation with n Step Sizes"
  [chromosome steps]
  (let [dim (count chromosome)
	bound 0.5
	tau (/ 1 (sqrt (* 2 (sqrt dim))))
	tauprime-rand (* (/ 1 (sqrt (* 2 dim))) (normal 0 1))
	new-steps-raw (amap steps i ret						
			    (* (aget steps i) 
			       (expt (Math/E) (+ tauprime-rand (* tau (normal 0 1))))))
	new-steps (amap new-steps-raw i ret
			(if (< (aget new-steps-raw i) bound) bound (aget new-steps-raw i)))
	new-pos (amap chromosome i ret
		      (+ (aget chromosome i) (* (aget new-steps i) (normal 0 1))))]
    (list new-pos new-steps)))
;; (def chromo (generate-chromosome 22 8))
;; (seq (second (adapted-mutation chromo (double-array (repeat 22 1.0)))))

(defn- do-offspring
  [acc parents-pair]
  (let [cross (crossover-inter (first parents-pair) (second parents-pair))
	mutated (adapted-mutation (first cross) (second cross))
	child {:tag 1 :chromosome (first mutated) :steps (second mutated) :fitness 0}
	old-pop-list (:poplist acc)]
    (assoc acc :poplist (conj old-pop-list child))))

(defn- generate-offspring
  [acc parents]
  (if DEBUG? (println "Generating Offspring..."))
  (reduce do-offspring acc parents))

(defn- survivor-selection
  [popu ftype popsize]
  (if DEBUG? (println "Survivor Selection..."))
  (let [survivors (take popsize (sort-by :fitness ftype 
					 (for [i (:poplist popu) :when (= (:tag i) 1)] i)))]
    (assoc popu :poplist survivors)))

(defn- parent-selection
  [popu popsize factor]
  (if DEBUG? (println "Parent Selection..."))
  (let [parents (scramble (take (* factor popsize) (cycle (:poplist popu))))
	splitted (split-at (* (/ factor 2) popsize) parents)]
    (map #(list %1 %2) (nth splitted 0) (nth splitted 1))))
;(def foo (init-population 10 22 8))
;(count (parent-selection foo 0.4 10))

(defn es
  "Starts the ES algorithm.
  Parameter: 
  fitness - defines the fitness function used. The only argument is the position of the particle (a double-vector).
  ftype - defines the type of optimization problem (minimize (<) or maximize (>)).
  dim - number of dimensions in solution.
  popsize - size of the population.
  offsp-factor - factor for the number of offspring created (usually 2).
  max-iterations - maximum number of iterations.
  "
  [fitness ftype dim popsize offsp-factor max-iterations]
  (let [popu (init-population popsize dim 8)
	popu-evaluated (evaluate-all-firstrun popu fitness)]

    (loop [runs max-iterations popu popu-evaluated]
      (if DEBUG? (println "best in generation:" 
			  (:fitness (first (sort-by :fitness ftype 
						    (:poplist popu))))))
      (if (zero? runs)
	(first (sort-by :fitness ftype (:poplist popu)))
	(let [parents (parent-selection popu popsize offsp-factor)
	      popu-with-children (generate-offspring popu parents)
	      popu-eva (evaluate-all popu-with-children fitness)
	      popu-surv (survivor-selection popu-eva ftype popsize)
	      new-pop (assoc popu-surv :poplist
			     (for [i (:poplist popu-surv)] (assoc i :tag 0)))]
	  (if DEBUG? (do (println "no. parent pairs:" (count parents))
			 (println "no. parents + children:" 
				  (count (:poplist popu-with-children)))
			 (println "no. survivors: " 
				  (count (:poplist popu-surv)))))
	  (recur (dec runs) new-pop))))))

(def DEBUG? false)

;; (seq (chromo-to-phenotype (:chromosome (es griewank < 10 30 2 100))))