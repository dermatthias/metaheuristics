(ns metaheuristics.iwo
  (:use [clojure.set :only (intersection difference)])
  (:use [clojure.contrib.math])
  (:use testfunctions))

;; rand from normal distribution
(defn- normal [mu sigma]
  (let [r (new java.util.Random)]
    (+ mu (* sigma (.nextGaussian r)))))

; invasive weed optimization
; --------------------------

(defstruct plant :seedlist :position :pfit :tfit)
(defstruct population :plantlist :gbest :gworst :minseed :maxseed :sigmaInit :sigmaFinal)

(defn- init-plant [nDimensions max]
  (agent (struct plant
		 (list)
		 (double-array (for [i (range nDimensions)] (rand max)))
		 (double 0.0)
		 (double 0.0))))

(defn- init-population [nPlants nDimensions nSeedMin nSeedMax sigmaInit sigmaFinal max]
  (let [plants (map (fn [_] (init-plant nDimensions max)) (range nPlants))
	gbest (init-plant nDimensions max)
	gworst (init-plant nDimensions max)]
    (agent (struct population plants gbest gworst (int nSeedMin) (int nSeedMax) sigmaInit sigmaFinal))))

(defn- set-bestworst [population]
  (let [plants (for [p (:plantlist population)] @p)
	sorted (sort-by :pfit < plants)]
    (assoc population :gbest (agent (first sorted)) :gworst (agent (last sorted)))))

;; eval seeds
(defn- eval-seed [seed fitness]
  (let [#^doubles pos (:position seed)
	fit (fitness pos)]
    (assoc seed :pfit fit)))

(defn- eval-plantseeds [plant fitness]
  (let [seeds (:seedlist @plant)]
    (doseq [s seeds] (send s eval-seed fitness))
    (apply await seeds)))

;; create new seeds
(defn- create-new-seed [pos sigma]
  (let [pos-offset (double-array (map (fn [_] (normal 0 sigma)) (range 22)))
	newpos (amap pos i ret
		     (+ (aget pos i) (aget pos-offset i)))]
    (agent (struct plant (list) newpos 0.0))))

(defn- generate-seeds [plant population maxIt modulation iteration]
  "Generates new seeds"
  (let [pfit (double (:pfit plant))
	#^doubles pos (:position plant)
	gbestfit (double (:pfit @(:gbest population)))
	gworstfit (double (:pfit @(:gworst population)))
	minSeed (int (:minseed population))
	maxSeed (int (:maxseed population))
	sigmaInit (double (:sigmaInit population))
	sigmaFinal (double (:sigmaFinal population))

	nSeeds (int (+ (* (/ (- pfit gworstfit) (- gbestfit gworstfit)) maxSeed)
		       (* (/ (- pfit gbestfit) (- gworstfit gbestfit)) minSeed)))	
	sigma (+ (* (/ (Math/pow (- maxIt iteration) modulation)
		       (Math/pow maxIt modulation))
		    (- sigmaInit sigmaFinal))
		 sigmaFinal)
	newSeeds (map (fn [_] (create-new-seed pos sigma)) (range nSeeds))]
    (assoc plant :seedlist newSeeds)))

(defn- competition [population limit style]
  (let [seeds (for [p (:plantlist population) s (:seedlist @p)] @s)
	plants (for [p (:plantlist population)] @p)       
	;; best up to limit
	newpop (take limit (sort-by :pfit < (concat seeds plants)))
	newagents (for [np newpop] (agent np))]
    (assoc population
      :plantlist newagents :gbest (agent (first newpop)) :gworst (agent (last newpop)))))

(defn- grow [fitness population maxIt nPlantsMax modulation iteration]
  ;; generate new seeds for each plant-agent
  (if DEBUG? (println "creating new seeds..."))
  (dorun (map #(send % generate-seeds @population maxIt modulation iteration) 
	      (:plantlist @population)))
  (apply await (:plantlist @population))
  (if DEBUG? (println "no. of seeds:" (count (for [p (:plantlist @population) seeds (:seedlist @p)]
				    seeds))))

  ;; update seed fitness
  (if DEBUG? (println "updating seed fitness..."))
  (dorun (pmap #(eval-plantseeds % fitness) (:plantlist @population)))

  ;; competitive exclusion
  (if DEBUG? (println "competition..."))
  (send population competition nPlantsMax 1)
  (await population)
  (if DEBUG? (println "best in generation" iteration":" (:pfit @(:gbest @population)) "\n---"))
  )

(defn iwo
  "Starts the IWO algorithm. Algorithm based on http://dx.doi.org/10.1016/j.ecoinf.2006.07.003
   Parameters:
   fitness - the fitness function to be used. only one parameter: the position of the particle (double-array).
   dim - number of dimensions in solution.
   nplants - number of initial plants.
   nplants-max - maximum number of plants in population.
   seed-min - minimum number of seeds generated per plant.
   seed-max - maximum number of seeds generated per plant.
   sigma-init - inital value for sigma (standard deviation).
   sigma-final - final value for sigma (standard deviation).
   modulation -  modulation index (usually 3). 
   max-iterations - maximum number of iterations.
   max-feat - maximum value for one feature.
   "
  [fitness dim nplants nplants-max seed-min seed-max sigma-init sigma-final modulation max-iterations max-feat]
  (let [population (init-population nplants dim seed-min seed-max sigma-init sigma-final max-feat)]
    ; one time eval of initial plants
    (dorun (map #(send % eval-seed fitness) (:plantlist @population)))
    (apply await (:plantlist @population))
    (send population set-bestworst)
    (await population)
    ; start IWO
    (dorun (map (fn [i] (grow fitness population max-iterations nplants-max modulation i))
		(range max-iterations)))
    ; return best solution
    @(:gbest @population)))

(def DEBUG? false)

;; testing
;; (seq (:position (iwo griewank 10 10 30 0 5 10 0.1 3 100 10.0)))