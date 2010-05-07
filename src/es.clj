(ns es
  (:use [clojure.set :only (intersection difference)])
  (:use [clojure.contrib.math])
  (:use test))

; helper functions
; ----------------
(defn replace-item
  "Gibt eine Liste zurück mit dem n-ten Element von l durch v ersetzt."
  [l n v]
  (concat (take n l) (list v) (drop (inc n) l)))

(defn bit-get
  "Liefert bit bei Index n"
  [x n]
  (if (true? (bit-test x n)) 1 0))

(defn scramble
  "Vermischt die übergebene Liste zufällig."
  [l]
  (let [items (java.util.ArrayList. l)]
	scramble (do (java.util.Collections/shuffle items) (seq items))))

(defn now []
  (java.sql.Timestamp. (.getTime (java.util.Date.))))
(defn timestamp []
  (.getTime (java.util.Date.)))

(defn rand-int [max number]
  (map (fn [_] (* max (rand))) (range number)))
;(rand-int 10 1)

;; rand from normal distribution
(defn normal [mu sigma]
  (let [r (new java.util.Random)]
    (+ mu (* sigma (.nextGaussian r)))))


; genetic algorithm
; -----------------
(defstruct individual :tag :chromosome :steps :fitness)
(defstruct population :poplist)

(defn generate-chromosome
  "Generiert ein zufälliges Chromosom der Länge n mit je 'bits' Bits pro Gen."
  [n bits]
  (let [max (int (- (Math/pow 2 bits) 1))]
    (double-array (rand-int max n))))
;(generate-chromosome 22 8)

(defn init-individual
  "Erzeugt ein Individuum. Tag 0 für Eltern, Tag 1 für Kinder."
  [n bits]
  (struct individual 0 (generate-chromosome n bits) 
	  (double-array (map #(* 10 %) (take 22 (repeatedly rand)))) 
	  0))
;(init-individual 22 8)

(defn init-population
  "Generiert die initiale Population mit n Individuen und dim Dimensionen."
  [n dim bits]
  (let [poplist (map (fn [_] (init-individual dim bits)) (range n))]
    (struct population poplist)))

(defn chromo-to-phenotype
  "Wandelt das Chromosom in seine Phenotyprepräsentation um. Genotyp: Werte zwischen 0 und 2^bits-1. Phenotyp die entsprechend normalisierten Werte zwischen 0 und 1, in der Summe gleich 1. Die 18 Genres werden zuvor einem Multiplikator reduziert."
  [chromosome gfactor]
  (let [genre-reduced (amap chromosome i ret
		     (if (> i 3)
		       (* gfactor (aget chromosome i))
		       (aget chromosome i)))
	sum (areduce genre-reduced i ret 0
		     (+ ret (aget genre-reduced i)))]
    (amap genre-reduced i ret
	  (/ (aget genre-reduced i) sum))))
;(chromo-to-phenotype (generate-chromosome 22 8) 0.2)

(defn evaluate-individual
  "Evaluiert ein übergebenes Individuum. Hier geschieht die eigentliche Arbeit, die Funktionen generate-neighbourhood und fitness-score lässt den kompletten Recommender für das aktuelle Ind. durchlaufen."
  [ind active others data gfactor]
  (let [w-geno (:chromosome ind)
	w-pheno (chromo-to-phenotype w-geno gfactor)

	allusers (generate-neighbourhood active others w-pheno data)
	mean (/ (reduce + (for [n allusers] (second n))) (- (:users data-info) 1))
	;;upperhalf (+ mean (/ (- (second (first allusers)) mean) 2.0))
	nhood (remove nil? (for [n allusers] (if (> (second n) mean) n)))
	fitval (fitness-score active nhood (keys (data active)) data data)]
    (assoc ind :chromosome w-geno :fitness (first fitval))))

(defn evaluate-all
  "Evaluiert die neu erzeugten Kinder der Population, jeweils per Aufruf von evaluate-individual."
  [active others popu data gfactor]
    (println "Evaluating population...")
  (let [agentlist (for [ind (:poplist popu) :when (= (:tag ind) 1)] (agent ind))]

    (dorun (map #(send %1 evaluate-individual active others data gfactor) agentlist))
    (apply await agentlist)

    (let [children (for [agent agentlist] @agent)
	  parents (for [p (:poplist popu) :when (= (:tag p) 0)] p)]
      (assoc popu :poplist (concat parents children)))))

(defn evaluate-all-firstrun
  "Evaluiert alle Individuen der Population, jeweils per Aufruf von evaluate-individual"
  [active others popu data gfactor]
  (let [agentlist (for [ind (:poplist popu)] (agent ind))]
    (dorun (map #(send %1 evaluate-individual active others data gfactor) agentlist))
    (apply await agentlist)
    (assoc popu :poplist (for [agent agentlist] @agent))))

(defn crossover-inter
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

(defn crossover-discrete
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


(defn adapted-mutation
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

(defn do-offspring
  "Erstellt aus jedem Elternpaar 1 Kind"
  [acc parents-pair]
  (let [cross (crossover-inter (first parents-pair) (second parents-pair))
	mutated (adapted-mutation (first cross) (second cross))
	child {:tag 1 :chromosome (first mutated) :steps (second mutated) :fitness 0}
	old-pop-list (:poplist acc)]
    (assoc acc :poplist (conj old-pop-list child))))

(defn generate-offspring
  "Ruft do-offspring mit der aktuellen Population so oft auf, bis alle Elternpaare abgearbeitet sind. Liefert die komplette Population (Kinder und Eltern) zurück."
  [acc parents]
  (println "Generating Offspring...")
  (reduce do-offspring acc parents))

(defn survivor-selection
  "Bestimmt, welche Individuen in die nächste Generation übernommen werden. 
   Nur die besten Nachkommen werden für die neue Generation selektiert, entspricht der typischen ES Strategie."
  [popu popsize]
  (println "Survivor Selection...")
  (let [survivors (take popsize (sort-by :fitness < (for [i (:poplist popu) :when (= (:tag i) 1)] i)))]
    (assoc popu :poplist survivors)))

(defn parent-selection
  "Erstellt popsize*factor viele Elternpaare. Die Fitness spielt keine Rolle, die Eltern werden zufällig gezogen."
  [popu popsize factor]
  (println "Parent Selection...")
  (let [parents (scramble (take (* factor popsize) (cycle (:poplist popu))))
	splitted (split-at (* (/ factor 2) popsize) parents)]
    (map #(list %1 %2) (nth splitted 0) (nth splitted 1))))
;(def foo (init-population 10 22 8))
;(count (parent-selection foo 0.4 10))

(defn es
  "Startet den GA mit einer Population von popsize und x-vielen Runs."
  [active others data gaargs filename]
  (println "esrec started." gaargs (timestamp) "\n First evaluation...")
  (let [popsize (nth gaargs 0)
	gfactor (nth gaargs 1)
	offspFactor (nth gaargs 2)
	iterations (last gaargs)
	popu (init-population popsize 22 8)
	popu-evaluated (evaluate-all-firstrun active others popu data gfactor)]

    (loop [runs iterations popu popu-evaluated]
      (println "best in generation:" (:fitness (first (sort-by :fitness < (:poplist popu)))))
      (if (zero? runs)
	(:chromosome (first (sort-by :fitness < (:poplist popu))))
	(let [parents (parent-selection popu popsize offspFactor)
	      popu-with-children (generate-offspring popu parents)
	      popu-eva (evaluate-all active others popu-with-children data gfactor)
	      popu-surv (survivor-selection popu-eva popsize)
	      ;setze alle survivor auf :tag 0
	      new-pop (assoc popu-surv :poplist
			     (for [i (:poplist popu-surv)] (assoc i :tag 0)))]
	  (do (println "no. parent pairs:" (count parents))
	      (println "no. parents + children:" (count (:poplist popu-with-children)))
	      (println "no. survivors: " (count (:poplist popu-surv))))
	  (recur (dec runs) new-pop))))))