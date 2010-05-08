(ns metaheuristics.pso
  (:use [clojure.contrib.math])
  (:use metaheuristics.testfunctions))

(defn- scramble
  [l]
  (let [items (java.util.ArrayList. l)]
	scramble (do (java.util.Collections/shuffle items) (seq items))))

; swarm
(defstruct particle :position :velocity :pbestpos :pbestfit)
(defstruct swarm :particlelist :gbest :vmax :hfit)

(defn- init-particle [nDimensions max]
  (agent (struct particle
		 (double-array (for [i (range nDimensions)] (rand max)))
		 (double-array (for [i (range nDimensions)] (rand max)))
		 (double-array (replicate nDimensions 1.0))
		 Double/MAX_VALUE)))

(defn- init-swarm [nParticles nDimensions vmaxDelta nHistoryFitness max]
  (let [particles (map (fn [_] (init-particle nDimensions max)) (range nParticles))
	best (init-particle nDimensions max)
	vmax (agent (double-array (repeat nDimensions (* vmaxDelta max))))
	hfit (agent (double-array (replicate nHistoryFitness Double/MAX_VALUE)))]
    (struct swarm particles best vmax hfit)))

(defn- inertia-weight 
  []  	 
  (double (+ 0.5 (/ (rand) 2.0))))

(defn- update-particle
  [particle fitness ftype gbest vmaxVec iteration]
  (let [#^doubles pos (:position particle)
	#^doubles vel (:velocity particle)
	#^doubles pbestpos (:pbestpos particle)
	pbestfit (double (:pbestfit particle))
	#^doubles gbestpos (:pbestpos @gbest)
	gbestfit (double (:pbestfit @gbest))

	; acceleration constants
	c1 1.494
	c2 1.494
	; random elements
	phi1 (double (rand))
	phi2 (double (rand))
	; inertia weight
	w (inertia-weight)

        ; update velocity
	newvel (amap vel idx ret
		     (+ (* c1 phi1 (- (aget gbestpos idx) (aget pos idx)))
			(* c2 phi2 (- (aget pbestpos idx) (aget pos idx)))
			(* (aget vel idx))))
	;velocity clamping
	vmax @vmaxVec
	clampvel (amap newvel i ret
			 (if (< (abs (aget newvel i)) (aget vmax i))
			   (aget newvel i)
			   (* (/ (aget vmax i) (abs (aget newvel i))) (aget newvel i))))

        ; update position
	newpos (amap pos idx ret
		     (+ (aget pos idx) (aget clampvel idx)))
	; update pbest
	fit (fitness newpos)
	newpbestfit (if (ftype fit pbestfit) fit pbestfit)
	newpbestpos (if (ftype fit pbestfit) newpos pbestpos)]
	; finally assign updates to agent
    (assoc particle :position newpos :velocity newvel 
	   :pbestpos newpbestpos :pbestfit newpbestfit)))

(defn- update-gbest [gbest ftype particles]
  (let [newbest (first (sort-by :pbestfit ftype (map #(deref %) particles)))
	#^doubles newpbestpos (:pbestpos newbest)
	newpbestfit (double (:pbestfit newbest))]
    (assoc gbest :pbestpos newpbestpos :pbestfit newpbestfit)))

(defn- update-hfit [hfit gbest]
  (let [gbestfit (double (:pbestfit @gbest))]
    (amap hfit i ret
	  (if (= i 0)
	    gbestfit
	    (aget hfit (- i 1))))))

(defn- update-vmax [vmax swarm iteration]
  (let [#^doubles nlastfit @(:hfit swarm)
	gbestfit (double (:pbestfit @(:gbest swarm)))
	prebeta (- 1.0 (* iteration 0.0001))
	beta (if (<= prebeta 0.0001) 0.0001 prebeta)]
    (amap vmax i ret
	  (if (every? (fn [x] (if (>= gbestfit x) true false)) nlastfit)
	    (* beta (aget vmax i))
	    (aget vmax i)))))

(defn- reset-vmax [vmax vmaxDelta max-feat]
  (double-array (repeat (count vmax) (* vmaxDelta max-feat))))

(defn- reset-particle [particle]
  (let [nDimensions (count (:position particle))
	newpos (double-array (for [i (range nDimensions)] (rand 255)))
	;;newvel (double-array (repeat nDimensions 0.0))
	newvel (double-array (for [i (range nDimensions)] (rand 255)))
	#^doubles newpbest (:pbestpos particle)]
    (assoc particle :position newpos :velocity newvel :pbestpos newpbest)))

(defn- fly [fitness ftype swarm iteration max-feat]
  ;;reset x random particles every 4 iteration  
  ;; (if (= (rem iteration 25) 0)
  ;;   (do
  ;;     (send (:vmax swarm) reset-vmax 1.0 max-feat)
  ;;     (await (:vmax swarm))
  ;;     (let [plist (map #(nth (:particlelist swarm) %)
  ;; 		       (take 15 (scramble (range 0 (count (:particlelist swarm))))))]
  ;; 	(dorun (map #(send %1 reset-particle) plist))
  ;; 	(apply await plist))))

  ; update particles
  (dorun (map #(send % update-particle fitness ftype (:gbest swarm) (:vmax swarm) iteration)
	      (:particlelist swarm)))
  (apply await (:particlelist swarm))

  ; update gbest
  (send (:gbest swarm) update-gbest ftype (:particlelist swarm))
  (await (:gbest swarm))
  (print (:pbestfit @(:gbest swarm)) " ")

  ;update  history fitness
  (send (:hfit swarm) update-hfit (:gbest swarm))
  (await (:hfit swarm))

  ;update vmax
  (send (:vmax swarm) update-vmax swarm iteration)
  (await (:vmax swarm))
  (println (aget @(:vmax swarm) 0) " ")
  )

(defn pso
  "Starts the PSO algorithm.
  Parameter: 
  fitness - defines the fitness function used in the PSO. The only argument is the position of the particle (a double-vector)
  ftype - defines the type of optimization problem (minimize (<) or maximize (>)).
  dim - number of dimensions in solution.
  nparticles - number of particles in swarm.
  vc-init - initial value for velocity clamping (usually 1.0)
  vc-hist - number of history value for velocity clamping (usually 0.8 * max-iterations).
  max-iterations - maximum number of iterations.
  max-feat - maximum value for one feature.
  "
  [fitness ftype dim nparticles vc-init vc-hist max-iterations max-feat]
  (let [swarm (init-swarm nparticles dim vc-init vc-hist max-feat)]
    (dorun (map (fn [i] (fly fitness ftype swarm i max-feat)) (range max-iterations)))
    @(:gbest swarm)))


;; (pso griewank < 10 30 10 1.0 (* 500 0.8) 10.0)

