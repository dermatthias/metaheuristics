(ns test)
(defn rastrigin [position]
  (let [sumvec (amap position i ret
	(+ (- (Math/pow (aget position i) 2)
	      (* 10 (Math/cos (* 2 (Math/PI) (aget position i)))))
	   10))]
    (areduce sumvec i ret 0
	     (+ ret (aget sumvec i)))))

(defn griewank [position]
  (let [prod (areduce position i ret 0
			 (* ret (Math/cos (/ (aget position i) (Math/sqrt (+ i 1))))))
	sum (/ (areduce position i ret 0
			(+ ret (Math/pow (aget position i) 2)))
	       4000)]
    (double (- (+ sum 1) prod))))
