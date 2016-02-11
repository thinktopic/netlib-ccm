(ns netlib-ccm.core-test
  (:require [clojure.test :refer :all]
            [netlib-ccm.core :refer :all]
            [clojure.core.matrix :as m]
            [clojure.core.matrix.protocols :as mp]
            [clojure.tools.namespace.repl :as r]))



(deftest basic-vector
  (testing "creating and operating on a dense vector"
    (let [vec1 (new-dense-vector 10)]
      (is (= (map double (repeat 10 0)))
          (m/eseq vec1)))
    (let [vec2 (new-dense-vector-from-data (range 10))
          sub1 (m/subvector vec2 5 5)]
      (is (= (map double (range 5 10))
             (m/eseq sub1))))
    (let [vec1 (new-dense-vector 10)
          vec2 (new-dense-vector-from-data (range 10))]
      (m/assign! vec1 vec2)
      (is (= (map double (range 10))
             (m/eseq vec1))))))



(deftest basic-matrix
  (testing "creating and operating on a dense matrix"
    (let [mat1 (m/array :netlib (partition 3 (range 1 10)))
          submat (m/submatrix mat1 1 2 1 2)]
      [mat1 submat])))


(defn dot-perftest
  []
  (let [iterations 1000]
    (doseq [size [10 100 1000 10000]
            impl [:vectorz :netlib]]
      (println "size" size "impl" impl)
      (let [v1 (m/array impl (range 1 (+ size 1)))]
        (time (dotimes [iter iterations]
                (m/length v1)))))))


(defn add-scaled-perftest
  []
  (let [iterations 1000]
    (doseq [size [10 100 1000 10000]
            impl [:vectorz :netlib]]
      (println "size" size "impl" impl)
      (let [v1 (m/array impl (range 1 (+ size 1)))]
        (time (dotimes [iter iterations]
                (m/add-scaled v1 v1 10)))))))



(defn inner-product-perftest
  []
  (let [iterations 10]
    (doseq [size [10 100 500]
            impl [:vectorz :netlib]]
      (println "size" size "impl" impl)
      (let [a (m/array impl (repeat size (range 1 (+ size 1))))
            b (m/array impl (repeat size (repeat size 1)))]
        (time (dotimes [iter iterations]
                (m/inner-product a b)))))))


(defn inner-product
  [impl size]
  (let [a (m/array impl (repeat size (range 1 (+ size 1))))
        b (m/array impl (repeat size (repeat size 1)))]
    (m/inner-product a b)))

(defn inner-product-shapes
  [impl]
  (let [a (m/array impl (repeat 3 1))
        b (m/array impl (partition 3 (range 1 10)))
        c (m/array impl (partition 3 (repeat 9 1)))]
    (println "m" b "v" a "c" c)
    (println "m*v" (m/shape (m/inner-product b a)) (m/inner-product b a))
    (println "v*m" (m/shape (m/inner-product a b)) (m/inner-product a b))
    (println "m*c" (m/inner-product b c))
    (println "c*m" (m/inner-product c b))))


(defn pre-scale
  []
  (let [size 9
        m (m/array :netlib (range 1 10))
        a (m/array :netlib (repeat 9 1))]
    (mp/element-map (mp/convert-to-nested-vectors m)
                    (fn [v] (mp/pre-scale a v)))))


(defn outer-product
  [impl size]
  (let [a (m/array impl (range 1 (+ size 1)))
        b (m/array impl (repeat size 1))]
    (m/outer-product a b)))

(defn print-result
  [tag item]
  (println tag (m/shape item) item))

(defn add-scaled-product
  [impl]
  (let [a (m/array impl (repeat 3 1))
        b (m/array impl (partition 3 (range 1 10)))
        c (m/array impl (partition 3 (repeat 9 1)))]
    (println "m" b "v" a "c" c "s" 2.0 "f" 2.0)
    (print-result "m + m*v*f" (mp/add-scaled-product b b a 2.0))
    (print-result "v + v*m*f" (mp/add-scaled-product a a b 2.0))
    (print-result "v + v*v*f" (mp/add-scaled-product b b b 2.0))
    (print-result "m + m*c*f" (mp/add-scaled-product b b c 2.0))
    (print-result "m + m*s*f" (mp/add-scaled-product b b 2.0 2.0))
    (print-result "v + v*s*f" (mp/add-scaled-product a a 2.0 2.0))))
