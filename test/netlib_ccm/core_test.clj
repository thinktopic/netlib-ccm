(ns netlib-ccm.core-test
  (:require [clojure.test :refer :all]
            [netlib-ccm.core :refer :all]
            [clojure.core.matrix :as m]
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


(defn outer-product
  [impl size]
  (let [a (m/array impl (range 1 (+ size 1)))
        b (m/array impl (repeat size 1))]
    (m/outer-product a b)))
