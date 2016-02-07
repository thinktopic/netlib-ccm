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
    (let [mat1 (m/array :netlib (partition 3 (range 1 10)))])))
