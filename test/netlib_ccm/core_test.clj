(ns netlib-ccm.core-test
  (:require [clojure.test :refer :all]
            [netlib-ccm.core :as nc]
            [clojure.core.matrix :as m]
            [clojure.core.matrix.protocols :as mp]
            [clojure.tools.namespace.repl :as r])
  (:import [com.github.fommil.netlib BLAS]
           [netlib_ccm Ops]))



(deftest basic-vector
  (testing "creating and operating on a dense vector"
    (let [vec1 (nc/new-dense-vector 10)]
      (is (= (map double (repeat 10 0)))
          (m/eseq vec1)))
    (let [vec2 (nc/new-dense-vector-from-data (range 10))
          sub1 (m/subvector vec2 5 5)]
      (is (= (map double (range 5 10))
             (m/eseq sub1))))
    (let [vec1 (nc/new-dense-vector 10)
          vec2 (nc/new-dense-vector-from-data (range 10))]
      (m/assign! vec1 vec2)
      (is (= (map double (range 10))
             (m/eseq vec1))))))


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
      (let [a (m/pack (m/array impl (repeat size (range 1 (+ size 1)))))
            b (m/pack (m/array impl (repeat size (repeat size 1))))]
        (time (dotimes [iter iterations]
                (m/inner-product a b)))))))


(defn transpose-perftest
  []
  (let [iterations 10]
    (doseq [size [10 100 500]
            impl [:vectorz :netlib]]
      (println "size" size "impl" impl)
      (let [a (m/array impl (repeat size  (range 1 (+ size 1))))]
        (time (dotimes [iter iterations]
                (m/transpose a)))))))


(defn inner-product
  [impl size]
  (let [a (m/array impl (repeat size (range 1 (+ size 1))))
        b (m/array impl (repeat size (repeat size 1)))]
    (m/inner-product a b)))

(defn compare-items
  [a b]
  (is (= (m/shape a) (m/shape b)))
  (is (= (map double (m/eseq a)) (map double (m/eseq b)))))

(defn inner-product [a b impl]
  (m/inner-product (m/array impl a)
                   (m/array impl b)))


(defn inner-compare [name a b]
  (testing name (compare-items (inner-product a b :vectorz)
                               (inner-product a b :netlib))))

(deftest test-inner-products
  []
  (let [a (repeat 3 1)
        b (partition 3 (range 1 10))
        c (partition 3 (repeat 9 1))]

    (inner-compare "m*v" b a)
    (inner-compare "v*m" a b)
    (inner-compare "m*c" b c)
    (inner-compare "c*m" c b)))


(defn outer-product
  [impl size]
  (let [a (m/array impl (range 1 (+ size 1)))
        b (m/array impl (repeat size 1))]
    (m/outer-product a b)))

(defn print-result
  [tag item]
  (println tag (m/shape item) item))

(defn element-multiply
  [a b impl]
  (let [a (m/array impl a)
        b (m/array impl b)]
    (mp/element-multiply a b)))

(defn test-elem-mul
  [name a b]
  (testing name (compare-items (element-multiply a b :vectorz)
                               (element-multiply a b :netlib))))


(deftest test-element-multiply
  (let [v (repeat 3 2)
        m (partition 3 (range 1 10))]
    (test-elem-mul "m*v" m v)
    (test-elem-mul "m*m" m m)
    (test-elem-mul "v*m" v m)
    (test-elem-mul "v*v" v v)))


(defn element-multiply!
  [a b impl]
  (let [a (m/array impl a)
        b (m/array impl b)]
    (mp/element-multiply! a b)))

(defn test-elem-mul!
  [name a b]
  (testing name (compare-items (element-multiply a b :vectorz)
                               (element-multiply a b :netlib))))


(deftest test-element-multiply!
  (let [v (repeat 3 2)
        m (partition 3 (range 1 10))]
    (test-elem-mul! "m*v" (m/clone m) v)
    (test-elem-mul! "m*m" (m/clone m) m)
    (test-elem-mul! "v*v" v v)))

(defn to-array-unless-scalar
  [item impl]
  (if (m/scalar? item)
    item
    (m/array impl item)))

(defn add-scaled-product
  [m a b factor impl]
  (m/add-scaled-product (to-array-unless-scalar m impl)
                        (to-array-unless-scalar a impl)
                        (to-array-unless-scalar b impl)
                        factor))

(defn do-test-add-scaled-product
  [name m a b factor]
  (testing name (compare-items (add-scaled-product m a b factor :vectorz)
                               (add-scaled-product m a b factor :netlib))))


(deftest test-add-scaled-product
  (let [a (repeat 3 1)
        b (partition 3 (range 1 10))
        c (partition 3 (repeat 9 1))]
    (do-test-add-scaled-product "m + m*v*f" b b a 2.0)
    ;;Not a valid test case; m needs to be the larger of either v m
    ;(print-result "v + v*m*f" (mp/add-scaled-product a a b 2.0))
    (do-test-add-scaled-product "v + v*v*f" b b b 2.0)
    (do-test-add-scaled-product "m + m*c*f" b b c 2.0)
    (do-test-add-scaled-product "m + m*s*f" b b 2.0 2.0)
    (do-test-add-scaled-product "v + v*s*f" a a 2.0 2.0)))


(defn perf-test-element-multiply
  []
  (doseq [elem-count [10 100 1000 10000 100000 1000000]]
    (println "elem-count" elem-count)
    (let [^netlib_ccm.core.DenseVector y (m/array :netlib (repeat elem-count 1))
          ^netlib_ccm.core.DenseVector x (m/array :netlib (repeat elem-count 2))
          ^netlib_ccm.core.DenseVector a (m/array :netlib (repeat elem-count 3))]
      (print "clojure:")
      (time (dotimes [iter 1000]
              (nc/dense-element-multiply! 1.0
                                          (.data a) 0
                                          (.data x) 0
                                          1.0
                                          (.data y) 0
                                          elem-count)))
      (print "java   :")
      (time (dotimes [iter 1000]
              (Ops/alphaAXPlusBetaY elem-count 1.0
                                    ^doubles (.data a) 0
                                    ^doubles (.data x) 0
                                    1.0
                                    ^doubles (.data y) 0))))))


(defn perf-test-axpy
  []
  (doseq [elem-count [10 50 100 150 200 1000 10000 100000 1000000 10000000]]
    (println "elem-count" elem-count)
    (let [^netlib_ccm.core.DenseVector y (m/array :netlib (repeat elem-count 1))
          ^netlib_ccm.core.DenseVector x (m/array :netlib (repeat elem-count 2))
          ^netlib_ccm.core.DenseVector a (m/array :netlib (repeat elem-count 3))
          op-len elem-count]
      (print "netlib:")
      (time (dotimes [iter 1000]
              (.daxpy (BLAS/getInstance) op-len 1.0
                      ^doubles (.data x) 0 1
                      ^doubles (.data y) 0 1)))
      (print "java   :")
      (time (dotimes [iter 1000]
              (Ops/axpy elem-count 1.0
                        ^doubles (.data x) 0
                        ^doubles (.data y) 0))))))


(defn perf-test-sqrt!
  []
  (doseq [elem-count [10 50 100 150 200 1000 10000 100000 1000000 10000000]]
    (println "elem-count" elem-count)
    (let [^netlib_ccm.core.DenseVector y (m/array :netlib (repeat elem-count 1))
          ^netlib_ccm.core.DenseVector x (m/array :netlib (repeat elem-count 2))
          ^netlib_ccm.core.DenseVector a (m/array :netlib (repeat elem-count 3))
          op-len elem-count]
      (print "library:")
      (time (dotimes [iter 1000]
              (m/sqrt! y)))
      (print "reify   :")
      (time (dotimes [iter 1000]
              (Ops/OpY elem-count (.data y) 0 (reify netlib_ccm.IUnaryOp
                                                (op [this input]
                                                  (Math/sqrt input))))))
      (print "macro   :")
      (time (dotimes [iter 1000]
              (nc/unary-op-macro! y (Math/sqrt lhs-value))))
      )))


(defn perf-test-binary-op!
  []
  (doseq [elem-count [10 50 100 150 200 1000 10000 100000 1000000]]
    (println "elem-count" elem-count)
    (let [^netlib_ccm.core.DenseVector y (m/array :netlib (repeat elem-count 1))
          ^netlib_ccm.core.DenseVector x (m/array :netlib (repeat elem-count 2))
          op-len elem-count
          y-view (.getStridedView y)
          x-view (.getStridedView x)]
      (print "closure:")
      (time (dotimes [iter 1000]
              (nc/strided-view-multiple-op! y-view
                                            x-view
                                            (fn [lhs-data lhs-offset rhs-data rhs-offset op-len]
                                              (Ops/OpXY op-len rhs-data rhs-offset
                                                        lhs-data lhs-offset
                                                        (reify netlib_ccm.IBinaryOp
                                                          (op [this rhs-val lhs-val]
                                                            (/ lhs-val rhs-val)))))
                                            (fn [^doubles lhs-data ^long lhs-offset
                                                 ^long op-len ^double rhs-val]
                                              (Ops/OpY op-len lhs-data lhs-offset
                                                       (reify netlib_ccm.IUnaryOp
                                                         (op [this lhs-val]
                                                           (/ lhs-val rhs-val))))))))
      (print "macro   :")
      (time (dotimes [iter 1000]
              (nc/strided-view-binary-java-op! y-view x-view (/ lhs-value rhs-value))))))
  )
