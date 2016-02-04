(ns netlib-ccm.core
  (:require [clojure.core.matrix :as ma]
            [clojure.core.matrix.protocols :as mp]
            [clojure.core.matrix.implementations :as mi]
            [clojure.reflect :as r])
  (:import [com.github.fommil.netlib BLAS]))


(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;;A view on the data that could be non-contiguous
(defrecord StridedView [^doubles data ^long offset ^long row-count ^long column-count ^long row-stride])

(defn strided-op
  "Perform an operation such as (assign! lhs rhs).
Op gets passed: lhs-double-array lhs-offset rhs-double-array rhs-offset op-amount"
  [^StridedView lhs ^StridedView rhs op]
  (let [rhs-data-stride (.column-count rhs)
        lhs-data-stride (.column-count lhs)]
      (loop [rhs-row 0]
        (when (< rhs-row (.row-count rhs))
          (loop [rhs-data-offset 0]
            (when (< rhs-data-offset rhs-data-stride)
              (let [rhs-data-copied (+ rhs-data-offset (* rhs-row (.column-count rhs)))
                    lhs-data-offset (rem rhs-data-copied lhs-data-stride)
                    lhs-row (quot rhs-data-copied lhs-data-stride)
                    op-amount (min (- lhs-data-stride lhs-data-offset)
                                   (- rhs-data-stride rhs-data-offset))
                    lhs-offset (+ (.offset lhs) (* lhs-row (.row-stride lhs)) lhs-data-offset)
                    rhs-offset (+ (.offset rhs) (* rhs-row (.row-stride rhs)) rhs-data-offset)]
                (op (.data lhs) lhs-offset (.data rhs) rhs-offset op-amount)
                (recur (+ rhs-data-offset op-amount)))))
          (recur (inc rhs-row))))))

(defn assign-strided-view!
  [^StridedView lhs ^StridedView rhs]
  (strided-op lhs rhs (fn [lhs-data lhs-offset rhs-data rhs-offset op-len]
                        (System/arraycopy ^doubles rhs-data ^long rhs-offset ^doubles lhs-data ^long lhs-offset
                                          ^long op-len))))

(defn clone-strided-view
  "Create a packed strided view from potentially non-dense view"
  ^StridedView [^StridedView source]
  (let [num-items (* (.row-count source) (.column-count source))
        ^doubles data (make-array ^doubles num-items)
        retval (->StridedView data 0 1 num-items num-items)]
    (assign-strided-view retval source)
    retval))


(definterface NetlibItem)

(definterface AbstractView
  (^StridedView getStridedView []))

(definterface AbstractVector
  (^double get [^long idx])
  (set [^long idx ^double val])
  (^long length [])
  (^netlib_ccm.core.AbstractVector clone []))

(definterface AbstractMatrix
  (^double get [^long row ^long column])
  (set [^long row ^long column ^double val])
  (^netlib_ccm.core.AbstractVector getRow [^long row])
  (setRow [^long row ^netlib_ccm.core.AbstractVector data])
  (^long rowCount [])
  (^long columnCount [])
  (^netlib_ccm.core.AbstractMatrix clone []))

(deftype DenseArray [^doubles data ^long offset ^long length]
  netlib_ccm.core.NetlibItem
  netlib_ccm.core.AbstractView
  (^StridedView getStridedView [] (->StridedView data offset 1 length length))

  netlib_ccm.core.AbstractVector
  (^double get [this ^long idx] (aget data (+ offset idx)))
  (set [this ^long idx ^double val] (aset data (+ offset idx) val))
  (^long length [this] length)
  (^netlib_ccm.core.AbstractVector clone [this] (let [retval (clone-strided-view (.getStridedView this))]
                                                  (->DenseArray (.data retval) 0 length)))

  clojure.lang.Seqable
  (seq [this] (map #(.get this %) (range length))))


(deftype DenseMatrix [^DenseArray data ^long row-count ^long column-count]
  netlib_ccm.core.NetlibItem
  netlib_ccm.core.AbstractView
  (^StridedView getStridedView [] (let [item-count (* row-count column-count)]
                                    (->StridedView (.data data) (.offset data) 1 item-count item-count)))
  netlib_ccm.core.AbstractMatrix
  (^double get [this ^long row ^long column] (aget ^doubles (.data data) (+ (.offset data)
                                                                            column
                                                                            (* row column-count))))
  (set [this ^long row ^long column ^double val] (aset ^doubles (.data data)
                                                       (+ (.offset data)
                                                          column
                                                          (* row column-count))
                                                       val))
  (^netlib_ccm.core.AbstractVector getRow [this ^long row] (->DenseArray ^doubles (.data data)
                                                                         (+ (.offset data)
                                                                            (* row column-count))
                                                                         column-count))
  (setRow [this ^long row ^netlib_ccm.core.AbstractVector data]
    (let [^AbstractView strided-data data
          ^AbstractView row-data (.getRow this row)]
      (assign-strided-view (.getStridedView row-data)
                           (.getStridedView strided-data))))
  (^long rowCount [this] row-count)
  (^long columnCount [this] column-count)
  (^netlib_ccm.core.AbstractMatrix clone [this] (deep-copy-matrix-view this)))

;;A strided view of data where you can have N contiguous elements separated in rows of length Y.
;;Think of a small submatrix in a larger one
(deftype StridedDenseMatrix [^DenseArray data ^long row-count ^long column-count ^long row-stride])
;;Interpret a submatrix as an array
(deftype StridedDenseArray [^StridedDenseMatrix data])

(defn deep-copy-array-view
  (^DenseArray [^DenseArray item ^long offset ^long length]
   (let [max-possible (min length (max 0 (- (.length item) offset)))
         new-data (make-array Double/TYPE max-possible)
         ^doubles data (.data item)]
     (System/arraycopy data (+ offset (.offset item)) new-data 0 max-possible)
     (->DenseArray new-data 0 max-possible)))
  (^DenseArray [^DenseArray item]
   (deep-copy-array-view item 0 (.length item))))


(defn deep-copy-matrix-view
  ^DenseMatrix [^DenseMatrix item]
  (let [new-data (deep-copy-array-view (.data item))]
    (->DenseMatrix new-data (.row-count item) (.column-count item))))


(defn in-range
  [^long item ^long min-eq ^long max]
  (and (>= item min-eq)
       (< max)))

(defn new-array-view
  [^long length]
  (let [new-data (make-array Double/TYPE length)]
    (->DenseArray new-data 0 length)))

(defn new-array-view-from-array
  [^doubles double-array]
  (let [len (count double-array)]
    (->DenseArray double-array 0 len)))

(defn double-array-from-data
  [data]
  (double-array (ma/eseq data)))

(defn new-array-view-from-data
  [data]
  (new-array-view-from-array (double-array-from-data data)))

(defn do-construct-matrix
  [data]
  (let [shape (ma/shape data)
        num-shape (count shape)]
    (case num-shape
      1 (new-array-view-from-array (double-array-from-data data))
      2 (->DenseMatrix (new-array-view-from-array
                        (double-array-from-data data))
                       (first shape)
                       (second shape)))))

(defn do-new-matrix-nd
  [shape]
  (case (count shape)
    1 (new-array-view (first shape))
    2 (->DenseMatrix (new-array-view (* (long (first shape)) (long (second shape)))
                                     (first shape) (second shape)))))

(eval
 `(extend-protocol mp/PImplementation
    ~@(mapcat (fn [sym]
                (cons sym
                      '(
                        (implementation-key [m] :netlib)
                        (supports-dimensionality? [m dims] (in-range dims 0 3))
                        (construct-matrix [m data] (do-construct-matrix data))
                        (new-vector [m length] (new-array-view length))
                        (new-matrix [m rows columns] (->DenseMatrix (new-array-view (* rows columns) rows columns)))
                        (new-matrix-nd [m shape] (do-new-matrix-nd shape)))))
        ['DenseArray 'DenseMatrix])))


(extend-protocol mp/PDimensionInfo
  DenseArray
  (dimensionality [m] 1)
  (get-shape [m] [(.length m)])
  (is-scalar? [m] false)
  (is-vector? [m] true)
  (dimension-count [m dim] (if (= 0 dim) (.length m) (throw (Exception. "Unsupported"))))
  DenseMatrix
  (dimensionality [m] 2)
  (get-shape [m] [(.row-count m) (.column-count m)])
  (is-scalar? [m] false)
  (is-vector? [m] true)
  (dimension-count [m dim]
    (case (long dim)
      0 (.row-count m)
      1 (.column-count m)
      (throw (Exception. "Unsupported")))))

(extend-protocol mp/PIndexedAccess
  DenseArray
  (get-1d [m row] (let [^DenseArray view m
                        view-data ^doubles (.data view)]
                    (aget view-data (+ (long row) (.offset view)))))
  (get-2d [m row column] (throw (Exception. "Unsupported")))
  (get-nd [m indexes] (if (= 1 (count indexes))
                        (mp/get-1d m (first indexes))
                        (throw (Exception. "Unsupported"))))
  DenseMatrix
  (get-1d [m row] (let [^DenseMatrix m m
                        row-stride (.column-count m)
                        row-offset (* row-stride row)
                        ^DenseArray m-data (.data m)]
                    (->DenseArray (.data m-data) (+ (.offset m-data) row-offset) row-stride)))
  (get-2d [m row column] (let [^DenseMatrix m m
                               row-stride (.column-count m)
                               row-offset (* row-stride row)
                               ^DenseArray m-data (.data m)
                               ^doubles ddata (.data m-data)]
                           (aget ddata (+ (.offset m-data) row-offset column))))
  (get-nd [m indexes] (let [idx-count (count indexes)]
                        (case idx-count
                          1 (mp/get-1d m (first indexes))
                          2 (mp/get-2d m (first indexes) (second indexes))
                          (throw (Exception. "Unsupported"))))))

(extend-protocol mp/PIndexedSettingMutable
  DenseArray
  (set-1d! [m row v] (let [^DenseArray m m
                           ^long row row
                           ^double v v]
                       (aset ^doubles (.data m) (+ (.offset m) row) v)))
  (set-2d! [m row column v] (throw (Exception. "Unsupported")))
  (set-nd [m indexes v]
    (case (count indexes)
      1 (mp/set-1d! m (first indexes) v)
      (throw (Exception. "Unsupported"))))
  DenseMatrix
  (set-1d! [m row v] (let [^DenseMatrix m m
                           ^DenseArray data (.data m)
                           ^DenseArray incoming v
                           row-stride (.column-count m)]
                       (System/arraycopy ^doubles (.data incoming)
                                         (.offset incoming)
                                         ^doubles (.data data)
                                         (+ (.offset data) (* row-stride row))
                                         row-stride)))
  (set-2d! [m row column v] (let [^DenseMatrix m m
                                  row-stride (.column-count m)
                                  m-offset (+ ^long column (* row-stride row))
                                  ^double v v
                                  ^DenseArray data (.data m)]
                              (aset ^doubles (.data data) (+ (.offset data) m-offset)) v))
  (set-nd! [m indexes v]
    (case (count indexes)
      1 (mp/set-1d! m (first indexes) v)
      2 (mp/set-2d! m (first indexes) (second indexes) v)
      (throw (Exception. "Unsupported")))))


(extend-protocol mp/PMatrixCloning
  DenseArray
  (clone [m] (deep-copy-array-view m))
  DenseMatrix
  (clone [m] (deep-copy-matrix-view m)))


(extend-protocol mp/PIndexedSetting
  DenseArray
  (set-1d [m row v] (let [retval (deep-copy-array-view m)]
                      (mp/set-1d retval row v)))
  (set-2d [m row column v] (throw (Exception. "Unsupported")))
  (set-nd [m indexes v]
    (case (count indexes)
      1 (mp/set-1d m indexes v)
      (throw (Exception. "Unsupported"))))
  (is-mutable? [m] true)

  DenseMatrix
  (set-1d [m row v] (let [retval (deep-copy-matrix-view m)]
                      (mp/set-1d! retval row v)))
  (set-2d [m row column v] (let [retval (deep-copy-matrix-view m)]
                             (mp/set-2d! retval row column v)))
  (set-nd [m indexes v]
    (case (count indexes)
      1 (mp/set-1d m (first indexes) v)
      2 (mp/set-2d m (first indexes) (second indexes) v)
      (throw (Exception. "Unsupported"))))
  (is-mutable? [m] true))


(extend-protocol mp/PTypeInfo
  DenseArray
  (element-type [m] Double/TYPE)
  DenseMatrix
  (element-type [m] Double/TYPE))


(extend-protocol mp/PValidateShape
  DenseArray
  (validate-shape [m] (let [^DenseArray m m
                            ^doubles data (.data m)
                            data-len (count data)
                            item-end (+ (.length m) (.offset m))]
                        (when (< (.offset m) 0)
                          (throw (Exception. "Array offset is less than zero")))
                        (when (> item-end data-len)
                          (throw (Exception. "Array end is past end of array data")))
                        [(.length m)]))
  DenseMatrix
  (validate-shape [m] (let [^DenseMatrix m m
                            ^DenseArray data (.data m)
                            ^doubles backing (.data data)
                            data-len (count data)
                            item-count (* (.row-count m) (.column-count m))
                            item-end (+ item-count (.offset data))]
                        (mp/validate-shape data)
                        (when (> item-end data-len)
                          (throw (Exception. "Matrix end is past end of data")))
                        [(.row-count m) (.column-count m)])))


;;TODO view or deep copy??
(extend-protocol mp/PRowColMatrix
  DenseArray
  (column-matrix [m data] (let [^DenseArray m m]
                            (->DenseMatrix m 1 (.length m))))
  (row-matrix [m data] (let [^DenseArray m m]
                         (->DenseMatrix m (.length m) 1))))


(extend-protocol mp/PMutableMatrixConstruction
  DenseArray
  (mutable-matrix [m] (deep-copy-array-view m))
  DenseMatrix
  (mutable-matrix [m] (deep-copy-matrix-view m)))


(extend-protocol mp/PMutableCoercion
  DenseArray
  (ensure-mutable [m] m)
  DenseMatrix
  (ensure-mutable [m] m))


(extend-protocol mp/PDense
  DenseArray
  (dense-coerce [m data] (new-array-view-from-data data))
  (dense [m] true)
  DenseMatrix
  (dense-coerce [m data] (new-array-view-from-data data))
  (dense [m] true))


(extend-protocol mp/PConversion
  DenseArray
  (convert-to-nested-vectors [m] (vec (.data ^DenseArray m)))
  DenseMatrix
  (convert-to-nested-vectors [m] (vec (.data ^DenseArray (.data ^DenseMatrix m)))))


(extend-protocol mp/PReshaping
  DenseArray
  (reshape [m shape] (let [^DenseArray m m
                           num-desired (apply * shape)
                           data-len (.length m)]
                       (when (> num-desired data-len)
                         (throw (Exception. "Attempt to reshape to a larger backing storage")))
                       (case (count shape)
                         1 (deep-copy-array-view m 0 num-desired)
                         2 (->DenseMatrix (deep-copy-array-view m 0 num-desired) (first shape) (second shape))
                         (throw (Exception. "Unsupported")))))
  DenseMatrix
  (reshape [m shape] (mp/reshape (.data ^DenseMatrix m) shape)))
