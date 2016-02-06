(ns netlib-ccm.core
  (:require [clojure.core.matrix :as ma]
            [clojure.core.matrix.protocols :as mp]
            [clojure.core.matrix.implementations :as mi]
            [clojure.reflect :as r])
  (:import [com.github.fommil.netlib BLAS]))


(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;;A view on the data that could be non-contiguous
;;Initial column count means the number of items in the
;;first row of the strided view.
;;Thus the strided view data actually starts at offset
;;  (+ offset (- column-count initial-column-count))
;;Furthermore initial-column-count and last-column-count
;;must be <= column-count
;;a. If the view has 1 row, then last-column-count must be zero.
;;b. (>= row-stride column-count)
;;There is slight ambiguity because given last-column-count if there are
;;any rows at all past the initial row you can add 1 to the row count
;;and set last-column-count to 0 and have the same object meaning.
(defrecord StridedView [^doubles data ^long offset ^long row-count ^long column-count
                        ^long row-stride ^long initial-column-count ^long last-column-count])

(defn new-strided-view
  (^StridedView [data offset row-count column-count
                 row-stride initial-column-count last-column-count]
   (when (and (= row-count 1)
              (or
               (not= 0 last-column-count)
               (= 0 initial-column-count)))
     (throw (Exception. "Invalid strided view format a.")))

   (when (< row-stride column-count)
     (throw (Exception. "Row stride is less than column count")))

   (->StridedView data offset row-count column-count
                  row-stride initial-column-count last-column-count))

  (^StridedView [data offset row-count column-count row-stride]
   (->StridedView data offset row-count column-count row-stride column-count 0))
  (^StridedView [^doubles data ^long offset ^long row-count ^long column-count]
   (->StridedView data offset row-count column-count column-count column-count 0))

  (^StridedView [^doubles data ^long offset ^long length]
   (->StridedView data offset 1 length length length 0)))


(defn check-strided-view-row-index
  [^StridedView view, ^long row-idx]
  (when (>= row-idx (.row-count view))
    (throw (Exception. "Attempt to access past end of view")))
  (when (< row-idx 0)
    (throw (Exception. "Attempt to access before beginning of view"))))

(defn get-strided-view-column-count
  ^long [^StridedView view, ^long row]
  (cond
    (= row 0) (.initial-column-count view)
    (= row (- (.row-count view) 1)) (.last-column-count view)
    :else (.column-count view)))

(defn get-strided-view-data-length
  ^long [^StridedView view]
  (+ (* (max 0 (- (.row-count view) 2)) (.column-count view))
     (.initial-column-count view)
     (.last-column-count view)))

(defn check-strided-view-data-offset
  [^StridedView view ^long data-offset]
  (when (< data-offset 0)
    (throw (Exception. "Attempt to access before beginning of view")))
  (when (>= data-offset (get-strided-view-data-length view))
    (throw (Exception. "Attempt to access past end of view"))))

(defn get-strided-view-row-from-data-offset
  ^long [^StridedView view ^long data-offset]
  (check-strided-view-data-offset view data-offset)
  (if (< data-offset (.initial-column-count view))
    0
    (let [data-offset (- data-offset (.initial-column-count view))]
      (+ 1 (quot data-offset (.column-count view))))))

(defn get-strided-view-row-length-from-data-offset
  ^long [^StridedView view ^long data-offset]
  (check-strided-view-data-offset view data-offset)
  (if (< data-offset (.initial-column-count view))
    (- (.initial-column-count view) data-offset)
    (let [rest-offset (- data-offset (.initial-column-count view))
          row-idx (+ 1 (quot rest-offset (.column-count view)))
          rest-leftover (rem rest-offset (.column-count view))]
      (if (= row-idx (- (.row-count view) 1))
        (- rest-leftover (.last-column-count view))
        (- rest-leftover (.column-count view))))))


(defn get-strided-view-total-offset
  ^long [^StridedView view ^long data-offset]
  (check-strided-view-data-offset view data-offset)
  (let [data-len (get-strided-view-data-length view)
        body-len (* (.column-count view) (max 0 (- (.row-count view) 2)))]
    (if (< data-offset (.initial-column-count view))
      (+ (+ (.offset view) (- (.column-count view) (.initial-column-count view)))
         data-offset)
      (let [data-rest (- data-offset (.initial-column-count view))]
        (+ (.offset view)
           (* (.row-stride view) (+ 1 (quot data-rest (.column-count view))))
           (rem data-rest (.column-count view)))))))


(defn create-sub-strided-view
  [^StridedView view ^long offset ^long length]
  (check-strided-view-data-offset view offset)
  (check-strided-view-data-offset view (+ offset (max 0 (- length 1))))
  (let [data-len (get-strided-view-data-length view)]
    (let [start-submat-offset (- (get-strided-view-total-offset view offset) (.offset view))
          ini-col-offset (rem start-submat-offset (.column-count view))
          ini-col-len (- (.column-count view) ini-col-offset)
          end-col-len (rem (- length ini-col-len) (.column-count view))
          body-num-rows (/ (- length (+ ini-col-len end-col-len))
                           (.row-stride view))
          ini-row (long (if-not (= 0 ini-col-len) 1 0))
          end-row (long (if-not (= 0 end-col-len) 1 0))]
      (new-strided-view (.data view) (.offset view)
                        (+ body-num-rows ini-row end-row)
                        (.column-count view)
                        (.row-stride view)
                        ini-col-len
                        end-col-len))))

(defn strided-op
  "Perform an operation such as (assign! lhs rhs).
Op gets passed: lhs-double-array lhs-offset rhs-double-array rhs-offset op-amount"
  [^StridedView lhs ^StridedView rhs op]
  (loop [rhs-row 0]
    (when (< rhs-row (.row-count rhs))
      (let [data-offset (+ (* (.row-stride rhs) (max 0 (- (.row-count rhs) 2)))
                           (.initial-column-count rhs))
            rhs-row-len (get-strided-view-column-count rhs rhs-row)]
       (loop [rhs-row-offset 0]
         (when (< rhs-row-offset rhs-row-len)
           (let [data-offset (+ data-offset rhs-row-offset)
                 rhs-total-offset (get-strided-view-total-offset data-offset)
                 rhs-row-len (- rhs-row-len rhs-row-offset)
                 lhs-row (get-strided-view-row-from-data-offset lhs data-offset)
                 lhs-row-len (get-strided-view-row-length-from-data-offset lhs data-offset)
                 lhs-total-offset (get-strided-view-total-offset data-offset)
                 op-amount (min lhs-row-len rhs-row-len)]
             (op (.data lhs) lhs-total-offset
                 (.data rhs) rhs-total-offset
                 op-amount)
             (recur (+ rhs-row-offset op-amount))))))
      (recur (inc rhs-row)))))


(defn assign-strided-view!
  "lhs must be at least as large as rhs"
  [^StridedView lhs ^StridedView rhs]
  (strided-op lhs rhs (fn [lhs-data lhs-offset rhs-data rhs-offset op-len]
                        (System/arraycopy ^doubles rhs-data ^long rhs-offset
                                          ^doubles lhs-data ^long lhs-offset
                                          ^long op-len))))

(defn clone-strided-view
  "Create a packed strided view from potentially non-dense view"
  ^StridedView [^StridedView source]
  (let [num-items (* (.row-count source) (.column-count source))
        ^doubles data (make-array Double/TYPE num-items)
        retval (new-strided-view data 0 1 num-items num-items)]
    (assign-strided-view! retval source)
    retval))

(defn validate-strided-view-shape [^StridedView data]
  (let [data-len (count (.data data))
        item-end (+ (.offset data)
                    (* (- (.row-count data) 1) (.row-stride data))
                    (.column-count data))]
    (when (> item-end data-len)
      (throw (Exception. "Invalid strided view")))))


(definterface NetlibItem)

(definterface AbstractView
  (^netlib_ccm.core.StridedView getStridedView []))

(definterface AbstractVector
  (^double get [^long idx])
  (set [^long idx ^double val])
  (^long length [])
  (^netlib_ccm.core.AbstractVector clone [])
  (validateShape []))

(definterface AbstractMatrix
  (^double get [^long row ^long column])
  (set [^long row ^long column ^double val])
  (^netlib_ccm.core.AbstractVector getRow [^long row])
  (setRow [^long row ^netlib_ccm.core.AbstractVector data])
  (^long rowCount [])
  (^long columnCount [])
  (^netlib_ccm.core.AbstractMatrix clone [])
  (validateShape []))

(declare ->DenseVector ->DenseMatrix)

(deftype DenseVector [^doubles data ^long offset ^long length]
  netlib_ccm.core.NetlibItem
  netlib_ccm.core.AbstractView
  (^StridedView getStridedView [this] (new-strided-view data offset 1 length length))

  netlib_ccm.core.AbstractVector
  (^double get [this ^long idx] (aget data (+ offset idx)))
  (set [this ^long idx ^double val] (aset data (+ offset idx) val))
  (^long length [this] length)
  (^netlib_ccm.core.AbstractVector clone [this]
   (let [retval (clone-strided-view (.getStridedView this))]
     (->DenseVector (.data retval) 0 length)))
  (validateShape [this] (validate-strided-view-shape (.getStridedView this)))

  clojure.lang.Seqable
  (seq [this] (map #(.get this %) (range length))))


(deftype DenseMatrix [^DenseVector data ^long row-count ^long column-count]
  netlib_ccm.core.NetlibItem
  netlib_ccm.core.AbstractView
  (^StridedView getStridedView [this]
   (let [item-count (* row-count column-count)]
     (new-strided-view (.data data) (.offset data)
                    1 item-count item-count)))
  netlib_ccm.core.AbstractMatrix
  (^double get [this ^long row ^long column]
   (aget ^doubles (.data data)
         (+ (.offset data)
            column
            (* row column-count))))
  (set [this ^long row ^long column ^double val]
    (aset ^doubles (.data data)
          (+ (.offset data)
             column
             (* row column-count))
          val))
  (^netlib_ccm.core.AbstractVector getRow [this ^long row]
   (->DenseVector ^doubles (.data data)
                 (+ (.offset data)
                    (* row column-count))
                 column-count))
  (setRow [this ^long row ^netlib_ccm.core.AbstractVector data]
    (let [^AbstractView strided-data data
          ^AbstractView row-data (.getRow this row)]
      (assign-strided-view! (.getStridedView row-data)
                            (.getStridedView strided-data))))
  (^long rowCount [this] row-count)
  (^long columnCount [this] column-count)
  (^netlib_ccm.core.AbstractMatrix clone [this]
   (let [retval (clone-strided-view (.getStridedView this))
         ary (->DenseVector (.data retval) 0 (* row-count column-count))]
     (->DenseMatrix ary row-count column-count)))
  (validateShape [this] (validate-strided-view-shape (.getStridedView this)))
  clojure.lang.Seqable
  (seq [this] (map #(.getRow this %) (range (.rowCount this)))))

(declare strided-view-to-vector)

;;A strided view of data where you can have N contiguous elements separated in rows of length Y.
;;Think of a small submatrix in a larger one
(deftype StridedMatrix [^StridedView data row-count column-count]
  netlib_ccm.core.NetlibItem
  netlib_ccm.core.AbstractView
  (^StridedView getStridedView [this] data)

  netlib_ccm.core.AbstractMatrix
  (^double get [this ^long row ^long column]
   (aget ^doubles (.data data)
         (get-strided-view-total-offset data (+ (* column-count row) column))))

  (set [this ^long row ^long column ^double val]
    (aset ^doubles (.data data)
          (get-strided-view-total-offset data (+ (* column-count row) column))
          val))

  (^netlib_ccm.core.AbstractVector getRow [this ^long row]
   (let [sub-view (create-sub-strided-view data (* column-count row) column-count)]
     (strided-view-to-vector sub-view)))

  (setRow [this ^long row ^netlib_ccm.core.AbstractVector data]
    (let [^AbstractView strided-data data
          ^AbstractView row-data (.getRow this row)]
      (assign-strided-view! (.getStridedView row-data)
                            (.getStridedView strided-data))))

  (^long rowCount [this] row-count)
  (^long columnCount [this] column-count)
  (^netlib_ccm.core.AbstractMatrix clone [this]
   (let [retval (clone-strided-view (.getStridedView this))
         ary (->DenseVector (.data retval) 0 (* (.row-count data) (.column-count data)))]
     (->DenseMatrix ary (.row-count data) (.column-count data))))
  (validateShape [this] (validate-strided-view-shape (.getStridedView this)))

  clojure.lang.Seqable
  (seq [this] (map #(.getRow this %) (range (.rowCount this)))))

;;Interpret a submatrix as an array
(deftype StridedVector [^StridedView data]
  netlib_ccm.core.NetlibItem
  netlib_ccm.core.AbstractView
  (^StridedView getStridedView [this] data)

  netlib_ccm.core.AbstractVector
  (^double get [this ^long idx]
   (aget ^doubles (.data data) (get-strided-view-total-offset data idx)))

  (set [this ^long idx ^double val]
    (aset ^doubles (.data data) (get-strided-view-total-offset data idx) val))

  (^long length [this] (* (.row-count data) (.column-count data)))

  (^netlib_ccm.core.AbstractVector clone [this]
   (strided-view-to-vector (clone-strided-view (.getStridedView this))))

  (validateShape [this] (validate-strided-view-shape (.getStridedView this)))

  clojure.lang.Seqable
  (seq [this] (map #(.get this %) (range (.length this)))))


(defn get-column
  ^StridedVector [^AbstractMatrix mat ^long column]
  (let [^AbstractView mat mat
        ^StridedView view (.getStridedView mat)]
    (when-not (= (.column-count view) (.initial-column-count view))
      (throw (Exception. "Unsupported")))
    (->StridedVector (new-strided-view (.data view)
                                      (+ (.offset view) column)
                                      (.row-count view)
                                      1
                                      (.row-stride view)))))

(defn strided-view-to-vector
  ^AbstractVector [^StridedView view]
  (if (or (= 1 (.row-count view))
          (= (.row-stride view) (.column-count view)))
    (->DenseVector (.data view) (get-strided-view-total-offset view 0)
                   (get-strided-view-data-length view))
    (->StridedVector view)))


(defn strided-view-to-matrix
  "A strided view can be a matrix when it addresses a block of memory
with the same number of columns in each row"
  ^AbstractMatrix [^StridedView view ^long num-rows ^long num-cols]
  (when-not (and (= (.column-count view) (.initial-column-count view))
                 (or (= 0 (.last-column-count view))
                     (= (.column-count view) (.last-column-count view))))
    (throw (Exception. "Cannot make matrix out of offset view")))
  (if (= (.column-count view) (.row-stride view))
    (->DenseMatrix (->DenseVector (.data view) (.offset view)
                                  (* num-rows num-cols))
                   num-rows
                   num-cols)
    (->StridedMatrix view num-rows num-cols)))


(defn set-column
  [^AbstractMatrix mat ^long column ^AbstractVector data]
  (let [^AbstractView col-view (get-column mat column)
        ^StridedView column (.getStridedView col-view)
        ^StridedView data-view (.getStridedView ^AbstractView data)]
    (assign-strided-view! column data-view)))


(defn new-dense-vector-from-strided-view
  ^DenseVector [^StridedView view]
  (let [new-view (clone-strided-view view)]
    (->DenseVector (.data new-view) 0 (* (.row-count new-view) (.column-count new-view)))))


(defn new-dense-matrix-from-strided-view
  ^DenseMatrix [^StridedView view row-count column-count]
  (let [^DenseVector new-vec (new-dense-vector-from-strided-view view)
        ^DenseMatrix retval (->DenseMatrix new-vec row-count column-count)]
    (.validateShape retval)
    retval))


(defn in-range
  [^long item ^long min-eq ^long max]
  (and (>= item min-eq)
       (< max)))

(defn new-dense-vector
  ^DenseVector [^long length]
  (let [new-data (make-array Double/TYPE length)]
    (->DenseVector new-data 0 length)))

(defn new-dense-vector-from-array
  ^DenseVector [^doubles double-array]
  (let [len (count double-array)]
    (->DenseVector double-array 0 len)))

(defn double-array-from-data
  ^doubles [data]
  (double-array (ma/eseq data)))

(defn new-dense-vector-from-data
  ^DenseVector [data]
  (new-array-view-from-array (double-array-from-data data)))

(defn do-construct-matrix
  [data]
  (let [shape (ma/shape data)
        num-shape (count shape)]
    (case num-shape
      1 (new-dense-vector-from-array (double-array-from-data data))
      2 (->DenseMatrix (new-dense-vector-from-array
                        (double-array-from-data data))
                       (first shape)
                       (second shape)))))

(defn do-new-matrix-nd
  [shape]
  (case (count shape)
    1 (new-dense-vector (first shape))
    2 (->DenseMatrix (new-dense-vector (* (long (first shape)) (long (second shape)))
                                      (first shape) (second shape)))))


(extend-protocol mp/PImplementation
  NetlibItem
  (implementation-key [m] :netlib)
  (supports-dimensionality? [m dims] (in-range dims 0 3))
  (construct-matrix [m data] (do-construct-matrix data))
  (new-vector [m length] (new-dense-vector length))
  (new-matrix [m rows columns] (->DenseMatrix (new-dense-vector (* rows columns)) rows columns))
  (new-matrix-nd [m shape] (do-new-matrix-nd shape)))


(extend-protocol mp/PDimensionInfo
  AbstractVector
  (dimensionality [m] 1)
  (get-shape [m] [(.length m)])
  (is-scalar? [m] false)
  (is-vector? [m] true)
  (dimension-count [m dim] (if (= 0 dim) (.length m) (throw (Exception. "Unsupported"))))
  AbstractMatrix
  (dimensionality [m] 2)
  (get-shape [m] [(.rowCount m) (.columnCount m)])
  (is-scalar? [m] false)
  (is-vector? [m] false)
  (dimension-count [m dim]
    (case (long dim)
      0 (.rowCount m)
      1 (.columnCount m)
      (throw (Exception. "Unsupported")))))

(extend-protocol mp/PIndexedAccess
  AbstractVector
  (get-1d [m row] (.get m row))
  (get-2d [m row column] (throw (Exception. "Unsupported")))
  (get-nd [m indexes] (if (= 1 (count indexes))
                        (mp/get-1d m (first indexes))
                        (throw (Exception. "Unsupported"))))
  AbstractMatrix
  (get-1d [m row] (.getRow m row))
  (get-2d [m row column] (.get m row column))
  (get-nd [m indexes] (let [idx-count (count indexes)]
                        (case idx-count
                          1 (mp/get-1d m (first indexes))
                          2 (mp/get-2d m (first indexes) (second indexes))
                          (throw (Exception. "Unsupported"))))))

(extend-protocol mp/PIndexedSettingMutable
  AbstractVector
  (set-1d! [m row v] (.set m row v))
  (set-2d! [m row column v] (throw (Exception. "Unsupported")))
  (set-nd [m indexes v]
    (case (count indexes)
      1 (mp/set-1d! m (first indexes) v)
      (throw (Exception. "Unsupported"))))
  AbstractMatrix
  (set-1d! [m row v] (.setRow m row v))
  (set-2d! [m row column v] (.set m row column v))
  (set-nd! [m indexes v]
    (case (count indexes)
      1 (mp/set-1d! m (first indexes) v)
      2 (mp/set-2d! m (first indexes) (second indexes) v)
      (throw (Exception. "Unsupported")))))


(extend-protocol mp/PMatrixCloning
  AbstractVector
  (clone [m] (.clone m))
  AbstractMatrix
  (clone [m] (.clone m)))


(extend-protocol mp/PIndexedSetting
  AbstractVector
  (set-1d [m row v] (let [retval (mp/clone m)]
                      (mp/set-1d retval row v)))
  (set-2d [m row column v] (throw (Exception. "Unsupported")))
  (set-nd [m indexes v]
    (case (count indexes)
      1 (mp/set-1d m indexes v)
      (throw (Exception. "Unsupported"))))
  (is-mutable? [m] true)

  AbstractMatrix
  (set-1d [m row v] (let [retval (mp/clone m)]
                      (mp/set-1d! retval row v)))
  (set-2d [m row column v] (let [retval (mp/clone m)]
                             (mp/set-2d! retval row column v)))
  (set-nd [m indexes v]
    (case (count indexes)
      1 (mp/set-1d m (first indexes) v)
      2 (mp/set-2d m (first indexes) (second indexes) v)
      (throw (Exception. "Unsupported"))))
  (is-mutable? [m] true))


(extend-protocol mp/PTypeInfo
  NetlibItem
  (element-type [m] Double/TYPE))


(extend-protocol mp/PValidateShape
  AbstractVector
  (validate-shape [m] (let [^AbstractView v m]
                        (validate-strided-view-shape (.getStridedView v))
                        [(.length m)]))
  AbstractMatrix
  (validate-shape [m] (let [^AbstractView v m]
                        (validate-strided-view-shape (.getStridedView v))
                        [(.rowCount m) (.columnCount m)])))


(extend-protocol mp/PMutableMatrixConstruction
  AbstractVector
  (mutable-matrix [m] (.clone ^AbstractVector m))
  AbstractMatrix
  (mutable-matrix [m] (.clone ^AbstractMatrix m)))


(extend-protocol mp/PMutableCoercion
  NetlibItem
  (ensure-mutable [m] m))


(defn dense-coerce-vec
  ^DenseVector [^AbstractVector item]
  (.clone item))


(defn dense-coerce-mat
  ^DenseMatrix [^AbstractMatrix item]
  (.clone item))


(extend-protocol mp/PDense
  NetlibItem
  (dense-coerce [m data] (new-dense-vector-from-data data))

  DenseVector
  (dense [m] m)
  DenseMatrix
  (dense [m] m)
  AbstractView
  (dense [m] (new-dense-vector-from-strided-view (.getStridedView ^AbstractView m))))


(extend-protocol mp/PConversion
  AbstractVector
  (convert-to-nested-vectors [m] (vec (seq m)))
  AbstractMatrix
  (convert-to-nested-vectors [m] (mapv mp/convert-to-nested-vectors (seq m))))


(extend-protocol mp/PReshaping
  AbstractView
  (reshape [m shape] (let [^DenseVector m (new-dense-vector-from-strided-view
                                           (.getStridedView ^AbstractView m))
                           num-desired (long (apply * shape))
                           data-len (.length m)
                           m (->DenseVector (.data m) 0 num-desired)]

                       (when (> num-desired data-len)
                         (throw (Exception. "Attempt to reshape to a larger backing storage")))
                       (case (count shape)
                         1 m
                         2 (->DenseMatrix m (first shape) (second shape))
                         (throw (Exception. "Unsupported"))))))


(extend-protocol mp/PPack
  DenseVector
  (pack [m] m)
  DenseMatrix
  (pack [m] m)
  AbstractVector
  (pack [m] (.clone ^AbstractVector m))
  AbstractMatrix
  (pack [m] (.clone ^AbstractMatrix m)))


(extend-protocol PMatrixSlices
  AbstractMatrix
  (get-row [m i] (.getRow ^AbstractMatrix m i))
  (get-column [m i] (get-column ^AbstractMatrix m i))
  (get-major-slice [m i] (.getRow ^AbstractMatrix m i))
  (get-slice [m dimension i]
    (case (long dimension)
      0 (mp/get-row m i)
      1 (mp/get-column m i)
      (throw (Exception. "Unsupported")))))


(extend-protocol mp/PMatrixRows
  AbstractMatrix
  (get-rows [m] (seq m)))


(extend-protocol mp/PMatrixColumns
  AbstractMatrix
  (get-columns [m] (map #(get-column m %) (range (.columnCount ^AbstractMatrix m)))))


(extend-protocol mp/PSubVector
  AbstractView
  (subvector [m start length]
    (let [^AbstractView view m
          ^StridedView data (.getStridedView view)]
      (strided-view-to-vector (create-sub-strided-view view start length)))))


(extend-protocol mp/PSubMatrix
  AbstractMatrix
  (submatrix [m dim-ranges]
    (let [^AbstractMatrix m m
          ^AbstractView mat-view m
          ^StridedView view (.getStridedView mat-view)
          num-dims (count dim-ranges)
          num-data-items (get-strided-view-data-length view)]
      (when-not (= 2 num-dims)
        (throw (Exception. "Number of dim ranges must be 2")))
      (when-not (and (= (.column-count view) (.initial-column-count view))
                     (= 0 (rem num-data-items (.columnCount m))))
        (throw (Exception. "Submatric views on offset matrixes are not supported")))
      (let [[[start-row num-rows] [start-col num-cols]] dim-ranges
            data-start-offset (+ (* (.columnCount m) start-row) start-col)]
        (when (or (> (+ start-row num-rows)
                     (.rowCount m))
                  (> (+ start-col num-cols)
                     (.columnCount m)))
          (throw (Exception. "Attempt to access outside of matrix")))
        (let [new-view (new-strided-view (.data view)
                                         (get-strided-view-total-offset view data-start-offset)
                                         num-rows
                                         num-cols
                                         (.row-stride view)
                                         num-cols
                                         0)]
          (strided-view-to-matrix new-view num-rows num-cols))))))
