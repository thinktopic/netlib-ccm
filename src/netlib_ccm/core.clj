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

(defn get-strided-view-row-length-from-row
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


(defn get-strided-view-data-offset-from-row
  ^long [^StridedView view ^long row-idx]
  (check-strided-view-row-index view row-idx)
  (if (= 0 row-idx)
    0
    (+ (.initial-column-count view) (* (.column-count view) (- row-idx 1)))))


(defn get-strided-view-row-length-from-data-offset
  ^long [^StridedView view ^long data-offset]
  (check-strided-view-data-offset view data-offset)
  (if (< data-offset (.initial-column-count view))
    (- (.initial-column-count view) data-offset)
    (let [rest-offset (- data-offset (.initial-column-count view))
          row-idx (+ 1 (quot rest-offset (.column-count view)))
          rest-leftover (rem rest-offset (.column-count view))]
      (if (= row-idx (- (.row-count view) 1))
        (- (.last-column-count view) rest-leftover)
        (- (.column-count view) rest-leftover)))))


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
    (let [first-row-len (if (< offset (.initial-column-count view))
                          (- (.initial-column-count view) offset)
                          (- (.column-count view)
                             (rem (- offset (.initial-column-count view))
                                  (.column-count view))))
          first-row-offset (- (.column-count view) first-row-len)]
      (if (<= length first-row-len)
        (new-strided-view (.data view) (.get-strided-view-total-offset view offset) length)
        (let [row-idx (get-strided-view-row-from-data-offset view offset)
              last-row-len (rem (+ first-row-offset length) (.column-count view))
              ary-offset (+ (.offset view) (* (.row-stride view) row-idx))
              body-num-rows (quot (- length (+ first-row-len last-row-len))
                                  (.column-count view))
              ini-row (long (if-not (= 0 first-row-len) 1 0))
              end-row (long (if-not (= 0 last-row-len) 1 0))]
          (println "row-idx" row-idx first-row-len first-row-offset last-row-len)
          (new-strided-view (.data view) ary-offset
                            (+ body-num-rows ini-row end-row)
                            (.column-count view)
                            (.row-stride view)
                            first-row-len
                            last-row-len))))))

(defn strided-op
  "Perform an operation such as (assign! lhs rhs).
  Op gets passed: lhs-double-array lhs-offset rhs-double-array rhs-offset op-amount"
  ([^StridedView lhs ^StridedView rhs op]
   (loop [rhs-row 0]
     (when (< rhs-row (.row-count rhs))
       (let [data-offset (get-strided-view-data-offset-from-row rhs rhs-row)
             rhs-row-len (get-strided-view-row-length-from-row rhs rhs-row)]
         (loop [rhs-row-offset 0]
           (when (< rhs-row-offset rhs-row-len)
             (let [data-offset (+ data-offset rhs-row-offset)
                   rhs-total-offset (get-strided-view-total-offset rhs data-offset)
                   rhs-row-len (- rhs-row-len rhs-row-offset)
                   lhs-row (get-strided-view-row-from-data-offset lhs data-offset)
                   lhs-row-len (get-strided-view-row-length-from-data-offset lhs data-offset)
                   lhs-total-offset (get-strided-view-total-offset lhs data-offset)
                   _ (println data-offset rhs-row-len lhs-row-len lhs-row rhs-row)
                   op-amount (min lhs-row-len rhs-row-len)]
               (op (.data lhs) lhs-total-offset
                   (.data rhs) rhs-total-offset
                   op-amount)
               (recur (+ rhs-row-offset op-amount))))))
       (recur (inc rhs-row)))))
  ([^StridedView lhs op]
   (loop [lhs-row 0]
     (when (< lhs-row (.row-count lhs))
       (let [row-len (get-strided-view-row-length-from-row lhs lhs-row)]
         (when (> row-len 0)
           (let [data-offset (get-strided-view-data-offset-from-row lhs lhs-row)
                 lhs-total-offset (get-strided-view-total-offset lhs data-offset)]
             (op (.data lhs) lhs-total-offset row-len))))
       (recur (inc lhs-row))))))



(defn assign-strided-view!
  "lhs must be at least as large as rhs.  Rhs.length must be an even multiple
of lhs.length"
  [^StridedView lhs ^StridedView rhs]
  (let [lhs-len (get-strided-view-data-length lhs)
        rhs-len (get-strided-view-data-length rhs)]
    ;;no op
    (when-not (= 0 rhs-len)
     (when-not (= 0 (rem lhs-len rhs-len))
       (throw (Exception. "Strided assignment: rhs-len not even multiple of lhs-len")))
     (if (= 1 rhs-len)
       (let [rhs-val (aget ^doubles (.data rhs) (get-strided-view-total-offset rhs 0))]
         (strided-op lhs (fn [^doubles data ^long offset ^long len]
                           (java.util.Arrays/fill data offset (+ offset len) rhs-val))))
       (let [num-reps (quot lhs-len rhs-len)
             copy-op (fn [lhs-data lhs-offset rhs-data rhs-offset op-len]
                       (println lhs-offset rhs-offset op-len)
                       (System/arraycopy ^doubles rhs-data ^long rhs-offset
                                         ^doubles lhs-data ^long lhs-offset
                                         ^long op-len))]
         (if (= 1 num-reps)
           (strided-op lhs rhs copy-op)
           (loop [rep 0]
             (when (< rep num-reps)
               (strided-op (create-sub-strided-view lhs (* rep rhs-len) rhs-len)
                           rhs copy-op)
               (recur (inc rep))))))))))


(defn clone-strided-view
  "Create a packed strided view from potentially non-dense view"
  ^StridedView [^StridedView source]
  (let [num-items (get-strided-view-data-length source)
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

(defn vector-iterator
  [^AbstractVector vec]
  (let [idx (atom 0)
        len (.length vec)]
    (reify java.util.Iterator
      (hasNext [this] (< @idx len))
      (next [this] (let [current @idx]
                     (swap! idx inc)
                     (.get vec current)))
      (remove [this] (throw (Exception. "Unsupported"))))))

(defn matrix-iterator
  [^AbstractMatrix mat]
  (let [idx (atom 0)
        len (.rowCount mat)]
    (reify java.util.Iterator
      (hasNext [this] (< @idx len))
      (next [this] (let [current @idx]
                     (swap! idx inc)
                     (.getRow mat current)))
      (remove [this] (throw (Exception. "Unsupported"))))))

(declare ->DenseVector ->DenseMatrix)

(defn vec-to-persistent-vec
  [^AbstractVector avec]
  (vec (seq avec)))

(defn mat-to-persistent-vec
  [^AbstractMatrix amat]
  (mapv (comp vec seq) amat))

(defn print-vector
  [^AbstractVector avec]
  (if (< (.length avec) 1000)
    (pr-str (vec-to-persistent-vec avec))
    (format "<vector of length %d>" (.length avec))))

(defn print-matrix
  [^AbstractMatrix amat]
  (if (< (* (.rowCount amat) (.columnCount amat)) 1000)
    (clojure.string/trim (with-out-str (clojure.pprint/pprint (mat-to-persistent-vec amat))))
    (format "<matrix of %dx%d>" (.rowCount amat) (.columnCount amat))))

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
  (seq [this] (map #(.get this %) (range length)))

  java.lang.Object
  (toString [this] (print-vector this))

  java.lang.Iterable
  (iterator [this] (vector-iterator this)))


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
  (seq [this] (map #(.getRow this %) (range (.rowCount this))))

  java.lang.Object
  (toString [this] (print-matrix this))

  java.lang.Iterable
  (iterator [this] (matrix-iterator this)))

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
         item-count (* row-count column-count)
         ary (->DenseVector (.data retval) 0 item-count)]
     (->DenseMatrix ary row-count column-count)))
  (validateShape [this] (validate-strided-view-shape (.getStridedView this)))

  clojure.lang.Seqable
  (seq [this] (map #(.getRow this %) (range (.rowCount this))))

  java.lang.Object
  (toString [this] (print-matrix this))

  java.lang.Iterable
  (iterator [this] (matrix-iterator this))
  )

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

  (^long length [this] (get-strided-view-data-length data))

  (^netlib_ccm.core.AbstractVector clone [this]
   (strided-view-to-vector (clone-strided-view (.getStridedView this))))

  (validateShape [this] (validate-strided-view-shape (.getStridedView this)))

  clojure.lang.Seqable
  (seq [this] (map #(.get this %) (range (.length this))))

  java.lang.Object
  (toString [this] (print-vector this))

  java.lang.Iterable
  (iterator [this] (matrix-iterator this)))


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
          (and (= (.row-stride view) (.column-count view))
               (= (.initial-column-count view) (.column-count view))
               (or (= (.last-column-count view) (.column-count view))
                   (= (.last-column-count view) 0))))
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
  (new-dense-vector-from-array (double-array-from-data data)))

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


(extend-protocol mp/PMatrixSlices
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
      (strided-view-to-vector (create-sub-strided-view data start length)))))


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
        (let [row-stride (if (= 1 (.row-count view))
                           (.columnCount m)
                           (.row-stride view))
              new-view (new-strided-view (.data view)
                                         (get-strided-view-total-offset view data-start-offset)
                                         ;;Add one row to account for the 0 last-column-count
                                         (+ num-rows 1)
                                         num-cols
                                         row-stride
                                         num-cols
                                         0)]
          (strided-view-to-matrix new-view num-rows num-cols))))))

;;double dispatch on type of source
(defprotocol PAbstractViewAssign
  (assign-source-to-view! [source m]))


(extend-protocol PAbstractViewAssign
  (Class/forName "[D")
  (assign-source-to-view! [source m]
    (let [^doubles source source
          ^AbstractView m m]
      (assign-strided-view! (.getStridedView m)
                            (new-strided-view source 0 (count source)))))

  AbstractView
  (assign-source-to-view! [source m]
    (let [^AbstractView source source
          ^AbstractView m m]
      (assign-strided-view! (.getStridedView m) (.getStridedView source))))

  clojure.lang.PersistentVector
  (assign-source-to-view! [source m]
    (let [^AbstractView m m
          ^doubles ddata (double-array-from-data source)]
      (assign-source-to-view! ddata m)))

  clojure.lang.ISeq
  (assign-source-to-view! [source m]
    (let [^AbstractView m m
          ^doubles ddata (double-array-from-data source)]
      (assign-source-to-view! ddata m)))

  Double
  (assign-source-to-view! [source m]
    (let [^double source source
          ^AbstractView m m]
      (strided-op (.getStridedView m)
                  (fn [^doubles data ^long offset ^long len]
                    (java.util.Arrays/fill data offset (+ offset len) source))))))


(extend-protocol mp/PAssignment
  AbstractView
  (assign! [m source] (assign-source-to-view! source m))
  (assign-array! [m arr] (mp/assign-array! m arr 0 (count arr)))
  (assign-array! [m arr offset len]
    (let [^AbstractView m m]
      (assign-strided-view! (.getStridedView m)
                            (new-strided-view arr offset len)))))


(extend-protocol mp/PMutableFill
  AbstractView
  (fill! [m value] (assign-source-to-view! value m)))


(defn is-strided-view-dense?
  "A dense strided view is one that has no offset and completely
fills its backing store"
  [^StridedView view]
  (and (= 0 (.offset view))
       (= (.column-count view) (.row-stride view))
       (= (.column-count view) (.initial-column-count view))
       (or (= (.column-count view) (.last-column-count view))
           (= 0 (.last-column-count view)))
       (= (count (.data view))
          (* (.row-count view) (.column-count view)))))


(extend-protocol mp/PDoubleArrayOutput
  AbstractView
  (to-double-array [m]
    (let [^AbstractView m m
          ^StridedView view (.getStridedView m)]
      (if (is-strided-view-dense? view)
        (.data view)
        (.data (clone-strided-view view)))))
  (as-double-array [m]
    (let [^AbstractView m m
          ^StridedView view (.getStridedView m)]
      (when (is-strided-view-dense? view)
        (.data view)))))


(extend-protocol mp/PVectorView
  AbstractView
  (as-vector [m]
    (let [^AbstractView m m]
      (strided-view-to-vector (.getStridedView m)))))


(def empty-vec (new-dense-vector 0))
(mi/register-implementation empty-vec)
