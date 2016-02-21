(ns netlib-ccm.core
  (:require [clojure.core.matrix :as ma]
            [clojure.core.matrix.protocols :as mp]
            [clojure.core.matrix.impl.mathsops :as mathsops]
            [clojure.core.matrix.implementations :as mi]
            [clojure.pprint])
  (:import [com.github.fommil.netlib BLAS]
           [netlib_ccm Ops IUnaryOp IBinaryOp
            StridedBuffer IStridedUnaryOp IStridedBinaryOp
            IStridedTernaryOp]))


(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(definterface NetlibItem)

(definterface AbstractView
  (^netlib_ccm.StridedBuffer getStridedBuffer []))

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
      (hasNext [this] (< ^long @idx len))
      (next [this] (let [current @idx]
                     (swap! idx inc)
                     (.get vec current)))
      (remove [this] (throw (Exception. "Unsupported"))))))

(defn matrix-iterator
  [^AbstractMatrix mat]
  (let [idx (atom 0)
        len (.rowCount mat)]
    (reify java.util.Iterator
      (hasNext [this] (< ^long @idx len))
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

(defn validate-strided-buffer-shape [^StridedBuffer data]
  (let [data-len (count (.data data))
        item-end (+ (.offset data)
                    (* (- (.row-count data) 1) (.row-stride data))
                    (.column-count data))]
    (when (> item-end data-len)
      (throw (Exception. "Invalid strided view")))))


(deftype DenseVector [^doubles data ^long offset ^long length]
  netlib_ccm.core.NetlibItem
  netlib_ccm.core.AbstractView
  (^StridedBuffer getStridedBuffer [this] (StridedBuffer. data offset length))

  netlib_ccm.core.AbstractVector
  (^double get [this ^long idx] (aget data (+ offset idx)))
  (set [this ^long idx ^double val] (aset data (+ offset idx) val))
  (^long length [this] length)
  (^netlib_ccm.core.AbstractVector clone [this]
   (let [^StridedBuffer retval (.clone (.getStridedBuffer this))]
     (->DenseVector (.data retval) 0 length)))
  (validateShape [this] (validate-strided-buffer-shape (.getStridedBuffer this)))

  clojure.lang.Seqable
  (seq [this] (map #(.get this %) (range length)))

  java.lang.Object
  (toString [this] (print-vector this))

  java.lang.Iterable
  (iterator [this] (vector-iterator this)))


(deftype DenseMatrix [^DenseVector data ^long row-count ^long column-count]
  netlib_ccm.core.NetlibItem
  netlib_ccm.core.AbstractView
  (^StridedBuffer getStridedBuffer [this]
   (let [item-count (* row-count column-count)]
     (StridedBuffer. (.data data) (.offset data) item-count)))
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
      (.assign (.getStridedBuffer row-data)
                            (.getStridedBuffer strided-data))))
  (^long rowCount [this] row-count)
  (^long columnCount [this] column-count)
  (^netlib_ccm.core.AbstractMatrix clone [this]
   (let [retval (.clone (.getStridedBuffer this))
         ary (->DenseVector (.data retval) 0 (* row-count column-count))]
     (->DenseMatrix ary row-count column-count)))
  (validateShape [this] (validate-strided-buffer-shape (.getStridedBuffer this)))
  clojure.lang.Seqable
  (seq [this] (map #(.getRow this %) (range (.rowCount this))))

  java.lang.Object
  (toString [this] (print-matrix this))

  java.lang.Iterable
  (iterator [this] (matrix-iterator this)))

(declare strided-view-to-vector)

;;A strided view of data where you can have N contiguous elements separated in rows of length Y.
;;Think of a small submatrix in a larger one
(deftype StridedMatrix [^StridedBuffer data ^long row-count ^long column-count]
  netlib_ccm.core.NetlibItem
  netlib_ccm.core.AbstractView
  (^StridedBuffer getStridedBuffer [this] data)

  netlib_ccm.core.AbstractMatrix
  (^double get [this ^long row ^long column]
   (aget ^doubles (.data data)
         (StridedBuffer/getTotalOffset data (+ (* column-count row) column))))

  (set [this ^long row ^long column ^double val]
    (aset ^doubles (.data data)
          (StridedBuffer/getTotalOffset data (+ (* column-count row) column))
          val))

  (^netlib_ccm.core.AbstractVector getRow [this ^long row]
   (let [sub-view (.createSubStridedBuffer data (* column-count row) column-count)]
     (strided-view-to-vector sub-view)))

  (setRow [this ^long row ^netlib_ccm.core.AbstractVector data]
    (let [^AbstractView strided-data data
          ^AbstractView row-data (.getRow this row)]
      (.assign (.getStridedBuffer row-data)
                            (.getStridedBuffer strided-data))))

  (^long rowCount [this] row-count)
  (^long columnCount [this] column-count)
  (^netlib_ccm.core.AbstractMatrix clone [this]
   (let [retval (.clone (.getStridedBuffer this))
         item-count (* row-count column-count)
         ary (->DenseVector (.data retval) 0 item-count)]
     (->DenseMatrix ary row-count column-count)))
  (validateShape [this] (validate-strided-buffer-shape (.getStridedBuffer this)))

  clojure.lang.Seqable
  (seq [this] (map #(.getRow this %) (range (.rowCount this))))

  java.lang.Object
  (toString [this] (print-matrix this))

  java.lang.Iterable
  (iterator [this] (matrix-iterator this))
  )

;;Interpret a submatrix as an array
(deftype StridedVector [^StridedBuffer data]
  netlib_ccm.core.NetlibItem
  netlib_ccm.core.AbstractView
  (^StridedBuffer getStridedBuffer [this] data)

  netlib_ccm.core.AbstractVector
  (^double get [this ^long idx]
   (aget ^doubles (.data data) (StridedBuffer/getTotalOffset data idx)))

  (set [this ^long idx ^double val]
    (aset ^doubles (.data data) (StridedBuffer/getTotalOffset data idx) val))

  (^long length [this] (StridedBuffer/getDataLength data))

  (^netlib_ccm.core.AbstractVector clone [this]
   (strided-view-to-vector (.clone (.getStridedBuffer this))))

  (validateShape [this] (validate-strided-buffer-shape (.getStridedBuffer this)))

  clojure.lang.Seqable
  (seq [this] (map #(.get this %) (range (.length this))))

  java.lang.Object
  (toString [this] (print-vector this))

  java.lang.Iterable
  (iterator [this] (matrix-iterator this)))

(defn strided-view-to-vector
  ^AbstractVector [^StridedBuffer view]
  (if (StridedBuffer/isDense view)
    (->DenseVector (.data view) (StridedBuffer/getTotalOffset view 0)
                   (StridedBuffer/getDataLength view))
    (->StridedVector view)))

(defn strided-view-to-matrix
  "A strided view can be a matrix when it addresses a block of memory
with the same number of columns in each row"
  ^AbstractMatrix [^StridedBuffer view ^long num-rows ^long num-cols]
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


(defn get-submatrix
  [^AbstractMatrix m start-row num-rows start-col num-cols]
  (let [^AbstractView mat-view m
        ^StridedBuffer view (.getStridedBuffer mat-view)
        num-data-items (StridedBuffer/getDataLength view)]
    (when-not (and (= (.column-count view) (.initial-column-count view))
                   (= 0 (rem num-data-items (.columnCount m))))
      (throw (Exception. "Submatric views on offset matrixes are not supported")))
    (let [^long start-row start-row
          ^long start-col start-col
          ^long num-rows num-rows
          ^long num-cols num-cols
          data-start-offset (+ (* (.columnCount m) start-row) start-col)]
        (when (or (> (+ start-row num-rows)
                     (.rowCount m))
                  (> (+ start-col num-cols)
                     (.columnCount m)))
          (throw (Exception. "Attempt to access outside of matrix")))
        (let [row-stride (if (= 1 (.row-count view))
                           (.columnCount m)
                           (.row-stride view))
              new-view (StridedBuffer. (.data view)
                                         (StridedBuffer/getTotalOffset view data-start-offset)
                                         ;;Add one row to account for the 0 last-column-count
                                         (+ num-rows 1)
                                         num-cols
                                         row-stride
                                         num-cols
                                         0)]
          (strided-view-to-matrix new-view num-rows num-cols)))))


(defn get-column
  ^StridedVector [^AbstractMatrix mat ^long column]
  (mp/as-vector (get-submatrix mat 0 (.rowCount mat) column 1)))


(defn set-column
  [^AbstractMatrix mat ^long column ^AbstractVector data]
  (let [^AbstractView col-view (get-column mat column)
        ^StridedBuffer column (.getStridedBuffer col-view)
        ^StridedBuffer data-view (.getStridedBuffer ^AbstractView data)]
    (.assign column data-view)))


(defn new-dense-vector-from-strided-view
  ^DenseVector [^StridedBuffer view]
  (let [new-view (.clone view)]
    (->DenseVector (.data new-view) 0 (* (.row-count new-view) (.column-count new-view)))))


(defn new-dense-matrix-from-strided-view
  ^DenseMatrix [^StridedBuffer view row-count column-count]
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


(defn new-dense-matrix
  ^DenseMatrix [^long row-count ^long column-count]
  (->DenseMatrix (new-dense-vector (* row-count column-count))
                 row-count column-count))

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
  (if (ma/scalar? shape)
    (new-dense-vector shape)
    (case (count shape)
      0 (new-dense-vector 1)
      1 (new-dense-vector (first shape))
      2 (new-dense-matrix (first shape) (second shape)))))

(defn view-as-dense-matrix
  ^DenseMatrix [^AbstractView v ^long num-rows ^long num-columns]
  (let [^StridedBuffer view (.getStridedBuffer v)]
    (when-not (StridedBuffer/isDense view)
      (throw (Exception. "Unsupported for non dense views")))
    (let [item-count (* num-rows num-columns)
          data-len (StridedBuffer/getDataLength view)]
      (when-not (= item-count data-len)
        (throw (Exception. (format "Unsupported as item-count data-len mismatch %d vs %d"
                                   item-count data-len))))
      (->DenseMatrix (strided-view-to-vector view) num-rows num-columns))))


(defprotocol ToNetlibType
  (to-netlib [item]))

(extend-protocol ToNetlibType
  (Class/forName "[D")
  (to-netlib [item] (->DenseVector item 0 (count item)))
  clojure.lang.PersistentVector
  (to-netlib [item] (do-construct-matrix item))
  AbstractView
  (to-netlib [item] item)
  clojure.lang.ISeq
  (to-netlib [item] (do-construct-matrix item))
  Number
  (to-netlib [item] (double item)))


(extend-protocol mp/PImplementation
  NetlibItem
  (implementation-key [m] :netlib)
  (supports-dimensionality? [m dims] (in-range dims 0 3))
  (construct-matrix [m data] (do-construct-matrix data))
  (new-vector [m length] (new-dense-vector length))
  (new-matrix [m ^long rows ^long columns]
    (->DenseMatrix (new-dense-vector (* rows columns)) rows columns))
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

(defn set-zero-d
  [^AbstractView view scalar]
  (mp/assign! view scalar))

(extend-protocol mp/PIndexedSettingMutable
  AbstractVector
  (set-1d! [m row v] (.set m row v))
  (set-2d! [m row column v] (throw (Exception. "Unsupported")))
  (set-nd! [m indexes v]
    (case (count indexes)
      0 (set-zero-d m v)
      1 (mp/set-1d! m (first indexes) v)
      (throw (Exception. "Unsupported"))))

  AbstractMatrix
  (set-1d! [m row v] (.setRow m row v))
  (set-2d! [m row column v] (.set m row column v))
  (set-nd! [m indexes v]
    (case (count indexes)
      0 (set-zero-d m v)
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
                        (validate-strided-buffer-shape (.getStridedBuffer v))
                        [(.length m)]))
  AbstractMatrix
  (validate-shape [m] (let [^AbstractView v m]
                        (validate-strided-buffer-shape (.getStridedBuffer v))
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
  AbstractView
  (dense [m]
    (let [^StridedBuffer new-data (.clone (.getStridedBuffer ^AbstractView m))
          dense-vec (->DenseVector (.data new-data) 0 (StridedBuffer/getDataLength new-data))]
      (if (instance? AbstractVector m)
        dense-vec
        (let [^AbstractMatrix m m]
          (->DenseMatrix dense-vec (.rowCount m) (.columnCount m))))))

  (dense-coerce [m data] (new-dense-vector-from-data data)))


(extend-protocol mp/PConversion
  AbstractVector
  (convert-to-nested-vectors [m] (vec (seq m)))
  AbstractMatrix
  (convert-to-nested-vectors [m] (mapv mp/convert-to-nested-vectors (seq m))))


(extend-protocol mp/PReshaping
  AbstractView
  (reshape [m shape] (let [^DenseVector m (new-dense-vector-from-strided-view
                                           (.getStridedBuffer ^AbstractView m))
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
          ^StridedBuffer data (.getStridedBuffer view)]
      (strided-view-to-vector (.createSubStridedBuffer data start length)))))


(extend-protocol mp/PSubMatrix
  AbstractMatrix
  (submatrix [m dim-ranges]
    (when-not (= 2 (count dim-ranges))
      (throw (Exception. "Must have exactly 2 ranges for submatrx")))
    (let [[[ start-row num-rows] [start-col num-cols]] dim-ranges]
      (get-submatrix m start-row num-rows start-col num-cols))))


;;double dispatch on type of source
(defprotocol PAbstractViewAssign
  (assign-source-to-view! [source m]))


(extend-protocol PAbstractViewAssign
  (Class/forName "[D")
  (assign-source-to-view! [source m]
    (let [^doubles source source
          ^AbstractView m m]
      (.assign (.getStridedBuffer m)
                            (StridedBuffer. source 0 (count source)))
      m))

  AbstractView
  (assign-source-to-view! [source m]
    (let [^AbstractView source source
          ^AbstractView m m]
      (when-not (identical? source m)
        (.assign (.getStridedBuffer m) (.getStridedBuffer source)))
      m))

  clojure.lang.PersistentVector
  (assign-source-to-view! [source m]
    (let [^AbstractView m m
          ^doubles ddata (double-array-from-data source)]
      (assign-source-to-view! ddata m)
      m))

  clojure.lang.ISeq
  (assign-source-to-view! [source m]
    (let [^AbstractView m m
          ^doubles ddata (double-array-from-data source)]
      (assign-source-to-view! ddata m)
      m))

  Number
  (assign-source-to-view! [source m]
    (let [source (double source)
          ^AbstractView m m]
      (if (= 0.0 source)
        (mp/scale! m 0.0)
        (StridedBuffer/unaryOperation (.getStridedBuffer m)
                    (fn [^doubles data ^long offset ^long len]
                      (java.util.Arrays/fill data offset (+ offset len) source))))
      m)))


(extend-protocol mp/PAssignment
  AbstractView
  (assign! [m source] (assign-source-to-view! source m))
  (assign-array! [m arr] (mp/assign-array! m arr 0 (count arr)))
  (assign-array! [m arr offset len]
    (let [^AbstractView m m]
      (.assign (.getStridedBuffer m)
               (StridedBuffer. arr offset len)))))


(extend-protocol mp/PMutableFill
  AbstractView
  (fill! [m value] (assign-source-to-view! value m)))


(extend-protocol mp/PDoubleArrayOutput
  AbstractView
  (to-double-array [m]
    (let [^AbstractView m m
          ^StridedBuffer view (.getStridedBuffer m)]
      (if (StridedBuffer/isCompleteDense view)
        (.data view)
        (.data (.clone view)))))
  (as-double-array [m]
    (let [^AbstractView m m
          ^StridedBuffer view (.getStridedBuffer m)]
      (when (StridedBuffer/isCompleteDense view)
        (.data view)))))


(extend-protocol mp/PVectorView
  AbstractView
  (as-vector [m]
    (let [^AbstractView m m]
      (strided-view-to-vector (.getStridedBuffer m)))))

(defn clone-abstract-view
  [^AbstractView view]
  (if (instance? AbstractMatrix view)
    (.clone ^AbstractMatrix view)
    (.clone ^AbstractVector view)))


(defn new-abstract-view
  [^AbstractView view]
  (if (instance? AbstractMatrix view)
    (let [^AbstractMatrix view view]
      (new-dense-matrix (.rowCount view) (.columnCount view)))
    (new-dense-vector (.length ^AbstractVector view))))


(defmacro unary-op-macro!
  [view & body]
  `(let [^"StridedBuffer" view# (.getStridedBuffer ^AbstractView ~view)
         ^"IUnaryOp" local-op# (reify netlib_ccm.IUnaryOp
                                 (op[this# ~'lhs-value]
                                   ~@body))]
     (StridedBuffer/unaryOperationOp view# local-op#)))

(defn unary-op!
  [^AbstractView view op]
  (unary-op-macro! view (op lhs-value)))


;;All unary operators
(defn unary-immutable-op
  [^AbstractView view op]
  (let [new-view (clone-abstract-view view)]
    (unary-op! new-view op)
    new-view))


(defn binary-op!
  "lhs gets (op lhs rhs)"
  [^AbstractView lhs ^AbstractView rhs bin-op]
  (let [^StridedBuffer rhs-view (.getStridedBuffer rhs)
        ^StridedBuffer lhs-view (.getStridedBuffer lhs)]
    (StridedBuffer/binaryOperationOp lhs-view rhs-view
                                      (reify IBinaryOp
                                                   (op [this lhs-val rhs-val]
                                                     (bin-op lhs-val rhs-val))))))


(defn binary-immutable-op
  [^AbstractView lhs ^AbstractView rhs op]
  (let [^AbstractView lhs-clone (clone-abstract-view lhs)]
    (binary-op! lhs-clone rhs op)
    lhs-clone))

(defprotocol SourceViewOp
  (source-view-op [source view op])
  (source-view-op! [source view op]))

(extend-protocol SourceViewOp
  (Class/forName "[D")
  (source-view-op [source view op]
    (let [^AbstractView view view]
      (binary-immutable-op view (->DenseVector source 0 (count source)))))
  (source-view-op! [source view op]
    (let [^AbstractView view view]
      (binary-op! view (->DenseVector source 0 (count source)) op)))

  AbstractView
  (source-view-op [source view op] (binary-immutable-op view source op))
  (source-view-op! [source view op] (binary-op! view source op))

  Number
  (source-view-op [source view op] (let [^AbstractView retval (clone-abstract-view view)]
                                     (unary-op! retval #(op % (double source)))
                                     retval))
  (source-view-op! [source view op] (unary-op! view #(op % (double source)))))


(extend-protocol mp/PFunctionalOperations
  AbstractView
  (element-seq [m] (let [view (.getStridedBuffer ^AbstractView m)]
                     (map #(aget ^doubles
                                 (.data view)
                                 (StridedBuffer/getTotalOffset view %))
                          (range (StridedBuffer/getDataLength view)))))

  (element-map!
    ([m f] (do (unary-op! m f) m))
    ([m f a] (do (source-view-op! (to-netlib a) m f) m))
    ([m f a more] (let [^StridedBuffer view (.getStridedBuffer ^AbstractView m)
                        num-items (StridedBuffer/getDataLength view)
                        rest-list (mapv to-netlib (concat a more))]
                    (loop [idx 0]
                      (when (< idx num-items)
                        (let [m-offset (StridedBuffer/getTotalOffset view idx)
                              double-val (double (apply f (aget ^doubles (.data view) m-offset)
                                                        (map #(mp/get-1d % idx) rest-list)))]
                          (aset ^doubles (.data view) m-offset double-val))
                        (recur (inc idx))))
                    m)))

  (element-map
    ([m f] (let [clone-view (clone-abstract-view m)]
             (mp/element-map! clone-view f)))
    ([m f a] (let [clone-view (clone-abstract-view m)]
               (mp/element-map! clone-view f a)))
    ([m f a more] (let [clone-view (clone-abstract-view m)]
                    (mp/element-map! clone-view f a more))))


  (element-reduce
    ([m f] (reduce f (mp/element-seq m)))
    ([m f init] (reduce f init (mp/element-seq m)))))


(defn dot-abstract-views
  ^double [^AbstractView lhs ^AbstractView rhs]
  (let [^StridedBuffer lhs-view (.getStridedBuffer lhs)
        ^StridedBuffer rhs-view (.getStridedBuffer rhs)
        ^doubles accum (make-array Double/TYPE 1)]
    (StridedBuffer/binaryOperation
     lhs-view rhs-view
     (reify IStridedBinaryOp
       (op [this lhs-data lhs-offset rhs-data rhs-offset op-len]
         (let [accum-val (aget accum 0)
               dot-val (.ddot (BLAS/getInstance) op-len
                              ^doubles lhs-data lhs-offset
                              ^doubles rhs-data rhs-offset)]
           (aset accum 0 (+ dot-val accum-val))))))
    (aget accum 0)))


(defprotocol AbstractViewVectorOps
  (abstract-view-dot [b a]))


(extend-protocol AbstractViewVectorOps
  (Class/forName "[D")
  (abstract-view-dot [b a] (dot-abstract-views (->DenseVector ^doubles b 0 (count b)) a))

  AbstractView
  (abstract-view-dot [b a] (dot-abstract-views a b))

  Object
  (abstract-view-dot [b a] (mp/vector-dot b a)))


(extend-protocol mp/PVectorOps
  AbstractView
  (vector-dot [a b] (abstract-view-dot b a))
  (length [a] (Math/sqrt (abstract-view-dot a a)))
  (length-squared [a] (abstract-view-dot a a))
  (normalize [a] (ma/div a (Math/sqrt (abstract-view-dot a a)))))

(defn get-arg-mat-type
  [item]
  (cond
    (instance? AbstractMatrix item) :matrix
    (instance? AbstractVector item) :vector
    (instance? Number item) :double
    :else
    (throw (Exception. (str "Unsupported type:" (type item))))))


(defn blas-axpy!
  "y gets ax + y.  x, y are abstract views and alpha is a double"
  [^double alpha x ^AbstractView y]
  (let [^AbstractView x (to-netlib x)
        alpha (double alpha)
        ^AbstractView y y
        ^StridedBuffer x-view (.getStridedBuffer x)
        ^StridedBuffer y-view (.getStridedBuffer y)
        x-len (StridedBuffer/getDataLength x-view)
        y-len (StridedBuffer/getDataLength y-view)]
    (when-not (= 0 (rem y-len x-len))
      (throw (Exception. "y-len must be even multiple of x-len")))
    (if (= 1 x-len)
      (let [x-val (aget ^doubles (.data x-view)
                        (StridedBuffer/getTotalOffset x-view 0))]
        (mp/matrix-add! y (* alpha x-val)))
      (let [axpy-op (reify IStridedBinaryOp
                      (op [this y-data y-offset x-data x-offset op-len]
                        ;;Found through some perf testing.  Probably different ratios
                        ;;on other architectures.
                        (if (< op-len 150)
                          (Ops/axpy op-len alpha
                                    x-data x-offset
                                    y-data y-offset)
                          (.daxpy (BLAS/getInstance) op-len alpha
                                  x-data x-offset 1
                                  y-data y-offset 1))))]
        (StridedBuffer/binaryOperation y-view x-view axpy-op))))
  y)



(extend-protocol mp/PAddScaledMutable
  AbstractView
  (add-scaled! [m a factor]
    (if (ma/scalar? a)
      (mp/matrix-add! m (* ^double a ^double factor))
      (blas-axpy! factor a m))))


(extend-protocol mp/PAddScaled
  AbstractView
  (add-scaled [m a factor]
    (mp/add-scaled! (clone-abstract-view m) a factor)))

(defn make-ones-data
  [new-len]
  (let [^DenseVector new-data (new-dense-vector new-len)
        ^StridedBuffer data-buf (.getStridedBuffer new-data)]
    (.set data-buf 1.0)
    new-data))

(defonce ones-data (atom (make-ones-data 1000)))

(defn ones
  "Array of ones of length m.  Don't assign to this or you get what you get."
  [m]
  (let [ones-len (if (ma/scalar? m)
                   (long m)
                   (StridedBuffer/getDataLength (.getStridedBuffer ^AbstractView m)))]
    (loop [^DenseVector current-ones @ones-data]
      (let [current-len (.length current-ones)]
        (if (> ones-len current-len)
          (let [new-len (* 2 ones-len)
                new-data (make-ones-data new-len)]
            (compare-and-set! ones-data current-ones new-data)
            (recur @ones-data))
          (ma/subvector current-ones 0 ones-len))))))


(extend-protocol mp/PMatrixAddMutable
  AbstractView
  (matrix-add! [m a]
    (if (ma/scalar? a)
      (let [^AbstractView m m
            ^double a a]
        (blas-axpy! a (ones m) m))
      (blas-axpy! 1.0 a m)))
  (matrix-sub! [m a]
    (if (ma/scalar? a)
      (let [^AbstractView m m
            a (* -1.0 ^double a)]
        (blas-axpy! a (ones m) m))
      (blas-axpy! -1.0 a m))))

(extend-protocol mp/PMatrixAdd
  AbstractView
  (matrix-add [m a]
    (mp/matrix-add! (clone-abstract-view m) a))
  (matrix-sub [m a]
    (mp/matrix-sub! (clone-abstract-view m) a)))


(extend-protocol mp/PMatrixMutableScaling
  AbstractView
  (scale! [m constant]
    (let [^StridedBuffer m-view (.getStridedBuffer ^AbstractView m)
          constant (double constant)]
      (if (StridedBuffer/isDense m-view)
        (let [op-len (StridedBuffer/getDataLength m-view)]
          (when (> op-len 0)
            (.dscal (BLAS/getInstance) op-len constant
                    (.data m-view) (StridedBuffer/getTotalOffset m-view 0) 1)))
        (StridedBuffer/unaryOperation m-view (fn [^doubles data ^long offset ^long len]
                                               (.dscal (BLAS/getInstance) len constant
                                                       data offset 1)))))
    m)
  (pre-scale [m constant] (mp/scale! m constant)))


(extend-protocol mp/PMatrixScaling
  AbstractView
  (scale [m constant] (mp/scale! (clone-abstract-view m) constant))
  (pre-scale [m constant] (mp/scale! (clone-abstract-view m) constant)))


(defn make-dense-vector
  ^DenseVector [^AbstractView m]
  (if (instance? DenseVector m)
    m
    (let [^StridedBuffer data-view (.getStridedBuffer ^AbstractView m)]
      (if (StridedBuffer/isDense data-view)
        (->DenseVector (.data data-view) (StridedBuffer/getTotalOffset data-view 0)
                       (StridedBuffer/getDataLength data-view))
        (let [new-view (.clone data-view)]
          (->DenseVector (.data new-view) 0 (StridedBuffer/getDataLength new-view)))))))


(defn make-dense-matrix
  (^DenseMatrix [^AbstractMatrix m]
   (if (instance? DenseMatrix m)
     m
     (let [^StridedBuffer data-view (.getStridedBuffer ^AbstractView m)]
       (if (StridedBuffer/isDense data-view)
         (->DenseMatrix (->DenseVector (.data data-view) (StridedBuffer/getTotalOffset data-view 0)
                                       (StridedBuffer/getDataLength data-view))
                        (.rowCount m)
                        (.columnCount m))
         ;;probably less expensive than operating on non-dense data
         (.clone m)))))

  (^DenseMatrix [m vector-to-matrix-conversion]
   (if (instance? AbstractMatrix m)
     (make-dense-matrix m)
     (let [^DenseVector dense-vec (make-dense-vector m)
           vec-len (.length dense-vec)
           column-count (if (= :row vector-to-matrix-conversion)
                          vec-len
                          1)
           row-count (if (= :row vector-to-matrix-conversion)
                       1
                       vec-len)]
       (->DenseMatrix dense-vec row-count column-count)))))


(defn matrix-or-vector?
  [item]
  (or (instance? AbstractMatrix item)
      (instance? AbstractVector item)))


(defn are-abstract-views-overlapping?
  [^AbstractView lhs ^AbstractView rhs]
  (StridedBuffer/areOverlapping (.getStridedBuffer lhs) (.getStridedBuffer rhs)))


(defn blas-gemm!
  "Defaults to a being a row-matrix if a vector and b being a column-matrix if a vector"
  ([trans-a? trans-b? alpha a b beta c]
   (when-not (and (matrix-or-vector? a)
                  (matrix-or-vector? b))
     (throw (Exception. "Unsupported")))
   (let [factor (double alpha)
         beta (double beta)
         ^DenseMatrix a (make-dense-matrix (to-netlib a) :row)
         ^DenseMatrix b (make-dense-matrix (to-netlib b) :column)
         ^AbstractView c c
         a-dims (if trans-a? [(.columnCount a) (.rowCount a)] [(.rowCount a) (.columnCount a)])
         b-dims (if trans-b? [(.columnCount b) (.rowCount b)] [(.rowCount b) (.columnCount b)])
         M (first a-dims)
         N (second b-dims)
         K (first b-dims)
         overlapping? (or (are-abstract-views-overlapping? a c)
                          (are-abstract-views-overlapping? b c))
         dest (if overlapping?
                (clone-abstract-view c)
                c)
         ;;C is a MxN matrix.
         ^DenseMatrix m (cond
                          (= 1 M) (make-dense-matrix dest :row)
                          (= 1 N) (make-dense-matrix dest :column)
                          :else
                          (make-dense-matrix dest))
         ^DenseVector a-data (.data a)
         ^DenseVector b-data (.data b)
         ^DenseVector m-data (.data m)
         trans-command-a (if trans-a? "t" "n")
         trans-command-b (if trans-b? "t" "n")]
     (when-not (and (= K (second a-dims))
                    (= M (.rowCount m))
                    (= N (.columnCount m)))
       (throw (Exception. (format "Incompatible matrix sizes: a %s b %s m %s"
                                  (str a-dims)
                                  (str b-dims)
                                  (str (ma/shape m))))))
     (.dgemm (BLAS/getInstance) trans-command-b trans-command-a N M K alpha
             (.data b-data) (.offset b-data) (.columnCount b)
             (.data a-data) (.offset a-data) (.columnCount a)
             beta
             (.data m-data) (.offset m-data) (.columnCount m))
     (assign-source-to-view! m c)))
  ([alpha a b beta c]
   (blas-gemm! false false alpha a b beta c)))


(defmulti typed-blas-gemm! (fn [trans-a? trans-b? alpha a b beta c]
                             [(get-arg-mat-type a)
                              (get-arg-mat-type b)]))

(defmethod typed-blas-gemm! [:matrix :matrix]
  [trans-a? trans-b? alpha a b beta c]
  (blas-gemm! trans-a? trans-b? alpha a b beta c))

(defmethod typed-blas-gemm! [:matrix :vector]
  [trans-a? trans-b? alpha a b beta c]
  (blas-gemm! trans-a? trans-b? alpha a b beta c))

(defmethod typed-blas-gemm! [:vector :matrix]
  [trans-a? trans-b? alpha a b beta c]
  (blas-gemm! trans-a? trans-b? alpha a b beta c))


(defmethod typed-blas-gemm! [:vector :vector]
  [trans-a? trans-b? alpha a b beta c]
  (let [outer-product? (or trans-a? trans-b?)]
    (if (outer-product?)
      (blas-gemm! trans-a? trans-b? alpha a b beta c)
      (+ (* ^double alpha ^double (dot-abstract-views a b))
         (* ^double beta ^double c)))))


(defn typed-blas-gemm-view-scalar!
  "A is a view, alpha b beta are scalars.  C is a view"
  [alpha a b beta c]
  (let [^double beta beta]
    (if (= 0.0 beta)
      (do
        (mp/assign! c a)
        (mp/scale! c b))
      (do
        (mp/scale! c beta)
        (blas-axpy! b a c)))
    c))


(defmethod typed-blas-gemm! [:vector :scalar]
  [trans-a? trans-b? alpha a b beta c]
  (typed-blas-gemm-view-scalar! alpha a b beta c))


(defmethod typed-blas-gemm! [:matrix :scalar]
  [trans-a? trans-b? alpha a b beta c]
  (typed-blas-gemm-view-scalar! alpha a b beta c))


;;General protocol for things of the shape:
;;c = alpha*A*b + beta*c
;;if c is a vector then we know that b is a vector
;;if c is a matrix, then we know that b is a matrix
;;A must always be a matrix
(extend-protocol mp/PAddInnerProductMutable
  AbstractView
  (add-inner-product! [m a b] (mp/add-inner-product! m a b 1.0))
  (add-inner-product! [m a b factor]
    (typed-blas-gemm! false false factor (to-netlib a) (to-netlib b) 1.0 m)))


(defmulti create-result-for-inner-product (fn [trans-a? trans-b? a b]
                                            [(get-arg-mat-type a)
                                             (get-arg-mat-type b)]))

(defn get-item-vector-shape-result
  (^long [trans? ^long vector-length rows-or-columns]
   (if-not trans?
     1
     vector-length)))

(defn get-item-matrix-shape-result
  (^long [trans? shape rows-or-columns]
   (let [^long row-count (shape 0)
         ^long column-count (shape 1)]
     (if trans?
       (if (= :row-count rows-or-columns)
         column-count
         row-count)
       (if (= :row-count rows-or-columns)
         row-count
         column-count)))))


(defn get-item-shape
  "Defaults to row vectors"
  ^long [trans? item rows-or-columns]
  (let [mat-type (get-arg-mat-type item)]
    (cond
      (= mat-type :matrix) (get-item-matrix-shape-result trans? (ma/shape item)
                                                  rows-or-columns)
      (= mat-type :vector) (get-item-vector-shape-result trans? (.length ^AbstractVector item)
                                                         rows-or-columns)
      :else 1)))


(defn shape-to-inner-product-result
  [^long row-count ^long column-count]
  (if (and (= 1 row-count)
           (= 1 column-count))
    0.0
    (let [^DenseVector backing (new-dense-vector (* row-count column-count))]
      (if (or (= 1 row-count)
              (= 1 column-count))
        backing
        (->DenseMatrix backing row-count column-count)))))


(defn create-result-for-gemm
  [trans-a? trans-b? a b]
  (let [scalar-a? (ma/scalar? a)
        scalar-b? (ma/scalar? b)]
    (if (or scalar-a?
            scalar-b?)
      (if (and scalar-a?
               scalar-b?)
        0.0
        (if scalar-a?
          (new-abstract-view b)
          (new-abstract-view a)))
      ;;Else neither are scalars and off we go...
      (let [row-count (get-item-shape trans-a? a :row-count)
            column-count (get-item-shape trans-b? b :column-count)]
        (shape-to-inner-product-result row-count column-count)))))


(defn blas-gemm
  [trans-a? trans-b? alpha a b]
  (let [a (to-netlib a)
        b (to-netlib b)]
    (typed-blas-gemm! trans-a? trans-b? alpha a b 0.0
                      (create-result-for-gemm trans-a? trans-b? a b))))


(extend-protocol mp/PMatrixProducts
  AbstractView
  ;;If these are both vectors, then this turns into dotproduct.
  (inner-product [a b]
    (blas-gemm false false 1.0 a b))

  (outer-product [m a]
    (if (mp/is-scalar? m)
      (mp/pre-scale a m)
      (mp/element-map (mp/convert-to-nested-vectors m)
                      (fn [v] (mp/pre-scale a v))))))


  ;;took a bit to find this one...
(defn blas-element-multiply!
  "y = alpha * (elem-mul a x) + b*y"
  [alpha a x beta y]
  (let [^DenseVector a (make-dense-vector a)
        ^DenseVector x (make-dense-vector x)
        ^AbstractView y y
        overlapping? (or (StridedBuffer/areOverlapping (.getStridedBuffer a)
                                                         (.getStridedBuffer y))
                         (StridedBuffer/areOverlapping (.getStridedBuffer x)
                                                         (.getStridedBuffer y)))
        dest (if overlapping?
               (clone-abstract-view y)
               y)
        ^DenseVector yd (make-dense-vector dest)
        alpha (double alpha)
        beta (double beta)]
    ;(println a x y dest yd)
    (when-not (and (= (.length a) (.length x))
                   (= (.length a) (.length yd)))
      (throw (Exception. "Element multiple of different length vectors unsupported")))
    (.dsbmv (BLAS/getInstance) "U" (.length a) 0
            alpha (.data a) (.offset a) 1
            (.data x) (.offset x) 1
            beta (.data yd) (.offset yd) 1)
    (assign-source-to-view! yd y)
    y))

(defn dense-element-multiply!
  [alpha a-data a-offset x-data x-offset beta y-data y-offset op-len]
  (let [data-len (long op-len)
        a-offset (long a-offset)
        x-offset (long x-offset)
        y-offset (long y-offset)
        alpha (double alpha)
        beta (double beta)
        ^doubles a-data a-data
        ^doubles x-data x-data
        ^doubles y-data y-data]
    (dotimes [idx data-len]
      (aset y-data (+ idx y-offset)
            (+ (* beta (aget y-data (+ idx y-offset)))
               (* alpha (aget x-data (+ idx x-offset))
                  (aget a-data (+ idx a-offset))))))))

(defn non-blas-element-multiply!
  [alpha a x beta y]
  (let [^StridedBuffer a-view (.getStridedBuffer ^AbstractView (to-netlib a))
        ^StridedBuffer x-view (.getStridedBuffer ^AbstractView (to-netlib x))
        ^StridedBuffer y-view (.getStridedBuffer ^AbstractView y)
        alpha (double alpha)
        beta (double beta)
        a-len (StridedBuffer/getDataLength a-view)
        x-len (StridedBuffer/getDataLength x-view)]
    (if (or (= 1 a-len)
            (= 1 x-len))
      (let [reduce-val (double
                        (if (= 1 a-len)
                          (aget (.data a-view) (StridedBuffer/getTotalOffset a-view 0))
                          (aget (.data x-view) (StridedBuffer/getTotalOffset x-view 0))))
            alpha (* alpha reduce-val)
            new-x (if (= 1 a-len)
                    x
                    a)]
        (blas-axpy! alpha new-x beta y))
      ;;Else we may need to actually consider all three objects
      (if (= 0.0 beta)
        (do
          ;;Assign to the overlapping views.  This takes care of identity
          ;;And non pathological cases.
          (if (StridedBuffer/areOverlapping y-view a-view)
            (.assign y-view a-view)
            (.assign y-view x-view))

          (let [src-view (if (StridedBuffer/areOverlapping y-view a-view)
                           x-view
                           a-view)]
            (StridedBuffer/binaryOperation
             y-view src-view
             (reify IStridedBinaryOp
               (op [this y-data y-offset x-data x-offset op-len]
                 (Ops/alphaAY op-len alpha x-data x-offset y-data y-offset))))))
        (StridedBuffer/ternaryOperation
         y-view a-view x-view
         (reify IStridedTernaryOp
           (op [this lhs-data lhs-offset mid-data mid-offset rhs-data rhs-offset op-len]
             (Ops/alphaAXPlusBetaY op-len alpha
                                   rhs-data rhs-offset
                                   mid-data mid-offset
                                   beta
                                   lhs-data lhs-offset))))))
    y))



(defmulti typed-element-multiply! (fn [m a result alpha beta]
                                    [(get-arg-mat-type m)
                                     (get-arg-mat-type a)]))

(defmethod typed-element-multiply! [:matrix :matrix]
  [m a result alpha beta]
  (non-blas-element-multiply! alpha m a beta result))

(defmethod typed-element-multiply! [:matrix :vector]
  [m a result alpha beta]
  (doall (map (fn [m-row res-row]
                (non-blas-element-multiply! alpha m-row a beta res-row))
              m
              result))
  result)


(defmethod typed-element-multiply! [:vector :vector]
  [a b result alpha beta]
  (non-blas-element-multiply! alpha b a beta result))


(defmethod typed-element-multiply! [:vector :matrix]
  [a b result alpha beta]
  (typed-element-multiply! b a result alpha beta))


(defn clone-larger-view
  [a b]
  (let [view-a (.getStridedBuffer ^AbstractView a)
        view-b (.getStridedBuffer ^AbstractView b)]
    (if (> (StridedBuffer/getDataLength view-a)
           (StridedBuffer/getDataLength view-b))
      (clone-abstract-view a)
      (clone-abstract-view b))))


(extend-protocol mp/PMatrixMultiplyMutable
  AbstractView
  (matrix-multiply! [m a]
    (if (ma/scalar? a)
      (ma/scale! m a)
      (assign-source-to-view! (mp/inner-product m a) m)))

  (element-multiply! [m a]
    (if (ma/scalar? a)
      (ma/scale! m a)
      (assign-source-to-view! (typed-element-multiply! (to-netlib a) m m
                                                       1.0 0.0)
                              m))))

(extend-protocol mp/PMatrixMultiply
  AbstractView
  (matrix-multiply [m a]
    (if (ma/scalar? a)
      (ma/scale m a)
      (mp/inner-product m a)))

  (element-multiply [m a]
    (if (ma/scalar? a)
      (ma/scale! m a)
      (typed-element-multiply! m (to-netlib a) (clone-larger-view m a) 1.0 0.0))))


(extend-protocol mp/PMatrixDivideMutable
  AbstractView
  (element-divide!
    ([m] (mp/element-map! m /))
    ([m a] (if (ma/scalar? a)
             (ma/scale! m (/ 1.0 (double a)))
             (let [m-view (.getStridedBuffer ^AbstractView m)
                   a-view (.getStridedBuffer ^AbstractView a)]
               (StridedBuffer/binaryOperationOp
                m-view a-view
                (reify IBinaryOp
                  (op [this lhs-val rhs-val]
                    (/ lhs-val rhs-val)))))))))


(extend-protocol mp/PMatrixDivide
  AbstractView
  (element-divide
    ([m] (mp/element-divide! (clone-abstract-view m)))
    ([m a] (mp/element-divide! (clone-abstract-view m) a))))


(extend-protocol mp/PAddScaledProductMutable
  AbstractView
  (add-scaled-product! [m a b ^double factor]
    (if (or (ma/scalar? a) (ma/scalar? b))
      (let [scalar-a? (ma/scalar? a)
            a (if scalar-a? b a)
            b (if scalar-a? a b)]
        (mp/add-scaled! m a (* factor ^double b)))
      (typed-element-multiply! (to-netlib a) (to-netlib b) m factor 1.0))))


(extend-protocol mp/PAddScaledProduct
  AbstractView
  (add-scaled-product [m a b factor]
    (mp/add-scaled-product! (clone-abstract-view m) (to-netlib a) (to-netlib b) factor)))


;;Turns out this is a common operation...
(extend-protocol mp/PTranspose
  AbstractView
  (transpose [m]
    (if (instance? AbstractVector m)
      m
      (let [^AbstractMatrix m m
            ^DenseMatrix retval (new-dense-matrix (.columnCount m) (.rowCount m))]
        (doall (map (fn [source-row dest-column]
                      (assign-source-to-view! source-row dest-column))
                    (ma/rows m)
                    (ma/columns retval)))
        retval))))


(extend-protocol mp/PElementCount
  AbstractView
  (element-count [m] (StridedBuffer/getDataLength (.getStridedBuffer ^AbstractView m))))



(eval
 `(extend-protocol mp/PMathsFunctionsMutable
    AbstractView
    ~@(map
        (fn [[fname op]]
          (let [fname (symbol (str fname "!"))]
            `(~fname [~'m] (unary-op-macro! ~'m (~op ~'lhs-value)))))
        mathsops/maths-ops)))


(defmacro maths-ops
  []
  `(extend-protocol mp/PMathsFunctions
     AbstractView
     ~@(map
         (fn [[fname op]]
           `(~fname [m#] (let [^AbstractView m# (clone-abstract-view m#)]
                            (unary-op-macro! m# (~op ~'lhs-value)))))
         mathsops/maths-ops)))

(maths-ops)

(def empty-vec (new-dense-vector 0))
7(mi/register-implementation empty-vec)
