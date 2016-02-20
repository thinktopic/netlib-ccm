(ns netlib-ccm.strided-view
  (:import [netlib_ccm Ops]))


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

   (when (< ^long row-stride ^long column-count)
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
  ;;special case 0 so that everyone that wants to get an offset for pos 0
  ;;doesn't have to check if it is within length of the vector
  (if (= 0 data-offset)
    (+ (.offset view) (max 0 (- (.column-count view) (.initial-column-count view))))
    (do
      (check-strided-view-data-offset view data-offset)
      (let [data-len (get-strided-view-data-length view)
            body-len (* (.column-count view) (max 0 (- (.row-count view) 2)))]
        (if (< data-offset (.initial-column-count view))
          (+ (+ (.offset view) (- (.column-count view) (.initial-column-count view)))
             data-offset)
          (let [data-rest (- data-offset (.initial-column-count view))]
            (+ (.offset view)
               (* (.row-stride view) (+ 1 (quot data-rest (.column-count view))))
               (rem data-rest (.column-count view)))))))))


(defn create-sub-strided-view
  [^StridedView view ^long offset ^long length]
  ;;special case for (= 0 length)
  (if (= 0 length)
    (new-strided-view (.data view) (get-strided-view-total-offset view 0) 0)
    (if (= length (get-strided-view-data-length view))
      view
      (do
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
              (new-strided-view (.data view) (get-strided-view-total-offset view offset) length)
              (let [row-idx (get-strided-view-row-from-data-offset view offset)
                    last-row-len (rem (+ first-row-offset length) (.column-count view))
                    ary-offset (+ (.offset view) (* (.row-stride view) row-idx))
                    body-num-rows (quot (- length (+ first-row-len last-row-len))
                                        (.column-count view))
                    ini-row (long (if-not (= 0 first-row-len) 1 0))
                    end-row 1]
                (new-strided-view (.data view) ary-offset
                                  (+ body-num-rows ini-row end-row)
                                  (.column-count view)
                                  (.row-stride view)
                                  first-row-len
                                  last-row-len)))))))))

(defn is-strided-view-dense?
  "A dense strided view is one has no gaps in its data"
  [^StridedView view]
  (or (= 1 (.row-count view))
      (= (.column-count view) (.row-stride view))))

(defn is-strided-view-complete-dense?
  "Complete dense means the strided view completely fills its backing store
and has no offset meaning it is accurately and completely represented by
  its backing store"
  [^StridedView view]
  (and (is-strided-view-dense? view)
       (= 0 (.offset view))
       (= (count (.data view))
          (get-strided-view-data-length view))))



(defn strided-op
  "Perform an operation such as (assign! lhs rhs).
  Op gets passed: lhs-double-array lhs-offset rhs-double-array rhs-offset op-amount"
  ([^StridedView lhs ^StridedView rhs op]
   (if (and (is-strided-view-dense? lhs)
            (is-strided-view-dense? rhs))
     (let [op-len (min (get-strided-view-data-length lhs)
                       (get-strided-view-data-length rhs))]
       (op (.data lhs) (get-strided-view-total-offset lhs 0)
           (.data rhs) (get-strided-view-total-offset rhs 0)
           op-len))
     ;;Less efficient case where one or both views are not dense and we have to
     ;;loop over rows.
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
                     op-amount (min lhs-row-len rhs-row-len)]
                 (op (.data lhs) lhs-total-offset
                     (.data rhs) rhs-total-offset
                     op-amount)
                 (recur (+ rhs-row-offset op-amount))))))
         (recur (inc rhs-row)))))
   lhs)
  ([^StridedView lhs op]
   (if (is-strided-view-dense? lhs)
     (let [op-len (get-strided-view-data-length lhs)]
       (op (.data lhs) (get-strided-view-total-offset lhs 0) op-len))
     (loop [lhs-row 0]
       (when (< lhs-row (.row-count lhs))
         (let [row-len (get-strided-view-row-length-from-row lhs lhs-row)]
           (when (> row-len 0)
             (let [data-offset (get-strided-view-data-offset-from-row lhs lhs-row)
                   lhs-total-offset (get-strided-view-total-offset lhs data-offset)]
               (op (.data lhs) lhs-total-offset row-len))))
         (recur (inc lhs-row)))))
   lhs))

(defmacro strided-view-binary-java-op!
  [lhs rhs & body]
  `(let [lhs-len# (get-strided-view-data-length ~lhs)
         rhs-len# (get-strided-view-data-length ~rhs)]
     ;;no op
     (when-not (= 0 rhs-len#)
       (when-not (= 0 (rem lhs-len# rhs-len#))
         (throw (Exception. "Strided assignment: rhs-len not even multiple of lhs-len")))
       (if (= 1 rhs-len#)
         (let [^"doubles" rhs-data# (.data ~rhs)
               ~'rhs-value (aget rhs-data#
                                 (get-strided-view-total-offset ~rhs 0))
               single-op# (reify netlib_ccm.IUnaryOp
                            (op [this# ^"double" ~'lhs-value]
                              ~@body))]
           (strided-op ~lhs (fn [data# offset# len#]
                              (Ops/OpY len# data# offset# single-op#))))
         (let [num-reps# (quot lhs-len# rhs-len#)
               multi-op# (reify netlib_ccm.IBinaryOp
                           (op [this# ^"double" ~'rhs-value ^"double" ~'lhs-value]
                             ~@body))]
           (loop [rep# 0]
             (when (< rep# num-reps#)
               (strided-op (create-sub-strided-view ~lhs (* rep# rhs-len#) rhs-len#)
                           ~rhs (fn [lhs-data# lhs-offset# rhs-data# rhs-offset# op-len#]
                                  (Ops/OpXY op-len# rhs-data# rhs-offset# lhs-data# lhs-offset#
                                            multi-op#)))
               (recur (inc rep#)))))))
     ~lhs))


(defn strided-view-multiple-op!
  [^StridedView lhs ^StridedView rhs multiple-val-op single-val-op]
  (let [lhs-len (get-strided-view-data-length lhs)
        rhs-len (get-strided-view-data-length rhs)]
    ;;no op
    (when-not (= 0 rhs-len)
     (when-not (= 0 (rem lhs-len rhs-len))
       (throw (Exception. "Strided assignment: rhs-len not even multiple of lhs-len")))
     (if (= 1 rhs-len)
       (let [rhs-val (aget ^doubles (.data rhs) (get-strided-view-total-offset rhs 0))]
         (strided-op lhs (fn [^doubles data ^long offset ^long len]
                           (single-val-op data offset len rhs-val))))
       (let [num-reps (quot lhs-len rhs-len)]
         (if (= 1 num-reps)
           (strided-op lhs rhs multiple-val-op)
           (loop [rep 0]
             (when (< rep num-reps)
               (strided-op (create-sub-strided-view lhs (* rep rhs-len) rhs-len)
                           rhs multiple-val-op)
               (recur (inc rep))))))))))


(defn strided-view-reference-same-data?
  [^StridedView lhs ^StridedView rhs]
  (or (identical? lhs rhs)
      (and (identical? (.data lhs) (.data rhs))
           (= (get-strided-view-data-length lhs)
              (get-strided-view-data-length rhs))
           (or (and (is-strided-view-dense? lhs)
                    (is-strided-view-dense? rhs))
               (= 0 (get-strided-view-data-length lhs))
               (and (= (.row-stride lhs) (.row-stride rhs))
                    (= (get-strided-view-total-offset lhs 0)
                       (get-strided-view-total-offset rhs 0)))))))


(defn assign-strided-view!
  "lhs must be at least as large as rhs.  Lhs.length must be a multiple
of Rhs.length"
  [^StridedView lhs ^StridedView rhs]
  (when-not (strided-view-reference-same-data? lhs rhs)
    (when (and (> (get-strided-view-data-length lhs) 0)
               (> (get-strided-view-data-length rhs) 0))
      (let [single-val-op (fn [^doubles data ^long offset ^long len ^double rhs-val]
                            (java.util.Arrays/fill data offset (+ offset len) rhs-val))
            multiple-val-op (fn [lhs-data lhs-offset rhs-data rhs-offset op-len]
                              (System/arraycopy ^doubles rhs-data ^long rhs-offset
                                                ^doubles lhs-data ^long lhs-offset
                                                op-len)
                              ;; Testing shows arraycopy to be faster than .dcopy in every
                              ;; case
                              ;; (.dcopy (BLAS/getInstance) (int op-len)
                              ;;         ^doubles rhs-data (int rhs-offset) 1
                              ;;         ^doubles lhs-data (int lhs-offset) 1)
                              )]
        (strided-view-multiple-op! lhs rhs multiple-val-op single-val-op)))))


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


(defn is-strided-view-data-in-range?
  "Assuming these items share a buffer and have nonzero length,
return true if they overlap"
  [^StridedView lhs ^StridedView rhs]
  (let [lhs-start (get-strided-view-total-offset lhs 0)
        rhs-start (get-strided-view-total-offset rhs 0)
        lhs-end (get-strided-view-total-offset lhs (- (get-strided-view-data-length lhs) 1))
        rhs-end (get-strided-view-total-offset rhs (- (get-strided-view-data-length rhs) 1))]
    (or (and (>= rhs-start lhs-start)
             (<= rhs-start lhs-end))
        (and (>= lhs-start rhs-start)
             (<= lhs-start rhs-end)))))


;;Also, blas doesn't like overlapping read/write memory
(defn are-strided-views-overlapping?
  [^StridedView lhs ^StridedView rhs]
  (and (identical? (.data lhs) (.data rhs))
       (and (> (get-strided-view-data-length lhs) 0)
            (> (get-strided-view-data-length rhs) 0))
       (is-strided-view-data-in-range? lhs rhs)))
