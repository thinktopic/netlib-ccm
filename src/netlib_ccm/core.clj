(ns netlib-ccm.core
  (:require [clojure.core.matrix :as ma]
            [clojure.core.matrix.protocols :as mp]
            [clojure.core.matrix.implementations :as mi]
            [clojure.reflect :as r])
  (:import [com.github.fommil.netlib BLAS]))


(defrecord array-view [^doubles data ^long offset ^long length])
(defrecord matrix-view [^array-view data ^long row-count ^long column-count])
