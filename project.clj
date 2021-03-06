(defproject thinktopic/netlib-ccm "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.6.0"]
                 [net.mikera/core.matrix "0.49.0"]
                 [com.github.fommil.netlib/all "1.1.2" :extension "pom"]]

  :java-source-paths ["java"]

  :profiles { :dev {:dependencies [[org.clojure/tools.namespace "0.2.11"]
                                   [net.mikera/vectorz-clj "0.43.0"]]
                    :plugins [[lein-nodisassemble "0.1.3"]]}}
  )
