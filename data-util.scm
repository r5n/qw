;;;; data-util.scm
;;;; Utility functions for working with datasets

(declare (usual-integrations))

(load "mat-mit")

(define (random-perm data)
  (error "Not implmented"))

(define (train-test-split data ratio)
  ;; labels should be the last column of the data
  (let* ((train-length (* ratio (length data))))
    (cons (take data (floor->exact train-length))
	  (drop data (ceiling->exact train-length)))))

(define (train-validate-test-split data train-perc v-perc)
  (error "Not implemented"))
