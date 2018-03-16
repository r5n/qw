;;;; K Nearest Neighbors

(load "mat-mit")

(define sample-data-set
  (list->mat '((1 1.1)
	       (1 1)
	       (0 0)
	       (0 0.1))))

(define labels
  (list->vector '(A A B B)))

(define (calculate-distances input data-set)
  (let* ((len (vector-length data-set))
	 (sq-diff-mat (mat-map (lambda (i) (square i))
			       (mat+mat (repvec-rows input len)
					(negate data-set))))
	 (sq-distances (make-initialized-vector
			len (lambda (i) (vector-fold + 0 (row-ref sq-diff-mat i)))))
	 (distances (vector-map (lambda (i) (sqrt i)) sq-distances)))
    distances))

(define (knn-classify input data-set labels k)
  ;; input : input vector to classify
  ;; data-set : matrix of training samples (n x m)
  ;; labels : a vector of labels (n x 1)
  ;; k : # of nearest neighbors
  (let ((d (calculate-distances input data-set)))
    d))

