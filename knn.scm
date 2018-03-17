;;;; K Nearest Neighbors

(load "mat-mit")

(define (vector-zip v1 v2)
  (make-initialized-vector
   (min (vector-length v1) (vector-length v2))
   (lambda (i)
     (cons (vector-ref v1 i) (vector-ref v2 i)))))

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
  (let* ((d (calculate-distances input data-set))
	 (dist-labl (vector-zip d labels))
	 (sorted-dist-labl (sort! dist-labl
				  (lambda (a b)
				    (< (car a) (car b)))))
	 (k-neighbors (vector-head sorted-dist-labl k))
	 (votes (make-initialized-vector (vector-length k-neighbors)
				       (lambda (i) (cdr (vector-ref k-neighbors i))))))
    (most-freq votes)))

;; quite ineffective. works for now.
(define (most-freq v)
  (let ((h (make-hash-table)))
    (for (i from 0 to (vector-length v))
      (hash-table/modify! h
			  (vector-ref v i)
			  0
			  (lambda (e) (+ e 1))))
    (caar (sort (hash-table->alist h)
		(lambda (a b) (> (cdr a) (cdr b)))))))
