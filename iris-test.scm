(load "iris-data")
(load "data-util")
(load "mat-mit") ;; need to change name of file

(define iris-dataset
  (let ((d (list->mat iris-data)))
    (for (i from 0 to (rows d))
      (let ((j (random (rows d))))
	(row-swp! d i j)))
    d))

(define split
  (train-test-split (mat->list iris-dataset) 0.5))

(define iris-trainset (list->mat (car split)))
(define iris-testset  (list->mat (cdr split)))

(define iris-train-labels
  (col-ref iris-trainset 4))
(define iris-test-labels
  (col-ref iris-testset 4))

(define iris-train
  (make-initialized-vector
   (rows iris-trainset)
   (lambda (i) (subvector (row-ref iris-trainset i) 0 4))))

(define iris-test
  (make-initialized-vector
   (rows iris-testset)
   (lambda (i) (subvector (row-ref iris-testset i) 0 4))))


;;;; Test KNN classifier
(load "knn")

(define (iris-knn-accuracy k)
  (let ((correct 0)
	(wrong 0))
    (begin
      (display "IRIS DATASET : KNN CLASSIFIER") (newline)
      (display "-----------------------------") (newline)
      (display "Training set size : ")
      (display (rows iris-train)) (newline)
      (display "Testing set size  : ")
      (display (rows iris-test)) (newline)
      (display "Training ...")
      (for (index from 0 to (rows iris-test))
	(let* ((sample (row-ref iris-test index))
	       (prediction (knn-classify sample iris-train iris-train-labels k)))
	  (if (equal? prediction (vector-ref iris-test-labels index))
	      (set! correct (+ correct 1))
	      (set! wrong (+ wrong 1)))))
      (newline)
      (display "Correct  : ") (display correct) (newline)
      (display "Wrong    : ") (display wrong) (newline)
      (display "Accuracy : ") (display (exact->inexact (/ correct (+ correct wrong))))
      (newline)
      (display "-----------------------------") (newline))
    #t))
