;;; A lot of these functions are written in an
;;; imperative style.

;;;; TODO
;; improve the way cofactor is written
;; improve the way mat-add is written

;;;; FEATURES TO ADD
;; inverse
;; gaussian elimination
;; solving linear equations
;; singular value decomposition

;;;; Macros
(define-syntax for
  (syntax-rules (from to downto by)
    ((_ (i from start to end) b1 ...)
     (do ((i start (+ i 1)))
	 ((= i end))
       b1 ...))
    ((_ (i from start to end by step) b1 ...)
     (do ((i start (+ i step)))
	 ((= i end))
       b1 ...))
    ((_ (i from start downto end) b1 ...)
     (do ((i start (- i 1)))
	 ((= i end))
       b1 ...))
    ((_ (i from start downto end by step) b1 ...)
     (do ((i start (- i step)))
	 ((= i end))
       b1 ...))))

(define-syntax mat-traverse
  (syntax-rules ()
    ((_ mat (i j) b1 ...)
     (for (i from 0 to (rows mat))
       (for (j from 0 to (cols mat))
	 b1 ...)))))

;;;; Matrices

(define (make-mat rows cols)
  (vector-map (lambda (x) (make-vector cols 0))
	      (make-vector rows 0)))

(define (dim mat)
  (cons (vector-length mat)
	(vector-length (vector-ref mat 0))))

(define (mat-ref mat i j)
  (vector-ref (vector-ref mat i) j))

(define (mat-set! mat i j val)
  (vector-set! (vector-ref mat i) j val))

(define (row-ref mat n)
  (vector-ref mat n))

(define (row-set! mat r val)
  (vector-set! mat r val))

(define (col-ref mat n)
  (vector-map (lambda (row) (vector-ref row n)) mat))

(define (diag-ref mat)
  (let ((i -1))
    (vector-map (lambda (x)
		  (begin (set! i (+ i 1))
			 (vector-ref x i))) mat)))
(define (rows mat)
  (car (dim mat)))

(define (cols mat)
  (cdr (dim mat)))

(define (is-square? mat)
  (eq? (rows mat) (cols mat)))

(define (mat->list mat)
  (vector->list
   (vector-map vector->list mat)))

(define (list->mat ls)
  (apply vector
	 (map (lambda (r) (apply vector r)) ls)))

(define (transpose mat)
  (apply vector (apply map vector (mat->list mat))))

;; simple algorithm. O(n^3)
(define (mat-mul A B)
  (if (= (cols A) (rows B))
      (let ((C (make-mat (rows A) (cols B))))
	(for (i from 0 to (rows A))
	  (for (j from 0 to (cols B))
	    (let ((sum 0))
	      (for (k from 0 to (cols A))
		(set! sum (+ sum (* (mat-ref A i k)
				    (mat-ref B k j)))))
	      (mat-set! C i j sum))))
	C)
      (error "Invalid dimensions")))

;; can make better using vector-map +
(define (mat-add A B)
  (if (equal? (dim A) (dim B))
      (let ((C (make-mat (rows A) (cols B))))
	(for (i from 0 to (rows A))
	  (for (j from 0 to (cols B))
	    (mat-set! C i j
		      (+ (mat-ref A i j)
			 (mat-ref B i j)))))
	C)
      (error "Invalid dimensions")))

(define (trace mat)
  (apply + (vector->list (diag-ref mat))))

;; special matrices
(define (iden n)
  (let ((id (make-mat n n)))
    (for (i from 0 to n)
      (mat-set! id i i 1))
    id))

(define (random-mat n r c)
  (let ((M (make-mat r c)))
    (for (i from 0 to r)
      (for (j from 0 to c)
	(mat-set! M i j (random n))))
    M))

;; row operations
(define (row-scl! mat r p)
  (for (i from 0 to (cols mat))
    (mat-set! mat r i (* p (mat-ref mat r i)))))

(define (row-swp! mat r1 r2)
  (let ((tmp (row-ref mat r1)))
    (row-set! mat r1 (row-ref mat r2))
    (row-set! mat r2 tmp)))

(define (row-add! mat r1 r2)
  (row-set! mat r1
	    (vector-map + (row-ref mat r1) (row-ref mat r2))))

;; determinants
(define (det-2x2 mat)
  (- (* (mat-ref mat 0 0) (mat-ref mat 1 1))
     (* (mat-ref mat 0 1) (mat-ref mat 1 0))))

(define (cofactor mat p q)
  (let* ((r (rows mat))
	 (c (cols mat))
	 (cof (make-mat (- r 1) (- c 1)))
	 (a 0) (b 0)) ;; a, b used to fill cof
    (for (i from 0 to r) ;; i, j to traverse mat
      (for (j from 0 to c)
	(if (and (not (= i p)) (not (= j q)))
	    (begin
	      (mat-set! cof a b (mat-ref mat i j))
	      (set! b (+ b 1))
	      (if (= b (- r 1))
		  (begin (set! b 0)
			 (set! a (+ a 1))))))))
    cof))

;;; Doolittle's LU Decomposition
;;
;; for each i = 0,1,2,...,n-1:
;; U(i,k) = A(i,k) - SUM(j=0,i)(L(i,j)*U(j,k))
;; for k = i,i+1,...,n-1 produces the kth row
;; of U
;;
;; L(i,k) = [A(i,k) - SUM(j=0,i)(L(i,j)*U(j,k))]/U(k,k)
;; for i=k+1,k+2,...,n-1 and L(i,i) = 1 produces
;; the kth column of L
;;
;; retuns a pair L U such that mat = L * U
(define (lu-decomp A)
  (let* ((n (rows A))
	 (L (make-mat n n))
	 (U (make-mat n n)))
    (for (i from 0 to n)
      (for (k from i to n)
	(let ((sum 0))
	  (for (j from 0 to i)
	    (set! sum (+ sum (* (mat-ref L i j)
				(mat-ref U j k)))))
	  (mat-set! U i k (- (mat-ref A i k) sum))))
      (for (k from i to n)
	(if (= i k)
	    (mat-set! L i i 1)
	    (let ((sum 0))
	      (for (j from 0 to i)
		(set! sum (+ sum (* (mat-ref L k j)
				    (mat-ref U j i)))))
	      (mat-set! L k i (/ (- (mat-ref A k i) sum)
			     (mat-ref U i i)))))))
    (cons L U)))

;;; Determinant (uses LU decomposition)
;; det(A) = det(LU) = det(L)*det(U)
;; Determinant of an upper or lower triangular
;; matrix is the product of the diagonal elements
;; Since L's diagonals are 1's, det(A) = det(U)
(define (det mat)
  (let* ((U (cdr (lu-decomp mat)))
	 (diags (diag-ref U))
	 (prod 1))
    (for (i from 0 to (vector-length diags))
      (set! prod (* prod (vector-ref diags i))))
    prod))

;;;; Test Matrices
(define m1 (list->mat '((2 0) (-3 -1))))
(define m2 (list->mat '((2 3 -7) (1 0 1))))
(define m3 (list->mat '((1 2 3) (4 5 6) (7 8 9))))
(define m4 (list->mat '((1 0) (0 1))))
(define m5 (list->mat '((24 50 12)
			(51 23 90)
			(-17 3 23))))
(define m6 (list->mat '((71 8 33)
			(8 76 24)
			(7 33 32))))
(define m7 (list->mat '((80 40 42 83)
			(2 70 63 28)
			(25 49 40 45)
			(86 30 93 34))))

;;;; Test functions
(define (test-lu-decomp mat)
  (let* ((decomp (lu-decomp mat))
	 (L (car decomp))
	 (U (cdr decomp)))
    (equal? mat (mat-mul L U))))
