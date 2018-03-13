;;;; MIT SCHEME specific version of `mat.scm`
(declare (usual-integrations))

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

;;;; Helper functions

(define (vector-fold f init vec)
  (if (zero? (vector-length vec))
      init
      (vector-fold f
		   (f (vector-ref vec 0) init)
		   (vector-tail vec 1))))

(define (mat-map f A)
  (vector-map (lambda (row) (vector-map (lambda (elt) (f elt)) row)) A))

(define (print-mat mat)
  (let ((max-strlen
	 (fold-left (lambda (a x) (max a
				       (string-length (number->string x))))
		    0 (concatenate (mat->list mat)))))
    (for (i from 0 to (rows mat))
      (if (or (= i 0) (= i (- (rows mat) 1)))
	  (display "[")
	  (display "|"))
      (for (j from 0 to (cols mat))
	(display (string-append
		  (make-string (- max-strlen
				  (string-length (number->string (mat-ref mat i j))))
			       #\space)
		  (number->string (mat-ref mat i j))
		  (if (= j (- (cols mat) 1))
		      ""
		      "  "))))
      (if (or (= i 0) (= i (- (rows mat) 1)))
	  (display "]")
	  (display "|"))
      (display "\n"))))

;;;; Matrices

(define (make-mat rows cols)
  (vector-map (lambda (x) (make-vector cols 0))
	      (make-vector rows 0)))

(define (dim mat)
  (cons (vector-length mat)
	(vector-length (vector-ref mat 0))))

(define (rows mat)
  (car (dim mat)))

(define (cols mat)
  (cdr (dim mat)))

(define (mat-ref mat i j)
  (vector-ref (vector-ref mat i) j))

(define (mat-set! mat i j val)
  (vector-set! (vector-ref mat i) j val))

(define (row-ref mat n)
  (vector-ref mat n))

(define (row-set! mat r val)
  (vector-set! mat r val))

(define (col-ref mat n)
  (make-initialized-vector (rows mat)
			   (lambda (i) (mat-ref mat i n))))

(define (diag-ref mat)
  (make-initialized-vector (rows mat)
			   (lambda (i) (mat-ref mat i i))))

(define (mat-copy mat)
  (make-initialized-vector (rows mat)
			   (lambda (i) (vector-copy (row-ref mat i)))))

(define (mat->list mat)
  (vector->list
   (vector-map vector->list mat)))

(define (list->mat ls)
  (apply vector
	 (map (lambda (row) (apply vector row)) ls)))

(define (transpose mat)
  (make-initialized-vector (cols mat)
			   (lambda (i) (col-ref mat i))))

(define (dot v1 v2)
  (if (= (vector-length v1) (vector-length v2))
      (vector-fold + 0 (vector-map * v1 v2))
      (error "Invalid dimensions -- DOT")))

;; simple algorithm, O(n^3)
(define (mat*mat A B)
  (let ((ca (cols A))
	(cb (cols B))
	(ra (rows A))
	(rb (rows B)))
    (if (= ca rb)
	(let ((C (make-mat ra cb)))
	  (for (i from 0 to ra)
	    (for (j from 0 to cb)
	      (let ((sum 0))
		(for (k from 0 to ca)
		  (set! sum (+ sum (* (mat-ref A i k)
				      (mat-ref B k j)))))
		(mat-set! C i j sum))))
	  C)
	(error "Invalid dimensions -- MAT*MAT"))))

;; b interpreted as a column vector
(define (mat*vec A b)
  (if (= (cols A) (vector-length b))
      (let ((C (make-vector (rows A) 0)))
	(for (i from 0 to (rows A))
	  (vector-set! C i (dot (row-ref A i) b)))
	C)
      (error "Invalid dimensions")))


(define (mat+mat A B)
  (vector-map (lambda (x y) (vector-map + x y)) A B))

(define (mat*num A n)
  (mat-map (lambda (x) (* x n)) A))

;; negate every element in matrix
(define (negate mat)
  (mat-map (lambda (x) (- x)) mat))

(define (mat+num A n)
  (mat-map (lambda (x) (+ x n)) A))

(define (trace mat)
  (vector-fold + 0 (diag-ref mat)))

;; row operations
(define (row-scl! mat r p)
  (for (i from 0 to (cols mat))
    (mat-set! mat r i (* p (mat-ref mat r i)))))

(define (row-swp! mat r1 r2)
  (vector-swap! mat r1 r2))

(define (row-swp! mat r1 r2)
  (let ((tmp (row-ref mat r1)))
    (row-set! mat r1 (row-ref mat r2))
    (row-set! mat r2 tmp)))

(define (row-add! mat r1 r2)
  (row-set! mat r1
	    (vector-map + (row-ref mat r1) (row-ref mat r2))))

;;; augment a matrix with single vector
(define (aug-mat mat vec)
  (let ((R (make-mat (rows mat) (+ 1 (cols mat)))))
    (for (i from 0 to (rows R))
      (row-set! R i (vector-append (row-ref mat i)
				   (vector (vector-ref vec i)))))
    R))


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
(define (doolittle-decomp A)
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

;;; Determinant (uses Doolittle LU decomposition)
;; det(A) = det(LU) = det(L)*det(U)
;; Determinant of an upper or lower triangular
;; matrix is the product of the diagonal elements
;; Since L's diagonals are 1's, det(A) = det(U)
(define (det-doolittle mat)
  (let* ((U (cdr (doolittle-decomp mat)))
	 (diags (diag-ref U))
	 (prod 1))
    (for (i from 0 to (vector-length diags))
      (set! prod (* prod (vector-ref diags i))))
    prod))

;;; LUP decompose
;;; code translated from https://en.wikipedia.org/wiki/LU_decomposition
;;; Matrix is changed. Contains both matrices L-E and U
;;; as A=(L-E)+U such that P*A = L*U.
;;; P is not stored as a matrix, but in an int vector of size N+1
;;; containing column indexes where P has 1. Last element P[N] = S+N
;;; where S is the number of row exchanges needed for determinant
;;; computation det(P) = (-1)^S
(define (lup-decompose! A tol)
  (let* ((N (rows A))
	 (P (make-initialized-vector (+ N 1) (lambda (e) e))))
    (for (i from 0 to N)
      (let ((maxA 0.0)
	    (imax i))
	(for (k from i to N)
	  (if (> (abs (mat-ref A k i)) maxA)
	      (begin (set! maxA (abs (mat-ref A k i)))
		     (set! imax k))))
	(if (< maxA tol)
	    (error "Matrix is degenerate"))
	(if (not (equal? imax i))
	    (let ((j (vector-ref P i)))
	      (vector-set! P i (vector-ref P imax))
	      (vector-set! P imax j)
	      (row-swp! A i imax)
	      (vector-set! P N (+ 1 (vector-ref P N)))))
	(for (j from (+ i 1) to N)
	  (mat-set! A j i (/ (mat-ref A j i)
			     (mat-ref A i i)))
	  (for (k from (+ i 1) to N)
	    (mat-set! A j k (- (mat-ref A j k)
			       (* (mat-ref A j i)
				  (mat-ref A i k))))))))
    (cons A P)))

;;; solve A*x = b using lup-decompose
(define (solve-Ax-b A b)
  (let* ((n (rows A))
	 (tmp (mat-copy A))
	 (lup (lup-decompose! tmp 1e-16))
	 (mat (car lup))
	 (perm (cdr lup))
	 (x (make-vector n 0)))
    (for (i from 0 to n)
      (vector-set! x i (vector-ref b (vector-ref perm i)))
      (for (k from 0 to i)
	(vector-set! x i (- (vector-ref x i)
			    (* (mat-ref mat i k)
			       (vector-ref x k))))))
    (for (i from (- n 1) downto -1)
      (for (k from (+ i 1) to n)
	(vector-set! x i (- (vector-ref x i)
			    (* (mat-ref mat i k)
			       (vector-ref x k)))))
      (vector-set! x i (/ (vector-ref x i)
			  (mat-ref mat i i))))
    x))

;;; Invert a matrix using lup-decompose
(define (invert A)
  (let* ((n (rows A))
	 (tmp (mat-copy A))
	 (lup (lup-decompose! tmp 1e-16))
	 (mat (car lup))
	 (perm (cdr lup))
	 (inv (make-mat n n)))
    (for (j from 0 to n)
      (for (i from 0 to n)
	(if (= j (vector-ref perm i))
	    (mat-set! inv i j 1.0)
	    (mat-set! inv i j 0.0))
	(for (k from 0 to i)
	  (mat-set! inv i j (- (mat-ref inv i j)
			       (* (mat-ref mat i k)
				  (mat-ref inv k j))))))
      (for (i from (- n 1) downto -1)
	(for (k from (+ i 1) to n)
	  (mat-set! inv i j (- (mat-ref inv i j)
			       (* (mat-ref mat i k)
				  (mat-ref inv k j)))))
	(mat-set! inv i j (/ (mat-ref inv i j)
			     (mat-ref mat i i)))))
    inv))

;;; Determinant of a matrix using lup-decompose
(define (det A)
  (let* ((n (rows A))
	 (tmp (mat-copy A))
	 (lup (lup-decompose! tmp 1e-16))
	 (mat (car lup))
	 (perm (cdr lup))
	 (d (mat-ref mat 0 0)))
    (for (i from 1 to n)
      (set! d (* d (mat-ref mat i i))))
    (if (zero? (remainder (- (vector-ref perm n) n) 2))
	d
	(- d))))

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
