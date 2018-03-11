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

;;; Chicken specific
(require-extension srfi-133) ; vector operations
(require-extension extras)   ; for random

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

;;; works only for square matrices as of now
(define (diag-ref mat)
  (vector-unfold (lambda (i)
		   (mat-ref mat i i)) (cols mat)))

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
  (vector-fold + 0 (diag-ref mat)))

(define (print-mat mat)
  (for (i from 0 to (rows mat))
    (display (row-ref mat i))
    (newline)))

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
  (vector-swap! mat r1 r2))

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

;;; reduce matrix to row echelon form
(define (forward-elim mat)
  (for (k from 0 to (rows mat))
    (let* ((i-max k)
	   (v-max (mat-ref mat i-max k)))
      (for (i from (+ k 1) to (rows mat))
	(if (> (abs (mat-ref mat i k)) v-max)
	    (set! v-max (mat-ref mat i k))
	    (set! i-max i)))
      (if (zero? (mat-ref mat k i-max))
	  k)
      (if (not (equal? i-max k))
	  (row-swp! mat k i-max))
      (for (i from (+ k 1) to (rows mat))
	(let ((f (/ (mat-ref mat i k) (mat-ref mat k k))))
	  (for (j from (+ k 1) to (rows mat))
	    (mat-set! mat i j (- (mat-ref mat i j)
			     (* (mat-ref mat k j) f))))
	  (mat-set! mat i k 0)))))
  -1)

;;; calculate value of unknowns
(define (back-sub mat)
  (let ((res (make-vector (rows mat) 0)))
    (for (i from (- (rows mat) 1) downto -1)
      (vector-set! res i (mat-ref mat i (rows mat)))
      (for (j from (+ i 1) to (rows mat))
	(vector-set! res i (- (vector-ref res i)
			      (* (mat-ref mat i j)
				 (vector-ref res j)))))
      (vector-set! res i (/ (vector-ref res i)
			    (mat-ref mat i i))))
    res))

;;; Gaussian elimination
;;; solve A*x = b
(define (gaussian-elim A b)
  (let* ((aug (aug-mat A b))
	 (singular-flag (forward-elim aug)))
    (if (equal? singular-flag -1)
	(back-sub aug)
	(begin (display "Singular Matrix") (newline)))))

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

