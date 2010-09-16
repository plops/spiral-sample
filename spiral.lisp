(defmacro with-arrays (arrays &body body)
  "Provides a corresponding accessor for each array as a local macro,
so that (ARRAY ...) corresponds to (AREF ARRAY ...)."
  `(macrolet ,(mapcar (lambda (array)
                        `(,array (&rest indices) `(aref ,',array ,@indices)))
                      arrays)
     ,@body))

(defmacro do-region ((indices end &optional (start '(0 0 0))) &body body)
  "Write intertwined loops to traverse a vector, an image or a volume."
  (unless (and (= (length indices)
                  (length end)))
    (error "Number of indices and interval-ends are not equal."))
  (labels ((rec (ind end start acc) ;; several loops
             (if (null end)
                 acc
                 (rec (cdr ind) (cdr end) (cdr start)
                      `((loop for ,(car ind) from ,(car start) 
                           below ,(car end) do ,@acc))))))
    (first (rec (reverse indices) ;; first index is outermost loop
                (reverse end)
                (reverse start) body))))

(defun write-pgm (filename img)
  (declare (string filename)
           ((array (unsigned-byte 8) 2) img)
           (values))
  (destructuring-bind (h w)
      (array-dimensions img)
    (declare ((integer 0 65535) w h))
    (with-open-file (s filename
                       :direction :output
                       :if-exists :supersede
                       :if-does-not-exist :create)
      (declare (stream s))
      (format s "P5~%~D ~D~%255~%" w h))
    (with-open-file (s filename 
                       :element-type '(unsigned-byte 8)
                       :direction :output
                       :if-exists :append)
      (let ((img1 (sb-ext:array-storage-vector img)))
        (write-sequence img1 s)))
    (values)))

(defun arc-length (a phi)
  (declare (single-float a phi)
	   (values single-float &optional))
  (let ((q (* phi (sqrt (1+ (* phi phi))))))
    (* .5 a (+ q (log q)))))

;; arc length archimedes spiral s(t)=a/2(p t+ln(p t)) with p=sqrt(1+t^2)
;; f(x) = a/2 (p x + ln(p x)) - s
;; f'(x) = a p
;; newton iteration to invert: x_{n+1}=x_n-f(x_n)/f'(x_n)

(defun find-zero (s a &optional (x 1.0))
  (declare (single-float s a x)
	   (values single-float &optional))
  (dotimes (i 12)
    (let* ((x2 (* x x))
	   (p (sqrt (1+ x2)))
	   (q (* p x))
	   (r (+ q (log q)))
	   (a/2 (* .5 a))
	   (f (- (* a/2 r) s))
	   (df (* a/2 (/ (* (1+ (* 2 x2)) (1+ q)) (* x p p)))))
      (incf x (- (/ f df)))))
  x)

#+nil
(find-zero 40.0 1.0 1.0)

;; solve 2d non-linear equations f(x,y)=0 g(x,y)=0 by newton method
;; jacobian J=((fx fy)(gx gy))
;; matrix inverse ((a b)(c d))^-1 = ((d -b)(-c a))/(ad-bc)
;; step P1=P0-J(P0)^-1 F(P0) with F=(f g), P0=(x y)
;; FindRoot[{a x == 100, a/2 (Sqrt[1+x^2]x+Log[Sqrt[1+x^2]x]) == 200 * 2 pi + 100/a}, {{x, 1}, {a, 1}}]
;; example n=200, bfp-radius=100 -> theta=25.3654, a=3.94238

(defun find-zero2 (&key (n 200) (bfp-radius 100.0))
  (declare (fixnum n)
	   (single-float bfp-radius)
	   (values single-float single-float &optional))
  (let* ((a 10.0)
	 (theta 15.0)
	 (pif #.(coerce pi 'single-float)))
    (declare ((single-float 0.0) theta a))
    (dotimes (i 14)
      (let* ((theta^2 (* theta theta))
	     (p (sqrt (1+ theta^2)))
	     (q (* p theta))
	     (r (+ q (log q)))
	     (f1 (- (* .5 a r) (* 2.0 pif n) (/ bfp-radius a)))
	     (f2 (- (* a theta) bfp-radius))
	     (f1t (/ (* .5 a (1+ (* 2 theta^2)) (1+ q)) 
		     (* p q)))
	     (f1a (+ (/ bfp-radius (* a a))
		     (* .5 r)))
	     (1/det (/ (- (* theta f1t) (* a f1a)))))
	(declare (single-float p))
	(incf theta (* 1/det (- (* f1a f2) (* theta f1))))
	(incf a (* 1/det (- (* a f1) (* f1t f2))))
	(format t "~a~%" (list a theta))))
    (values a theta)))

#+nil
(find-zero2 :n 10 :bfp-radius 7.0)

(defun clamp (x)
  (declare (fixnum x)
	   (values (unsigned-byte 8) &optional))
  (cond ((< x 0) 0)
	((< 255 x) 255)
	(t x)))

(defun draw-spiral (&key (h 256) (w 256) (points 100) (bfp-radius 100.0))
  (declare (fixnum h w points)
	   (single-float bfp-radius)
	   (values (simple-array (unsigned-byte 8) 2) &optional))
  (let* ((img (make-array (list h w) :element-type '(unsigned-byte 8)))
	 (wh (floor w 2))
	 (hh (floor h 2))
	 (a (find-zero2 :n points :bfp-radius bfp-radius))
	 (ds (* a 2.0 (coerce pi 'single-float))))
    (with-arrays (img)
      (dotimes (p points)
	(let* ((phi (find-zero (* p ds) a))
	       (z (* a phi (exp (complex 0.0 phi))))
	       (x (realpart z))
	       (y (imagpart z))
	       (i (clamp (floor (+ x wh))))
	       (j (clamp (floor (+ y hh)))))
	  (setf (img j i) 255))))
    img))
#+nil
(write-pgm "/dev/shm/o.pgm" (draw-spiral :points 600 :bfp-radius 80.0))

