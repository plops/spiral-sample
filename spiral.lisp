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

;; arc length archimedes spiral s(t)=a/2(p t+ln(p t)) with p=sqrt(1+t^2)
;; f(x) = a/2 (p x + ln(p x)) - s
;; f'(x) = a p
;; newton iteration to invert: x_{n+1}=x_n-f(x_n)/f'(x_n)
(defun f-over-df (s a x)
  (declare (single-float s a x)
	   (values single-float &optional))
  (let ((p (sqrt (1+ (* x x)))))
   (- (* .5 (+ x (/ (log (* p x)) p)))
      (/ s (* a p)))))

(defun find-zero (s a &optional (x 1.0))
  (declare (single-float s a x)
	   (values single-float &optional))
  (dotimes (i 10)
    (format t "~a~%" (incf x (- (f-over-df s a x)))))
  x)

(find-zero 13.0 1.0 1.0)

(defun draw-spiral (&key (h 256) (w 256) (points/2 100) (turns/2 4.0))
  (declare (fixnum h w points/2)
	   (single-float turns/2)
	   (values (simple-array (unsigned-byte 8) 2) &optional))
  (let ((img (make-array (list h w) :element-type '(unsigned-byte 8)))
	(wh (floor w 2))
	(hh (floor h 2))
	(a (* .01 w)))
    (with-arrays (img)
      (loop for p from (- points/2) upto points/2 do
	   (let* ((phi (* 2.0 pi (abs p) turns/2 (/ points/2)))
		  (z (* a phi (exp (complex 0.0 phi))))
		  (x (realpart z))
		  (y (imagpart z))
		  (i (floor (+ x wh)))
		  (j (floor (+ y hh))))
	     (setf (img j i) 255))))
    img))

(write-pgm "/dev/shm/o.pgm" (draw-spiral))