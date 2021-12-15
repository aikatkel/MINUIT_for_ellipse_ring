C     ===================
      Program Test_Minuit
C     ===================
      Implicit DoublePrecision (A-H,O-Z)

      Character*8 CVar
      
      External FCN
      External Futil

C     ================
C     Interactive Call
C     ================

      Call Minuit(FCN,Futil)

      Stop
      End

C     ===============================================
      Subroutine FCN(NPAR,GRAD,FVAL,XVAL,IFLAG,FUTIL)
C     ===============================================

      Implicit DoublePrecision (A-H,O-Z)

      External Futil

!     Parameter (Np=100)

      Dimension YY(512,512),new_img(512,512)
      Dimension XVal(*),Grad(*)


C     +------------------+
C     |  Initialization  |
C     +------------------+


      If (IFLAG.EQ.1) Then

C     =========
C     Read Data
C     =========

          Open (Unit=11,file='big.txt',Status='OLD')

          Do j = 1,512
             Read (11,*) (YY(i,j), i=1,512)
          Enddo

          Close (11)

          Write (*,*) 'Matrix successfully loaded'
          max_inp = 0.0d0
          Do i=1,512
              Do j=1,512
                  If (YY(i,j).GT.max_inp) Then
                  	max_inp = YY(i,j)
                  EndIf
              EndDo
          EndDo
          
          Do i=1,512
              Do j=1,512
                  YY(i,j) = YY(i,j)/max_inp
                  If (YY(i,j)>0.4d0) Then
                      YY(i,j) = 1.0d0
                  Else
                      YY(i,j) = 0.0d0
                  EndIf
              EndDo
          EndDo
          
      EndIf

C     +-------------------+
C     |  Calclulate Chi2  |
C     +-------------------+

      x0     = Xval(1)
      y0     = Xval(2)
      sx     = Xval(3)
      sy     = Xval(4)
      drx    = Xval(5)
      dry    = Xval(6)
      phi    = Xval(7)

      Chi2 = 0.0d0
      
      
      
      Do i=1,512
          Do j=1,512
              x=j
              y=i
              new_img(i,j) = 0.0d0
              rr = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0))
              the = datan2(y-y0,x-x0)
              part_a = sy*sy*cos(the-phi)*cos(the-phi)
              part_b = sx*sx*sin(the-phi)*sin(the-phi)
              srr = (sx*sy) / dsqrt(part_a + part_b)
              p_a = (sy+dry)*(sy+dry)*cos(the-phi)*cos(the-phi)
              p_b = (sx+drx)*(sx+drx)*sin(the-phi)*sin(the-phi)
              srdr = (sx+drx)*(sy+dry)/dsqrt(p_a + p_b)
              
              If ((rr.LE.srdr).AND.(rr.GE.srr)) Then
                  new_img(i,j) = 1.0d0
              Endif
          EndDo
      EndDo
      FVal = 0.0d0
      
      Do i=1,512
          Do j=1,512
              diff = (new_img(i,j)-YY(i,j))
              FVal = FVal + (diff*diff)
          EndDo
      EndDo
      
      
C     +--------------------------------+
C     |  Fit Ended - Write Parameters  |
C     +--------------------------------+    


      If (IFLAG.EQ.3) Then
      
         Write (*,*)
         
         Do i=1,NPAR
            Write (*,101) i,XVal(i)
         EndDo
         
         Write (*,*)
         
  101 Format (4x,I2,'. Parameter =',F12.4)
  
      EndIf
      
      Return
      End

C     ================
      Subroutine Futil
C     ================
      
      Return
      End
