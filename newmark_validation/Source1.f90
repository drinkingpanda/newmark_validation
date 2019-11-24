program  test_newmark
   use ODE_solver
    implicit none

    real*8  :: delt_tt ,t,tmp(3)
    real*8  :: M(5),C(5),K(5),b(3),X(3),X1(3),X2(3),Y(3),Y1(3),Y2(3)
    real*8  :: pi = acos(-1.0d0)
    real*8  :: w1,w2,w3,Q
    
    integer :: ia(4),ja(5)
    integer :: nnz,N,lb
    integer :: i ,P, N_max
    
    P = 0.5
   
    !construct analytical result: [u1,u2,u3]=sin(w1*t)[1,0,0] + sin(w2*t)*[1,1,1] + sin(w3*t)*[0,0,1]
    w1 = 1.0d0
    w2 = 2.0d0
    w3 = 30.0d0
    
    delt_tt = 2*pi/max(w1,w2,w3)/0.5
    
    
    ia = [1,3,4,6]
    ja = [1,2,2,2,3]
    
    C  = 0.0d0
    
    M  = [1.0d0,0.0d0,        &
      &         1.0d0,        &
      &         0.0d0,1.0d0   ]
    
    
    K  = [ w1**2,  w2**2 - w1**2 ,            &
      &                    w2**2 ,            &
      &            w2**2 - w3**2 ,     w3**2  ]
    
    b  = [0.0d0,0.0d0,0.0d0]
    
    N   = 3
    lb  = 1
    nnz = 5
    
    ! X is old value,Y is new value
    X = [  0.0d0 ,  0.0d0 ,  0.0d0 ]
    X1= [  w1+w2 ,     w2 ,  w2+w3 ]
    X2= [  0.0d0 ,  0.0d0 ,  0.0d0 ]

    !new value
    Y = 0.0d0
    Y1= 0.0d0
    Y2= 0.0d0
    
    N_max = 2000

    open (901,file = 'data.txt')
    
    write(901,106)'U1','U2','U3','U1_exact','U2_exact','U3_exact'
    
    do i=1,N_max
        
        print*,"iteration: ",i
        
        call newMark_bsr( N,lb,nnz,delt_tt ,ia,ja,M,C,K,b,X,X1,X2,Y,Y1,Y2)
         
        X  = Y
        X1 = Y1
        X2 = Y2
        
        t   = (i-1)*delt_tt
        
        tmp = 1.0d0*sin(w1*t)*[1.0d0,0.0d0,0.0d0]      &
         &   +1.0d0*sin(w2*t)*[1.0d0,1.0d0,1.0d0]      &
         &   +1.0d0*sin(w3*t)*[0.0d0,0.0d0,1.0d0]
        

        write(901,116)X(1:3),tmp(1:3)
 
 
        
    enddo
    
    close(901)
    
 106 format(  6( a10   ,2x)  )
 116 format( 6( f10.7 ,2x ) )
    
    
    
    
end program test_newmark