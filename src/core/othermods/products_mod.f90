MODULE M_products
!Contains functions to calculate the scalar and vector products

USE GLOBAL

IMPLICIT NONE

  PRIVATE
  PUBLIC :: CROSS, DOT
  
 CONTAINS
!-----------------------------------------------------------------------------!   
 FUNCTION CROSS(A,B)
  IMPLICIT NONE
 !+ vector product

  REAL(num),DIMENSION(3) :: A,B,CROSS

  CROSS(1) = A(2)*B(3) - A(3)*B(2)
  CROSS(2) = A(3)*B(1) - A(1)*B(3)
  CROSS(3) = A(1)*B(2) - A(2)*B(1)
 
 END FUNCTION CROSS
!-----------------------------------------!   
 FUNCTION DOT(A,B)
  IMPLICIT NONE
 !+scalar product
  REAL(num), DIMENSION(3) :: A,B
  REAL(num) :: DOT

  DOT = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)

 END FUNCTION DOT
!-----------------------------------------!   
END MODULE M_products
