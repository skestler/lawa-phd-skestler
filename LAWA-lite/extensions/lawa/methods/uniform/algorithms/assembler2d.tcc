/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Kristina Steih, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

namespace lawa{
    
template<typename T, typename Basis>
Assembler2D<T, Basis>::
Assembler2D(Basis& _basis)
    : basis(_basis)
{        
}


template<typename T, typename Basis>
template<typename BilinearForm>
flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >
Assembler2D<T, Basis>::
assembleStiffnessMatrix(BilinearForm& a, int J_x, int J_y, T tol)
{   
    typedef typename Basis::FirstBasisType FirstBasis;
    typedef typename Basis::SecondBasisType SecondBasis;
    FirstBasis b1 = basis.first;
    SecondBasis b2 = basis.second;
    
    UniformIndex2D<Basis> I(basis, J_x, J_y);
    
    flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > A(basis.dim(J_x,J_y), basis.dim(J_x,J_y));
                                                 
     /* ============  v = Scaling Fct x Scaling Fct ==========================*/
     //std::cout << "===== v = SF * SF =======" << std::endl;
     Range<int> Rvx = b1.mra.rangeI(b1.j0);
     Range<int> Rvy = b2.mra.rangeI(b2.j0);
     Range<int> Rux = b1.mra.rangeI(b1.j0);
     Range<int> Ruy = b2.mra.rangeI(b2.j0);
     for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
       for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){

         /* u = Scaling Fct x Scaling Fct */ 
         Rux = b1.mra.rangeI(b1.j0);
         Ruy = b2.mra.rangeI(b2.j0);    
         for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
           for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){    
                         
             T val = a(XBSpline, b1.j0, kvx, XBSpline, b2.j0, kvy,
                       XBSpline, b1.j0, kux, XBSpline, b2.j0, kuy);
             if(fabs(val) > tol){
                 A(I(XBSpline, b1.j0, kvx, XBSpline, b2.j0, kvy), 
                   I(XBSpline, b1.j0, kux, XBSpline, b2.j0, kuy)) = val;
             }
             
           }
         }

         /* u = Scaling Fct x Wavelet */ 
         Rux = b1.mra.rangeI(b1.j0);
         for(int juy = b2.j0; juy <= basis.J2_max(J_x, J_y, b1.j0-1) - 1; ++juy){
           Ruy = b2.rangeJ(juy); 
           for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
             for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                 
               T val = a(XBSpline, b1.j0, kvx, XBSpline,  b2.j0, kvy,
                         XBSpline, b1.j0, kux, XWavelet, juy, kuy);
               if(fabs(val) > tol){
                   A(I(XBSpline, b1.j0, kvx, XBSpline, b2.j0, kvy),
                     I(XBSpline, b1.j0, kux, XWavelet, juy, kuy)) = val;
               }
               
             }
           }  
         }

         /* u = Wavelet x Scaling Function */ 
         Ruy = b2.mra.rangeI(b2.j0);
         for(int jux = b1.j0; jux <= basis.J1_max(J_x, J_y, b2.j0-1) - 1; ++jux){
           Rux = b1.rangeJ(jux); 
           for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
             for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
               
               T val = a(XBSpline, b1.j0, kvx, XBSpline, b2.j0, kvy,
                         XWavelet,  jux, kux, XBSpline, b2.j0, kuy);
               if(fabs(val) > tol){
                   A(I(XBSpline, b1.j0, kvx, XBSpline, b2.j0, kvy),
                     I(XWavelet, jux, kux, XBSpline, b2.j0, kuy)) = val;
               }

             }
           }
         }

         /* u = Wavelet x Wavelet */ 
         for(int jux = b1.j0; jux <= basis.J1_max(J_x, J_y, b2.j0) - 1; ++jux){
           Rux = b1.rangeJ(jux); 
           for(int juy = b2.j0; juy <= basis.J2_max(J_x, J_y, jux) - 1; ++juy){    
             Ruy = b2.rangeJ(juy); 
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                 
                 T val = a(XBSpline,  b1.j0, kvx, XBSpline, b2.j0, kvy,
                           XWavelet, jux, kux,  XWavelet, juy, kuy);
                 if(fabs(val) > tol){
                     A(I(XBSpline, b1.j0, kvx, XBSpline, b2.j0, kvy),
                       I(XWavelet, jux, kux, XWavelet, juy, kuy)) = val;
                 }

               }
             }
           }            
         }   
       }
     }    

     /* ============  v = Scaling Fct x Wavelet ==========================*/
    //std::cout << "===== v = SF * W =======" << std::endl;

     Rvx = b1.mra.rangeI(b1.j0);
     for(int jvy = b2.j0; jvy <= basis.J2_max(J_x, J_y, b1.j0-1) - 1; jvy++){
       Rvy = b2.rangeJ(jvy);
       for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
         for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){

           /* u = Scaling Fct x Scaling Fct */ 
           Rux = b1.mra.rangeI(b1.j0);
           Ruy = b2.mra.rangeI(b2.j0);    
           for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
             for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
               
               T val = a(XBSpline, b1.j0, kvx, XWavelet, jvy, kvy,
                         XBSpline, b1.j0, kux, XBSpline,  b2.j0, kuy);
               if(fabs(val) > tol){
                   A(I(XBSpline, b1.j0, kvx, XWavelet, jvy, kvy),
                     I(XBSpline, b1.j0, kux, XBSpline, b2.j0, kuy)) = val;
               }

             }
           }

           /* u = Scaling Fct x Wavelet */ 
           Rux = b1.mra.rangeI(b1.j0);
           for(int juy = b2.j0; juy <= basis.J2_max(J_x, J_y, b1.j0-1) - 1; ++juy){
             Ruy = b2.rangeJ(juy); 
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                
                 T val = a(XBSpline, b1.j0, kvx, XWavelet, jvy, kvy,
                           XBSpline, b1.j0, kux, XWavelet, juy, kuy);
                 if(fabs(val) > tol){  
                     A(I(XBSpline, b1.j0, kvx, XWavelet, jvy, kvy),
                       I(XBSpline, b1.j0, kux, XWavelet, juy, kuy)) = val;
                 }                
                
               }
             }  
           }

           /* u = Wavelet x Scaling Function */ 
           Ruy = b2.mra.rangeI(b2.j0);
           for(int jux = b1.j0; jux <= basis.J1_max(J_x, J_y, b2.j0-1) - 1; ++jux){
             Rux = b1.rangeJ(jux); 
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                   
                 T val = a(XBSpline,  b1.j0, kvx, XWavelet, jvy, kvy,
                           XWavelet, jux, kux,   XBSpline, b2.j0, kuy);
                 if(fabs(val) > tol){
                     A(I(XBSpline, b1.j0, kvx, XWavelet, jvy, kvy), 
                       I(XWavelet, jux, kux, XBSpline, b2.j0, kuy)) = val;
                 }
                 
               }
             }
           }

           /* u = Wavelet x Wavelet */ 
           for(int jux = b1.j0; jux <= basis.J1_max(J_x, J_y, b2.j0) -1; ++jux){
             Rux = b1.rangeJ(jux);
             for(int juy = b2.j0; juy <= basis.J2_max(J_x, J_y, jux) -1; ++juy){
               Ruy = b2.rangeJ(juy);
               for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
                 for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){

                   T val = a(XBSpline,  b1.j0, kvx, XWavelet, jvy, kvy,
                             XWavelet, jux, kux,   XWavelet, juy, kuy);
                   if(fabs(val) > tol){
                       A(I(XBSpline, b1.j0, kvx, XWavelet, jvy, kvy), 
                         I(XWavelet, jux, kux, XWavelet, juy, kuy)) = val;
                   }

                 }
               } 
             }
           }
         }
       }
     }

     /* ============  v = Wavelet x Scaling Fct ==========================*/
     //cout << "===== v = W * SF =======" << endl;

     Rvy = b2.mra.rangeI(b2.j0);
     for(int jvx = b1.j0; jvx <= basis.J1_max(J_x, J_y, b2.j0-1) - 1; jvx++){
       Rvx = b1.rangeJ(jvx);
       for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
         for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){

           /* u = Scaling Fct x Scaling Fct */ 
           Rux = b1.mra.rangeI(b1.j0);
           Ruy = b2.mra.rangeI(b2.j0);  
           for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
             for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
               
               T val = a(XWavelet, jvx, kvx,   XBSpline, b2.j0, kvy,
                         XBSpline,  b1.j0, kux, XBSpline, b2.j0, kuy);
               if(fabs(val) > tol){
                   A(I(XWavelet, jvx, kvx, XBSpline, b2.j0, kvy),
                     I(XBSpline, b1.j0, kux, XBSpline, b2.j0, kuy)) = val;
               }
               
             }
           }

           /* u = Scaling Fct x Wavelet */ 
           Rux = b1.mra.rangeI(b1.j0);
           for(int juy = b2.j0; juy <= basis.J2_max(J_x, J_y, b1.j0-1) - 1; ++juy){
             Ruy = b2.rangeJ(juy); 
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                 
                 T val = a(XWavelet, jvx, kvx,   XBSpline, b2.j0, kvy,
                           XBSpline,  b1.j0, kux, XWavelet, juy, kuy);
                 if(fabs(val) > tol){
                     A(I(XWavelet, jvx, kvx, XBSpline, b2.j0, kvy),
                       I(XBSpline, b1.j0, kux, XWavelet, juy, kuy)) = val;
                 }

               }
             }  
           }

           /* u = Wavelet x Scaling Function */ 
           Ruy = b2.mra.rangeI(b2.j0);
           for(int jux = b1.j0; jux <= basis.J1_max(J_x, J_y, b2.j0-1) - 1; ++jux){
             Rux = b1.rangeJ(jux); 
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                 
                 T val = a(XWavelet, jvx, kvx, XBSpline, b2.j0, kvy,
                           XWavelet, jux, kux, XBSpline, b2.j0, kuy);
                 if(fabs(val) > tol){
                     A(I(XWavelet, jvx, kvx, XBSpline, b2.j0, kvy),
                       I(XWavelet, jux, kux, XBSpline, b2.j0, kuy)) = val;
                 }
                 
               }
             }
           }

           /* u = Wavelet x Wavelet */ 
           for(int jux = b1.j0; jux <= basis.J1_max(J_x, J_y, b2.j0) -1; ++jux){
             Rux = b1.rangeJ(jux);
             for(int juy = b2.j0; juy <= basis.J2_max(J_x, J_y, jux) -1; ++juy){
               Ruy = b2.rangeJ(juy);
               for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
                 for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                   
                   T val = a(XWavelet, jvx, kvx, XBSpline,  b2.j0, kvy,
                             XWavelet, jux, kux, XWavelet, juy, kuy);
                   if(fabs(val) > tol){
                       A(I(XWavelet, jvx, kvx, XBSpline, b2.j0, kvy),
                         I(XWavelet, jux, kux, XWavelet, juy, kuy)) = val;
                   }

                 }
               } 
             }
           }
         }
       }
     }   

     /* ============  v = Wavelet x Wavelet ==========================*/
     //cout << "===== v = W * W =======" << endl;  
     for(int jvx = b1.j0; jvx <= basis.J1_max(J_x, J_y, b2.j0) -1; ++jvx){
       Rvx = b1.rangeJ(jvx);
       for(int jvy = b2.j0; jvy <= basis.J2_max(J_x, J_y, jvx) -1; ++jvy){
         Rvy = b2.rangeJ(jvy);
         for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
           for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){

             /* u = Scaling Fct x Scaling Fct */  
             Rux = b1.mra.rangeI(b1.j0);
             Ruy = b2.mra.rangeI(b2.j0);
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){

                 T val = a(XWavelet, jvx, kvx,   XWavelet, jvy, kvy,
                           XBSpline,  b1.j0, kux, XBSpline,  b2.j0, kuy);
                 if(fabs(val) > tol){
                     A(I(XWavelet, jvx, kvx, XWavelet, jvy, kvy),
                       I(XBSpline,  b1.j0, kux, XBSpline, b2.j0, kuy)) = val;
                 }

               }
             }

             /* u = Scaling Fct x Wavelet */ 
             Rux = b1.mra.rangeI(b1.j0);
             for(int juy = b2.j0; juy <= basis.J2_max(J_x, J_y, b1.j0-1) - 1; ++juy){
               Ruy = b2.rangeJ(juy); 
               for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
                 for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                   
                   T val = a(XWavelet, jvx, kvx,  XWavelet, jvy, kvy,
                             XBSpline, b1.j0, kux, XWavelet, juy, kuy);
                   if(fabs(val) > tol){
                       A(I(XWavelet, jvx, kvx, XWavelet, jvy, kvy),
                         I(XBSpline, b1.j0, kux, XWavelet, juy, kuy)) = val;
                   }

                 }
               }  
             }

             /* u = Wavelet x Scaling Function */ 
             Ruy = b2.mra.rangeI(b2.j0);
             for(int jux = b1.j0; jux <= basis.J1_max(J_x, J_y, b2.j0-1) - 1; ++jux){
               Rux = b1.rangeJ(jux); 
               for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
                 for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                 
                   T val = a(XWavelet, jvx, kvx, XWavelet, jvy, kvy,
                             XWavelet, jux, kux, XBSpline,  b2.j0, kuy);
                   if(fabs(val) > tol){
                       A(I(XWavelet, jvx, kvx, XWavelet, jvy, kvy),
                         I(XWavelet, jux, kux, XBSpline, b2.j0, kuy)) = val;
                   }  

                 }
               }
             }

             /* u = Wavelet x Wavelet */ 
             for(int jux = b1.j0; jux <= basis.J1_max(J_x, J_y, b2.j0) -1; ++jux){
               Rux = b1.rangeJ(jux);
               for(int juy = b2.j0; juy <= basis.J2_max(J_x, J_y, jux) -1; ++juy){
                 Ruy = b2.rangeJ(juy);
                 for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
                   for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){

                     T val = a(XWavelet, jvx, kvx, XWavelet, jvy, kvy,
                               XWavelet, jux, kux, XWavelet, juy, kuy);
                     if(fabs(val) > tol){
                         A(I(XWavelet, jvx, kvx, XWavelet, jvy, kvy),
                           I(XWavelet, jux, kux, XWavelet, juy, kuy)) = val;
                     }

                   }
                 } 
               }
             }           
           }
         } 
       }
     }
    
    A.finalize();
    return A;
 
}  


template<typename T, typename Basis>
template<typename RHSIntegral>
flens::DenseVector<flens::Array<T> >
Assembler2D<T, Basis>::assembleRHS(RHSIntegral& rhs, int J_x, int J_y)
{
    typedef typename Basis::FirstBasisType FirstBasis;
    typedef typename Basis::SecondBasisType SecondBasis;
    FirstBasis b1 = basis.first;
    SecondBasis b2 = basis.second;
    
    UniformIndex2D<Basis> I(basis, J_x, J_y);
     
    flens::DenseVector<flens::Array<T> > f(basis.dim(J_x, J_y));

    /*  ============  v = Scaling Fct x Scaling Fct ==========================*/
    //std::cout << "SF x SF : " << std::endl;
    Range<int> Rvx = b1.mra.rangeI(b1.j0);
    Range<int> Rvy = b2.mra.rangeI(b2.j0);
    for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
      for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){
          
          f(I(XBSpline, b1.j0, kvx, XBSpline, b2.j0, kvy)) 
              = rhs(XBSpline, b1.j0, kvx, XBSpline, b2.j0, kvy);

      }
    }

    /* ============  v = Scaling Fct x Wavelet ==========================*/
    //std::cout << "SF x W : " << std::endl;
    Rvx = b1.mra.rangeI(b1.j0);
    for(int jvy = b2.j0; jvy <= basis.J2_max(J_x, J_y, b1.j0-1) - 1; ++jvy){
      Rvy = b2.rangeJ(jvy);
      for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
        for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){
        
          f(I(XBSpline, b1.j0, kvx, XWavelet, jvy, kvy)) 
              = rhs(XBSpline, b1.j0, kvx, XWavelet, jvy, kvy);
                     
        }
      }  
    }

    /* ============  v = Wavelet x Scaling Fct ==========================*/
    //std::cout << "W x SF : " << std::endl;
    Rvy = b2.mra.rangeI(b2.j0);
    for(int jvx = b1.j0; jvx <= basis.J1_max(J_x, J_y, b2.j0-1) - 1; ++jvx){
      Rvx = b1.rangeJ(jvx);
      for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
        for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){
            
          f(I(XWavelet, jvx, kvx, XBSpline, b2.j0, kvy)) 
              = rhs(XWavelet, jvx, kvx, XBSpline, b2.j0, kvy);
          
        }
      }  
    }

    /*  ============  v = Wavelet x Wavelet ==========================*/
    //std::cout << "W x W : " << std::endl;
    for(int jvx = b1.j0; jvx <= basis.J1_max(J_x, J_y, b2.j0) -1; ++jvx){
      Rvx = b1.rangeJ(jvx);
      for(int jvy = b2.j0; jvy <= basis.J2_max(J_x, J_y, jvx) - 1; ++jvy){
        Rvy = b2.rangeJ(jvy);
        for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
          for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){
            
            f(I(XWavelet, jvx, kvx, XWavelet, jvy, kvy)) 
                = rhs(XWavelet, jvx, kvx, XWavelet, jvy, kvy);

          }
        }
      }
    }
    
    return f;
}


template<typename T, typename Basis>
template<typename Preconditioner>
flens::DiagonalMatrix<T>    
Assembler2D<T, Basis>::
assemblePreconditioner(Preconditioner& P, int J_x, int J_y)
{
    typedef typename Basis::FirstBasisType FirstBasis;
    typedef typename Basis::SecondBasisType SecondBasis;
    FirstBasis b1 = basis.first;
    SecondBasis b2 = basis.second;
    
    UniformIndex2D<Basis> I(basis, J_x, J_y);
    
    flens::DenseVector<flens::Array<T> > D(basis.dim(J_x, J_y));

    /* SF x SF */
    Range<int> Rx = b1.mra.rangeI(b1.j0);
    Range<int> Ry = b2.mra.rangeI(b2.j0);   
    for(int kx = Rx.firstIndex(); kx <= Rx.lastIndex(); ++kx){
        for(int ky = Ry.firstIndex(); ky <= Ry.lastIndex(); ++ky){
            D(I(XBSpline, b1.j0, kx, XBSpline, b2.j0, ky)) = P(XBSpline, b1.j0, kx, XBSpline, b2.j0, ky);
            
        }
    }
    
    /* SF x W */
    Rx = b1.mra.rangeI(b1.j0);
    for(int jy = b2.j0; jy <= basis.J2_max(J_x, J_y, b1.j0-1)-1; ++jy){
        Ry = b2.rangeJ(jy);
        for(int kx = Rx.firstIndex(); kx <= Rx.lastIndex(); ++kx){
            for(int ky = Ry.firstIndex(); ky <= Ry.lastIndex(); ++ky){
                D(I(XBSpline, b1.j0, kx, XWavelet, jy, ky)) = P(XBSpline, b1.j0, kx, XWavelet, jy, ky);                
            }
        }
    }
    
    /* W x SF */
    Ry = b2.mra.rangeI(b2.j0);
    for(int jx = b1.j0; jx <= basis.J1_max(J_x, J_y, b2.j0-1) -1; ++jx){
        Rx = b1.rangeJ(jx);
        for(int kx = Rx.firstIndex(); kx <= Rx.lastIndex(); ++kx){
            for(int ky = Ry.firstIndex(); ky <= Ry.lastIndex(); ++ky){
                D(I(XWavelet, jx, kx, XBSpline, b2.j0, ky)) = P(XWavelet, jx, kx, XBSpline, b2.j0, ky);            
            }
        }      
    }
    
    /* W x W */
    for(int jx = b1.j0; jx <= basis.J1_max(J_x, J_y, b2.j0) - 1; ++jx){
        Rx = b1.rangeJ(jx);
        for(int jy = b2.j0; jy <= basis.J2_max(J_x, J_y, jx) - 1; ++jy){
            Ry = b2.rangeJ(jy);
            for(int kx = Rx.firstIndex(); kx <= Rx.lastIndex(); ++kx){
                for(int ky = Ry.firstIndex(); ky <= Ry.lastIndex(); ++ky){
                    D(I(XWavelet, jx, kx, XWavelet, jy, ky)) = P(XWavelet, jx, kx, XWavelet, jy, ky);                
                }
            } 
        }
    }
    
    return flens::DiagonalMatrix<T>(D);
    
}
    
} // namespace lawa

