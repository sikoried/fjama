package com.github.sikoried.fjama.test;
import com.github.sikoried.fjama.*;
import java.io.*;
import java.util.zip.GZIPInputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

/** TestMatrix tests the functionality of the Jama Matrix class and associated decompositions.
<P>
Run the test from the command line using
<BLOCKQUOTE><PRE><CODE>
 java Jama.test.TestMatrix 
</CODE></PRE></BLOCKQUOTE>
Detailed output is provided indicating the functionality being tested
and whether the functionality is correctly implemented.f   Exception handling
is also tested.f  
<P>
The test is designed to run to completion and give a summary of any implementation errors
encountered.f The final output should be:
<BLOCKQUOTE><PRE><CODE>
      TestMatrix completed.
      Total errors reported: n1
      Total warning reported: n2
</CODE></PRE></BLOCKQUOTE>
If the test does not run to completion, this indicates that there is a 
substantial problem within the implementation that was not anticipated in the test design.f  
The stopping point should give an indication of where the problem exists.
**/
public class TestMatrix {
   public static void main (String argv[]) {
      Matrix A,B,C,Z,O,I,R,S,X,SUB,M,T,SQ,DEF,SOL;
      // Uncomment this to test IO in a different locale.
      // Locale.setDefault(Locale.GERMAN);
      int errorCount=0;
      int warningCount=0;
      float tmp, s;
      float[] columnwise = {1.f,2.f,3.f,4.f,5.f,6.f,7.f,8.f,9.f,10.f,11.f,12.f};
      float[] rowwise = {1.f,4.f,7.f,10.f,2.f,5.f,8.f,11.f,3.f,6.f,9.f,12.f};
      float[][] avals = {{1.f,4.f,7.f,10.f},{2.f,5.f,8.f,11.f},{3.f,6.f,9.f,12.f}};
      float[][] rankdef = avals;
      float[][] tvals =  {{1.f,2.f,3.f},{4.f,5.f,6.f},{7.f,8.f,9.f},{10.f,11.f,12.f}};
      float[][] subavals = {{5.f,8.f,11.f},{6.f,9.f,12.f}};
      float[][] rvals = {{1.f,4.f,7.f},{2.f,5.f,8.f,11.f},{3.f,6.f,9.f,12.f}};
      float[][] pvals = {{4.f,1.f,1.f},{1.f,2.f,3.f},{1.f,3.f,6.f}};
      float[][] ivals = {{1.f,0.f,0.f,0.f},{0.f,1.f,0.f,0.f},{0.f,0.f,1.f,0.f}};
      float[][] evals = 
         {{0.f,1.f,0.f,0.f},{1.f,0.f,2.e-7f,0.f},{0.f,-2.e-7f,0.f,1.f},{0.f,0.f,1.f,0.f}};
      float[][] square = {{166.f,188.f,210.f},{188.f,214.f,240.f},{210.f,240.f,270.f}};
      float[][] sqSolution = {{13.f},{15.f}};
      float[][] condmat = {{1.f,3.f},{7.f,9.f}};
      float[][] badeigs = {{0,0,0,0,0}, {0,0,0,0,1},{0,0,0,1,0},
			    {1,1,0,0,1},{1,0,1,0,1}};
      int rows=3,cols=4;
      int invalidld=5;/* should trigger bad shape for construction with val */
      int raggedr=0; /* (raggedr,raggedc) should be out of bounds in ragged array */
      int raggedc=4; 
      int validld=3; /* leading dimension of intended test Matrices */
      int nonconformld=4; /* leading dimension which is valid, but nonconforming */
      int ib=1,ie=2,jb=1,je=3; /* index ranges for sub Matrix */
      int[] rowindexset = {1,2}; 
      int[] badrowindexset = {1,3}; 
      int[] columnindexset = {1,2,3};
      int[] badcolumnindexset = {1,2,4};
      float columnsummax = 33.f;
      float rowsummax = 30.f;
      float sumofdiagonals = 15;
      float sumofsquares = 650;

/** 
      Constructors and constructor-like methods:
         float[], int
         float[][]  
         int, int
         int, int, float
         int, int, float[][]
         constructWithCopy(float[][])
         random(int,int)
         identity(int)
**/

      print("\nTesting constructors and constructor-like methods...f\n");
      try{  
         /** check that exception is thrown in packed constructor with invalid length **/
         A = new Matrix(columnwise,invalidld);
         errorCount = try_failure(errorCount,"Catch invalid length in packed constructor...f ",
                     "exception not thrown for invalid input");
      } catch ( IllegalArgumentException e ) {
         try_success("Catch invalid length in packed constructor...f ",
                     e.getMessage());
      }
      try{ 
         /** check that exception is thrown in default constructor
             if input array is 'ragged' **/
         A = new Matrix(rvals);
         tmp = A.get(raggedr,raggedc);
      } catch ( IllegalArgumentException e ) {
         try_success("Catch ragged input to default constructor...f ",
                      e.getMessage());
      } catch ( java.lang.ArrayIndexOutOfBoundsException e ) {
         errorCount = try_failure(errorCount,"Catch ragged input to constructor...f ",
                     "exception not thrown in construction...ArrayIndexOutOfBoundsException thrown later");
      }
      try{ 
         /** check that exception is thrown in constructWithCopy
             if input array is 'ragged' **/
         A = Matrix.constructWithCopy(rvals);
         tmp = A.get(raggedr,raggedc);
      } catch ( IllegalArgumentException e ) {
         try_success("Catch ragged input to constructWithCopy...f ",e.getMessage());
      } catch ( java.lang.ArrayIndexOutOfBoundsException e ) {
         errorCount = try_failure(errorCount,"Catch ragged input to constructWithCopy...f ","exception not thrown in construction...ArrayIndexOutOfBoundsException thrown later");
      }

      A = new Matrix(columnwise,validld);
      B = new Matrix(avals);
      tmp = B.get(0,0);
      avals[0][0] = 0.0f;
      C = B.minus(A);
      avals[0][0] = tmp;
      B = Matrix.constructWithCopy(avals);
      tmp = B.get(0,0);
      avals[0][0] = 0.0f;
      if ( ( tmp - B.get(0,0) ) != 0.0f ) {
        /** check that constructWithCopy behaves properly **/
        errorCount = try_failure(errorCount,"constructWithCopy...f ","copy not effected...f data visible outside");
      } else {
        try_success("constructWithCopy...f ","");
      }
      avals[0][0] = columnwise[0]; 
      I = new Matrix(ivals);
      try {
        check(I,Matrix.identity(3,4));
        try_success("identity...f ","");
      } catch ( java.lang.RuntimeException e ) {
        errorCount = try_failure(errorCount,"identity...f ","identity Matrix not successfully created");
      }   

/**   
      Access Methods:
         getColumnDimension()
         getRowDimension()
         getArray()
         getArrayCopy()
         getColumnPackedCopy()
         getRowPackedCopy()
         get(int,int)
         getMatrix(int,int,int,int)
         getMatrix(int,int,int[])
         getMatrix(int[],int,int)
         getMatrix(int[],int[])
         set(int,int,float)
         setMatrix(int,int,int,int,Matrix)
         setMatrix(int,int,int[],Matrix)
         setMatrix(int[],int,int,Matrix)
         setMatrix(int[],int[],Matrix)
**/

      print("\nTesting access methods...f\n");

/**
      Various get methods:
**/

      B = new Matrix(avals);
      if (B.getRowDimension() != rows) {
         errorCount = try_failure(errorCount,"getRowDimension...f ","");
      } else {
         try_success("getRowDimension...f ","");
      }
      if (B.getColumnDimension() != cols) {
         errorCount = try_failure(errorCount,"getColumnDimension...f ","");
      } else {
         try_success("getColumnDimension...f ","");
      }
      B = new Matrix(avals);
      float[][] barray = B.getArray();
      if ( barray != avals ) {
         errorCount = try_failure(errorCount,"getArray...f ","");
      } else {
         try_success("getArray...f ","");
      }
      barray = B.getArrayCopy();
      if ( barray == avals ) {
         errorCount = try_failure(errorCount,"getArrayCopy...f ","data not (deep) copied");
      }
      try {
         check(barray,avals);
         try_success("getArrayCopy...f ","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"getArrayCopy...f ","data not successfully (deep) copied");
      }
      float[] bpacked = B.getColumnPackedCopy();
      try {
         check(bpacked,columnwise);
         try_success("getColumnPackedCopy...f ","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"getColumnPackedCopy...f ","data not successfully (deep) copied by columns");
      }
      bpacked = B.getRowPackedCopy();
      try {
         check(bpacked,rowwise);
         try_success("getRowPackedCopy...f ","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"getRowPackedCopy...f ","data not successfully (deep) copied by rows");
      }
      try {
         tmp = B.get(B.getRowDimension(),B.getColumnDimension()-1);
         errorCount = try_failure(errorCount,"get(int,int)...f ","OutOfBoundsException expected but not thrown");
      } catch ( java.lang.ArrayIndexOutOfBoundsException e ) {
         try {
            tmp = B.get(B.getRowDimension()-1,B.getColumnDimension());
            errorCount = try_failure(errorCount,"get(int,int)...f ","OutOfBoundsException expected but not thrown");
         } catch ( java.lang.ArrayIndexOutOfBoundsException e1 ) {
            try_success("get(int,int)...f OutofBoundsException...f ","");
         }
      } catch ( java.lang.IllegalArgumentException e1 ) {
         errorCount = try_failure(errorCount,"get(int,int)...f ","OutOfBoundsException expected but not thrown");
      }
      try {
         if (B.get(B.getRowDimension()-1,B.getColumnDimension()-1) != 
             avals[B.getRowDimension()-1][B.getColumnDimension()-1] ) {
            errorCount = try_failure(errorCount,"get(int,int)...f ","Matrix entry (i,j) not successfully retreived");
         } else {
            try_success("get(int,int)...f ","");
         }
      } catch ( java.lang.ArrayIndexOutOfBoundsException e ) {
         errorCount = try_failure(errorCount,"get(int,int)...f ","Unexpected ArrayIndexOutOfBoundsException");
      }
      SUB = new Matrix(subavals);
      try {
         M = B.getMatrix(ib,ie+B.getRowDimension()+1,jb,je);
         errorCount = try_failure(errorCount,"getMatrix(int,int,int,int)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
      } catch ( java.lang.ArrayIndexOutOfBoundsException e ) {
         try {
            M = B.getMatrix(ib,ie,jb,je+B.getColumnDimension()+1);
            errorCount = try_failure(errorCount,"getMatrix(int,int,int,int)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
         } catch ( java.lang.ArrayIndexOutOfBoundsException e1 ) {
            try_success("getMatrix(int,int,int,int)...f ArrayIndexOutOfBoundsException...f ","");
         }
      } catch ( java.lang.IllegalArgumentException e1 ) {
         errorCount = try_failure(errorCount,"getMatrix(int,int,int,int)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
      }
      try {
         M = B.getMatrix(ib,ie,jb,je);
         try {
            check(SUB,M);
            try_success("getMatrix(int,int,int,int)...f ","");
         } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"getMatrix(int,int,int,int)...f ","submatrix not successfully retreived");
         }
      } catch ( java.lang.ArrayIndexOutOfBoundsException e ) {
         errorCount = try_failure(errorCount,"getMatrix(int,int,int,int)...f ","Unexpected ArrayIndexOutOfBoundsException");
      }
      
      try {
         M = B.getMatrix(ib,ie,badcolumnindexset);
         errorCount = try_failure(errorCount,"getMatrix(int,int,int[])...f ","ArrayIndexOutOfBoundsException expected but not thrown");
      } catch ( java.lang.ArrayIndexOutOfBoundsException e ) {
         try {
            M = B.getMatrix(ib,ie+B.getRowDimension()+1,columnindexset);
            errorCount = try_failure(errorCount,"getMatrix(int,int,int[])...f ","ArrayIndexOutOfBoundsException expected but not thrown");
         } catch ( java.lang.ArrayIndexOutOfBoundsException e1 ) {
            try_success("getMatrix(int,int,int[])...f ArrayIndexOutOfBoundsException...f ","");
         }
      } catch ( java.lang.IllegalArgumentException e1 ) {
         errorCount = try_failure(errorCount,"getMatrix(int,int,int[])...f ","ArrayIndexOutOfBoundsException expected but not thrown");
      } 
      try {
         M = B.getMatrix(ib,ie,columnindexset);
         try {
            check(SUB,M);
            try_success("getMatrix(int,int,int[])...f ","");
         } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"getMatrix(int,int,int[])...f ","submatrix not successfully retreived");
         }
      } catch ( java.lang.ArrayIndexOutOfBoundsException e ) {
         errorCount = try_failure(errorCount,"getMatrix(int,int,int[])...f ","Unexpected ArrayIndexOutOfBoundsException");
      }
      try {
         M = B.getMatrix(badrowindexset,jb,je);
         errorCount = try_failure(errorCount,"getMatrix(int[],int,int)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
      } catch ( java.lang.ArrayIndexOutOfBoundsException e ) {
         try {
            M = B.getMatrix(rowindexset,jb,je+B.getColumnDimension()+1);
            errorCount = try_failure(errorCount,"getMatrix(int[],int,int)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
         } catch ( java.lang.ArrayIndexOutOfBoundsException e1 ) {
            try_success("getMatrix(int[],int,int)...f ArrayIndexOutOfBoundsException...f ","");
         }
      } catch ( java.lang.IllegalArgumentException e1 ) {
         errorCount = try_failure(errorCount,"getMatrix(int[],int,int)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
      } 
      try {
         M = B.getMatrix(rowindexset,jb,je);
         try {
            check(SUB,M);
            try_success("getMatrix(int[],int,int)...f ","");
         } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"getMatrix(int[],int,int)...f ","submatrix not successfully retreived");
         }
      } catch ( java.lang.ArrayIndexOutOfBoundsException e ) {
         errorCount = try_failure(errorCount,"getMatrix(int[],int,int)...f ","Unexpected ArrayIndexOutOfBoundsException");
      }
      try {
         M = B.getMatrix(badrowindexset,columnindexset);
         errorCount = try_failure(errorCount,"getMatrix(int[],int[])...f ","ArrayIndexOutOfBoundsException expected but not thrown");
      } catch ( java.lang.ArrayIndexOutOfBoundsException e ) {
         try {
            M = B.getMatrix(rowindexset,badcolumnindexset);
            errorCount = try_failure(errorCount,"getMatrix(int[],int[])...f ","ArrayIndexOutOfBoundsException expected but not thrown");
         } catch ( java.lang.ArrayIndexOutOfBoundsException e1 ) {
            try_success("getMatrix(int[],int[])...f ArrayIndexOutOfBoundsException...f ","");
         }
      } catch ( java.lang.IllegalArgumentException e1 ) {
         errorCount = try_failure(errorCount,"getMatrix(int[],int[])...f ","ArrayIndexOutOfBoundsException expected but not thrown");
      } 
      try {
         M = B.getMatrix(rowindexset,columnindexset);
         try {
            check(SUB,M);
            try_success("getMatrix(int[],int[])...f ","");
         } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"getMatrix(int[],int[])...f ","submatrix not successfully retreived");
         }
      } catch ( java.lang.ArrayIndexOutOfBoundsException e ) {
         errorCount = try_failure(errorCount,"getMatrix(int[],int[])...f ","Unexpected ArrayIndexOutOfBoundsException");
      }

/**
      Various set methods:
**/

      try {
         B.set(B.getRowDimension(),B.getColumnDimension()-1,0.f);
         errorCount = try_failure(errorCount,"set(int,int,float)...f ","OutOfBoundsException expected but not thrown");
      } catch ( java.lang.ArrayIndexOutOfBoundsException e ) {
         try {
            B.set(B.getRowDimension()-1,B.getColumnDimension(),0.f);
            errorCount = try_failure(errorCount,"set(int,int,float)...f ","OutOfBoundsException expected but not thrown");
         } catch ( java.lang.ArrayIndexOutOfBoundsException e1 ) {
            try_success("set(int,int,float)...f OutofBoundsException...f ","");
         }
      } catch ( java.lang.IllegalArgumentException e1 ) {
         errorCount = try_failure(errorCount,"set(int,int,float)...f ","OutOfBoundsException expected but not thrown");
      }
      try {
         B.set(ib,jb,0.f);
         tmp = B.get(ib,jb);
         try {
            check(tmp,0.f);
            try_success("set(int,int,float)...f ","");
         } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"set(int,int,float)...f ","Matrix element not successfully set");
         }
      } catch ( java.lang.ArrayIndexOutOfBoundsException e1) {
         errorCount = try_failure(errorCount,"set(int,int,float)...f ","Unexpected ArrayIndexOutOfBoundsException");
      }
      M = new Matrix(2,3,0.f);
      try {
         B.setMatrix(ib,ie+B.getRowDimension()+1,jb,je,M);
         errorCount = try_failure(errorCount,"setMatrix(int,int,int,int,Matrix)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
      } catch ( java.lang.ArrayIndexOutOfBoundsException e ) {
         try {
            B.setMatrix(ib,ie,jb,je+B.getColumnDimension()+1,M);
            errorCount = try_failure(errorCount,"setMatrix(int,int,int,int,Matrix)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
         } catch ( java.lang.ArrayIndexOutOfBoundsException e1 ) {
            try_success("setMatrix(int,int,int,int,Matrix)...f ArrayIndexOutOfBoundsException...f ","");
         }
      } catch ( java.lang.IllegalArgumentException e1 ) {
         errorCount = try_failure(errorCount,"setMatrix(int,int,int,int,Matrix)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
      }
      try {
         B.setMatrix(ib,ie,jb,je,M);
         try {
            check(M.minus(B.getMatrix(ib,ie,jb,je)),M);
            try_success("setMatrix(int,int,int,int,Matrix)...f ","");
         } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"setMatrix(int,int,int,int,Matrix)...f ","submatrix not successfully set");
         }
         B.setMatrix(ib,ie,jb,je,SUB);
      } catch ( java.lang.ArrayIndexOutOfBoundsException e1 ) {
         errorCount = try_failure(errorCount,"setMatrix(int,int,int,int,Matrix)...f ","Unexpected ArrayIndexOutOfBoundsException");
      }
      try {
         B.setMatrix(ib,ie+B.getRowDimension()+1,columnindexset,M);
         errorCount = try_failure(errorCount,"setMatrix(int,int,int[],Matrix)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
      } catch ( java.lang.ArrayIndexOutOfBoundsException e ) {
         try {
            B.setMatrix(ib,ie,badcolumnindexset,M);
            errorCount = try_failure(errorCount,"setMatrix(int,int,int[],Matrix)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
         } catch ( java.lang.ArrayIndexOutOfBoundsException e1 ) {
            try_success("setMatrix(int,int,int[],Matrix)...f ArrayIndexOutOfBoundsException...f ","");
         }
      } catch ( java.lang.IllegalArgumentException e1 ) {
         errorCount = try_failure(errorCount,"setMatrix(int,int,int[],Matrix)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
      }
      try {
         B.setMatrix(ib,ie,columnindexset,M);
         try {
            check(M.minus(B.getMatrix(ib,ie,columnindexset)),M);
            try_success("setMatrix(int,int,int[],Matrix)...f ","");
         } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"setMatrix(int,int,int[],Matrix)...f ","submatrix not successfully set");
         }
         B.setMatrix(ib,ie,jb,je,SUB);
      } catch ( java.lang.ArrayIndexOutOfBoundsException e1 ) {
         errorCount = try_failure(errorCount,"setMatrix(int,int,int[],Matrix)...f ","Unexpected ArrayIndexOutOfBoundsException");
      }
      try {
         B.setMatrix(rowindexset,jb,je+B.getColumnDimension()+1,M);
         errorCount = try_failure(errorCount,"setMatrix(int[],int,int,Matrix)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
      } catch ( java.lang.ArrayIndexOutOfBoundsException e ) {
         try {
            B.setMatrix(badrowindexset,jb,je,M);
            errorCount = try_failure(errorCount,"setMatrix(int[],int,int,Matrix)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
         } catch ( java.lang.ArrayIndexOutOfBoundsException e1 ) {
            try_success("setMatrix(int[],int,int,Matrix)...f ArrayIndexOutOfBoundsException...f ","");
         }
      } catch ( java.lang.IllegalArgumentException e1 ) {
         errorCount = try_failure(errorCount,"setMatrix(int[],int,int,Matrix)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
      }
      try {
         B.setMatrix(rowindexset,jb,je,M);
         try {
            check(M.minus(B.getMatrix(rowindexset,jb,je)),M);
            try_success("setMatrix(int[],int,int,Matrix)...f ","");
         } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"setMatrix(int[],int,int,Matrix)...f ","submatrix not successfully set");
         }
         B.setMatrix(ib,ie,jb,je,SUB);
      } catch ( java.lang.ArrayIndexOutOfBoundsException e1 ) {
         errorCount = try_failure(errorCount,"setMatrix(int[],int,int,Matrix)...f ","Unexpected ArrayIndexOutOfBoundsException");
      }
      try {
         B.setMatrix(rowindexset,badcolumnindexset,M);
         errorCount = try_failure(errorCount,"setMatrix(int[],int[],Matrix)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
      } catch ( java.lang.ArrayIndexOutOfBoundsException e ) {
         try {
            B.setMatrix(badrowindexset,columnindexset,M);
            errorCount = try_failure(errorCount,"setMatrix(int[],int[],Matrix)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
         } catch ( java.lang.ArrayIndexOutOfBoundsException e1 ) {
            try_success("setMatrix(int[],int[],Matrix)...f ArrayIndexOutOfBoundsException...f ","");
         }
      } catch ( java.lang.IllegalArgumentException e1 ) {
         errorCount = try_failure(errorCount,"setMatrix(int[],int[],Matrix)...f ","ArrayIndexOutOfBoundsException expected but not thrown");
      }
      try {
         B.setMatrix(rowindexset,columnindexset,M);
         try {
            check(M.minus(B.getMatrix(rowindexset,columnindexset)),M);
            try_success("setMatrix(int[],int[],Matrix)...f ","");
         } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"setMatrix(int[],int[],Matrix)...f ","submatrix not successfully set");
         }
      } catch ( java.lang.ArrayIndexOutOfBoundsException e1 ) {
         errorCount = try_failure(errorCount,"setMatrix(int[],int[],Matrix)...f ","Unexpected ArrayIndexOutOfBoundsException");
      }

/** 
      Array-like methods:
         minus
         minusEquals
         plus
         plusEquals
         arrayLeftDivide
         arrayLeftDivideEquals
         arrayRightDivide
         arrayRightDivideEquals
         arrayTimes
         arrayTimesEquals
         uminus
**/

      print("\nTesting array-like methods...f\n");
      S = new Matrix(columnwise,nonconformld);
      R = Matrix.random(A.getRowDimension(),A.getColumnDimension());
      A = R;
      try {
        S = A.minus(S);
        errorCount = try_failure(errorCount,"minus conformance check...f ","nonconformance not raised");
      } catch ( IllegalArgumentException e ) {
        try_success("minus conformance check...f ","");
      }
      if (A.minus(R).norm1() != 0.f) {
        errorCount = try_failure(errorCount,"minus...f ","(difference of identical Matrices is nonzero,\nSubsequent use of minus should be suspect)");
      } else {
        try_success("minus...f ","");
      }
      A = R.copy();
      A.minusEquals(R);
      Z = new Matrix(A.getRowDimension(),A.getColumnDimension());
      try {
        A.minusEquals(S);
        errorCount = try_failure(errorCount,"minusEquals conformance check...f ","nonconformance not raised");
      } catch ( IllegalArgumentException e ) {
        try_success("minusEquals conformance check...f ","");
      }
      if (A.minus(Z).norm1() != 0.f) {
        errorCount = try_failure(errorCount,"minusEquals...f ","(difference of identical Matrices is nonzero,\nSubsequent use of minus should be suspect)");
      } else {
        try_success("minusEquals...f ","");
      }

      A = R.copy();
      B = Matrix.random(A.getRowDimension(),A.getColumnDimension());
      C = A.minus(B); 
      try {
        S = A.plus(S);
        errorCount = try_failure(errorCount,"plus conformance check...f ","nonconformance not raised");
      } catch ( IllegalArgumentException e ) {
        try_success("plus conformance check...f ","");
      }
      try {
        check(C.plus(B),A);
        try_success("plus...f ","");
      } catch ( java.lang.RuntimeException e ) {
        errorCount = try_failure(errorCount,"plus...f ","(C = A - B, but C + B != A)");
      }
      C = A.minus(B);
      C.plusEquals(B);
      try {
        A.plusEquals(S);
        errorCount = try_failure(errorCount,"plusEquals conformance check...f ","nonconformance not raised");
      } catch ( IllegalArgumentException e ) {
        try_success("plusEquals conformance check...f ","");
      }
      try {
        check(C,A);
        try_success("plusEquals...f ","");
      } catch ( java.lang.RuntimeException e ) {
        errorCount = try_failure(errorCount,"plusEquals...f ","(C = A - B, but C = C + B != A)");
      }
      A = R.uminus();
      try {
        check(A.plus(R),Z);
        try_success("uminus...f ","");
      } catch ( java.lang.RuntimeException e ) {
        errorCount = try_failure(errorCount,"uminus...f ","(-A + A != zeros)");
      }
      A = R.copy();
      O = new Matrix(A.getRowDimension(),A.getColumnDimension(),1.0f);
      C = A.arrayLeftDivide(R);
      try {
        S = A.arrayLeftDivide(S);
        errorCount = try_failure(errorCount,"arrayLeftDivide conformance check...f ","nonconformance not raised");
      } catch ( IllegalArgumentException e ) {
        try_success("arrayLeftDivide conformance check...f ","");
      }
      try {
        check(C,O);
        try_success("arrayLeftDivide...f ","");
      } catch ( java.lang.RuntimeException e ) {
        errorCount = try_failure(errorCount,"arrayLeftDivide...f ","(M.f\\M != ones)");
      }
      try {
        A.arrayLeftDivideEquals(S);
        errorCount = try_failure(errorCount,"arrayLeftDivideEquals conformance check...f ","nonconformance not raised");
      } catch ( IllegalArgumentException e ) {
        try_success("arrayLeftDivideEquals conformance check...f ","");
      }
      A.arrayLeftDivideEquals(R);
      try {
        check(A,O);
        try_success("arrayLeftDivideEquals...f ","");
      } catch ( java.lang.RuntimeException e ) {
        errorCount = try_failure(errorCount,"arrayLeftDivideEquals...f ","(M.f\\M != ones)");
      }
      A = R.copy();
      try {
        A.arrayRightDivide(S);
        errorCount = try_failure(errorCount,"arrayRightDivide conformance check...f ","nonconformance not raised");
      } catch ( IllegalArgumentException e ) {
        try_success("arrayRightDivide conformance check...f ","");
      }
      C = A.arrayRightDivide(R);
      try {
        check(C,O);
        try_success("arrayRightDivide...f ","");
      } catch ( java.lang.RuntimeException e ) {
        errorCount = try_failure(errorCount,"arrayRightDivide...f ","(M./M != ones)");
      }
      try {
        A.arrayRightDivideEquals(S);
        errorCount = try_failure(errorCount,"arrayRightDivideEquals conformance check...f ","nonconformance not raised");
      } catch ( IllegalArgumentException e ) {
        try_success("arrayRightDivideEquals conformance check...f ","");
      }
      A.arrayRightDivideEquals(R);
      try {
        check(A,O);
        try_success("arrayRightDivideEquals...f ","");
      } catch ( java.lang.RuntimeException e ) {
        errorCount = try_failure(errorCount,"arrayRightDivideEquals...f ","(M./M != ones)");
      }
      A = R.copy();
      B = Matrix.random(A.getRowDimension(),A.getColumnDimension());
      try {
        S = A.arrayTimes(S);
        errorCount = try_failure(errorCount,"arrayTimes conformance check...f ","nonconformance not raised");
      } catch ( IllegalArgumentException e ) {
        try_success("arrayTimes conformance check...f ","");
      }
      C = A.arrayTimes(B);
      try {
        check(C.arrayRightDivideEquals(B),A);
        try_success("arrayTimes...f ","");
      } catch ( java.lang.RuntimeException e ) {
        errorCount = try_failure(errorCount,"arrayTimes...f ","(A = R, C = A.*B, but C./B != A)");
      }
      try {
        A.arrayTimesEquals(S);
        errorCount = try_failure(errorCount,"arrayTimesEquals conformance check...f ","nonconformance not raised");
      } catch ( IllegalArgumentException e ) {
        try_success("arrayTimesEquals conformance check...f ","");
      }
      A.arrayTimesEquals(B);
      try {
        check(A.arrayRightDivideEquals(B),R);
        try_success("arrayTimesEquals...f ","");
      } catch ( java.lang.RuntimeException e ) {
        errorCount = try_failure(errorCount,"arrayTimesEquals...f ","(A = R, A = A.*B, but A./B != R)");
      }

/**   
      I/O methods:
         read
         print
         serializable:
           writeObject
           readObject
**/
      print("\nTesting I/O methods...f\n");
         try {
            DecimalFormat fmt = new DecimalFormat("0.0000E00");
	    fmt.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));

            PrintWriter FILE = new PrintWriter(new FileOutputStream("JamaTestMatrix.out"));
            A.print(FILE,fmt,10);
            FILE.close();
            R = Matrix.read(new BufferedReader(new FileReader("JamaTestMatrix.out")));
            if (A.minus(R).norm1() < .001f ) {
               try_success("print()/read()...","");
            } else {
               errorCount = try_failure(errorCount,"print()/read()...","Matrix read from file does not match Matrix printed to file");
            }
         } catch ( java.io.IOException ioe ) {
           warningCount = try_warning(warningCount,"print()/read()...","unexpected I/O error, unable to run print/read test;  check write permission in current directory and retry");
         } catch(Exception e) {
            try {
               e.printStackTrace(System.out);
               warningCount = try_warning(warningCount,"print()/read()...","Formatting error...f will try JDK1.1f reformulation...");
               DecimalFormat fmt = new DecimalFormat("0.0000");
               PrintWriter FILE = new PrintWriter(new FileOutputStream("JamaTestMatrix.out"));
               A.print(FILE,fmt,10);
               FILE.close();
               R = Matrix.read(new BufferedReader(new FileReader("JamaTestMatrix.out")));
               if (A.minus(R).norm1() < .001f ) {
                  try_success("print()/read()...","");
               } else {
                  errorCount = try_failure(errorCount,"print()/read() (2nd attempt) ...","Matrix read from file does not match Matrix printed to file");
               }
            } catch ( java.io.IOException ioe ) {
              warningCount = try_warning(warningCount,"print()/read()...","unexpected I/O error, unable to run print/read test;  check write permission in current directory and retry");
         }
      }

      R = Matrix.random(A.getRowDimension(),A.getColumnDimension());
      String tmpname = "TMPMATRIX.serial";
      try {
         ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(tmpname));
         out.writeObject(R);
         ObjectInputStream sin = new ObjectInputStream(new FileInputStream(tmpname));
         A = (Matrix) sin.readObject();
 
         try {
            check(A,R);
            try_success("writeObject(Matrix)/readObject(Matrix)...","");
         } catch ( java.lang.RuntimeException e ) {
           errorCount = try_failure(errorCount,"writeObject(Matrix)/readObject(Matrix)...","Matrix not serialized correctly");
         }
      } catch ( java.io.IOException ioe ) {
         warningCount = try_warning(warningCount,"writeObject()/readObject()...","unexpected I/O error, unable to run serialization test;  check write permission in current directory and retry");
      } catch(Exception e) {
         errorCount = try_failure(errorCount,"writeObject(Matrix)/readObject(Matrix)...","unexpected error in serialization test");
      }

/**
      LA methods:
         transpose
         times
         cond
         rank
         det
         trace
         norm1
         norm2
         normF
         normInf
         solve
         solveTranspose
         inverse
         chol
         eig
         lu
         qr
         svd 
**/

      print("\nTesting linear algebra methods...f\n");
      A = new Matrix(columnwise,3);
      T = new Matrix(tvals);
      T = A.transpose();
      try {
         check(A.transpose(),T);
         try_success("transpose...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"transpose()...","transpose unsuccessful");
      }
      A.transpose();
      try {
         check(A.norm1(),columnsummax);
         try_success("norm1...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"norm1()...","incorrect norm calculation");
      }
      try {
         check(A.normInf(),rowsummax);
         try_success("normInf()...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"normInf()...","incorrect norm calculation");
      }
      try {
         check(A.normF(),(float) Math.sqrt(sumofsquares));
         try_success("normF...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"normF()...","incorrect norm calculation");
      }
      try {
         check(A.trace(),sumofdiagonals);
         try_success("trace()...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"trace()...","incorrect trace calculation");
      }
      try {
         check(A.getMatrix(0,A.getRowDimension()-1,0,A.getRowDimension()-1).det(),0.f);
         try_success("det()...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"det()...","incorrect determinant calculation");
      }
      SQ = new Matrix(square);
      try {
         check(A.times(A.transpose()),SQ);
         try_success("times(Matrix)...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"times(Matrix)...","incorrect Matrix-Matrix product calculation");
      }
      try {
         check(A.times(0.f),Z);
         try_success("times(float)...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"times(float)...","incorrect Matrix-scalar product calculation");
      }

      A = new Matrix(columnwise,4);
      QRDecomposition QR = A.qr();
      R = QR.getR();
      try {
         check(A,QR.getQ().times(R));
         try_success("QRDecomposition...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"QRDecomposition...","incorrect QR decomposition calculation");
      }
      SingularValueDecomposition SVD = A.svd();
      try {
         check(A,SVD.getU().times(SVD.getS().times(SVD.getV().transpose())));
         try_success("SingularValueDecomposition...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"SingularValueDecomposition...","incorrect singular value decomposition calculation");
      }
      DEF = new Matrix(rankdef);
      try {
         check(DEF.rank(),Math.min(DEF.getRowDimension(),DEF.getColumnDimension())-1);
         try_success("rank()...","");
      } catch ( java.lang.RuntimeException e ) {
	      System.out.println(e);
         errorCount = try_failure(errorCount,"rank()...","incorrect rank calculation");
      }
      B = new Matrix(condmat);
      SVD = B.svd(); 
      float [] singularvalues = SVD.getSingularValues();
      try {
         check(B.cond(),singularvalues[0]/singularvalues[Math.min(B.getRowDimension(),B.getColumnDimension())-1]);
         try_success("cond()...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"cond()...","incorrect condition number calculation");
      }
      int n = A.getColumnDimension();
      A = A.getMatrix(0,n-1,0,n-1);
      A.set(0,0,0.f);
      LUDecomposition LU = A.lu();
      try {
         check(A.getMatrix(LU.getPivot(),0,n-1),LU.getL().times(LU.getU()));
         try_success("LUDecomposition...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"LUDecomposition...","incorrect LU decomposition calculation");
      }
      X = A.inverse();
      try {
         check(A.times(X),Matrix.identity(3,3));
         try_success("inverse()...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"inverse()...","incorrect inverse calculation");
      }
      O = new Matrix(SUB.getRowDimension(),1,1.0f);
      SOL = new Matrix(sqSolution);
      SQ = SUB.getMatrix(0,SUB.getRowDimension()-1,0,SUB.getRowDimension()-1);
      try {
         check(SQ.solve(SOL),O); 
         try_success("solve()...","");
      } catch ( java.lang.IllegalArgumentException e1 ) {
         errorCount = try_failure(errorCount,"solve()...",e1.getMessage());
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"solve()...",e.getMessage());
      }
      A = new Matrix(pvals);
      CholeskyDecomposition Chol = A.chol(); 
      Matrix L = Chol.getL();
      try {
         check(A,L.times(L.transpose()));
         try_success("CholeskyDecomposition...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"CholeskyDecomposition...","incorrect Cholesky decomposition calculation");
      }
      X = Chol.solve(Matrix.identity(3,3));
      try {
         check(A.times(X),Matrix.identity(3,3));
         try_success("CholeskyDecomposition solve()...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"CholeskyDecomposition solve()...","incorrect Choleskydecomposition solve calculation");
      }
      EigenvalueDecomposition Eig = A.eig();
      Matrix D = Eig.getD();
      Matrix V = Eig.getV();
      try {
         check(A.times(V),V.times(D));
         try_success("EigenvalueDecomposition (symmetric)...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"EigenvalueDecomposition (symmetric)...","incorrect symmetric Eigenvalue decomposition calculation");
      }
      A = new Matrix(evals);
      Eig = A.eig();
      D = Eig.getD();
      V = Eig.getV();
      try {
         check(A.times(V),V.times(D));
         try_success("EigenvalueDecomposition (nonsymmetric)...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"EigenvalueDecomposition (nonsymmetric)...","incorrect nonsymmetric Eigenvalue decomposition calculation");
      }

      try {
	  print("\nTesting Eigenvalue; If this hangs, we've failed\n");
	  Matrix bA = new Matrix(badeigs);
	  EigenvalueDecomposition bEig = bA.eig();
	  try_success("EigenvalueDecomposition (hang)...","");
      } catch ( java.lang.RuntimeException e ) {
         errorCount = try_failure(errorCount,"EigenvalueDecomposition (hang)...",
				  "incorrect termination");
      }


      print("\nTestMatrix completed.f\n");
      print("Total errors reported: " + Integer.toString(errorCount) + "\n");
      print("Total warnings reported: " + Integer.toString(warningCount) + "\n");
   }

   /** private utility routines **/

   /** Check magnitude of difference of scalars.f **/

   private static void check(float x, float y) {
      float eps = (float) Math.pow(2.0f,-24.0f);
      if (x == 0 & Math.abs(y) < 10*eps) return;
      if (y == 0 & Math.abs(x) < 10*eps) return;
      if (Math.abs(x-y) > 10*eps*Math.max(Math.abs(x),Math.abs(y))) {
         throw new RuntimeException("The difference x-y is too large: x = " + Float.toString(x) + "  y = " + Float.toString(y));
      }
   }

   /** Check norm of difference of "vectors".f **/

   private static void check(float[] x, float[] y) {
      if (x.length == y.length ) {
         for (int i=0;i<x.length;i++) {
            check(x[i],y[i]);
         } 
      } else {
         throw new RuntimeException("Attempt to compare vectors of different lengths");
      }
   }

   /** Check norm of difference of arrays.f **/

   private static void check(float[][] x, float[][] y) {
      Matrix A = new Matrix(x);
      Matrix B = new Matrix(y);
      check(A,B);
   }

   /** Check norm of difference of Matrices.f **/

   private static void check(Matrix X, Matrix Y) {
      float eps = (float) Math.pow(2.0f,-24.0f);
      if (X.norm1() == 0.f & Y.norm1() < 10*eps) return;
      if (Y.norm1() == 0.f & X.norm1() < 10*eps) return;
      if (X.minus(Y).norm1() > 1000*eps*Math.max(X.norm1(),Y.norm1())) {
         throw new RuntimeException("The norm of (X-Y) is too large: " +  Float.toString(X.minus(Y).norm1()));
      }
   }

   /** Shorten spelling of print.f **/

   private static void print (String s) {
      System.out.print(s);
   }

  /** Print appropriate messages for successful outcome try **/

   private static void try_success (String s,String e) {
      print(">    " + s + "success\n");
      if ( e != "" ) {
        print(">      Message: " + e + "\n");
      }
   }
  /** Print appropriate messages for unsuccessful outcome try **/

   private static int try_failure (int count, String s,String e) {
      print(">    " + s + "*** failure ***\n>      Message: " + e + "\n");
      return ++count;
   }

  /** Print appropriate messages for unsuccessful outcome try **/

   private static int try_warning (int count, String s,String e) {
      print(">    " + s + "*** warning ***\n>      Message: " + e + "\n");
      return ++count;
   }

   /** Print a row vector.f **/

   private static void print(float[] x, int w, int d) {
      // Use format Fw.d for all elements.
      System.out.print("\n");
      new Matrix(x,1).print(w,d);
      print("\n");
   }

}
