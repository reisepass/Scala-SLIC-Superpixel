import ij._
import ij.io.Opener
import breeze.linalg._
import breeze.stats.DescriptiveStats._
import breeze.stats._
import breeze.numerics._
import ij.plugin.Duplicator


object scratch {
   def main(args: Array[String]): Unit = {

    val tft0 = System.currentTimeMillis()
    // val examplePath = "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/testing_groundtruth.tif"
    //val examplePath = "/home/mort/workspace/back_ScalaSLIC/data/training.tif"
    val examplePath = "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/training_80_84_165.tif"
    // val examplePath = "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/training_48_40_40.tif"
    //val examplePath = "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/training_9_9_10.tif"
    // val examplePath = "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/training_saveAsimgJ.tif"
    val examplePa3 = "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/testing_groundtruth_first5.tif"
    val examplePath2 = "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/training.tif"
   val S = 20 // Superpixel center grid interval 

    val opener = new Opener();
    val imp2 = opener.openImage(examplePath);

    val tft1 = System.currentTimeMillis()
    println((tft1 - tft0) + " OpenImage() time")
    //imp2.show() 

    val aStack = imp2.getStack

    val v1 = aStack.getVoxel(1, 1, 1)
    val v2 = aStack.getVoxel(1, 2, 1)
    val v3 = aStack.getVoxel(1, 1, 2)

    val zsize = aStack.getSize //Zsize 
    val allTheD = imp2.getDimensions
    val xDim = allTheD(0)
    val yDim = allTheD(1)
    val zDim = allTheD(3)
    println("(" + xDim + "," + yDim + "," + zDim + ") image Dims")
    assert(xDim > 0 & yDim > 0 & zDim > 0)
    
    
    /*
    val img3 = opener.openImage("/home/mort/workspace/dissolve-struct/data/generated/MSRC_ObjCategImageDatabase_v2/Images1_10_s.bmp")
    val stack3 = img3.getStack()
    val x2 = stack3.getWidth
    val y2 = stack3.getHeight
    val z2 = stack3.getSize
    */
    
    println("fun times")
    
    
    val tFu = System.currentTimeMillis()
    val copyImg = Array.fill(xDim,yDim,zDim)(0.0)
    for(x <- 0 until xDim; y <- 0 until yDim; z <- 0 until zDim){
      copyImg(x)(y)(z)=aStack.getVoxel(x,y,z)
    }
    println(" copying image over took "+(System.currentTimeMillis()-tFu))
    
    /*
     * 
     * Test Some More 
     * 
     * 
     */
  
         
      def printSuperPixels(clusterAssign:Array[Array[Array[Int]]],imp2:ImagePlus, borderCol: Double = 300.0, label: String = "NoName") {
      val imp2_pix = new Duplicator().run(imp2);
      val aStack_supPix = imp2_pix.getStack

      for (vX <- 0 until xDim) {
        for (vY <- 0 until yDim) {
          var lastLabel = -1
          for (vZ <- 0 until zDim) {
            if (clusterAssign(vX)(vY)(vZ) != lastLabel) {
              aStack_supPix.setVoxel(vX, vY, vZ, borderCol)
            }
            lastLabel = clusterAssign(vX)(vY)(vZ)
          }
        }
      }

      for (vX <- 0 until xDim) {
        for (vZ <- 0 until zDim) {
          var lastLabel = -1
          for (vY <- 0 until yDim) {

            if (clusterAssign(vX)(vY)(vZ) != lastLabel) {
              aStack_supPix.setVoxel(vX, vY, vZ, borderCol)
            }
            lastLabel = clusterAssign(vX)(vY)(vZ)
          }
        }
      }

      for (vZ <- 0 until zDim) {

        for (vY <- 0 until yDim) {
          var lastLabel = -1
          for (vX <- 0 until xDim) {
            if (clusterAssign(vX)(vY)(vZ) != lastLabel) {
              aStack_supPix.setVoxel(vX, vY, vZ, borderCol)
            }
            lastLabel = clusterAssign(vX)(vY)(vZ)
          }
        }
      }
      imp2_pix.setStack(aStack_supPix)
      //imp2_pix.show()
      IJ.saveAs(imp2_pix, "tif", "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/image" + label + ".tif");
      //TODO free this memory after saving 
    }
    
    
    
    val distFn = (a:Double,b:Double)=>sqrt(Math.pow(a-b,2))
    val sumFn  = (a:Double,b:Double)=>(a+b)
    val normFn = (a:Double,b:Int)   =>(a/b)
    
     
    
    
    
    val allIm =  new SLIC[Double](distFn,sumFn,normFn,copyImg,20,15,minChangePerIter=0.002,connectivityOption="Imperative",debug=false)
     printSuperPixels(allIm.calcSuperPixels(),imp2,300.0,"Imperative")
     
     val allFn =  new SLIC[Double](distFn,sumFn,normFn,copyImg,20,15,connectivityOption="Functional",debug=false)
     printSuperPixels(allFn.calcSuperPixels(),imp2,300.0,"_Functional")
     
     val allno = new SLIC[Double](distFn,sumFn,normFn,copyImg,20,15,connectivityOption="None",debug=false)
     printSuperPixels(allno.calcSuperPixels(),imp2,300.0,label="_ConnectivityNotEnforced")

     
      val alld=  new SLIC[Double](distFn,sumFn,normFn,copyImg,20,15,minChangePerIter=0.002,connectivityOption="Functional2",debug=true)
     printSuperPixels(alld.calcSuperPixels(),imp2,300.0,"debugfn2")
     
   }
}