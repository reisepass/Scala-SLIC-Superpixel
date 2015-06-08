import ij._
import ij.io.Opener

object scratch {
   def main(args: Array[String]): Unit = {

    val tft0 = System.currentTimeMillis()
    // val examplePath = "/home/mort/workspace/ScalaSLIC/data/testing_groundtruth.tif"
    val examplePath = "/home/mort/workspace/ScalaSLIC/data/training.tif"
    //val examplePath = "/home/mort/workspace/ScalaSLIC/data/training_80_84_165.tif"
    // val examplePath = "/home/mort/workspace/ScalaSLIC/data/training_48_40_40.tif"
    //val examplePath = "/home/mort/workspace/ScalaSLIC/data/training_9_9_10.tif"
    // val examplePath = "/home/mort/workspace/ScalaSLIC/data/training_saveAsimgJ.tif"
    val examplePa3 = "/home/mort/workspace/ScalaSLIC/data/testing_groundtruth_first5.tif"
    val examplePath2 = "/home/mort/workspace/ScalaSLIC/data/training.tif"
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
    
    
    
    val img3 = opener.openImage("/home/mort/workspace/dissolve-struct/data/generated/MSRC_ObjCategImageDatabase_v2/Images1_10_s.bmp")
    val stack3 = img3.getStack()
    val x2 = stack3.getWidth
    val y2 = stack3.getHeight
    val z2 = stack3.getSize
    
    println("fun times")
   }
}