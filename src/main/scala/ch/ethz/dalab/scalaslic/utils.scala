package ch.ethz.dalab.scalaslic

import javax.imageio.ImageIO
import java.awt.image.BufferedImage
import java.io.File
import ij.IJ
import ij.ImagePlus
import ij.plugin.Duplicator
import ij.io.Opener

/**
 * @author mort
 */
object utils {

  val workingDir = new File(".")

 sys.process.Process(Seq("mkdir", "printout"), workingDir).! 

  
   def fileReadling(fileLoc: String): Array[Array[Array[Int]]] = {
    val opener = new Opener();
    val img = opener.openImage(fileLoc);
    val aStack = img.getStack
    //Note: The properties  aStack.getBitDepth() and aStack.getColorModel()  are useful for non gray-scale data 
    val xDim = aStack.getWidth
    val yDim = aStack.getHeight
    val zDim = aStack.getSize

    val copyImage = Array.fill(xDim, yDim, zDim) { 0 }
    for (x <- 0 until xDim; y <- 0 until yDim; z <- 0 until zDim) {
      copyImage(x)(y)(z) = aStack.getVoxel(x, y, z).asInstanceOf[Int] //Depending on color model this may cause inaccuracies, but here we prefer simple data types. 
    }

    copyImage

  }
  
  
  def printSuperPixels(clusterAssign: Array[Array[Array[Int]]], imp2: ImagePlus, borderCol: Double = Int.MaxValue.asInstanceOf[Double], label: String = "NoName") {

    val imp2_pix = new Duplicator().run(imp2);
    val aStack_supPix = imp2_pix.getStack

    val xDim = aStack_supPix.getWidth
    val yDim = aStack_supPix.getHeight
    val zDim = aStack_supPix.getSize

    if (zDim > 5) {
      for (vX <- 0 until xDim) {
        for (vY <- 0 until yDim) {
          var lastLabel = clusterAssign(vX)(vY)(0)
          for (vZ <- 0 until zDim) {
            if (clusterAssign(vX)(vY)(vZ) != lastLabel) {
              aStack_supPix.setVoxel(vX, vY, vZ, borderCol)
            }
            lastLabel = clusterAssign(vX)(vY)(vZ)
          }
        }
      }
    }

    for (vX <- 0 until xDim) {
      for (vZ <- 0 until zDim) {
        var lastLabel = clusterAssign(vX)(0)(vZ)
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
        var lastLabel = clusterAssign(0)(vY)(vZ)
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
    IJ.saveAs(imp2_pix, "tif", "printout/" + label + "_seg.tif"); //TODO this should not be hardcoded 
    //TODO free this memory after saving 
  }

}