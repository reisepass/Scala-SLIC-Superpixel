import ij._
import ij.io.Opener
import breeze.linalg._
import breeze.stats.DescriptiveStats._
import breeze.stats._
import breeze.numerics._
import breeze.linalg._
import ij.plugin.Duplicator
import scala.util.matching.Regex
import java.io.File
import ch.ethz.dalab.scalaslic.SLIC
import ch.ethz.dalab.scalaslic.DatumCord

object scratch2 {

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
    IJ.saveAs(imp2_pix, "bmp", "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/image" + label + ".bmp");
    //TODO free this memory after saving 
  }

  def main(args: Array[String]): Unit = {

    val zGrid = (floor(1 / 2) until 1 by 3).toList
    val xGrid = (round(1 / 2) until 1 by 3).toList

    val path = "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/original"
    val lostofBMP = Option(new File(path).list).map(_.filter(_.endsWith(".bmp"))).get

    lostofBMP.foreach { fName =>
      {

        val examplePath = "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/original/" + fName
        val opener = new Opener();
        val img = opener.openImage(examplePath);
        val aStack = img.getStack
        val bitDep = aStack.getBitDepth()
        val colMod = aStack.getColorModel()
        val xDim = aStack.getWidth
        val yDim = aStack.getHeight
        val zDim = aStack.getSize

        val copyImage = Array.fill(xDim, yDim, zDim) { (0, 0, 0) }
        for (x <- 0 until xDim; y <- 0 until yDim; z <- 0 until zDim) {
          val curV = aStack.getVoxel(x, y, z).asInstanceOf[Int]
          copyImage(x)(y)(z) = (colMod.getRed(curV), colMod.getGreen(curV), colMod.getBlue(curV))
        }

        val distFn = (a: (Int, Int, Int), b: (Int, Int, Int)) => sqrt(Math.pow(a._1 - b._1, 2) + Math.pow(a._2 - b._2, 2) + Math.pow(a._3 - b._3, 2))
        val sumFn = (a: (Int, Int, Int), b: (Int, Int, Int)) => ((a._1 + b._1, a._2 + b._2, a._3 + a._3))
        val normFn = (a: (Int, Int, Int), n: Int) => ((a._1 / n, a._2 / n, a._3 / n))

        val allGr = new SLIC[(Int, Int, Int)](distFn, sumFn, normFn, copyImage, 30, 15, minChangePerIter = 0.002, connectivityOption = "Imperative", debug = false)
        val mask = allGr.calcSuperPixels()

        val histBinsPerCol = 3
        val histWidt = 255 / histBinsPerCol
        val featureFn = (data: List[DatumCord[(Int, Int, Int)]]) => {
          val redHist = Array.fill(histBinsPerCol) { 0 }
          val greenHist = (Array.fill(histBinsPerCol) { 0 })
          val blueHist = (Array.fill(histBinsPerCol) { 0 })
          data.foreach(a => {
            redHist(min(histBinsPerCol - 1, (a.cont._1 / histWidt))) += 1
            greenHist(min(histBinsPerCol - 1, (a.cont._2 / histWidt))) += 1
            blueHist(min(histBinsPerCol - 1, (a.cont._3 / histWidt))) += 1

          })
          val all = List(redHist, greenHist, blueHist).flatten
          val mySum = all.sum
          Vector(all.map(a => (a.asInstanceOf[Double] / mySum)).toArray)
        }
        val (feat, ma) = allGr.prepareGraph(mask, featureFn)

        val allNo = new SLIC[(Int, Int, Int)](distFn, sumFn, normFn, copyImage, 30, 15, minChangePerIter = 0.002, connectivityOption = "None", debug = false)
        printSuperPixels(allNo.calcSuperPixels(), img, label = (fName.substring(0, fName.length - 4) + "None"))

        val allFn = new SLIC[(Int, Int, Int)](distFn, sumFn, normFn, copyImage, 30, 15, minChangePerIter = 0.002, connectivityOption = "Functional", debug = false)
        printSuperPixels(allFn.calcSuperPixels(), img, label = (fName.substring(0, fName.length - 4) + "Functional"))

        val allIm = new SLIC[(Int, Int, Int)](distFn, sumFn, normFn, copyImage, 30, 15, minChangePerIter = 0.002, connectivityOption = "Imperative", debug = false)
        printSuperPixels(allIm.calcSuperPixels(), img, label = (fName.substring(0, fName.length - 4) + "Imperative"))

      }
    }
  }
}